#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <limits>

#ifndef HAVE_BAM_H
#error You do not have bam library
#endif
#include <bam.h>

//using namespace std;
#include <recrec.hxx>
#include <tktools.hxx>
#include <distsnp.hxx>
using std::string;
using std::vector;
using std::map;
using std::sort;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::stringstream;
using std::exception;
using std::invalid_argument;
using std::runtime_error;
using std::pair;
using std::atoi;
using std::hex;
using std::dec;
using std::endl;
using namespace tkbio;

using tktools::util::get_argument_string;
using tktools::util::get_argument_integer;
using tktools::util::get_argument_float;
using tktools::util::has_option;
using tktools::io::file_exists;
using tktools::io::get_file_size;
using tktools::split_items;
using tktools::bio::convert_chromosome_to_code;
using tktools::bio::convert_code_to_chromosome;

const char FILE_SEPARATOR = '/';

namespace {
    bool cache_position_comparator(const cache_position* lhs, const cache_position* rhs) {
        if (lhs->chromosome != rhs->chromosome) {
            return lhs->chromosome < rhs->chromosome;
        } else {
            return lhs->position < rhs->position;
        }
    }
}


//#define TEST 1
dbsnp_locus::dbsnp_locus(size_t position, const string& reference, const string& alternative) {//, int num_strains) {
    _position = position;
    _reference = reference;
    _alternative = alternative;
    //num_strains = 0;
    _num_strains = 0;//num_strains;
    _strains = NULL;
//     _strains = new unsigned char[_num_strains];
//     for (int i = 0; i < _num_strains; i++) {
//         _strains[i] = 0x00;
//     }
}

dbsnp_locus::dbsnp_locus(size_t position, const string& rsid, const string& reference, const string& alternative, int num_strains) {
    _position = position;
    _reference = reference;
    _alternative = alternative;
    //num_strains = 0;
    _rsid = rsid;
    _num_strains = num_strains;
    _strains = new unsigned char[_num_strains];
    for (int i = 0; i < _num_strains; i++) {
        _strains[i] = 0x00;
    }
}

dbsnp_locus::~dbsnp_locus() {
    delete[] _strains;
}

unsigned char dbsnp_locus::get_genotype(int strain_index) const {
    return _strains[strain_index];
}

void dbsnp_locus::set_genotype(int index, char const* info) {
    if (index < 0 || index > _num_strains) {
        throw std::out_of_range("cannot set genotype at the index");
    }
    if (_strains == NULL) {
        _strains = new unsigned char[_num_strains];
        for (int i = 0; i < _num_strains; i++) {
            _strains[i] = (unsigned char)0x00;
        }
    }
    if (info[0] == '.') {
        _strains[index] = 0x00;
    }
    int col = 0;
    const char* ptr = info;
    unsigned char genotype = 0x00;
    for (;;) {
        char c = *ptr;
        if (c == '\0') {
            break;
        } else if (c == ':') {
            break;
        } else if (c == '/') {
            col = 1;
            //if (++col >= 2) break;
        } else if (c >= '0' && c <= '9') {
            int an = atoi(ptr);
            unsigned char allele;
            if (an < 0) {
                allele = (unsigned char)0;
            } else if (an >= 15) {
                allele = (unsigned char)15;
            } else {
                allele = an + 1;
            }
            //cout << col << " ;; " << ptr << " => " << (int)allele << endl;
            if (col == 0) {
                genotype = allele << 4;
                col = -1;
            } else if (col == 1) {
                genotype |= allele;
                col = -1;
            }
        }
        ptr++;
    }
    //cout << info << "\t" << hex << (int)genotype << dec << endl;
  
    _strains[index] = genotype;
}

string dbsnp_locus::to_string() const {
    stringstream ss;
    ss << _position << "\t" << _reference << "\t" << _alternative;
    for (int i = 0; i < _num_strains; i++) {
        unsigned char g = _strains[i];
        if ((g & 240) == 0 || (g & 15) == 0) {
            ss << "\t" << ".";
        } else {
            ss << "\t" << (int)((g >> 4) - 1) << "/" << (int)((g & 15) - 1);
        }
//        ss << "\t" << (int)(g >> 4) << "/" << (int)(g & 15);
//        ss << "\t" << (int)g << ";" << (int)(g >> 4) << "/" << (int)(g & 15);
    }
    return ss.str();
}

string dbsnp_locus::to_string(const string& chromosome) const {
    stringstream ss;
    ss << chromosome << "\t" << _position << "\t";
    ss << _reference << "\t" << _alternative;
    for (int i = 0; i < _num_strains; i++) {
        unsigned char g = _strains[i];
        if ((g & 240) == 0 || (g & 15) == 0) {
            ss << "\t" << ".";
        } else {
            ss << "\t" << (int)((g >> 4) - 1) << "/" << (int)((g & 15) - 1);
        }
    }
    return ss.str();
}


const int dbsnp_file::CACHE_RETAIN_TOLERANCE = 10000000;
const size_t dbsnp_file::CACHE_LINE_INTERVAL = 100000;

dbsnp_file::dbsnp_file(const char* filename, const vector<string>& strains) {
    _strains = strains;
    _filename = filename;
    _cache_chromosome = "not_loaded";
    _cache_range_start = _cache_range_end = 0;
}

string dbsnp_file::get_cache_filename(const char* filename) {
    int slen = strlen(filename);
    string dirname;
    string basename;
    for (int i = slen - 1; i >= 0; i--) {
        if (filename[i] == FILE_SEPARATOR) {
            dirname = string(filename).substr(0, i + 1);
            basename = string(filename).substr(i + 1);
            break;
        }
    }
    return dirname + "." + basename + ".snpcache";
}

namespace {
    const unsigned int MAGIC = 0x129B6F1;
    const unsigned int MAGIC2 = 0x129B6F2;
    int header_size = 1024;
    const unsigned int CACHE_VERSION = 100;
}

void dbsnp_file::save_cache(const char* filename) const throw (exception) {
    ofstream fo(filename);
    if (!fo.is_open()) {
        throw invalid_argument("cannot save cache");
    }

    // assign probes to each chromosome
    map<string,vector<cache_position const*> > probes;
    //map<string,int> probe_size;
    for (int i = 0; i < (int)_indicators.size(); i++) {
        const cache_position* cp = _indicators[i];
        map<string,vector<cache_position const*> >::iterator it = probes.find(cp->chromosome);
        if (it == probes.end()) {//find(cp->chromosome)) {
            vector<cache_position const*> pos;
            pos.push_back(cp);
            probes[cp->chromosome] = pos;
        } else {
            it->second.push_back(cp);
        }
    }

    // magic number
    //fo.write(reinterpret_cast<char const*>(&MAGIC), sizeof(unsigned int));
    fo.write(reinterpret_cast<char const*>(&MAGIC2), sizeof(unsigned int));
    fo.write(reinterpret_cast<char const*>(&CACHE_VERSION), sizeof(unsigned int));

    // the number of chromosomes
    int num_chromosomes = probes.size();
    fo.write(reinterpret_cast<char const*>(&num_chromosomes), sizeof(int));
    
    // filename
    int len = _filename.size() + 1;
    fo.write(reinterpret_cast<char const*>(&len), sizeof(int));
    fo.write(_filename.c_str(), sizeof(char) * len);

    // strains
    int num_strains = _strains.size();
    fo.write(reinterpret_cast<char const*>(&num_strains), sizeof(int));
    for (int i = 0; i < num_strains; i++) {
        int len = _strains[i].size() + 1;
        const char* ptr = _strains[i].c_str();
        fo.write(reinterpret_cast<char const*>(&len), sizeof(int));
        fo.write(reinterpret_cast<char const*>(ptr), sizeof(char) * len);
    }

    // names of chromosome, length, name (ends with 0), address, number of cache probes
    int location = header_size;

    // chromosomes
    //map<string,int> pointers;
    vector<pair<string,int> > pointers;
    for (map<string,vector<cache_position const*> >::const_iterator it = probes.begin(); it != probes.end(); it++) {
        int num_probes = it->second.size();
        int buffer_size = sizeof(size_t) * 2 * num_probes;
        int slen = it->first.size() + 1;
        fo.write(reinterpret_cast<char const*>(&slen), sizeof(int));
        fo.write(it->first.c_str(), sizeof(char) * slen);
        fo.write(reinterpret_cast<char const*>(&location), sizeof(int));
        fo.write(reinterpret_cast<char*>(&num_probes), sizeof(int));

        int minpos = 0;
        int maxpos = 0;
        if (num_probes > 0) {
            minpos = (it->second)[0]->position;
            maxpos = (it->second)[it->second.size() - 1]->position;
            cerr << it->first << "\t" << minpos << "-" << maxpos << endl;
        }
        fo.write(reinterpret_cast<char*>(&minpos), sizeof(int));
        fo.write(reinterpret_cast<char*>(&maxpos), sizeof(int));

        pointers.push_back(make_pair(it->first, location));
        location += buffer_size;
    }
    
    for (vector<pair<string,int> >::const_iterator it = pointers.begin(); it != pointers.end(); it++) {
        string chromosome = it->first;
        fo.seekp(it->second);
        const vector<cache_position const*>& points = probes[chromosome];
        for (int i = 0; i < (int)points.size(); i++) {
            cache_position const* p = points[i];
            fo.write(reinterpret_cast<const char*>(&(p->position)), sizeof(size_t));
            fo.write(reinterpret_cast<const char*>(&(p->file_position)), sizeof(size_t));
        }
    }
    fo.close();
}

dbsnp_file* dbsnp_file::load_cache(const char* filename) throw (exception) {
    dbsnp_file* snpfile = NULL;
    ifstream fi(filename);
    if (fi.is_open() == false) {
        throw runtime_error("cannot open cache file");
    }

    // magic number
    unsigned int magic;
    unsigned int version = 0;
    fi.read(reinterpret_cast<char*>(&magic), sizeof(unsigned int));
    if (magic != MAGIC) {
        if (magic == MAGIC2) {
            fi.read(reinterpret_cast<char*>(&version), sizeof(unsigned int));
        } else {
            fi.close();
            throw runtime_error("magic number inconsistent");
        }
    }

    // the number of chromosomes
    int num_chromosomes;
    fi.read(reinterpret_cast<char*>(&num_chromosomes), sizeof(int));

    // filename
    int string_length;
    fi.read(reinterpret_cast<char*>(&string_length), sizeof(int));
    char* filename_vcf = new char[string_length];
    fi.read(filename_vcf, sizeof(char) * string_length);

    // strains
    int num_strains;
    fi.read(reinterpret_cast<char*>(&num_strains), sizeof(int));
    vector<string> strains;
    for (int i = 0; i < num_strains; i++) {
        int len;
        fi.read(reinterpret_cast<char*>(&len), sizeof(int));
        char* nbuf = new char[len];
        fi.read(nbuf, sizeof(char) * len);
        strains.push_back(string(nbuf));
        delete[] nbuf;
    }
    snpfile = new dbsnp_file(filename_vcf, strains);
    delete[] filename_vcf;

    // chromosomes
    // 0-4: length of name
    // 4-x: chromosome name
    // x-x+4:address
    // x+4-x+8:number of probes
    // x+8-x+16:minimum position, maximum_position (version >= 100)
    map<string,pair<int,int> > probe_range;
    for (int i = 0; i < num_chromosomes; i++) {
        int slen;
        char* cbuf;
        int num_probes, location;
        int minpos, maxpos;
        fi.read(reinterpret_cast<char*>(&slen), sizeof(int));
        cbuf = new char[slen];
        fi.read(cbuf, sizeof(char) * slen);
        fi.read(reinterpret_cast<char*>(&location), sizeof(int));
        fi.read(reinterpret_cast<char*>(&num_probes), sizeof(int));
        if (version >= 100) {
            fi.read(reinterpret_cast<char*>(&minpos), sizeof(int));
            fi.read(reinterpret_cast<char*>(&maxpos), sizeof(int));
            probe_range[cbuf] = std::make_pair(minpos, maxpos);
        }
        //cout << cbuf << "\t" << num_probes << " probes from " << hex << location << dec << endl;
        size_t ptr = fi.tellg();
        fi.seekg(location);
        for (int j = 0; j < num_probes; j++) {
            size_t cpos, fpos;
            fi.read(reinterpret_cast<char*>(&cpos), sizeof(size_t));
            fi.read(reinterpret_cast<char*>(&fpos), sizeof(size_t));
            snpfile->_indicators.push_back(new cache_position(cbuf, cpos, fpos));
            if (version < 100) {
                if (j == 0) {
                    
                }
            }
        }
        fi.seekg(ptr);
        delete[] cbuf;
    }

    sort(snpfile->_indicators.begin(), snpfile->_indicators.end(), cache_position_comparator);
    if (version >= 100) {
        snpfile->_cover_range = probe_range;
    }

    fi.close();
    return snpfile;
}

dbsnp_file::~dbsnp_file() {
    for (int i = 0; i < (int)_cache.size(); i++) {
        delete _cache[i];
    }
    for (int i = 0; i < (int)_indicators.size(); i++) {
        delete _indicators[i];
    }
}

void dbsnp_file::load_snps(const string& chromosome, int start, int end) {
    //cout << "CACHED RANGE:" << _cache_range_start << ", " << _cache_range_end << endl;
    if (_cache_chromosome == chromosome) {
        if (_cache_range_start <= start && end <= _cache_range_end) {
            return;
        }
    } else {
        _cache_chromosome = chromosome;
        _cache_range_start = _cache_range_end = 0;
    }
    //cout << "start ? " << (_cache_range_start <= start) << ", end ? " << (end <= _cache_range_end) << endl;
    for (int i = 0; i < (int)_cache.size(); i++) {
        delete _cache[i];
    }
    _cache.erase(_cache.begin(), _cache.end());

    if (_indicators.size() == 0) {
        return;
    }
    int left = 0;
    int right = _indicators.size();
    int center = (left + right) / 2;
    while (true) {
        const cache_position* probe = _indicators[center];
        int p = (int)probe->position;
        //cout << center << "\t" << left << "," << right << "\t" << (start) << "-" << (end) << "\t" << chromosome << "/" << probe->chromosome << ":" << (p) << "\t" << (p - start) << ", " << (end - p) << endl;
        if (probe->chromosome < chromosome) {
            left = center + 1;
        } else if (probe->chromosome > chromosome) {
            right = center; 
        } else if (p < start) {//(int)probe->position < start) {
            left = center + 1;
        } else if (p > end) {//(int)probe->position > end) {
            right = center;
        } else {
            break;
        }
        if (left >= right) {
            break;
        }
        center = (left + right) / 2;
    }

    int index = center;
    size_t left_limit = _indicators[0]->file_position;//header_size;
    size_t right_limit = _indicators[_indicators.size() - 1]->file_position;//0;//(int)_indicators.size();
    int left_index = 0;//center;
    int right_index = _indicators.size() - 1;//center;
    while (index >= 0) {
        const cache_position* probe = _indicators[index];
        if (probe->chromosome != chromosome || (int)probe->position < start) {
            left_limit = probe->file_position;
            left_index = index;
            //cout << index << ":" << left_limit << endl;
            break;
        }
        index--;
    }
    // if (index == 0 || _indicators[index]->chromosome != chromosome) {
    //     _cache_range_start = 0;
    // }

    index = center;
    while (index < (int)_indicators.size()) {
        const cache_position* probe = _indicators[index];
        //cout << probe->position << ", " << probe->file_position << endl;
        if (probe->chromosome != chromosome || (int)probe->position > end) {
            right_limit = probe->file_position;
            right_index = index;
            //cout << index << ":" << right_limit << endl;
            break;
        }
        index++;
    }
    //cout << left_index << " , " << right_index << endl;

    if (_indicators[center]->chromosome == chromosome) {
        _cache_range_start = std::min((int)_indicators[left_index]->position + 1, start);
        _cache_range_end = std::max((int)_indicators[right_index]->position - 1, end);
        if (right_index >= (int)_indicators.size() || _indicators[right_index]->chromosome != chromosome) {
            _cache_range_end = std::numeric_limits<int>::max();
        }
    } else {
        if (index >= (int)_indicators.size() - 1 || _indicators[index + 1]->chromosome > chromosome) {
            _cache_range_end = std::numeric_limits<int>::max();
        }
        if (index == 0 || _indicators[index - 1]->chromosome < chromosome) {
            _cache_range_start = 0;
        }
    }
    if (left_limit >= right_limit) {
        return;
    }

    ifstream fi(_filename.c_str());
    fi.seekg(left_limit);
    int col_strain = 9;
    for (;;) {
        string line;
        getline(fi, line);
        vector<string> items = split_items(line, '\t');
        //cout << items[0] << ":" << items[1] << "\t" << items[2] << endl;
        if (items[0] != chromosome) break;
        if (items.size() >= _strains.size() + col_strain) {
            int pos = atoi(items[1].c_str());
            dbsnp_locus* snp = new dbsnp_locus(pos, items[2], items[3], items[4], _strains.size());
            for (int i = 0; i < (int)_strains.size(); i++) {
                snp->set_genotype(i, items[i + col_strain].c_str());
            }
            _cache.push_back(snp);
        }
        if (fi.eof() || (right_limit > 0 && fi.tellg() >= (long long)right_limit)) break;
    }
    fi.close();
}

// get snps in the given range
// if snps out of range were selected, snps will be loaded via load_snps method
vector<dbsnp_locus const*> dbsnp_file::get_snps(string chromosome, int start, int end) const throw (exception) {
    vector<dbsnp_locus const*> snps;
    string cname;
    if (chromosome.find("chr") == 0) {
        cname = chromosome.substr(3, chromosome.size());
    } else {
        cname = chromosome;
    }
    if (_cover_range.size() > 0) {
        map<string,pair<int,int> >::const_iterator it = _cover_range.find(cname);
        if (it != _cover_range.end()) {
            pair<int,int> crange = it->second;
            if (start < crange.first) {
                start = crange.first;
            }
            if (end > crange.second) {
                end = crange.second;
            }
        } else {
            return snps;
        }
    }
    if (start > end) {
        return snps;
    }
    const_cast<dbsnp_file*>(this)->load_snps(cname, start, end);
    for (int i = 0; i < (int)(_cache.size()); i++) {
        dbsnp_locus const* cp = _cache[i];
        if (start <= (int)cp->_position && (int)cp->_position <= end) {
            snps.push_back(cp);
        }
    }
    return snps;
}

//
dbsnp_file* dbsnp_file::load_dbsnp(const char* filename, bool force_uncached) throw (exception) {
    string cache = get_cache_filename(filename);
    if (!force_uncached) {
        //cout << cache << endl;
        if (file_exists(cache.c_str())) {
            try {
                dbsnp_file* sf = load_cache(cache.c_str());
                return sf;
            } catch (exception& e) {
                cerr << e.what() << endl;
            }
        }
    }
    ifstream fi(filename);
    if (fi.is_open() == false) {
        throw invalid_argument("cannot open VCF file");
    }
    dbsnp_file* snpfile = NULL;
    size_t next_position = 0;
    size_t interval = CACHE_LINE_INTERVAL;
    string prev_chrom = "";
    size_t num_lines = 0;
    int cover_minimum = 0;
    int prev_pos = 0;
    //int num_strains = 0;
    unsigned int minimum_columns = std::numeric_limits<unsigned int>::max();
    while (!fi.eof()) {
        string line;
        size_t fpos = fi.tellg();
        getline(fi, line);
        if (snpfile == NULL) {
            if (line.find("#CHROM") == 0) {
                vector<string> strains;
                vector<string> items = split_items(line, '\t');
                if (items.size() > 9) {
                    minimum_columns = items.size();
                } else {
                    minimum_columns = items.size();
                }
                for (int i = 9; i < (int)items.size(); i++) {
                    strains.push_back(items[i]);
                }
                snpfile = new dbsnp_file(filename, strains);
            }
        } else {
            vector<string> items = split_items(line, '\t');
            if (items.size() >= minimum_columns) {
                if (++num_lines % 1000 == 0) {
                    cerr << " " << (num_lines / 1000) << " " << items[0] << ":" << items[1] << "        \r";
                }
                        
                string chrom = items[0];
                size_t position = std::atoi(items[1].c_str());
                if (chrom != prev_chrom) {
                    if (prev_chrom != "") {
                        snpfile->_cover_range[prev_chrom] = std::make_pair(cover_minimum, prev_pos);
                    }
                    next_position = 0;
                    prev_chrom = chrom;
                    cover_minimum = (int)position;
                }
                if (position >= next_position) {
                    snpfile->add_cache(chrom, position, fpos);
                    next_position += interval;
                }
                prev_pos = (int)position;
            }
        }
    }
    if (prev_chrom != "" && prev_pos > 0) {
        snpfile->_cover_range[prev_chrom] = std::make_pair(cover_minimum, prev_pos);
    } else {
        throw runtime_error("invalid VCF file : no SNPs");
    }
    snpfile->save_cache(cache.c_str());
    return snpfile;
}

void dbsnp_file::add_cache(const string& chromosome, size_t pos, size_t fpos) {
    _indicators.push_back(new cache_position(chromosome, pos, fpos));
}

namespace {
    unsigned char get_genotype_code(const string& genotype) {
        if (genotype.size() < 3) {
            return 0x00;
        }
        int left = (int)(genotype.c_str()[0] - '0') + 1;
        int right = (int)(genotype.c_str()[2] - '0') + 1;
        if (left > 0 && right > 0 && left <= 9 && right <= 9) {
            return (unsigned char)(((left << 4) & 240) | (right & 15));
        } else {
            return 0x00;
        }
    }

    void display_snp_distances(const dbsnp_file* dbsnp,//const vector<dbsnp_locus*>& snps,
                               const string& chromosome, int start, int end,
                               const map<int,unsigned char>& genotypes,
                               ostream& ost) {
        if (genotypes.size() == 0) {
            return;
        }
        vector<dbsnp_locus const*> snps = dbsnp->get_snps(chromosome, start, end);
        int num = dbsnp->strain_size();
        int* matches = new int[num + 1];
        int* unmatches = new int[num + 1];
        int* half = new int[num + 1];
        int* others = new int[num + 1];
        bool flag_xy = chromosome.find("X") != string::npos && chromosome.find("Y") != string::npos;
        for (int i = 0; i <= num; i++) {
            matches[i] = unmatches[i] = half[i] = others[i] = 0;
        }
        for (int i = 0; i < (int)snps.size(); i++) {
            dbsnp_locus const* snp = snps[i];
            int pos = snp->position();
            map<int,unsigned char>::const_iterator it = genotypes.find(pos);
            if (it != genotypes.end()) {
                unsigned char gt = it->second;
                if (flag_xy && ((gt >> 4) != (gt & 15))) {
                    continue;
                }
                if (gt == 0x11) {
                    matches[0] ++;
                } else if ((gt & 240) == 0x10 || (gt & 15) == 0x1) {
                    half[0] ++;
                } else {
                    unmatches[0]++;
                }

                for (int j = 0; j < num; j++) {
                    unsigned char gs;
                    gs = snp->get_genotype(j);
                    if (gs != 0) {
                        if (gs == 0x11) {
                            others[j + 1] ++;
                        } else if (gt == gs) {
                            matches[j + 1]++;
                        } else {
                            unsigned char g1 = gt >> 4;
                            unsigned char g2 = gt & 15;
                            unsigned char g3 = gs >> 4;
                            unsigned char g4 = gs & 15;
                            if (g1 == g3 || g1 == g4 || g2 == g3 || g2 == g4) {
                                half[j + 1]++;
                            } else {
                                unmatches[j + 1]++;
                            }
                        }
                    }
                }
            }
        }
        ost << chromosome << "\t" << start << "\t" << end << "\t" << snps.size();
        ost << std::setprecision(4);
        for (int i = 0; i <= num; i++) {
            int s0 = matches[i];
            int s1 = unmatches[i];
            int s2 = half[i];
            int s3 = others[i];
            int total = s0 + s1 + s2;
            double score;
            if (total == 0) {
                score = 0;
            } else {
                score = (double)(s1 + s2 * 0.5) * 100.0 / total;
            }
            ost << "\t" << score;
            ost << ";" << s0 << "/" << s2 << "/" << s1 << "/" << s3;
        }
        ost << "\n" << std::flush;
        delete[] matches;
        delete[] unmatches;
        delete[] half;
        delete[] others;
    }
}

void snp_distance::main(int argc, char** argv) throw (exception) {
    int window_size = get_argument_integer(argc, argv, "w", 1000000);
    int step = get_argument_integer(argc, argv, "s", window_size);
    const char* filename_input = get_argument_string(argc, argv, "i", NULL);
    const char* filename_output = get_argument_string(argc, argv, "o", NULL);
    const char* filename_snps = get_argument_string(argc, argv, "V", NULL);
    bool verbose = has_option(argc, argv, "verbose");
    ostream* ost = &cout;

    if (verbose) {
        cerr << "Input file         : " << filename_input << endl;
        cerr << "Output             : " << (filename_output == NULL ? "stdout" : filename_output) << endl;
        cerr << "Snps               : " << filename_snps << endl;
        cerr << "Window size        : " << window_size << endl;
        cerr << "Step               : " << step << endl;
    }


    ifstream fi(filename_input);
    if (fi.is_open() == false) {
        throw invalid_argument("cannot open input file");
    }
    if (filename_output != NULL) {
        ofstream* fo = new ofstream(filename_output);
        if (fo->is_open() == false) {
            delete fo;
            throw invalid_argument("cannot open output file");
        }
        ost = fo;
    }
  
    dbsnp_file* dbsnp = dbsnp_file::load_dbsnp(filename_snps);
  
    map<int,unsigned char> genotypes;
    int start = 0;
    int interval = step;//window_size;
    string chromosome;

    *ost << "#CHROM\tstart\tend\tsnps\tREF";
    for (int i = 0; i < (int)dbsnp->strain_size(); i++) {
        *ost << "\t" << dbsnp->get_strain(i);
    }
    *ost << "\n";
    while (!fi.eof()) {
        string line;
        getline(fi, line);
        //cout << line << endl;
        vector<string> items = split_items(line, '\t');
        vector<double> scores;
        //if (items[0] == chromosome && items.size() > 5) {
        if (items.size() > 5 && line.c_str()[0] != '#') {
            int pos = atoi(items[1].c_str());
            if (chromosome != items[0] || pos >= start + window_size) {
                if (genotypes.size() > 0) {
                    display_snp_distances(dbsnp, chromosome, start, start + window_size, genotypes, *ost);
                }
                //map<int,unsigned char> saved;
                if (chromosome != items[0]) {
                    genotypes.erase(genotypes.begin(), genotypes.end());
                    chromosome = items[0];
                    start = (pos / interval) * interval;
                } else {
		  start += interval;
                    for (map<int,unsigned char>::iterator it = genotypes.begin(); it != genotypes.end(); it++) {
		      if (it->first < start) {
                            genotypes.erase(it);
                        }
                    }
                }
            } 
            //const string& genotype = items[4];
            unsigned char gcode = get_genotype_code(items[4]);
            if (gcode > 0) {
                //cout << pos << " " << (int)gcode << endl;
                genotypes[pos] = gcode;
            }
        }
    }
    display_snp_distances(dbsnp, chromosome, start, start + window_size, genotypes, *ost);

    delete dbsnp;
}
