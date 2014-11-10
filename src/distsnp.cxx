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
#include <gtf.hxx>

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
dbsnp_locus::dbsnp_locus(size_t position, const string& reference, const string& alternative) {
    _position = position;
    _reference = reference;
    _alternative = alternative;
    _num_strains = 0;//num_strains;
    _strains = NULL;
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

int dbsnp_file::get_strain_index(const string& strain) const {
  for (int i = 0; i < (int)_strains.size(); i++) {
    if (strain == _strains[i]) {
      return i;
    }
  }
  return -1;
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
            //cerr << it->first << "\t" << minpos << "-" << maxpos << endl;
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
	  return;
        }
    }
    if (start > end) {
      return;
    }
    //cout << "CACHED RANGE:" << _cache_range_start << ", " << _cache_range_end << endl;
    if (_cache_chromosome == cname) {//hromosome) {
        if (_cache_range_start <= start && end <= _cache_range_end) {
            return;
        }
    } else {
      _cache_chromosome = cname;//chromosome;
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
        if (probe->chromosome < cname) {//chromosome) {
            left = center + 1;
        } else if (probe->chromosome > cname) {//chromosome) {
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
        if (probe->chromosome != cname || (int)probe->position < start) {
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
        if (probe->chromosome != cname || (int)probe->position > end) {
            right_limit = probe->file_position;
            right_index = index;
            //cout << index << ":" << right_limit << endl;
            break;
        }
        index++;
    }
    //cout << left_index << " , " << right_index << endl;

    if (_indicators[center]->chromosome == cname) {
        _cache_range_start = std::min((int)_indicators[left_index]->position + 1, start);
        _cache_range_end = std::max((int)_indicators[right_index]->position - 1, end);
        if (right_index >= (int)_indicators.size() || _indicators[right_index]->chromosome != cname) {//chromosome) {
            _cache_range_end = std::numeric_limits<int>::max();
        }
    } else {
      if (index >= (int)_indicators.size() - 1 || _indicators[index + 1]->chromosome > cname) {//chromosome) {
            _cache_range_end = std::numeric_limits<int>::max();
        }
        if (index == 0 || _indicators[index - 1]->chromosome < cname) {
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
        if (items[0] != cname) break;//chromosome) break;
        if (_strains.size() == 0 || items.size() >= _strains.size() + col_strain) {
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

dbsnp_locus const* dbsnp_file::get_snp(const string& chromosome, int position) const {
  const_cast<dbsnp_file*>(this)->load_snps(chromosome, position, position + 1);
    int left = 0;
    int right = _cache.size();
    for (;;) {
      if (left == right) break;
      int index = (left + right) >> 1;
      dbsnp_locus const* cp = _cache[index];
      if ((int)cp->_position < position) {
	left = index + 1;
      } else if ((int)cp->_position > position) {
	right = index;
      } else {
	return cp;
      }
    }
    return NULL;
}

// get snps in the given range
// if snps out of range were selected, snps will be loaded via load_snps method
vector<dbsnp_locus const*> dbsnp_file::get_snps(string chromosome, int start, int end) const throw (exception) {
    vector<dbsnp_locus const*> snps;
    const_cast<dbsnp_file*>(this)->load_snps(chromosome, start, end);
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

polymorphic_allele::polymorphic_allele(int chromosome, int position) {
    _chromosome = chromosome;
    _position = position;
    for (int i = 0; i < 8; i++) _bases[i] = 0;
    _reference = _alternative = '\0';
}

polymorphic_allele::polymorphic_allele(int chromosome, int position, int const* bases1, int const* bases2) {
    _chromosome = chromosome;
    _position = position;
    for (int i = 0; i < 4; i++) {
        _bases[i] = bases1[i];
        _bases[i + 4] = bases2[i];
    }
    _reference = _alternative = '\0';
}

void polymorphic_allele::set_variation(char ref, char alt) {
    _reference = ref;
    _alternative = alt;
}

void polymorphic_allele::resolve_variation() {
    int n1 = _bases[0] + _bases[4];
    int n2 = 0;
    int i1 = 0;
    int i2 = 0;
    for (int i = 1; i < 4; i++) {
        int n = _bases[i] + _bases[i + 4];
        if (n > n2) {
            if (n > n1) {
                n2 = n1; i2 = i1;
                n1 = n; i1 = i;
            } else {
                n2 = n; i2 = i;
            }
        }
    }
    _reference = "ACGT"[i1];
    _alternative = "ACGT"[i2];
}

char polymorphic_allele::reference() const {
    if (_reference != '\0') return _reference;
    const_cast<polymorphic_allele*>(this)->resolve_variation();
    return _reference;
}

char polymorphic_allele::alternative() const {
    if (_alternative != '\0') return _alternative;
    const_cast<polymorphic_allele*>(this)->resolve_variation();
    return _alternative;
}

void polymorphic_allele::set_bases(int slot, int a, int c, int g, int t) throw (out_of_range) {
    if (slot == 0) {
        _bases[0] = a; 
        _bases[1] = c; 
        _bases[2] = g; 
        _bases[3] = t; 
    } else if (slot == 1) {
        _bases[4] = a; 
        _bases[5] = c; 
        _bases[6] = g; 
        _bases[7] = t; 
    } else {
        throw out_of_range("slot must be 0 or 1");
    }
}

void polymorphic_allele::get_bases(int slot, int*& counts) const throw (out_of_range) {
    if (slot == 0) {
        for (int i = 0; i < 4; i++) counts[i] = _bases[i];
    } else if (slot == 1) {
        for (int i = 0; i < 4; i++) counts[i] = _bases[i + 4];
    } else {
        throw out_of_range("slot must be 0 or 1");
    }
}


string polymorphic_allele::to_string() const {
    stringstream ss;
    // CHROM POS ID REF ALT QUAL FILTER INFO FORMAT #1 #2
    string cname = tktools::bio::convert_code_to_chromosome(_chromosome);
    ss << cname << "\t" << (_position + 1) << "\t.\t" << reference() << "\t" << alternative() << "\t100.0\tPASS\t.\tGT:GP";
    int refid = 0;
    int altid = 0;
    switch (reference()) {
    case 'A': refid = 0; break;
    case 'C': refid = 1; break;
    case 'G': refid = 2; break;
    case 'T': refid = 3; break;
    }
    switch (alternative()) {
    case 'A': altid = 0; break;
    case 'C': altid = 1; break;
    case 'G': altid = 2; break;
    case 'T': altid = 3; break;
    }

    for (int i = 0; i < 2; i++) {
        int const* row = i == 0 ? _bases : _bases + 4;
        string gtype = "0/0"; // 0/1, 1/1
        int n1 = row[refid];
        int n2 = row[altid];
        if (n1 == 0 && n2 > 0) {
            gtype = "1/1";
        } else if (n1 > 0 && n2 > 0) {
            gtype = "0/1";
        } else if (n1 == 0 && n2 == 0) {
            gtype = "*";
        }
        ss << "\t" << gtype;
        for (int j = 0; j < 4; j++) {
            ss << ":" << row[j];
        }
    }
    return ss.str();
}

namespace tkbio {
    bool operator == (const polymorphic_allele& lhs, const polymorphic_allele& rhs) {
        return lhs._chromosome == rhs._chromosome && lhs._position == rhs._position;
    }
}

float polymorphic_allele::frequency(int slot) const throw (out_of_range) {
    int index;
    if (slot == 0) {
        index = 0;
    } else if (slot == 1) {
        index = 4;
    } else {
        throw out_of_range("slot must be 0 or 1");
    }
    int max_index = index;
    int max = _bases[index];
    int second_index = -1;
    int second = 0;
    for (int i = 1; i < 3; i++) {
        int cnt = _bases[i + index];
        if (cnt > second) {
            if (cnt > max) {
                second = max;
                second_index = max_index;
                max = cnt;
                max_index = i + index;
            } else {
                second = cnt;
                second_index = i + index;
            }
        }
    }
    if (second_index < 0) { // monomorphic
        return 0.0f;
    }
    int total = max + second;
    if (total == 0) {
        return 0.0f;
    }
    return (float)second / (total);
}

int polymorphic_allele::coverage(int slot) const throw (out_of_range) {
    if (slot == 0) {
        return _bases[0] + _bases[1] + _bases[2] + _bases[3];
    } else if (slot == 1) {
        return _bases[4] + _bases[5] + _bases[6] + _bases[7];
    } else {
        throw out_of_range("slot must be 0 or 1");
    }
}

polymorphic_allele polymorphic_allele::NO_DATA(-1,-1);

string allele_count::to_string() const {
    stringstream ss;
    ss << "H/H:" << this->homo_homo << " " << "H/h:" << this->homo_hetero << " "
       << "h/H:" << this->hetero_homo << " " << "h/h:" << this->hetero_hetero
       << " h/h*:" << this->hetero_hetero_diff
       << " other:" << this->uncountable;
    return ss.str();
}

double allele_count::homozygosity(bool whole) const {
    if (whole) {
        double h = homo_homo * 1.0 + hetero_homo * 0.0 + homo_hetero * 0.0 + homo_homo * 1.0 + hetero_hetero_diff;
        return h / (homo_homo + hetero_homo + homo_hetero + hetero_hetero + hetero_hetero_diff);
    } else {
        return (double)hetero_hetero / (hetero_hetero + hetero_hetero_diff);
    }
}

allele_count& operator += (allele_count& lhs, const allele_count& rhs) {
    lhs.homo_homo += rhs.homo_homo;
    lhs.homo_hetero += rhs.homo_hetero;
    lhs.hetero_homo += rhs.hetero_homo;
    lhs.hetero_hetero += rhs.hetero_hetero;
    lhs.hetero_hetero_diff += rhs.hetero_hetero_diff;
    lhs.uncountable += rhs.uncountable;
    return lhs;
}

denovo_snp::denovo_snp(const gtffile* gtf, int chromosome, int start, int stop) throw (logic_error) {
    if (start < 0) start = 0;
    _genes = gtf;
    _size = 0;
    _position = NULL;
    _start = start;
    _stop = stop;
    _count1 = _count2 = NULL;//new int*[0];
    _mapped1 = _mapped2 = 0;
    _reserved = 0;
    initialize_buffer(chromosome, start, stop);
}

denovo_snp::denovo_snp(int chromosome, int start, int stop) throw (logic_error) {
    if (start < 0) start = 0;
    _genes = NULL;
    _size = 0;
    _position = NULL;
    _start = start;
    _stop = stop;
    _count1 = _count2 = NULL;//new int*[0];
    _mapped1 = _mapped2 = 0;
    _reserved = 0;
    initialize_buffer_without_gtf(chromosome, start, stop);
}

denovo_snp::~denovo_snp() {
    release_buffer();
}

void denovo_snp::release_buffer() {
    //if (_count1 != NULL) {
    for (int i = 0; i < _reserved; i++) {
        delete[] _count1[i];
        delete[] _count2[i];
    }
    delete[] _count1;
    delete[] _count2;
    delete[] _position;
}

void denovo_snp::set_scope_without_gtf(int chromosome, int start, int stop) {
    if (chromosome != _chromosome || _stop < start || _start > stop || _reserved == 0) {
        release_buffer();
        initialize_buffer(chromosome, start, stop);
        return;
    }

    //int minpos = _start < start ? _start : start;
    //int maxpos = _stop > stop ? _stop : stop;
    //cerr << _reserved << ", " << minpos << "-" << maxpos << ", " << (maxpos - minpos) << endl;
    //cerr << _start << "-" << _stop << " => " << start << "-" << stop << endl;
    if (stop - start <= _reserved) {//maxpos - minpos <= _reserved) {
        for (int i = 0; i < stop - start; i++) {
            _position[i] = start + i;
        }
        if (_start < start) {
            // ------    current
            //   ------  next
            int offset = start - _start;
            for (int i = 0, j = offset; j < _stop - _start; i++, j++) {
                // for (int k = 0; k < 4; k++) {
                //     _count1[i][k] = _count1[j][k];
                //     _count2[i][k] = _count2[j][k];
                // }
                memcpy(_count1[i], _count1[j], sizeof(int) * 4);
                memcpy(_count2[i], _count2[j], sizeof(int) * 4);
            }
            //cerr << "clear " << _size - offset << ", " << (stop - start) << endl;
            for (int i = _size - offset; i < stop - start; i++) {
                //cerr << "clear at " << i << " / " << _size << " : " << hex << (void*)_count1[i] << ", " << (void*)_count2[i] << dec << endl;
                for (int j = 0; j < 4; j++) {
                    _count1[i][j] = _count2[i][j] = 0;
                }
            }
        } else {
            //   ------
            // ------
            int offset = _start - start;
            for (int i = stop - start - 1, j = (_stop - _start - 1) - offset; j >= 0; i--, j--) {
                memcpy(_count1[i], _count1[j], sizeof(int) * 4);
                memcpy(_count2[i], _count2[j], sizeof(int) * 4);
            }
            for (int i = 0; i < offset; i++) {
                for (int j = 0; j < 4; j++) {
                    _count1[i][j] = _count2[i][j] = 0;
                }
            }
        }
        _start = start;
        _stop = stop;
        _size = stop - start;
        _mapped1 = _mapped2 = 0;
        return;
    } else {
        //cerr << 
        delete[] _position;
        int** cbuf1 = _count1;
        int** cbuf2 = _count2;
        int offset = 0;
        int span = _size;
        //cout << "SIZE = " << _size << " => " << stop - start << endl;
        if (_start > start) {
            offset = _start - start;
        }
        int ps = _start;
        initialize_buffer(chromosome, start, stop);
        //cout << "REVISED SIZE = " << _size << endl;
        if (span > _size - offset) {
            span = _size - offset;
        }
        //cout << span << endl;
        for (int i = 0; i < span; i++) {//begin; i < tail; i++) {
            int pos = ps + i;
            if (_start <= pos && pos < _stop) {
                int index = pos - _start;
                //cerr << index << endl;
                memcpy(_count1 + index, cbuf1 + i, sizeof(int) * 4);
                memcpy(_count2 + index, cbuf2 + i, sizeof(int) * 4);
            }
            delete[] cbuf1[i];
            delete[] cbuf2[i];
        }
        //     int* c1 = _count1[i + offset];
        //     int* c2 = _count2[i + offset];
        //     for (int j = 0; j < 4; j++) {
        //         c1[j] = cbuf1[i][j];
        //         c2[j] = cbuf2[i][j];
        //     }
        //     delete[] cbuf1[i];
        //     delete[] cbuf2[i];
        // }
        delete[] cbuf1;
        delete[] cbuf2;
    }
}

void denovo_snp::set_scope(int chromosome, int start, int stop) {
    //cerr << __func__ << endl;
    if (_genes == NULL) {
        //cerr << __func__ << endl;
        set_scope_without_gtf(chromosome, start, stop);
        return;
    }
    if (_chromosome != chromosome || _size == 0) {
        release_buffer();
        initialize_buffer(chromosome, start, stop);
    } else {
        if (_position[_size - 1] >= start || _position[0] < stop) {
            vector<int> available;
            // int begin = _size;
            // for (int i = 0; i < _size; i++) {
            //     int pos = _position[i];
            //     if (pos >= start) {
            //         begin = i;
            //         break;
            //     }
            // }
            // int end = 0;
            // for (int i = _size - 1; i >= 0; i--) {
            //     int pos = _position[i];
            //     if (pos <= stop) {
            //         end = i + 1;
            //         break;
            //     }
            // }
            int old_size = _size;
            int* old_position = _position;
            int** old_count1 = _count1;
            int** old_count2 = _count2;

            initialize_buffer(chromosome, start, stop);
            int buffer_size = old_size + _size;
            int* pbuf = new int[buffer_size];
            int** cbuf1 = new int*[buffer_size];
            int** cbuf2 = new int*[buffer_size];
            for (int i = 0; i < buffer_size; i++) {
                cbuf1[i] = new int[4];
                cbuf2[i] = new int[4];
                for (int j = 0; j < 4; j++) {
                    cbuf1[i][j] = cbuf2[i][j] = 0;
                }
            }

            int i0 = 0, i1 = 0;
            int index = 0;
            for (;;) {
                int p0 = _position[i0];
                int p1 = old_position[i1];

                if (p0 < p1) {
                    pbuf[index++] = p0;
                    i0 ++;
                } else if (p0 > p1) {
                    pbuf[index++] = p1;
                    if (++i1 == old_size) {
                        for (int i = i0; i < _size; i++) {
                            pbuf[index++] = _position[i];
                        }
                        break;
                    }
                } else {
                    for (int j = 0; j < 4; j++) {
                        cbuf1[index][j] = old_count1[i1][j];
                        cbuf2[index][j] = old_count2[i1][j];
                    }
                    pbuf[index++] = p0;
                    i1++;
                    i0++;
                }
                if (i0 == _size) {
                    for (int i = i1; i < old_size; i++) {
                        memcpy(cbuf1[index], old_count1[i], sizeof(int) * 4);
                        memcpy(cbuf2[index], old_count2[i], sizeof(int) * 4);
                        pbuf[index++] = old_position[i];
                    }
                    break;
                } 
                if (i1 == _size) {
                    for (int i = i0; i < _size; i++) {
                        pbuf[index++] = _position[i];
                    }
                    break;
                }
            }

            release_buffer();
            _size = index;
            _reserved = buffer_size;
            _position = pbuf;
            _count1 = cbuf1;
            _count2 = cbuf2;
            _start = start;
            _stop = stop;
            _chromosome = chromosome;
            _mapped1 = _mapped2 = 0;
        }
    }
}

int denovo_snp::get_index(int position) const {
    int left = 0;
    int right = _size;
    for (;;) {
        if (left == right) {
            return -1;
        }
        int center = (left + right) >> 1;
        int p = _position[center];
        if (p < position) {
            left = center + 1;
        } else if (p > position) {
            right = center;
        } else {
            return center;
        }
    }
}

polymorphic_allele denovo_snp::get_allele(int position) const {
    int index = get_index(position);
    if (index < 0) {
        return polymorphic_allele::NO_DATA;
    }
    polymorphic_allele pa(_chromosome, position);
    int const* locus1 = _count1[index];
    int const* locus2 = _count2[index];
    pa.set_bases(0, locus1[0], locus1[1], locus1[2], locus1[3]);
    pa.set_bases(1, locus2[4], locus2[5], locus2[6], locus2[7]);
    return pa;
}

void denovo_snp::initialize_buffer_without_gtf(int chromosome, int start, int stop) throw (logic_error) {
    if (stop < start) {
        throw logic_error("start >= stop");
    }
    _genes = NULL;
    _chromosome = chromosome;
    _start = start;
    _stop = stop;
    _reserved = _size = stop - start;// + 1;
    _position = new int[_size];
    for (int i = 0; i < _size; i++) {
        _position[i] = i + start;
    }
    _count1 = new int*[_size];
    _count2 = new int*[_size];
    for (int i = 0; i < _size; i++) {
        _count1[i] = new int[4];
        _count2[i] = new int[4];
        for (int j = 0; j < 4; j++) {
            _count1[i][j] = _count2[i][j] = 0;
        }
    }
}

void denovo_snp::initialize_buffer(int chromosome, int start, int stop) throw (logic_error) {
    //cerr << "BUFFER " << hex << (void*) _genes << dec << "      " << endl;
    if (_genes == NULL) {
        initialize_buffer_without_gtf(chromosome, start, stop);
        return;
    }
    int* posbuffer = new int[stop - start];// + 1];
    for (int i = 0, span = stop - start; i < span; i++) {
        posbuffer[i] = 0;
    }
    _chromosome = chromosome;
    vector<const gtfgene*> genes = _genes->find_genes(_chromosome, start, stop);
    //cout << genes.size() << " genes found\n";
    for (int i = 0; i < (int)genes.size(); i++) {
        const vector<gtfexon>& exons = genes[i]->exons();
        //cout << "ADD " << genes[i]->to_string() << endl;
        for (int j = 0; j < (int)exons.size(); j++) {
            int p5 = exons[j].position5();
            int p3 = exons[j].position3();
            if (p5 < start) p5 = start;
            if (p3 >= stop) p3 = stop - 1;
            for (int p = p5; p <= p3; p++) {
                posbuffer[p - start] = 1;
            }
        }
    }
    // count 
    int length = 0;
    for (int i = 0, span = stop - start; i < span; i++) {
        if (posbuffer[i] != 0) {
            length++;
        }
    }
    //cout << "buffer size = " << length << endl;

    _size = length;
    _reserved = _size;
    _position = new int[length];
    _count1 = new int*[length];
    _count2 = new int*[length];
    for (int i = 0, span = stop - start, index = 0; i < span; i++) {
        if (posbuffer[i] != 0) {
            _position[index] = i + start;
            _count1[index] = new int[4];
            _count2[index] = new int[4];
            for (int j = 0; j < 4; j++) {
                _count1[j] = _count2[j] = 0;
            }
            index++;
        }
    }
    delete[] posbuffer;
}

void denovo_snp::set_quality(int qual) {
    _quality = qual;
}

void denovo_snp::add_read(int slot, bam1_t const* read) throw (out_of_range) {
    uint32_t const* cigar = bam1_cigar(read);
    uint8_t const* sequence = bam1_seq(read);
    uint8_t const* qptr = bam1_qual(read);
    int len = read->core.n_cigar;
    int rpos = read->core.pos; // reference position
    int qpos = 0; // query position
    int offset = 4;
    if (slot != 0 && slot != 1) {
        throw out_of_range("slot is 0 or 1");
    }
    int** count;
    if (slot == 0) {
        count = _count1;
        _mapped1++;
    } else if (slot == 1) {
        count = _count2;
        _mapped2++;
    } else {
        throw out_of_range("slot must be 0 or 1");
    }
//    int 
    bool flag = false;
    static int base2index[16] = {-1,0,1,-1,2,-1,-1,-1,3,-1,-1,-1,-1,-1,-1,-1};
    //static char base2nuc[16] = "0AC.G...T......";
    //string seq = "";
    //stringstream ss;
    unsigned char qthr = _quality;
    if (_quality >= 255) {
        qthr = (unsigned char)255;
    } else if (_quality < 0) {
        qthr = 0;
    }
    for (int i = 0; i < len; i++) {
        int op = bam_cigar_op(cigar[i]);
        int slen = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) { // M
            int index = -1;
            for (int j = 0; j < slen; j++) {
                if (index < 0) {
                    index = get_index(rpos);
                } else {
                    if (_position[++index] != rpos) {
                        index = -1;
                    }
                }
                if (index >= 0) {
                    uint8_t base = (sequence[qpos >> 1] >> offset) & 15;
                    //cout << qpos << "\t" << (int)base << endl;
                    int bi = base2index[base];
                    if (index >= 0) {
                        unsigned char qual = qptr[qpos];
                        if (qual > qthr) {
                            count[index][bi] ++;
                            flag = true;
                        }
                    }
                    //ss << base2nuc[base];
                }
                qpos++;
                rpos++;
                //offset = 4 - offset;
                offset ^= 4;
            }
        } else if (op == BAM_CINS) { // I
            qpos += slen;
            offset ^= ((slen & 1) << 2);
        } else if (op == BAM_CDEL) { // D
            rpos += slen;
        } else {
            break;
        }
    }
    //cout << ss.str() << endl;
    if (flag) {
        if (slot == 0) {
            _mapped1++;
        } else {
            _mapped2++;
        }
    }
}

namespace {
    void find_major_alleles(int const* bases, int& b1, int& n1, int& b2, int& n2) {
        n1 = bases[0];
        b1 = 0;
        n2 = 0;
        b2 = 1;
        for (int i = 1; i < 4; i++) {
            int n = bases[i];
            if (n > n2) {
                if (n > n1) {
                    n2 = n1; b2 = b1;
                    n1 = n; b1 = i;
                } else {
                    n2 = n; b2 = i;
                }
            }
            n = bases[i + 4];
        }
    }
    bool has_third_allele(int const* bases) {
        int n = 0;
        for (int i = 0; i < 4; i++) {
            n += (int)(bases[i] != 0 ? 1 : 0);
        }
        return n > 2;
    }
}

vector<polymorphic_allele> denovo_snp::get_polymorphism(int coverage, double hetero) const {
    vector<polymorphic_allele> alleles;
    for (int i = 0; i < _size; i++) {
        int b11, b12, b21, b22;
        int n11, n12, n21, n22;
        find_major_alleles(_count1[i], b11, n11, b12, n12);
        int c1 = n11 + n12;
        if (c1 < coverage) continue;
        find_major_alleles(_count2[i], b21, n21, b22, n22);
        int c2 = n21 + n22;
        if (c2 < coverage) continue;
        if ((c1 * hetero < n12 && n22 == 0) 
            || (c2 * hetero < n22 && n12 == 0)
            || (n22 == 0 && n12 == 0 && b11 != b21)) {
            alleles.push_back(polymorphic_allele(_chromosome, _position[i], _count1[i], _count2[i]));
        }
    }
    return alleles;
}

vector<polymorphic_allele> denovo_snp::get_polymorphism(int coverage, double heterogeneity, int start, int stop) const {
    vector<polymorphic_allele> alleles;
    for (int i = 0; i < _size; i++) {
        int p = _position[i];
        if (p < start) continue;
        if (p >= stop) break;
        int b11, b12, b21, b22;
        int n11, n12, n21, n22;
        find_major_alleles(_count1[i], b11, n11, b12, n12);
        int c1 = n11 + n12;
        if (c1 < coverage) continue;
        find_major_alleles(_count2[i], b21, n21, b22, n22);
        int c2 = n21 + n22;
        if (c2 < coverage) continue;
        if ((c1 * heterogeneity < n12 && n22 == 0) 
            || (c2 * heterogeneity < n22 && n12 == 0)) {
            //cout << p << ":" << n11 << "," << n12 << " // " << n21 << "," << n22 << endl;
            alleles.push_back(polymorphic_allele(_chromosome, _position[i], _count1[i], _count2[i]));
        }
    }
    return alleles;
}

allele_count denovo_snp::count_allele_types(int coverage, double hetero_threshold, int start, int stop, chromosome_seq const* chromosome) const {
    allele_count count;
    for (int i = 0; i < _size; i++) {
        int p = _position[i];
        if (p < start) continue;
        if (p >= stop || p >= chromosome->length()) break;
        int b11, b12, b21, b22;
        int n11, n12, n21, n22;
        find_major_alleles(_count1[i], b11, n11, b12, n12);
        int c1 = n11 + n12;
        if (c1 < coverage || has_third_allele(_count1[i])) continue;
        find_major_alleles(_count2[i], b21, n21, b22, n22);
        int c2 = n21 + n22;
        if (c2 < coverage || has_third_allele(_count2[i])) continue;
        // if (b11 != b21) {
        //     cerr << p << " " << chromosome->get_base(p+1) << " " << b11 << " " << b21 << endl;
        // }
        double r1 = c1 * hetero_threshold;
        double r2 = c2 * hetero_threshold;
        int h1, h2;
        if (r1 < n12) {        // 1 is hetero
            h1 = 1;
        } else if (n12 == 0) { // 1 is homo
            h1 = 0;
        } else {
            count.uncountable ++;
            continue;
        }
        if (r2 < n22) {        // 2 is hetero
            h2 = 1;
        } else if (n22 == 0) { // 2 is homoe
            h2 = 0;
        } else {
            count.uncountable++;
            continue;
        }

        if (h1 == 0 && h2 == 0) {
            if (b11 == b21) { // homo for both
                count.homo_homo++;
            } else {
                count.uncountable++;
            }
        } else if (h1 == 1 && h2 == 1) { // hetero hetero
            int b6 = chromosome->get_base_id(p + 1);
            if (b11 == b6) {
                if (b21 == b6) {
                    if (b12 == b22) {
                        count.hetero_hetero++;
                    } else {
                        count.hetero_hetero_diff++;
                    }
                } else if (b22 == b6) {
                    if (b12 == b21) {
                        count.hetero_hetero++;
                    } else {
                        count.hetero_hetero_diff++;
                    }
                }
            } else if (b12 == b6) {
                if (b21 == b6) {
                    if (b11 == b22) {
                        count.hetero_hetero++;
                    } else {
                        count.hetero_hetero_diff++;
                    }
                } else if (b22 == b6) {
                    if (b11 == b21) {
                        count.hetero_hetero++;
                    } else {
                        count.hetero_hetero_diff++;
                    }
                }
            } else {
                count.uncountable++;
            }
        } else if (h1 == 0 && h2 == 1) { // homo hetero
            if (b21 == b11 || b22 == b11) {
                count.homo_hetero++;
            } else {
                count.uncountable++;
            }                
        } else if (h1 == 1 && h2 == 0) { // hetero homo
            if (b11 == b21 || b12 == b11) {
                count.hetero_homo++;
            } else {
                count.uncountable++;
            }
        }
    }
    return count;
}

namespace {
    string get_sample_name(string filename, char delimiter='/') {
        int pos = 0;
        for (int i = filename.size() - 1; i >= 0; i--) {
            char c = filename.c_str()[i];
            if (c == delimiter) {
                pos = i + 1;
                break;
            }
        }
        int tail = filename.size();
        for (int i = pos; i < (int)filename.size(); i++) {
            char c = filename.c_str()[i];
            if (c == '.') {
                tail = i;
                break;
            }
        }
        return filename.substr(pos, tail - pos);
    }
}

void denovo_snp::enumerate_hetero(int argc, char** argv) throw (exception) {
    try {
        const char* filename1 = get_argument_string(argc, argv, "1", "/mnt/smb/tae/stap/shira/BAM6/Sample6.bam");
        const char* filename2 = get_argument_string(argc, argv, "2", "/mnt/smb/tae/stap/shira/BAM12/Sample12.bam");
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        const char* filename_gtf = get_argument_string(argc, argv, "G", NULL);
        const char* filename_chrm = get_argument_string(argc, argv, "s", NULL);
        int coverage = get_argument_integer(argc, argv, "c", 20);
        double heterozygosity = get_argument_float(argc, argv, "z", 0.33);
        int chunk_size = get_argument_integer(argc, argv, "w",  200000);
        int margin_size = get_argument_integer(argc, argv, "m", 100000);
        int num_files = 2;
        int maximum_reads_in_window = get_argument_integer(argc, argv, "x", 0);
        int quality_threshold = get_argument_integer(argc, argv, "q", 10);
        bool verbose = has_option(argc, argv, "verbose");
        bool global_homozygosity = has_option(argc, argv, "I");
        const char* target_chromosome = get_argument_string(argc, argv, "C", NULL);

        // if (filename_gtf == NULL || file_exists(filename_gtf) == false) {
        //     throw logic_error("cannot open GTF file");
        // }
        if (verbose) {
            cerr << "filename 1 : " << filename1 << endl;
            cerr << "filename 2 : " << filename2 << endl;
            cerr << "gtf        : " << (filename_gtf == NULL ? "none" : filename_gtf)  << endl;
            cerr << "heterozygosity : " << heterozygosity << endl;
            cerr << "coverage   : " << coverage << endl;
            cerr << "chunk size : " << chunk_size << endl;
            cerr << "margin     : " << margin_size << endl;
            cerr << "max reads  : " << maximum_reads_in_window << endl;
            cerr << "quality    : " << quality_threshold << endl;
            cerr << "output     : " << (filename_output == NULL ? "stdout" : filename_output) << endl;
            if (filename_chrm != NULL) {
                cerr << "genome     : " << filename_chrm << endl;
            }
        }
        if (file_exists(filename1) == false || file_exists(filename2) == false) {
            throw logic_error("two bam files required");
        }

        //cout << "GTF\n";
        gtffile* gtf = NULL;
        if (filename_gtf != NULL) {
            gtffile::load_gtf(filename_gtf);
            if (verbose) {
                cerr << "GTF file " << gtf->size() << " genes\n";
            }
        }

        // if (true) {
        //     //cout << "test\n";
        //     denovo_snp* snp = new denovo_snp(gtf, 1, 105000000, 119000000);
        //     snp->set_scope(1, 115000000, 120000000);
        //     delete snp;
        //     return;
        // }

        ostream* ost = &cout;
        bamFile* bamfiles = new bamFile[num_files];
        bam_header_t** headers = new bam_header_t*[num_files];
        bam1_t** reads = new bam1_t*[num_files];
        int* status = new int[num_files];

        if (filename_output != NULL) {
            ofstream* fo = new ofstream(filename_output);
            if (fo->is_open() == false) {
                throw invalid_argument("cannot open output file");
            }
            ost = fo;
        }

        // if (verbose) {
        //     cerr << "opening files\n";
        // }
        bamfiles[0] = bam_open(filename1, "rb");
        bamfiles[1] = bam_open(filename2, "rb");
        for (int i = 0; i < num_files; i++) {
            headers[i] = bam_header_read(bamfiles[i]);
            reads[i] = bam_init1();
            status[i] = 0;
        }

        // if (verbose) {
        //     cerr << "check integrity\n";
        // }
        if (check_header_consistency(num_files, headers) == false) {
            throw runtime_error("incompatible bam files");
        }

        map<int,int> bc2cc;
        {
            set<string> accepted;
            if (target_chromosome != NULL) {
                vector<string> items = split_items(target_chromosome, ',');
                for (int i = 0; i < (int)items.size(); i++) {
                    cerr << items[i] << endl;
                    accepted.insert(items[i]);
                }
            }
            bam_header_t* h = headers[0];
            for (int i = 0; i < h->n_targets; i++) {
                const char* hn = h->target_name[i];
                if (strncmp(hn, "chr", 3) == 0) {
                    hn += 3;
                } 
                int code = -1;
                if (accepted.size() == 0 || accepted.find(hn) != accepted.end()) {
                    //cout << "accepted " << hn << endl;
                    code = convert_chromosome_to_code(hn);
                }
                if (code < 0 || code > 128) {
                    code = -1;
                }
                if (code >= 0) {
                    bc2cc[i] = code;
                }
            }
        }
        //cout << bc2cc.size() << endl;
        //exit(0);
        int current_chromosome = -1;
        int next_chromosome = -1;
        int position = 0;
        int next_position = position + chunk_size;
        //int steps = 0;
        vector<chromosome_seq*> chromosomes;

        if (filename_chrm != NULL) {
            cerr << filename_chrm << endl;
            chromosomes = chromosome_seq::load_genome(filename_chrm);
        } else if (global_homozygosity) {
            throw logic_error("measuring homozygosity requires reference genome");
        }

        // if (verbose) {
        //     cerr << "set chromosome names\n";
        // }
        map<int,chromosome_seq const*> bam2chrm = chromosome_seq::map_chromosome(headers[0], chromosomes);

        denovo_snp* detector = NULL;
        set<int> touched_chromosome;
        allele_count snp_groups;
        chromosome_seq const* chromosome = NULL;
        
//        bool debug = false;
        if (!global_homozygosity) {
            *ost << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
            *ost << get_sample_name(filename1) << "\t" << get_sample_name(filename2) << endl;
        }
        for (;;) {
            bool chromosome_change = false;
            bool finished = false;
            for (int i = 0; i < num_files; i++) {
                //cerr << i << endl;
                for (;;) {
                    bam1_t const* r;
                    if (status[i] != 0) {
                        r = reads[i];
                    } else {
                        if (bam_read1(bamfiles[i], reads[i]) > 0) {
                            r = reads[i];
                        } else {
                            finished = true;
                            chromosome_change = false;
                            break;
                        }
                    }
                    //bam1_t const* r = reads[i];
                    int chrm = r->core.tid;
                    int pos = r->core.pos;
                    
                    if (chrm != current_chromosome) {
                        if (touched_chromosome.find(chrm) != touched_chromosome.end()) {
                            //cerr << "touched " << i << "\n";
                            status[i] = 0;
                            continue;
                        }
                        chromosome_change = true;
                        if (next_chromosome <= 0) {
                            next_chromosome = chrm;
                        }
                        status[i] = 1;
                        break;
                    }
                    if (pos > next_position) {
                        status[i] = 1;
                        break;
                    }
                    
                    if (detector == NULL) {
                        if (bc2cc.find(current_chromosome) != bc2cc.end()) {
                            int chrm = bc2cc[current_chromosome];
                            //cerr << "instanciate new detector      \n";
                            detector = new denovo_snp(gtf, chrm, position - margin_size, position + chunk_size + margin_size);
                            detector->set_quality(quality_threshold);
                        }
                    }
                    if (detector != NULL) {
                        if ((maximum_reads_in_window == 0 || detector->mapped(i) < maximum_reads_in_window) && detector->has_buffer()) {
                            detector->add_read(i, r);
                        }
                    }
                    status[i] = 0;
                }
            }
            
            // detect sample specific SNPs
            if (current_chromosome >= 0) {// && detector != NULL) {
                //cerr << "analyzing " << current_chromosome << "       " << endl;
                //cerr << headers[0]->target_name[0] << ", " << headers[0]->target_name[1] << ", " << endl;
                if (detector != NULL) {
                    if (global_homozygosity) {
                        if (chromosome != NULL) {
                            allele_count cnt = detector->count_allele_types(coverage, heterozygosity, position, position + chunk_size, chromosome);
                            snp_groups += cnt;
                            // snp_groups.homo_homo += cnt.homo_homo;
                            // snp_groups.hetero_homo += cnt.hetero_homo;
                            // snp_groups.homo_hetero += cnt.homo_hetero;
                            // snp_groups.hetero_hetero += cnt.hetero_hetero;
                            // snp_groups.uncountable += cnt.uncountable;
                            if (verbose) {
                                cerr << " " << headers[0]->target_name[current_chromosome] << ":" << position << "-" << next_position;// << "   " << __LINE__;
                                cerr <<" " <<  snp_groups.to_string() << "        \r";
                            }
                        }
                    } else {
                        if (verbose) {
                            cerr << " " << headers[0]->target_name[current_chromosome] << ":" << position << "-" << next_position;// << "   " << __LINE__;
                            //cerr << current_chromosome << ":" << position << "-" << next_position << " ";
                            if (detector != NULL) {
                                for (int i = 0; i < num_files; i++) {
                                    cerr << " " << i << ":" << detector->mapped(i);
                                }
                            }
                            cerr << "     \r";
                        }
                        vector<polymorphic_allele> alleles = detector->get_polymorphism(coverage, heterozygosity, position, position + chunk_size);
                        for (int i = 0; i < (int)alleles.size(); i++) {
                            *ost << alleles[i].to_string() << endl << flush;
                        }
                    }
                }
            }

            if (finished) {
                break;
            }

            if (chromosome_change) {
                touched_chromosome.insert(current_chromosome);
                int pos_min = numeric_limits<int>::max();
                if (verbose) {
                    if (global_homozygosity) {
                        cerr << snp_groups.to_string() << endl;
                        cerr << "Homozygosity (total)   : " << snp_groups.homozygosity() << endl;
                        cerr << "Homozygosity (partial) : " << snp_groups.homozygosity(false) << endl;
                    }
                    cerr << "change chromosome to " << headers[0]->target_name[next_chromosome] << "      \r";//, sweep " << total_reads << "reads from " << position << "\r";

                    // int total_reads = 0;
                    // for (int i = 0; i < num_files; i++) {
                    //     //total_reads += detectors[i]->size();
                    // }

                    
                    // cerr << "change chromosome to " << headers[0]->target_name[next_chromosome] << ", sweep " << total_reads << "reads from " << position << "\r";
                }
                if (!finished) {
                    for (int i = 0; i < num_files; i++) {
                        //detectors[i]->sweep(current_chromosome);
                        pair<int,int> region;// = detectors[i]->span();
                        pos_min = pos_min < region.first ? pos_min : region.first;
                    }
                    if (pos_min == numeric_limits<int>::max()) {
                        pos_min = 0;
                    }
                    current_chromosome = next_chromosome;
                    next_chromosome = -1;
                    position = pos_min;
                    next_position = position + chunk_size;
                    if (detector != NULL) {
                        detector->set_scope(bc2cc[current_chromosome], position, next_position + margin_size);
                    }
                    if (bam2chrm.find(current_chromosome) != bam2chrm.end()) {
                        chromosome = bam2chrm[current_chromosome];
                    } else {
                        chromosome = NULL;
                        cerr << "no chromosome named " << current_chromosome << endl;
                    }
                } else {
                    break;
                }
            } else { // next position
                // for (int i = 0; i < num_files; i++) {
                //     //detectors[i]->sweep(current_chromosome, 0, next_position);// - margin_size);
                // }
                int pos_min = numeric_limits<int>::max();
                for (int i = 0; i < num_files; i++) {
                    pair<int,int> span;// = detectors[i]->span();
                    if (span.second > 0 && pos_min > span.first) {
                        pos_min = span.first;
                    }
                    //detectors[i]->sweep(current_chromosome);
                }
                if (pos_min == numeric_limits<int>::max()) {
                    pos_min = next_position;
                } else if (pos_min > next_position) {
                    next_position = pos_min;
                }
                position = next_position;
                next_position = position + chunk_size;
                if (detector != NULL) {
                    detector->set_scope(detector->_chromosome, position, next_position + margin_size);
                }
            }
        }

        if (global_homozygosity) {
            *ost << snp_groups.to_string() << endl;
            *ost << "Homozygosity (total)   : " << snp_groups.homozygosity() << endl;
            *ost << "Homozygosity (partial) : " << snp_groups.homozygosity(false) << endl;
            
        }
        
        if (filename_output != NULL) {
            dynamic_cast<ofstream*>(ost)->close();
            delete ost;
            ost = &cout;
        }

        for (int i = 0; i < num_files; i++) {
            if (reads[i] != NULL) {
                bam_destroy1(reads[i]);
            }
            bam_header_destroy(headers[i]);
            bam_close(bamfiles[i]);
            //delete detectors[i];
        }
        delete detector;
        delete[] reads;
        delete[] headers;
        delete[] bamfiles;
        delete gtf;
        delete[] status;
        //return 0;

    } catch (exception& e) {
        throw;
    }
}

