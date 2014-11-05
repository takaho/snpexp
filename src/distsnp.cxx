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

using std::strncmp;
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
        //cout << items[0] << ":" << items[1] << "\t" << items[2] << endl;
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
}

polymorphic_allele::polymorphic_allele(int chromosome, int position, int const* bases1, int const* bases2) {
    _chromosome = chromosome;
    _position = position;
    for (int i = 0; i < 4; i++) {
        _bases[i] = bases1[i];
        _bases[i + 4] = bases2[i];
    }
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

denovo_snp::denovo_snp(const gtffile* gtf, int chromosome, int start, int stop) throw (logic_error) {
    if (start < 0) start = 0;
    _genes = gtf;
    _size = 0;
    _position = NULL;
    _start = start;
    _stop = stop;
    _count1 = _count2 = NULL;//new int*[0];
    _mapped = 0;
    _reserved = 0;
    initialize_buffer(gtf, chromosome, start, stop);
}

denovo_snp::~denovo_snp() {
    release_buffer();
}

void denovo_snp::release_buffer() {
    if (_count1 != NULL) {
        for (int i = 0; i < _reserved; i++) {
            delete[] _count1[i];
            delete[] _count2[i];
        }
    }
    delete[] _count1;
    delete[] _count2;
    delete[] _position;
}

void denovo_snp::set_scope(int chromosome, int start, int stop) {
    if (_chromosome != chromosome || _size == 0) {
        release_buffer();
        initialize_buffer(_genes, chromosome, start, stop);
    } else {
        if (_position[_size - 1] >= start || _position[0] < stop) {
            vector<int> available;
            int begin = _size;
            for (int i = 0; i < _size; i++) {
                int pos = _position[i];
                if (pos >= start) {
                    begin = i;
                    break;
                }
            }
            int end = 0;
            for (int i = _size - 1; i >= 0; i--) {
                int pos = _position[i];
                if (pos <= stop) {
                    end = i + 1;
                    break;
                }
            }
            int old_size = _size;
            int* old_position = _position;
            int** old_count1 = _count1;
            int** old_count2 = _count2;
            //cout << "SIZE=" << _size << endl;
            //cout << "OLD:" << hex << (void*)old_count1 << ", " << (void*)old_count2 << dec << endl;
            // for (int i = 0; i < _size; i++) {
            //     cout << i << ":OLD:" << hex << (void*)old_count1[i] << ", " << (void*)old_count2[i] << dec << endl;
            // }

            initialize_buffer(_genes, chromosome, start, stop);
            int buffer_size = old_size + _size;
            //cout << "prepared buffer size " << buffer_size << endl;
            //release_buffer();
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
            //cout << "OLD:" << hex << (void*)old_count1[0] << ", " << (void*)old_count2[0] << dec << endl;
            // i0 : index of current fields
            // i1 : index of old fields
            int i0 = 0, i1 = 0;
            int index = 0;
            for (;;) {
                int p0 = _position[i0];
                int p1 = old_position[i1];
                //cout << "INDEX:" << index << ", " << buffer_size - index << " = " << _size - i0 << " + " << old_size - i1 << endl;
                //cout << i0 << ":" << p0 << ", " << i1 << ":" << p1 << endl;
                //cout << "OLD:" << hex << (void*)old_count1[0] << ", " << (void*)old_count2[0] << dec << endl;

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
                    //cout << "OLD:" << hex << (void*)old_count1[0] << ", " << (void*)old_count2[0] << dec << endl;
                    //cout << index << " : pos=" << p0 << endl;
                    for (int j = 0; j < 4; j++) {
                        //cout << j << endl;
                        //                        cout << hex << (void*)cbuf1[index] << ", " << old_count1[i1] << dec << endl;
                        cbuf1[index][j] = old_count1[i1][j];
                        cbuf2[index][j] = old_count2[i1][j];
                    }
                    pbuf[index++] = p0;
                    i1++;
                    i0++;
                }
                if (i0 == _size) {
                    //cout << "fill tail\n";
                    for (int i = i1; i < old_size; i++) {
                        //cout << i << " / " << old_size << endl;
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
            //cout << "INDEX=" << index << endl;
            release_buffer();
            //cout << "released\n";
            _size = index;
            _reserved = buffer_size;
            _position = pbuf;
            _count1 = cbuf1;
            _count2 = cbuf2;
            _start = start;
            _stop = stop;
            _chromosome = chromosome;
            _mapped = 0;
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

void denovo_snp::initialize_buffer(const gtffile* gtf, int chromosome, int start, int end) throw (logic_error) {
    int* posbuffer = new int[end - start];
    for (int i = 0, span = end - start; i < span; i++) {
        posbuffer[i] = 0;
    }
    _chromosome = chromosome;
    vector<const gtfgene*> genes = gtf->find_genes(_chromosome, start, end);
    //cout << genes.size() << " genes found\n";
    for (int i = 0; i < (int)genes.size(); i++) {
        const vector<gtfexon>& exons = genes[i]->exons();
        //cout << "ADD " << genes[i]->to_string() << endl;
        for (int j = 0; j < (int)exons.size(); j++) {
            int p5 = exons[j].position5();
            int p3 = exons[j].position3();
            if (p5 < start) p5 = start;
            if (p3 >= end) p3 = end - 1;
            for (int p = p5; p <= p3; p++) {
                posbuffer[p - start] = 1;
            }
        }
    }
    // count 
    int length = 0;
    for (int i = 0, span = end - start; i < span; i++) {
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
    for (int i = 0, span = end - start, index = 0; i < span; i++) {
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

void denovo_snp::add_read(int slot, bam1_t const* read) throw (out_of_range) {
    const uint32_t* cigar = bam1_cigar(read);
    const uint8_t* sequence = bam1_seq(read);
    int len = read->core.n_cigar;
    int rpos = read->core.pos; // reference position
    int qpos = 0; // query position
    int offset = 4;
    if (slot != 0 && slot != 1) {
    }
    int** count;
    if (slot == 0) {
        count = _count1;
    } else if (slot == 1) {
        count = _count2;
    } else {
        throw out_of_range("slot must be 0 or 1");
    }
//    int 
    for (int i = 0; i < len; i++) {
        int op = bam_cigar_op(cigar[i]);
        int slen = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) { // M
//        if (op == BAM_CMATCH) {
            int index = -1;
            //int index = get_index(rpos);
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
                    if (base == 1) { // A
                        count[index][0] ++;
                    } else if (base == 2) { // C
                        count[index][1] ++;
                    } else if (base == 4) { // G
                        count[index][2] ++;
                    } else if (base == 8) { // T
                        count[index][3] ++;
                    } 
                }
                qpos++;
                offset ^= 4;
            }
//            ss << 'M';
        } else if (op == BAM_CINS) {
//            ss << 'I';
            qpos += slen;
            offset ^= ((slen & 1) << 2);
        } else if (op == BAM_CDEL) {
            rpos += slen;
//            ss << 'D';
        } else {
            break;
        }
    }
}

namespace {
    void find_major_alleles(int const* bases, int& b1, int& n1, int& b2, int& n2) {
    }
}

vector<polymorphic_allele> denovo_snp::get_polymorphism(int coverage, double heterogeneity) const {
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
        if ((c1 * heterogeneity < n12 && n22 == 0) 
            || (c2 * heterogeneity < n22 && n12 == 0)) {
            alleles.push_back(polymorphic_allele(_chromosome, _position[i], _count1[i], _count2[i]));
        }
    }
    return alleles;
}


void denovo_snp::enumerate_hetero(int argc, char** argv) throw (exception) {
    try {
        const char* filename1 = get_argument_string(argc, argv, "1", "/mnt/smb/tae/stap/shira/BAM6/Sample6.bam");
        const char* filename2 = get_argument_string(argc, argv, "2", "/mnt/smb/tae/stap/shira/BAM12/Sample12.bam");
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        const char* filename_gtf = get_argument_string(argc, argv, "G", NULL);
        const char* filename_chrm = get_argument_string(argc, argv, "s", NULL);
        int coverage = get_argument_integer(argc, argv, "c", 20);
        double heterozygosity = get_argument_float(argc, argv, "z", 0.5);
        int chunk_size = get_argument_integer(argc, argv, "w",  200000);
        int margin_size = get_argument_integer(argc, argv, "m", 100000);
        int num_files = 2;
        int maximum_reads_in_window = get_argument_integer(argc, argv, "x", 0);
        bool verbose = has_option(argc, argv, "verbose");

        if (verbose) {
            cerr << "filename 1 : " << filename1 << endl;
            cerr << "filename 2 : " << filename2 << endl;
            cerr << "gtf        : " << filename_gtf << endl;
            cerr << "heterozygosity : " << heterozygosity << endl;
            cerr << "coverage   : " << coverage << endl;
            cerr << "chunk size : " << chunk_size << endl;
            cerr << "margin     : " << margin_size << endl;
            cerr << "max reads  : " << maximum_reads_in_window << endl;
            cerr << "output     : " << (filename_output == NULL ? "stdout" : filename_output) << endl;
        }

        //cout << "GTF\n";
        gtffile* gtf = gtffile::load_gtf(filename_gtf);

        if (true) {
            //cout << "test\n";
            denovo_snp* snp = new denovo_snp(gtf, 1, 105000000, 119000000);
            snp->set_scope(1, 115000000, 120000000);
            delete snp;
            return;
        }

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
            bam_header_t* h = headers[0];
            for (int i = 0; i < h->n_targets; i++) {
                const char* hn = h->target_name[i];
                if (std::strncmp(hn, "chr", 3) == 0) {
                    hn += 3;
                } 
                int code = convert_chromosome_to_code(hn);
                bc2cc[i] = code;
            }
        }
        int current_chromosome = -1;
        int next_chromosome = -1;
        int position = 0;
        int next_position = position + chunk_size;
        int steps = 0;
        vector<chromosome_seq*> chromosomes = chromosome_seq::load_genome(filename_chrm);
        map<int,chromosome_seq const*> bam2chrm = chromosome_seq::map_chromosome(headers[0], chromosomes);

        denovo_snp* detector = NULL;

//        bool debug = false;
        
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
                        int chromosome_code = 0;
                        detector = new denovo_snp(gtf, chromosome_code, position - margin_size, position + chunk_size + margin_size);
                    }
                    
                    if ((maximum_reads_in_window == 0 || detector->mapped() < maximum_reads_in_window) && detector->has_buffer()) {
                        detector->add_read(i, r);
                    }
                    status[i] = 0;
                }
            }
            
            // detect gaps and insertions
            if (current_chromosome >= 0) {
                if (verbose) {
                    if (++steps % 1000 == 0) {
                        cerr << " " << headers[0]->target_name[current_chromosome] << ":" << position << "-" << next_position << " ";
                        //cerr << current_chromosome << ":" << position << "-" << next_position << " ";
                        for (int i = 0; i < num_files; i++) {
                            cerr << " " << i << ":";// << ;//detectors[i]->size();
                        }
                        cerr << "     \r";
                    }
                }
            }

            if (finished) {
                break;
            }

            if (chromosome_change) {
                int pos_min = numeric_limits<int>::max();
                if (verbose) {
                    int total_reads = 0;
                    for (int i = 0; i < num_files; i++) {
                        //total_reads += detectors[i]->size();
                    }
                    cerr << "change chromosome to " << headers[0]->target_name[next_chromosome] << ", sweep " << total_reads << "reads from " << position << "\r";
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
                } else {
                    break;
                }
            } else { // next position
                for (int i = 0; i < num_files; i++) {
                    //detectors[i]->sweep(current_chromosome, 0, next_position);// - margin_size);
                }
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

