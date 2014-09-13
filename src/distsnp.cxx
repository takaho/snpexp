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
const char FILE_SEPARATOR = '/';
dbsnp_locus::dbsnp_locus(size_t position, string reference, string alternative, int num_strains) {
    _position = position;
    _reference = reference;
    _alternative = alternative;
    _num_strains = num_strains;
    _strains = new unsigned char[_num_strains];
    for (int i = 0; i < _num_strains; i++) {
        _strains[i] = 0x00;
    }
}

dbsnp_locus::~dbsnp_locus() {
    delete[] _strains;
}

void dbsnp_locus::set_genotype(int index, char const* info) {
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


dbsnp_file::dbsnp_file(const char* filename, const vector<string>& strains) {
    _strains = strains;
    _filename = filename;
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
    int header_size = 1024;
}

void dbsnp_file::save_cache(const char* filename) const throw (exception) {
    ofstream fo(filename);
    if (!fo.is_open()) {
        throw invalid_argument("cannot save cache");
    }

    // assign probes to each chromosome
    map<string,vector<cache_position const*> > probes;
    map<string,int> probe_size;
    for (int i = 0; i < (int)_indicators.size(); i++) {
        const cache_position* cp = _indicators[i];
        map<string,vector<cache_position const*> >::iterator it = probes.find(cp->chromosome);
        if (it == probes.end()) {//find(cp->chromosome)) {
            vector<cache_position const*> pos;
            pos.push_back(cp);
            probes[cp->chromosome] = pos;//.insert(make_pair(cp->chromosome, cp));
            probe_size[cp->chromosome] = 1;
        } else {
            it->second.push_back(cp);
            probe_size[cp->chromosome]++;
        }
    }
    cout << " " << probes.size() << " chromosomes\n";

    // magic number
    fo.write(reinterpret_cast<char const*>(&MAGIC), sizeof(unsigned int));

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
    
    for (map<string,int>::const_iterator it = probe_size.begin(); it != probe_size.end(); it++) {
        int num_probes = it->second;
        int slen = it->first.size() + 1;
        fo.write(reinterpret_cast<char const*>(&slen), sizeof(int));
        fo.write(it->first.c_str(), sizeof(char) * slen);
        fo.write(reinterpret_cast<char const*>(&location), sizeof(int));
        fo.write(reinterpret_cast<char*>(&num_probes), sizeof(int));
        location += sizeof(size_t) * 2 * it->second + sizeof(int) * 2;
    }
    
    fo.seekp(header_size);
    for (map<string,vector<cache_position const*> >::const_iterator it = probes.begin(); it != probes.end(); it++) {
        string chromosome = it->first;
        //cout << chromosome << " : " << it->second.size() << " probes\n";
        const vector<cache_position const*>& points = it->second;
        for (int i = 0; i < (int)points.size(); i++) {//size; i++) {
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
    //size_t size;
    unsigned int magic;
    fi.read(reinterpret_cast<char*>(&magic), sizeof(unsigned int));
    if (magic != MAGIC) {
        fi.close();
        throw runtime_error("magic number inconsistent");
    }
    int num_chromosomes;
    fi.read(reinterpret_cast<char*>(&num_chromosomes), sizeof(int));
    int string_length;
    fi.read(reinterpret_cast<char*>(&string_length), sizeof(int));
    //string filename;
    char* filename_vcf = new char[string_length];
    fi.read(filename_vcf, sizeof(char) * string_length);
    //filename = buffer;
    //delete[] buffer;
    //snpfile->_filename = buffer;

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
    for (int i = 0; i < num_chromosomes; i++) {
        int slen;
        char* cbuf;
        int num_probes, location;
        //size_t cpos, fpos;
        fi.read(reinterpret_cast<char*>(&slen), sizeof(int));
        cbuf = new char[slen];
        fi.read(cbuf, sizeof(char) * slen);
        fi.read(reinterpret_cast<char*>(&location), sizeof(int));
        fi.read(reinterpret_cast<char*>(&num_probes), sizeof(int));
        //cerr << cbuf << "\t" << num_probes << " probes from " << hex << location << dec << endl;
        size_t ptr = fi.tellg();
        fi.seekg(location);
        for (int j = 0; j < num_probes; j++) {
            size_t cpos, fpos;
            fi.read(reinterpret_cast<char*>(&cpos), sizeof(size_t));
            fi.read(reinterpret_cast<char*>(&fpos), sizeof(size_t));
            snpfile->_indicators.push_back(new cache_position(cbuf, cpos, fpos));
        }
        fi.seekg(ptr);
        delete[] cbuf;
    }

    sort(snpfile->_indicators.begin(), snpfile->_indicators.end(), cache_position_comparator);

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
    if (_cache.size() > 0) {
      if ((int)_cache[0]->_position >= start && (int)_cache[_cache.size() - 1]->_position <= end) {
            return;
        }
    }
    int left = 0;
    int right = _indicators.size();
    int center = (left + right) / 2;
    while (true) {
        const cache_position* probe = _indicators[center];
        //cout << probe->chromosome << "\t" << probe->position << "\t" << probe->file_position << endl;
        if (probe->chromosome < chromosome) {
            left = center + 1;
        } else if (probe->chromosome > chromosome) {
            right = center; 
        } else if ((int)probe->position < start) {
            left = center + 1;
        } else if ((int)probe->position > end) {
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
    size_t left_limit = header_size;
    size_t right_limit = 0;//(int)_indicators.size();
    while (index >= 0) {
        const cache_position* probe = _indicators[index];
        if (probe->chromosome != chromosome || (int)probe->position < start) {
            left_limit = probe->file_position;
            //cout << index << ":" << left_limit << endl;
            break;
        }
        index--;
    }

    index = center + 1;
    while (index < (int)_indicators.size()) {
        const cache_position* probe = _indicators[index];
        //cout << probe->position << ", " << probe->file_position << endl;
        if (probe->chromosome != chromosome || (int)probe->position > end) {
            right_limit = probe->file_position;
            //cout << index << ":" << right_limit << endl;
            break;
        }
        index++;
    }
    for (int i = 0; i < (int)_cache.size(); i++) {
        delete _cache[i];
    }
    _cache.erase(_cache.begin(), _cache.end());
    ifstream fi(_filename.c_str());
    //if (right_limit = 
    fi.seekg(left_limit);
    int col_strain = 9;
    //cout << center << ":" << left_limit << ", " << right_limit << "       " << endl;
    for (;;) {
        string line;
        getline(fi, line);
        vector<string> items = split_items(line, '\t');
        if (items.size() >= _strains.size() + col_strain) {
            int pos = atoi(items[1].c_str());
            dbsnp_locus* snp = new dbsnp_locus(pos, items[3], items[4], _strains.size());
            for (int i = 0; i < (int)_strains.size(); i++) {
                snp->set_genotype(i, items[i + col_strain].c_str());
            }
            _cache.push_back(snp);
            //cout << snp->to_string(chromosome) << "\n";
        }
        if (fi.eof() || (right_limit > 0 && fi.tellg() >= (long long)right_limit)) break;
    }
    fi.close();
}

vector<dbsnp_locus const*> dbsnp_file::get_snps(string chromosome, int start, int end) const throw (exception) {
    string cname;
    if (chromosome.find("chr") == 0) {
        cname = chromosome.substr(3, chromosome.size());
    } else {
        cname = chromosome;
    }
    vector<dbsnp_locus const*> snps;
    const_cast<dbsnp_file*>(this)->load_snps(cname, start, end);
    for (int i = 0; i < (int)(_cache.size()); i++) {
        dbsnp_locus const* cp = _cache[i];
        if (start <= (int)cp->_position && (int)cp->_position <= end) {
            snps.push_back(cp);
        }
    }
    return snps;
}

dbsnp_file* dbsnp_file::load_dbsnp(const char* filename, bool force_uncached) throw (exception) {
    string cache = get_cache_filename(filename);
    if (!force_uncached) {
        //cout << cache << endl;
        if (file_exists(cache.c_str())) {
            dbsnp_file* sf = load_cache(cache.c_str());
            return sf;
        }
    }
    ifstream fi(filename);
    if (fi.is_open() == false) {
        throw invalid_argument("cannot open VCF file");
    }
    dbsnp_file* snpfile = NULL;//new dbsnp_file(filename);//NULL;
    size_t next_position = 0;
    size_t interval = 100000;
    string prev_chrom = "";
//#define DEBUG
//#ifdef DEBUG
    size_t num_lines = 0;
//#endif
    while (!fi.eof()) {
        string line;
        size_t fpos = fi.tellg();
        getline(fi, line);
        if (snpfile == NULL) {//strains.size() == 0) {
            if (line.find("#CHROM") == 0) {
                vector<string> strains;
                vector<string> items = split_items(line, '\t');
                //cout << line << endl;
                for (int i = 9; i < (int)items.size(); i++) {
                    strains.push_back(items[i]);
                    //cout << strains.size() << ":" << items[i] << endl;
                }
                if (strains.size() > 0) {
                    snpfile = new dbsnp_file(filename, strains);
                }
                //cout << line << endl;
            }
        } else {
            vector<string> items = split_items(line, '\t');
//            cout << items.size() << ", " << snpfile->strain_number() + 9 << endl;
            if ((int)items.size() >= snpfile->strain_number() + 9) {
//#ifdef DEBUG
                if (++num_lines % 1000 == 0) {
                    cerr << " " << (num_lines / 1000) << " " << items[0] << ":" << items[1] << "        \r";
                    // if (num_lines > 1000000) {
                    //     break;
                    // }
                }
//#endif
                        
                string chrom = items[0];
                if (chrom != prev_chrom) {
                    next_position = 0;
                    prev_chrom = chrom;
                }
                size_t position = std::atoi(items[1].c_str());
                if (position >= next_position) {
                    cout << snpfile->_indicators.size() << " add probe at " << fpos << " of " << position << " / " << next_position << "    \r";
                    snpfile->add_cache(chrom, position, fpos);//new cache_position(chrom, position, fpos));
                    next_position += interval;
                }
            } else {
                cout << "too few columns " << items.size() << " / " << snpfile->strain_number() << endl;
            }
        }
    }
    //cout << "save cache\n";
    snpfile->save_cache(cache.c_str());
    return snpfile;
}

void dbsnp_file::add_cache(const string& chromosome, size_t pos, size_t fpos) {
    _indicators.push_back(new cache_position(chromosome, pos, fpos));
}

    
// // int main(int argc, char** argv) {
// //     try {
// //         const char* filename_snps = get_argument_string(argc, argv, "V", NULL);

// //         if (true) {
// //             dbsnp_file* dbsnp = dbsnp_file::load_dbsnp(filename_snps, has_option(argc, argv, "f"));
// //             int start = 10000000;
// //             int end = start + 10000;
// //             string chromosome = "17";
// //             vector<dbsnp_locus const*> snps = dbsnp->get_snps("17", start, end);
// //             for (int i = 0; i < (int)snps.size(); i++) {
// //                 cout << snps[i]->to_string(chromosome) << endl;
// //             }
// //             //cout << snps.size() << endl;
            
// //             delete dbsnp;
// //             exit(0);
// //         }
// //         const char* filename = get_argument_string(argc, argv, "i", NULL);
// //         int coverage = get_argument_integer(argc, argv, "c", 10);
// //         double heterozygosity = get_argument_float(argc, argv, "H", 0.25);
// //         bool verbose = has_option(argc, argv, "verbose");
// //         int window_size = get_argument_integer(argc, argv, "w", 20000);
// //         int window_margin = get_argument_integer(argc, argv, "m", 2000);
// //         //const char* filename_genome = get_argument_string(argc, argv, "g", NULL);
// //         const char* filename_output = get_argument_string(argc, argv, "o", NULL);

// //         if (verbose) {
// //             cerr << "Filename         : " << filename << endl;
// //             cerr << "Minimum coverage : " << coverage << endl;
// //             cerr << "Heterozygosity   : " << heterozygosity << endl;
// //             if (filename_snps != NULL) {
// //                 cerr << "Variation file   : " << filename_snps << endl;
// //             }
// //             //cerr << "Genome file      : " << filename_genome << endl;
// //             cerr << "Coverage         : " << coverage << endl;
// //             cerr << "Window size      : " << window_size << endl;
// //         }

// //         ostream* ost = &cout;
// //         vector<recfragment*> fragments;
// //         map<string,const recfragment*> pairs;
// //         bam1_t* read;
// //         bam_header_t* header;
// //         bamFile bamfile;

// //         if (filename_output != NULL) {
// //             ofstream* fo = new ofstream(filename_output);
// //             if (fo->is_open() == false) {
// //                 delete fo;
// //                 throw invalid_argument("cannot open output stream");
// //             }
// //             ost = fo;
// //         }
// //         // if (verbose) {
// //         //     cerr << "loading genomic sequences " << flush;
// //         // }
// //         //vector<chromosome_seq*> fasta_files = chromosome_seq::load_genome(filename_genome);
// // //        map<int,pair<int,unsigned char*> >  genome = load_genome(filename_genome);
// //         // if (verbose) {
// //         //     cerr << fasta_files.size() << " chromosomes\n";
// //         // }

// //         // select standard snps only
// //         // if (filename_snps != NULL) {
// //         //     remove_nonestablished_bases(filename_snps, fasta_files);
// //         // }

// //         bamfile = bam_open(filename, "rb");
// //         header = bam_header_read(bamfile);
// //         read = bam_init1();

// // /*
// // first
// //      <--------------window_size----------------->
// //      |------|============================|------|
// //       margin         processing           margin
// // region_start                               region_end


// // second
// //     |==================|-----|
// //    region_start             region_end



// // retain_start
// // analysis_start
// // analysis_end
// // retain_end

// //  */
// //         int retain_start = 0;
// //         int retain_end = window_size;
// //         int analysis_start = retain_start + window_margin;
// //         int analysis_end = retain_end - window_margin;
// //         int current_bamchrm = -1;
// //         string chromosome;
// //         //chromosome_seq const* chromosome = NULL;
// //         int current_chromosome = -1;
// //         int next_chromosome = -1;
// //         //int chromosome_length = 0;
// //         string chromosome_name;
// //         size_t max_fragments = 10000;
// //         //const unsigned char* chromosome_sequence = NULL;

// //         while (bam_read1(bamfile, read) > 0) {
// //             bool processing = false;

// //             if (current_bamchrm != read->core.tid) {
// //                 // finish chromosome
// //                 processing = true;
// //                 const char* name = header->target_name[read->core.tid];
// //                 chromosome = name;
// //                 //next_chromosome = chromosome_seq::get_chromosome_code(strlen(name), name);
// //             } else if (read->core.pos > retain_end) {
// //                 // process current region
// //                 if (fragments.size() > 0) {
// //                     processing = true;
// //                 }
// //             }
// //             //cout << fragments.size() << endl;
// //             if (processing) {
// //                 // analysis
// //                 if (fragments.size() > coverage && chromosome != "") {//_current != NULL) {//chromosome != NULL) {
// //                     if (verbose) {
// //                         cerr << " " << chromosome << ":" << analysis_start << "-" << analysis_end << ":" << fragments.size() << "            \r";
// //                     }
// //                     if (fragments.size() >= max_fragments) {
// //                         cerr << "skip " << chromosome << ":" << analysis_start << "-" << analysis_end << " because of too much fragments " << fragments.size() << "       \n";
// //                     } else {
// //                         analyze_strain();
// //                         //detect_recombination(fragments, chromosome, analysis_start, analysis_end, coverage, heterozygosity, *ost);
// //                         *ost << flush;
// //                     }
// //                 }

// //                 // change chromosome
// //                 if (next_chromosome != current_chromosome) {
// //                     //cerr << "change chromosoe\n";
//                     chromosome = next_chromosome;
//                     // chromosome = "";//NULL;
//                     // for (int i = 0; i < (int)fasta_files.size(); i++) {
//                     //     if (fasta_files[i]->code() == next_chromosome) {
//                     //         chromosome = fasta_files[i];
//                     //         break;
//                     //     }
//                     // }
//                     // if (chromosome == NULL) {
//                     //     if (verbose) {
//                     //         cerr << "cannot find " << next_chromosome << endl;
//                     //     }
//                     //     current_chromosome = -1;
//                     // } else {
//                     //     current_chromosome = next_chromosome;
//                     // }
                    
//                     retain_start = 0;
//                     analysis_start = 0;
//                     analysis_end = window_size - window_margin;
//                     retain_end = window_size;

//                     // delete all fragments
//                     for (vector<recfragment*>::iterator it = fragments.begin(); it != fragments.end(); it++) {
//                         delete *it;
//                     }
//                     fragments.erase(fragments.begin(), fragments.end());
//                     current_bamchrm = read->core.tid;
//                     current_chromosome = next_chromosome;
//                 } else {
//                     // delete outsiders
//                     {
//                         for (int i = 0; i < (int)fragments.size(); i++) {
//                             if (fragments[i]->position3() < analysis_start) {
//                                 delete fragments[i];
//                                 fragments[i] = NULL;
//                             } else if (fragments[i]->position5() >= analysis_start) {
//                                 break;
//                             }
//                         }
//                         int index = 0;
//                         for (int i = 0; i < (int)fragments.size(); i++) {
//                             if (fragments[i] != NULL) {
//                                 fragments[index++] = fragments[i];
//                             }
//                         }
//                         fragments.erase(fragments.begin() + index, fragments.end());
//                     }

//                     // revise analyzing window
//                     analysis_start = analysis_end;
//                     if (fragments.size() > 0 && analysis_start < fragments[0]->position5()) {
//                         analysis_start = fragments[0]->position5();
//                     }
//                     retain_start = analysis_start - window_margin;
//                     retain_end = retain_start + window_size;
//                     analysis_end = retain_end - window_margin;
//                 }
//             }
//             if (current_chromosome > 0 && fragments.size() < max_fragments) {
//                 recfragment* frag = new recfragment(current_chromosome, read);
//                 fragments.push_back(frag);
//             }
//         }

//         // process remaining
//         if (fragments.size() > coverage && chromosome != "") {
//             analyze_strain();
//             //detect_recombination(fragments, chromosome, analysis_start, analysis_end,
//             //coverage, heterozygosity, *ost);
//         }


//         // close output file
//         if (filename_output != NULL) {
//             dynamic_cast<ofstream*>(ost)->close();
//             delete ost;
//         }

//         // clean up fragments
//         for (vector<recfragment*>::iterator it = fragments.begin(); it != fragments.end(); it++) {
//             delete *it;
//         }
//         // clean up sequences
//         // for (vector<chromosome_seq*>::iterator it = fasta_files.begin(); it != fasta_files.end(); it++) {
//         //     delete *it;
//         // }

//         // clean up BAM
//         bam_destroy1(read);
//         bam_header_destroy(header);
//         bam_close(bamfile);

//         return 0;
//     } catch (exception& e) {
//         cerr << e.what() << endl;
//         return -1;
//     }
// }
