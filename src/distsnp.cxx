#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <cstdlib>

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
using std::strlen;
using std::exception;
using std::invalid_argument;
using std::runtime_error;
using std::pair;
using std::atoi;

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
    cout << probes.size() << " chromosomes\n";
    fo.write(reinterpret_cast<char const*>(&MAGIC), sizeof(unsigned int));
    int num_chromosomes = probes.size();
    fo.write(reinterpret_cast<char const*>(&num_chromosomes), sizeof(int));
    size_t location = header_size;
    for (map<string,int>::const_iterator it = probe_size.begin(); it != probe_size.end(); it++) {
        int num_probes = it->second;
        int slen = it->first.size();
        fo.write(reinterpret_cast<char const*>(&slen), sizeof(int));
        fo.write(it->first.c_str(), sizeof(char) * (slen + 1));
        fo.write(reinterpret_cast<char const*>(&location), sizeof(size_t));
        //fo.write(it->
        location += sizeof(size_t) * 2;
    }
    
    int len = _filename.size() + 1;
    fo.write(reinterpret_cast<char const*>(&len), sizeof(int));
    fo.write(_filename.c_str(), sizeof(char) * len);
    int num_strains = _strains.size();
    fo.write(reinterpret_cast<char const*>(&num_strains), sizeof(int));
    for (int i = 0; i < num_strains; i++) {
        int len = _strains[i].size();
        const char* ptr = _strains[i].c_str();
        fo.write(reinterpret_cast<char const*>(&len), sizeof(int));
        fo.write(reinterpret_cast<char const*>(ptr), sizeof(char) * (_strains[i].size() + 1));
    }
    //char zero = '\0';
    fo.seekp(header_size);
    // for (int i = sizeof(int) + sizeof(unsigned int) ; i < header_size; i++) {
    //     fo.write(&zero, sizeof(char));
    // }
    // 0-4 length of chromosome label
    // 4-x chromosome
    // x-x+4 number of indicators
    // x+5-x+8 position of snp
    // x+8-x+12 position of file
    for (map<string,vector<cache_position const*> >::const_iterator it = probes.begin(); it != probes.end(); it++) {
        string chromosome = it->first;
        cout << chromosome << " : " << it->second.size() << " probes\n";
        const vector<cache_position const*>& pos = it->second;
        int size = chromosome.size() + 1;
        //fo.write(&size, sizeof(size_t));
        fo.write(reinterpret_cast<char*>(&size), sizeof(int));
        fo.write(chromosome.c_str(), sizeof(const char) * (chromosome.size() + 1));
        //size = chromosome.size()
        //fo.write((size_t)chromosome.size(), sizeof(char));
        const vector<cache_position const*>& points = it->second;
        size = pos.size();
        //size = pos.size() * 2 * sizeof(size_t);
        fo.write(reinterpret_cast<char*>(&size), sizeof(int));
        for (int i = 0; i < size; i++) {
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
    // skip header
    fi.seekg(header_size);
    for (int i = 0; i < num_chromosomes; i++) {
        int size_buffer;
        //size_t num_items, number;
        fi.read(reinterpret_cast<char*>(&size_buffer), sizeof(int));
        char* buffer = new char[size_buffer];
        fi.read(buffer, sizeof(char) * size_buffer);
        if (buffer[size_buffer - 1] != '\0') {
            buffer[size_buffer - 1] = '\0';
        }
        string chromosome(buffer);
        delete[] buffer;
        fi.read(reinterpret_cast<char*>(&size_buffer), sizeof(int));
        for (int i = 0; i < size_buffer; i++) {
            size_t pos, fpos;
            fi.read(reinterpret_cast<char*>(&pos), sizeof(size_t));
            fi.read(reinterpret_cast<char*>(&fpos), sizeof(size_t));
            snpfile->_indicators.push_back(new cache_position(chromosome, pos, fpos));
        }
    }
    fi.close();
    return snpfile;
}

dbsnp_file::~dbsnp_file() {
    for (int i = 0; i < (int)_cache.size(); i++) {
        delete _cache[i];
    }
}

vector<dbsnp_locus*> dbsnp_file::load_snps(const string& chromosome, int start, int end) throw (exception) {
    //vector<dbsnp_locus*>
    return _cache;
}

dbsnp_file* dbsnp_file::load_dbsnp(const char* filename) throw (exception) {
    string cache = get_cache_filename(filename);
    cout << cache << endl;
    if (file_exists(cache.c_str())) {
        dbsnp_file* sf = load_cache(cache.c_str());
        return sf;
    }
    ifstream fi(filename);
    if (fi.is_open() == false) {
        throw invalid_argument("cannot open VCF file");
    }
    dbsnp_file* snpfile = NULL;//new dbsnp_file(filename);//NULL;
    size_t next_position = 0;
    size_t interval = 100000;
    string prev_chrom = "";
#define DEBUG
#ifdef DEBUG
    size_t num_lines = 0;
#endif
    while (!fi.eof()) {
        string line;
        size_t fpos = fi.tellg();
        getline(fi, line);
        if (snpfile == NULL) {//strains.size() == 0) {
            if (line.find("#CHROM") == 0) {
                vector<string> strains;
                vector<string> items = split_items(line, '\t');
                cout << line << endl;
                for (int i = 9; i < (int)items.size(); i++) {
                    strains.push_back(items[i]);
                    cout << strains.size() << ":" << items[i] << endl;
                }
                if (strains.size() > 0) {
                    snpfile = new dbsnp_file(filename, strains);
                }
                cout << line << endl;
            }
        } else {
            vector<string> items = split_items(line, '\t');
//            cout << items.size() << ", " << snpfile->strain_number() + 9 << endl;
            if (items.size() >= snpfile->strain_number() + 9) {
#ifdef DEBUG
                if (++num_lines % 1000 == 0) {
                    cerr << " " << (num_lines / 1000) << " " << items[0] << ":" << items[1] << "        \r";
                    if (num_lines > 100000) {
                        break;
                    }
                }
#endif
                        
                string chrom = items[0];
                if (chrom != prev_chrom) {
                    next_position = 0;
                    prev_chrom = chrom;
                }
                size_t position = std::atoi(items[1].c_str());
                if (position >= next_position) {
                    cout << snpfile->_indicators.size() << " add probe at " << fpos << " of " << position << " / " << next_position << "    \n";
                    snpfile->add_cache(chrom, position, fpos);//new cache_position(chrom, position, fpos));
                    next_position += interval;
                }
            } else {
                cout << "too few columns " << items.size() << " / " << snpfile->strain_number() << endl;
            }
        }
    }
    cout << "save cache\n";
    snpfile->save_cache(cache.c_str());
    return snpfile;
}

void dbsnp_file::add_cache(const string& chromosome, size_t pos, size_t fpos) {
    _indicators.push_back(new cache_position(chromosome, pos, fpos));
}

namespace {
    void analyze_strain() {
    }
}
    
int main(int argc, char** argv) {
    try {
        const char* filename_snps = get_argument_string(argc, argv, "V", NULL);

        if (true) {
            dbsnp_file* dbsnp = dbsnp_file::load_dbsnp(filename_snps);
            delete dbsnp;
            exit(0);
        }
        const char* filename = get_argument_string(argc, argv, "i", NULL);
        int coverage = get_argument_integer(argc, argv, "c", 10);
        double heterozygosity = get_argument_float(argc, argv, "H", 0.25);
        bool verbose = has_option(argc, argv, "verbose");
        int window_size = get_argument_integer(argc, argv, "w", 20000);
        int window_margin = get_argument_integer(argc, argv, "m", 2000);
        //const char* filename_genome = get_argument_string(argc, argv, "g", NULL);
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);

        if (verbose) {
            cerr << "Filename         : " << filename << endl;
            cerr << "Minimum coverage : " << coverage << endl;
            cerr << "Heterozygosity   : " << heterozygosity << endl;
            if (filename_snps != NULL) {
                cerr << "Variation file   : " << filename_snps << endl;
            }
            //cerr << "Genome file      : " << filename_genome << endl;
            cerr << "Coverage         : " << coverage << endl;
            cerr << "Window size      : " << window_size << endl;
        }

        ostream* ost = &cout;
        vector<recfragment*> fragments;
        map<string,const recfragment*> pairs;
        bam1_t* read;
        bam_header_t* header;
        bamFile bamfile;

        if (filename_output != NULL) {
            ofstream* fo = new ofstream(filename_output);
            if (fo->is_open() == false) {
                delete fo;
                throw invalid_argument("cannot open output stream");
            }
            ost = fo;
        }
        // if (verbose) {
        //     cerr << "loading genomic sequences " << flush;
        // }
        //vector<chromosome_seq*> fasta_files = chromosome_seq::load_genome(filename_genome);
//        map<int,pair<int,unsigned char*> >  genome = load_genome(filename_genome);
        // if (verbose) {
        //     cerr << fasta_files.size() << " chromosomes\n";
        // }

        // select standard snps only
        // if (filename_snps != NULL) {
        //     remove_nonestablished_bases(filename_snps, fasta_files);
        // }

        bamfile = bam_open(filename, "rb");
        header = bam_header_read(bamfile);
        read = bam_init1();

/*
first
     <--------------window_size----------------->
     |------|============================|------|
      margin         processing           margin
region_start                               region_end


second
    |==================|-----|
   region_start             region_end



retain_start
analysis_start
analysis_end
retain_end

 */
        int retain_start = 0;
        int retain_end = window_size;
        int analysis_start = retain_start + window_margin;
        int analysis_end = retain_end - window_margin;
        int current_bamchrm = -1;
        string chromosome;
        //chromosome_seq const* chromosome = NULL;
        int current_chromosome = -1;
        int next_chromosome = -1;
        //int chromosome_length = 0;
        string chromosome_name;
        size_t max_fragments = 10000;
        //const unsigned char* chromosome_sequence = NULL;

        while (bam_read1(bamfile, read) > 0) {
            bool processing = false;

            if (current_bamchrm != read->core.tid) {
                // finish chromosome
                processing = true;
                const char* name = header->target_name[read->core.tid];
                chromosome = name;
                //next_chromosome = chromosome_seq::get_chromosome_code(strlen(name), name);
            } else if (read->core.pos > retain_end) {
                // process current region
                if (fragments.size() > 0) {
                    processing = true;
                }
            }
            //cout << fragments.size() << endl;
            if (processing) {
                // analysis
                if (fragments.size() > coverage && chromosome != "") {//_current != NULL) {//chromosome != NULL) {
                    if (verbose) {
                        cerr << " " << chromosome << ":" << analysis_start << "-" << analysis_end << ":" << fragments.size() << "            \r";
                    }
                    if (fragments.size() >= max_fragments) {
                        cerr << "skip " << chromosome << ":" << analysis_start << "-" << analysis_end << " because of too much fragments " << fragments.size() << "       \n";
                    } else {
                        analyze_strain();
                        //detect_recombination(fragments, chromosome, analysis_start, analysis_end, coverage, heterozygosity, *ost);
                        *ost << flush;
                    }
                }

                // change chromosome
                if (next_chromosome != current_chromosome) {
                    //cerr << "change chromosoe\n";
                    chromosome = next_chromosome;
                    // chromosome = "";//NULL;
                    // for (int i = 0; i < (int)fasta_files.size(); i++) {
                    //     if (fasta_files[i]->code() == next_chromosome) {
                    //         chromosome = fasta_files[i];
                    //         break;
                    //     }
                    // }
                    // if (chromosome == NULL) {
                    //     if (verbose) {
                    //         cerr << "cannot find " << next_chromosome << endl;
                    //     }
                    //     current_chromosome = -1;
                    // } else {
                    //     current_chromosome = next_chromosome;
                    // }
                    
                    retain_start = 0;
                    analysis_start = 0;
                    analysis_end = window_size - window_margin;
                    retain_end = window_size;

                    // delete all fragments
                    for (vector<recfragment*>::iterator it = fragments.begin(); it != fragments.end(); it++) {
                        delete *it;
                    }
                    fragments.erase(fragments.begin(), fragments.end());
                    current_bamchrm = read->core.tid;
                    current_chromosome = next_chromosome;
                } else {
                    // delete outsiders
                    {
                        for (int i = 0; i < (int)fragments.size(); i++) {
                            if (fragments[i]->position3() < analysis_start) {
                                delete fragments[i];
                                fragments[i] = NULL;
                            } else if (fragments[i]->position5() >= analysis_start) {
                                break;
                            }
                        }
                        int index = 0;
                        for (int i = 0; i < (int)fragments.size(); i++) {
                            if (fragments[i] != NULL) {
                                fragments[index++] = fragments[i];
                            }
                        }
                        fragments.erase(fragments.begin() + index, fragments.end());
                    }

                    // revise analyzing window
                    analysis_start = analysis_end;
                    if (fragments.size() > 0 && analysis_start < fragments[0]->position5()) {
                        analysis_start = fragments[0]->position5();
                    }
                    retain_start = analysis_start - window_margin;
                    retain_end = retain_start + window_size;
                    analysis_end = retain_end - window_margin;
                }
            }
            if (current_chromosome > 0 && fragments.size() < max_fragments) {
                recfragment* frag = new recfragment(current_chromosome, read);
                fragments.push_back(frag);
            }
        }

        // process remaining
        if (fragments.size() > coverage && chromosome != "") {
            analyze_strain();
            //detect_recombination(fragments, chromosome, analysis_start, analysis_end,
            //coverage, heterozygosity, *ost);
        }


        // close output file
        if (filename_output != NULL) {
            dynamic_cast<ofstream*>(ost)->close();
            delete ost;
        }

        // clean up fragments
        for (vector<recfragment*>::iterator it = fragments.begin(); it != fragments.end(); it++) {
            delete *it;
        }
        // clean up sequences
        // for (vector<chromosome_seq*>::iterator it = fasta_files.begin(); it != fasta_files.end(); it++) {
        //     delete *it;
        // }

        // clean up BAM
        bam_destroy1(read);
        bam_header_destroy(header);
        bam_close(bamfile);

        return 0;
    } catch (exception& e) {
        cerr << e.what() << endl;
        return -1;
    }
}
