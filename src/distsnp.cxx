#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <fstream>

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

string dbsnp_file::get_cache_filename(const char* filename) {
    int slen = strlen(filename);
    string dirname;
    string basename;
    for (int i = slen - 1; i >= 0; i--) {
        if (filename[i] == FILE_SEPARATOR) {
            dirname = string(filename).substr(0, i + 1);
            basename = string(filename).substr(i + 1);
        }
    }
    return dirname + "." + basename + ".snpcache";
}

void dbsnp_file::save_cache(const char* filename) const throw (exception) {
    ofstream fo(filename);
    if (!fo.is_open()) {
        throw invalid_argument("cannot save cache");
    }
    map<string,vector<cache_position const*> > indicators;
    for (int i = 0; i < (int)_indicators.size(); i++) {
        const cache_position* cp = _indicators[i];
        map<string,vector<cache_position const*> >::iterator it = indicators.find(cp->chromosome);
        if (it == indicators.find(cp->chromosome)) {
            vector<cache_position const*> pos;
            pos.push_back(cp);
            indicators.insert(make_pair(cp->chromosome, cp));
        } else {
            it->second.push_back(cp);
        }
    }
    for (map<string,vector<cache_position const*> >::const_iterator it = indicators.begin(); it != indicators.end(); it++) {
        string chromosome = it->first;
        const vector<cache_position const*>& pos = it->second;
        size_t size = chromosome.size() + 1;
        //fo.write(&size, sizeof(size_t));
        fo.write(reinterpret_cast<char*>(&size), sizeof(size_t));
        fo.write(chromosome.c_str(), sizeof(const char) * (chromosome.size() + 1));
        //size = chromosome.size()
        //fo.write((size_t)chromosome.size(), sizeof(char));
        size = pos.size() * 2 * sizeof(size_t);
        fo.write(reinterpret_cast<char*>(&size), sizeof(size_t));
        for (int i = 0; i < (int)pos.size(); i++) {
            fo.write(reinterpret_cast<const char*>(&pos[i]->position), sizeof(size_t));
            fo.write(reinterpret_cast<const char*>(&pos[i]->file_position), sizeof(size_t));
        }
    }
    fo.close();
}

dbsnp_file* dbsnp_file::load_cache(const char* filename) throw (exception) {
    dbsnp_file* snpfile = NULL;
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
    dbsnp_file* snpfile = new dbsnp_file();//NULL;
    size_t next_position = 0;
    size_t interval = 10000;
    string prev_chrom = "";
    while (!fi.eof()) {
        string line;
        size_t fpos = fi.tell();
        getline(fi, line);
        if (snpfile == NULL) {//strains.size() == 0) {
            if (line.find("#CHROM") == 0) {
                vector<string> strains;
                vector<string> items = split_items(line, '\t');
                for (int i = 9; i < (int)items.size(); i++) {
                    strains.push_back(items[i]);
                }
                if (strains.size() > 0) {
                    snpfile = new dbsnp_file(filename, strains);
                }
            }
        } else {
            vector<string> items = split_items(line, '\t');
            if (items.size() >= 9 + snpfile->strain_number()) {
                string chrom = items[0];
                if (chrom != prev_chrom) {
                    next_position = 0;
                }
                size_t position = atoi(items[1]);
                if (position >= next_position) {
                    snpfile->add_cache(new cache_position(chrom, position, fpos));
                    next_position += interval;
                }
            }
        }
    }
    snpfile->save_cache(cache.c_str());
    return snpfile;
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
