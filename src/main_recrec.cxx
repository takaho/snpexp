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
using std::string;
using std::vector;
using std::map;
using std::sort;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::stringstream;


#include <tktools.hxx>
#include <gtf.hxx>
#include <recrec.hxx>

using namespace tkbio;

using tktools::util::get_argument_string;
using tktools::util::get_argument_integer;
using tktools::util::get_argument_float;
using tktools::util::has_option;
using tktools::split_items;
using tktools::bio::convert_chromosome_to_code;
using tktools::bio::convert_code_to_chromosome;

namespace {
    vector<hetero_locus>
    scan_heterozygous_loci(vector<recfragment*>& fragments,
                           chromosome_seq const* chromosome,
                           int start, int end,
                           int minimum_coverage=10,
                           float hetero_threshold=0.2f) throw (exception) {
        int freq[5];
        vector<hetero_locus> candidates;
        bool gapped = false;
        for (int pos = start; pos < end; pos++) {
            int refid = chromosome->get_base_id(pos);
            if (refid < 0) continue;
            for (int i = 0; i < 5; i++) freq[i] = 0;
            int coverage = 0;
            for (int j = 0; j < (int)fragments.size(); j++) {
                const recfragment* fr = fragments[j];
                if (fr->position5() <= pos && pos < fr->position3()) {
                    int num = 0;
                    int index = fr->get_base_id(pos, num);
                    if (index >= 0) {
                        freq[index] += num;
                    }
                    //coverage += num;
                }
            }
            int altid = -1;
            int altnum = 0;
            int refnum = freq[refid];
            for (int i = 0; i < 5; i++) {
                if (i != refid && altnum < freq[i]) {
                    altid = i;
                    altnum = freq[i];
                }
            }
            if (altid == 4) {
                if (gapped) {
                    continue;
                } else {
                    gapped = true;
                }
            } else {
                gapped = false;
            }
            coverage = freq[refid] + altnum;
            if (coverage >= minimum_coverage && refnum > 1 && altnum > 1) {
                int threshold = (int)(coverage * hetero_threshold + 0.5);
                if (threshold <= altnum && threshold <= refnum) {
                    candidates.push_back(hetero_locus(chromosome->code(), pos, refid, refnum, altid, altnum));
                }
            }
        }
        return candidates;
    }

    void detect_recombination(vector<recfragment*>& fragments,
                              chromosome_seq const* chromosome,
                              int start, int end,
                              int minimum_coverage=15,
                              float hetero_threshold=0.3f,
                              ostream& ost=cout) throw (exception) {
        // bind pairs
        recfragment::bundle_pairs(fragments);

        vector<hetero_locus> hetero_loci = scan_heterozygous_loci(fragments, chromosome, start, end, minimum_coverage, hetero_threshold);
        float diff_degree = 0.8f;
        bool displayed = false;
        for (int i = 0; i < (int)fragments.size(); i++) {
            //cout << i << ":" << flush;
            pair<int,int> site = fragments[i]->get_recombination(hetero_loci, diff_degree);
            if (site.first < site.second) {
                if (!displayed) {
                    ost << "\n";
                    //cout << chromosome_name << ":" << start << "-" << end << endl;
                    displayed = true;
                }
                string pattern = fragments[i]->get_recombination_pattern(hetero_loci);
                ost << chromosome->name() << ":" << site.first << "-" << site.second 
                    << "\t" << chromosome->name() << ":" << fragments[i]->position5() << "-" << fragments[i]->position3()
                    << "\t" << pattern << "\t" << fragments[i]->name() << endl;
            // } else {
            //     cout << "no_recombination";
            }
            //cout << endl;
        }
    }

    void remove_nonestablished_bases(const char* vcffile, vector<chromosome_seq*>& chromosomes) throw (exception) {
        cout << vcffile << endl;
        ifstream fs(vcffile);
        chromosome_seq* current = NULL;
        string chrom_current;
        if (!fs.is_open()) {
            throw invalid_argument("cannot open vcf file");
        }
        //cerr << "open snps\n";
        size_t num_lines = 0;
        int last_position = 0;
        while (!fs.eof()) {
            string line;
            getline(fs, line);
            if (line.c_str()[0] == '#') {
                continue;
            }
            if (++num_lines % 1000 == 0) {
                cerr << " " << num_lines / 1000 << " " << chrom_current
                     << ":" << last_position
                     << "         \r";
            }
            int col = 0;
            int pos = 0;
            //cout << line << endl;
            char const* ptr = line.c_str();
            for (int i = 0; i < (int)line.size(); i++) {
                char c = ptr[i];
                //cout << i << ":" << c << "\t" << col << endl;
                if (c == '\t') {
                    if (col == 0) {
                        string chrom_line = line.substr(0, i);
                        if (chrom_line != chrom_current) {
                            int code = convert_chromosome_to_code(chrom_line.c_str());
                            current = NULL;
                            for (int j = 0; j < (int)chromosomes.size(); j++) {
                                if (code == chromosomes[j]->code()) {
                                    current = chromosomes[j];
                                    cerr << "importing SNPs from " << current->name() << "          \r";
                                    break;
                                }
                            }
                            chrom_current = chrom_line;
                            last_position = 0;
                            if (current == NULL) break;
                        } else if (current == NULL) {
                            break;
                        }
                    } else if (col == 1) {
                        int position = atoi(line.c_str() + pos) - 1;
                        vector<string> items = split_items(line, '\t');
                        //cout << current->name() << ":" << position << "\t" << current->get_base(position) << "\t" << items[3] << "\t" << items[4] << endl;
                        current->mask(last_position, position);
                        // for (int pos = last_position; pos < position; pos++) {
                        //     current->set_base(pos, 'N');
                        // }
                        last_position = position + 1;
                        break;
                    }
                    pos = i + 1;
                    col += 1;
                }
            }
        }
        fs.close();
    }     
}
    

int main(int argc, char** argv) {
    try {
        const char* filename = get_argument_string(argc, argv, "i", NULL);
        int coverage = get_argument_integer(argc, argv, "c", 10);
        double heterozygosity = get_argument_float(argc, argv, "H", 0.25);
        bool verbose = has_option(argc, argv, "verbose");
        int window_size = get_argument_integer(argc, argv, "w", 20000);
        int window_margin = get_argument_integer(argc, argv, "m", 2000);
        const char* filename_genome = get_argument_string(argc, argv, "g", NULL);
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        const char* filename_snps = get_argument_string(argc, argv, "V", NULL);

        if (verbose) {
            cerr << "Filename         : " << filename << endl;
            cerr << "Minimum coverage : " << coverage << endl;
            cerr << "Heterozygosity   : " << heterozygosity << endl;
            if (filename_snps != NULL) {
                cerr << "Variation file   : " << filename_snps << endl;
            }
            cerr << "Genome file      : " << filename_genome << endl;
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
        if (verbose) {
            cerr << "loading genomic sequences " << flush;
        }
        vector<chromosome_seq*> fasta_files = chromosome_seq::load_genome(filename_genome);
//        map<int,pair<int,unsigned char*> >  genome = load_genome(filename_genome);
        if (verbose) {
            cerr << fasta_files.size() << " chromosomes\n";
        }

        // select standard snps only
        if (filename_snps != NULL) {
            remove_nonestablished_bases(filename_snps, fasta_files);
        }

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
        chromosome_seq const* chromosome = NULL;
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
                next_chromosome = chromosome_seq::get_chromosome_code(strlen(name), name);
            } else if (read->core.pos > retain_end) {
                // process current region
                if (fragments.size() > 0) {
                    processing = true;
                }
            }
            //cout << fragments.size() << endl;
            if (processing) {
                // analysis
                if (fragments.size() > coverage && chromosome != NULL) {
                    if (verbose) {
                        cerr << " " << chromosome->name() << ":" << analysis_start << "-" << analysis_end << ":" << fragments.size() << "            \r";
                    }
                    if (fragments.size() >= max_fragments) {
                        cerr << "skip " << chromosome->name() << ":" << analysis_start << "-" << analysis_end << " because of too much fragments " << fragments.size() << "       \n";
                    } else {
                        detect_recombination(fragments, chromosome, analysis_start, analysis_end, coverage, heterozygosity, *ost);
                        *ost << flush;
                    }
                }

                // change chromosome
                if (next_chromosome != current_chromosome) {
                    //cerr << "change chromosoe\n";
                    chromosome = NULL;
                    for (int i = 0; i < (int)fasta_files.size(); i++) {
                        if (fasta_files[i]->code() == next_chromosome) {
                            chromosome = fasta_files[i];
                            break;
                        }
                    }
                    if (chromosome == NULL) {
                        if (verbose) {
                            cerr << "cannot find " << next_chromosome << endl;
                        }
                        current_chromosome = -1;
                    } else {
                        current_chromosome = next_chromosome;
                    }
                    
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
        if (fragments.size() > coverage && chromosome != NULL) {
            detect_recombination(fragments, chromosome, analysis_start, analysis_end,
                                 coverage, heterozygosity, *ost);
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
        for (vector<chromosome_seq*>::iterator it = fasta_files.begin(); it != fasta_files.end(); it++) {
            delete *it;
        }

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
