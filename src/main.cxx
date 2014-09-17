#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <set>

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
using std::endl;
using std::set;

#include <tktools.hxx>
#include <gtf.hxx>
#include <recrec.hxx>
#include <fragmentprocessor.hxx>
#include <distsnp.hxx>

using namespace tkbio;

using tktools::util::get_argument_string;
using tktools::util::get_argument_integer;
using tktools::util::get_argument_float;
using tktools::util::has_option;
using tktools::split_items;
using tktools::bio::convert_chromosome_to_code;
using tktools::bio::convert_code_to_chromosome;

namespace {
    bool check_header_consistency(int num, bam_header_t** headers) {
        for (int i = 0; i < num; i++) {
            bam_header_t* hi = headers[i];
            for (int j = 0; j < i; j++) {
                bam_header_t* hj = headers[j];
                if (hi->n_targets != hj->n_targets) {
                    return false;
                }
                for (int k = 0; k < hi->n_targets; k++) {
                    if (strcmp(hi->target_name[k], hj->target_name[k]) != 0) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
}    


namespace {
    void show_help() {
        cerr << "recrec [command] [options] bamfile1 bamfile2 ...\n";
        cerr << "OPTIONS\n";
        cerr << " -c <number>  : minimum coverage for analysis (10)`n";
        cerr << " -w <number>  : window size to analysis chunk, set more than double of fragment size (3000)\n";
        cerr << " -H <ratio>   : heterozygosity to determine SNPs gruop\n";
        cerr << " -g <filename> : genome sequence\n";
        cerr << " -o <filename> : output filename (stdout)\n";
        cerr << " -q <number>   : threshold of mapping quality (10)\n";
        cerr << " -D <number>   : display verbosity (0:all, 1:countable, 2:hetero)\n";
    }
}

int main(int argc, char** argv) {
    try {
        //const char* filename = get_argument_string(argc, argv, "i", NULL);
        string command = argv[1];
        int coverage = get_argument_integer(argc, argv, "c", 10);
        double heterozygosity = get_argument_float(argc, argv, "H", 0.25);
        bool verbose = has_option(argc, argv, "verbose");
        int window_size = get_argument_integer(argc, argv, "w", 3000);
        int window_margin = get_argument_integer(argc, argv, "m", 500);
        const char* filename_genome = get_argument_string(argc, argv, "g", NULL);
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        const char* filename_snps = get_argument_string(argc, argv, "V", NULL);
        size_t max_fragments = get_argument_integer(argc, argv, "x", 10000);
        int mapping_quality = get_argument_integer(argc, argv, "q", 10);
	int display_mode  = get_argument_integer(argc, argv, "D", 1);
        dbsnp_file* dbsnp = NULL;
        vector<string> filenames;

        int mode = 0;
        fragment_processor* processor = NULL;

        if (has_option(argc, argv, "h")) {
            show_help();
            return 0;
        }

        if (command == "rec" || command == "recombination") {
            mode = 0;
        } else if (command == "count") {
            mode = 1;
        } else if (command == "dist") {
            mode = 2;
        }

        for (int i = 2; i < argc; i++) {
            if (argv[i][0] != '-' && tktools::io::file_exists(argv[i]) && strcmp(argv[i] + strlen(argv[i]) - 4, ".bam") == 0) {
                cout << argv[i] << endl;
                filenames.push_back(argv[i]);
            }
        }

        int num_files = filenames.size();
        if (window_margin * 4 > window_size) {
            window_margin = window_size / 4;
        }

        if (verbose) {
            cerr << "Minimum coverage : " << coverage << endl;
            cerr << "Heterozygosity   : " << heterozygosity << endl;
            if (filename_snps != NULL) {
                cerr << "Variation file   : " << filename_snps << endl;
            }
            if (filename_genome != NULL) {
                cerr << "Genome file      : " << filename_genome << endl;
            }
            cerr << "Coverage         : " << coverage << endl;
            cerr << "Window size      : " << window_size << endl;
            cerr << "Window margin    : " << window_margin << endl;
            cerr << "Quality          : " << mapping_quality << endl;
	    if (mode == 1) {
	      cerr << "Display mode     : " << display_mode << endl;
	    }
            for (int i = 0; i < num_files; i++) {
                cerr << "Filename         : " << filenames[i] << endl;
            }
        }

        if (mode == 0) {
            processor = new recombination_detector(coverage, heterozygosity);
            processor->set_display_mode(display_mode);
        } else if (mode == 1) {
            snp_enumerator* _proc = new snp_enumerator(coverage, heterozygosity);
            _proc->set_display_mode(display_mode);
            processor = _proc;//new snp_enumerator(coverage, heterozygosity);
        } else if (mode == 2) {
        }
        if (processor == NULL) {
            throw invalid_argument("no effective mode");
        }

        if (filename_snps != NULL) {
            if (verbose) {
                cerr << "loading VCF\n";
            }
            dbsnp = dbsnp_file::load_dbsnp(filename_snps);
            processor->set_vcf(dbsnp);
        }
        processor->set_quality_threshold(mapping_quality);

        ostream* ost = &cout;
        vector<recfragment*> fragments;
        map<string,const recfragment*> pairs;
        //filenames.push_back(filename);

        //int num_files = 1;
        bamFile* bamfiles = new bamFile[num_files];
        bam_header_t** headers = new bam_header_t*[num_files];
        bam1_t** reads = new bam1_t*[num_files];
        int* status = new int[num_files];
        //bool filled = new bool[num_files];
        //bam_header_t* header;
        //bamFile bamfile;

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
        // if (filename_snps != NULL) {
        //     remove_nonestablished_bases(filename_snps, fasta_files);
        // }

        // bamfile = bam_open(filename, "rb");
        // header = bam_header_read(bamfile);
        // read = bam_init1();

        for (int i = 0; i < num_files; i++) {
            bamfiles[i] = bam_open(filenames[i].c_str(), "rb");
            headers[i] = bam_header_read(bamfiles[i]);
            reads[i] = bam_init1();
            status[i] = 0x00;
        }
        if (check_header_consistency(num_files, headers) == false) {
            throw runtime_error("incompatible bam files");
        }

/*
first
     <--------------window_size----------------->
     |------|============================|------|
      margin         processing           margin
region_start                               region_end


second
    |==================|-----|
   region_start             region_end

*/
        int retain_start = 0;
        int retain_end = window_size;
        int analysis_start = retain_start + window_margin;
        int analysis_end = retain_end - window_margin;
        int current_bamchrm = -1;
        chromosome_seq const* chromosome = NULL;
        int current_chromosome = -1;
        //string chromosome_name;

        // status
        // 0x00   : read next
        // 0x01   : read of nect section is stored
        // 0x02   : different chromosome
        // 0x04   : touched
        // 0x80   : reach end
        bool skipping = false;
        while (true) {
            int num_chr_finish = 0;
            vector<recfragment*> fragments;
            set<int> snp_positions;
            if (dbsnp != NULL && chromosome != NULL) {
                vector<dbsnp_locus const*> snps = dbsnp->get_snps(chromosome->name(), analysis_start, analysis_end);
                //cout << chromosome->name() << ":" << analysis_start << "-" << analysis_end << " // " << snps.size() << "\n";
                for (int i = 0; i < (int)snps.size(); i++) {
                    snp_positions.insert(snps[i]->position());
                }
            }
            for (int i = 0; i < num_files; i++) {
                if ((status[i] & 0x04) != 0) {
                    bam1_t const* r = reads[i];
                    // put the stored read in range of the chromosome 
                    if (r->core.tid == current_bamchrm) {
                        if (retain_start <= r->core.pos && r->core.pos < retain_end) {
                            recfragment* frag;
                            if (dbsnp != NULL) {
                                frag = new recfragment(current_chromosome, reads[i], i, snp_positions);
                            } else {
                                frag = new recfragment(current_chromosome, reads[i], i);
                            }
                            fragments.push_back(frag);
                            if (status[i] != 0xff) {
                                status[i] = 0x00;
                            }
                        }
                    } else {
                        num_chr_finish++;
                    }
                }
                if ((status[i] & 0x04) == 0) {
                    for (;;) {
                        // read until the end of the box
                        if (bam_read1(bamfiles[i], reads[i]) > 0) {
                            bam1_t const* r = reads[i];
                            //check_quality(reads[i]);
                            if (chromosome == NULL && !skipping) {
                                //cout << "determining chromosome\n";
                                int ccode = convert_chromosome_to_code(headers[i]->target_name[r->core.tid]);
                                for (int j = 0; j < (int)fasta_files.size(); j++) {
                                    if (fasta_files[j]->code() == ccode) {
                                        chromosome = fasta_files[j];
                                    }
                                }
                                //if (chromosome == NULL) { // no chromosome
                            }

                            //cout << i << ":" << r->core.tid << ":" << r->core.pos << "\t" << current_bamchrm << "/" << r->core.tid << "\t" << analysis_start << "," << analysis_end << endl;
                            if (current_bamchrm != r->core.tid) { // next chromosome
                                status[i] = 0x02 | 0x04;
                                num_chr_finish ++;
                                //cout << "change chromosome\n";
                                break;
                            } else if (r->core.pos > retain_end) { // next window
                                status[i] = 0x01 | 0x04;
                                //cout << "window filled in " << i << endl;
                                break;
                            } else {
                                if (fragments.size() <= max_fragments && chromosome != NULL) {
//                            cout << "PUT2" << endl;
                                    recfragment* frag;
                                    if (dbsnp != NULL) {
                                        frag = new recfragment(current_chromosome, reads[i], i, snp_positions);
                                    } else {
                                        frag = new recfragment(current_chromosome, reads[i], i);
                                    }
//                                    recfragment* frag = new recfragment(current_chromosome, reads[i], i);
                                    fragments.push_back(frag);
                                    status[i] |= 0x04;
                                }
                            }
                        } else {
                            status[i] = 0xff; // finished
                            num_chr_finish ++;
                        }
                    }
                }
            }

            //cout << "Fragment size = " << fragments.size() << endl;
            if ((int)fragments.size() > coverage && chromosome != NULL) {
                // integrate pairs
                recfragment::bundle_pairs(fragments);
                // process
                //process_fragments(fragments, chromosome, analysis_start, analysis_end);//, coverage, heterozygosity, *ost);
                //detect_recombination(fragments, chromosome, analysis_start, analysis_end, coverage, heterozygosity, *ost);

                //cerr << chromosome << ":" << analysis_start << "-" << analysis_end << "//" << snp_positions.size() << "       \r";
                processor->process_fragments(fragments, chromosome, analysis_start, analysis_end, *ost);
            }

            if (num_chr_finish == num_files) { // change chromosome
                int next_chr = -1;
                for (int i = 0; i < num_files; i++) {
                    bam1_t const* r = reads[i];
                    if (r->core.tid != current_bamchrm && (next_chr < 0 || r->core.tid < next_chr)) {
                        next_chr = r->core.tid;
                    }
                }
                if (next_chr < 0) {
                    break;
                }
                if (verbose) {
                    cerr << "change chromosome  " << headers[0]->target_name[next_chr] << endl;
                }
                current_bamchrm = next_chr;
                
                chromosome = NULL;
                //if (next_chr >= 0) {
                int ccode = convert_chromosome_to_code(headers[0]->target_name[next_chr]);
                for (int j = 0; j < (int)fasta_files.size(); j++) {
                    if (fasta_files[j]->code() == ccode) {
                        chromosome = fasta_files[j];
                    }
                }
                skipping = (chromosome == NULL);
//                } else {
//                }

                retain_start = 0;
                analysis_start = 0;
                retain_end = window_size;
                analysis_end = retain_end - window_margin;

                for (vector<recfragment*>::iterator it = fragments.begin(); it != fragments.end(); it++) {
                    delete *it;
                }
                fragments.erase(fragments.begin(), fragments.end());
                //chromosome = NULL;
            } else {
                analysis_start = analysis_end;
                if (fragments.size() > 0 && analysis_start < fragments[0]->position5()) {
                    analysis_start = fragments[0]->position5();
                }
                retain_start = analysis_start - window_margin;
                retain_end = retain_start + window_size;
                analysis_end = retain_end - window_margin;

                int tail = 0;
                for (int i = 0; i < (int)fragments.size(); i++) {
                    //for (vector<recfragment*>::iterator it = fragments.begin(); it != fragments.end(); it++) {
                    //delete *it;
                    const recfragment* f = fragments[i];
                    if (f->position3() > analysis_end) {
                        fragments[tail++] = fragments[i];
                    } else {
                        delete f;
                    }
                }
                fragments.erase(fragments.begin() + tail, fragments.end());
                //fragments.erase(fragments.begin(), fragments.end());
            }

        }

        if (verbose) {
            cerr << "finish processing\n";
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
        // bam_destroy1(read);
        // bam_header_destroy(header);
        // bam_close(bamfile);
        for (int i = 0; i < num_files; i++) {
            bam_destroy1(reads[i]);
            bam_header_destroy(headers[i]);
            bam_close(bamfiles[i]);
        }
        delete[] reads;
        delete[] headers;
        delete[] bamfiles;

        // dbsnp
        delete dbsnp;

        return 0;
    } catch (exception& e) {
        show_help();
        cerr << e.what() << endl;
        return -1;
    }
}
