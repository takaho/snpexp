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
#include <varstr.hxx>

using namespace tkbio;

using tktools::util::get_argument_string;
using tktools::util::get_argument_integer;
using tktools::util::get_argument_float;
using tktools::util::has_option;
using tktools::split_items;
using tktools::bio::convert_chromosome_to_code;
using tktools::bio::convert_code_to_chromosome;

// namespace {
//     bool check_header_consistency(int num, bam_header_t** headers) {
//         for (int i = 0; i < num; i++) {
//             bam_header_t* hi = headers[i];
//             for (int j = 0; j < i; j++) {
//                 bam_header_t* hj = headers[j];
//                 if (hi->n_targets != hj->n_targets) {
//                     return false;
//                 }
//                 for (int k = 0; k < hi->n_targets; k++) {
//                     if (strcmp(hi->target_name[k], hj->target_name[k]) != 0) {
//                         return false;
//                     }
//                 }
//             }
//         }
//         return true;
//     }
// }    

namespace {
  void filter_snps_by_strain(int argc, char** argv) throw (exception) {
    const char* filename_snps = get_argument_string(argc, argv, "V", NULL);
    const char* filename_count = get_argument_string(argc, argv, "i", NULL);
    const char* filename_output = get_argument_string(argc, argv, "o", NULL);
    const char* strain = get_argument_string(argc, argv, "S", "129S1");
    bool verbose = has_option(argc, argv, "verbose");

    if (verbose) {
      cerr << "SNPs      : " << filename_snps << endl;
      cerr << "Input     : " << filename_count << endl;
      cerr << "Output    : " << (filename_output == NULL ? "stdout" : filename_output) << endl;
      cerr << "Strain    : " << strain << endl;
    }

    ifstream fi(filename_count);
    dbsnp_file* dbsnp = dbsnp_file::load_dbsnp(filename_snps);
    int strain_index = dbsnp->get_strain_index(strain);
    ostream* ost = &cout;
    cerr << strain_index << endl;

    if (strain_index < 0) {
      delete dbsnp;
      throw invalid_argument("no strain");
    }

    if (filename_output != NULL) {
      ofstream* fo = new ofstream(filename_output);
      if (fo->is_open() == false) {
	delete fo;
	throw invalid_argument("cannot open output stream");
      }
      ost = fo;
    }
    size_t num_lines = 0;
    while (!fi.eof()) {
      string line;
      getline(fi, line);
      num_lines ++;
      vector<string> items = split_items(line, '\t');
      if (verbose && num_lines % 1000 == 0) {
	cerr << (num_lines / 1000) << "k " << items[0] << ":" << items[1] << "       \r";
      }
      dbsnp_locus const* snp = dbsnp->get_snp(items[0], atoi(items[1].c_str()));
      if (snp != NULL) {
	unsigned char ref = snp->get_genotype(0);
	unsigned char alt = snp->get_genotype(strain_index);
	if (ref != alt) {
	  *ost << line << endl;
	}
      }
    }
    fi.close();
    delete dbsnp;
    if (verbose) {
      cerr << "loading genomic sequences " << flush;
    }
    if (filename_output != NULL) {
      dynamic_cast<ofstream*>(ost)->close();
      delete ost;
    }
  }
}


namespace {
    void show_help() {
        cerr << "recrec [command] [options] bamfile1 bamfile2 ...\n";
        cerr << "COMMANDS\n";
        cerr << " rec/recombination   : detect recombination junctions\n";
        cerr << " count               : enumerate nucleotide frequency of SNP sites\n";
        cerr << " dist                : read <count> command output and estimate difference from strains\n";
        cerr << "OPTIONS\n";
        cerr << " -c <number>  : minimum coverage for analysis (10)`n";
        cerr << " -w <number>  : window size to analysis chunk, set more than double of fragment size (3000)\n";
        cerr << " -H <ratio>   : heterozygosity to determine SNPs gruop\n";
        cerr << " -s <filename> : genome sequence\n";
        cerr << " -o <filename> : output filename (stdout)\n";
        cerr << " -q <number>   : threshold of mapping quality (10)\n";
        cerr << " -D <number>   : display verbosity (0:all, 1:countable, 2:hetero)\n";
        cerr << " -G <filename> : annotation file in GTF format\n";
        cerr << " --pared       : paired\n";

        cerr << "OPTIONS for recombination\n";
        cerr << " -l <number>   : the number of SNPs to define haplotypes (default 2)\n";
        cerr << " -g <number>   : maximum gaps between snps for haplotype detection (default 0)\n";
        cerr << " -b <number>   : allele balance parameter (default 0.5)\n";
        cerr << "OPTIONS for dist command\n";
        cerr << " -i <filename> : the output filename of count command\n";
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
        bool paired = has_option(argc, argv, "-paired");
        int window_margin = get_argument_integer(argc, argv, "m", 0);
        const char* filename_genome = get_argument_string(argc, argv, "s", NULL);
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        const char* filename_snps = get_argument_string(argc, argv, "V", NULL);
        const char* filename_gtf = get_argument_string(argc, argv, "G", NULL);
        size_t max_fragments = get_argument_integer(argc, argv, "x", 10000);
        int mapping_quality = get_argument_integer(argc, argv, "q", 10);
	int display_mode  = get_argument_integer(argc, argv, "D", 1);
        dbsnp_file* dbsnp = NULL;
        vector<string> filenames;
        double allele_balance = get_argument_float(argc, argv, "b", 0.5);
        bool exit_on_y = has_option(argc, argv, "Y");
        int gap_tolerance = get_argument_integer(argc, argv, "g", 0);
        int snp_stretch = get_argument_integer(argc, argv, "l", 2);

        int mode = 0;
        fragment_processor* processor = NULL;
        gtffile* gtf = NULL;

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
        } else if (command == "filter") {
	  mode = 3;
        } else if (command == "str") {
             mode = 4;
        } else if (command == "cmp") {
            mode = 5;
        }
        cout << mode << endl;
        if (mode == 5) {
            try {
                denovo_snp::enumerate_hetero(argc, argv);
                return 0;
            } catch (exception& e) {
                throw;
            }
        }

        if (mode == 4) {
            try {
                str_collection::detect_str(argc, argv);
                return 0;
            } catch (exception& e) {
                throw;
            }
        }

	if (mode == 3) {
	  try {
	    filter_snps_by_strain(argc, argv);
	    return 0;
	  } catch (exception& e) {
	    cerr << "error in filtering\n";
	    throw;
	  }
	}

        if (filename_genome == NULL) {
            throw invalid_argument("no genomic sequence given");
        }

	// if (mode == 2) {
        //     try {
        //         snp_distance::main(argc, argv);
        //         return 0;
        //     } catch (exception& e) {
        //         cerr << "error while mode 2\n";
        //         throw;
        //     }
	// }

        for (int i = 2; i < argc; i++) {
            if (argv[i][0] != '-' && tktools::io::file_exists(argv[i]) && strcmp(argv[i] + strlen(argv[i]) - 4, ".bam") == 0) {
                //cout << argv[i] << endl;
                filenames.push_back(argv[i]);
            }
        }

        int num_files = filenames.size();
        if (window_margin * 6 > window_size) {
            window_margin = window_size / 6;
        } else if (window_margin < 100) {
            window_margin = 100;
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
            recombination_detector* p = new recombination_detector(coverage, heterozygosity);
            p->set_allele_balance(allele_balance);
            p->set_haplotype_parameters(snp_stretch, gap_tolerance);
            //processor = new recombination_detector(coverage, heterozygosity);
            processor = p;
            processor->set_display_mode(display_mode);
        } else if (mode == 1) {
            snp_enumerator* _proc = new snp_enumerator(coverage, heterozygosity);
            _proc->set_display_mode(display_mode);
            processor = _proc;//new snp_enumerator(coverage, heterozygosity);
        // } else if (mode == 4) {
        //     if (filename_snps != NULL) {
        //         throw invalid_argument("STR mode does not require SNPs");
        //     }
        //     str_detector* _proc = new str_detector(coverage, heterozygosity);
        //     _proc->set_display_mode(display_mode);
        //     processor = _proc;
        } else if (mode == 2) {
            strain_estimator* _proc = new strain_estimator();
            _proc->set_coverage(coverage);
            _proc->set_heterozygosity(heterozygosity);
            processor = _proc;
            dynamic_cast<strain_estimator*>(processor)->to_string();
        }

        if (processor == NULL) {
            throw invalid_argument("no effective mode");
        }
        if (filename_gtf != NULL) {
            if (verbose) cerr << "loading GTF file " << filename_gtf << "\n";
            gtf = gtffile::load_gtf(filename_gtf);
            processor->set_gtf(gtf);
        }

        if (filename_snps != NULL) {
            if (verbose) {
                cerr << "loading VCF\n";
            }
            dbsnp = dbsnp_file::load_dbsnp(filename_snps);

            // vector<dbsnp_locus const*> snps;
            // snps = dbsnp->get_snps("8",  10000000, 12000000);
            // cout << snps.size() << endl;
            // snps = dbsnp->get_snps("11",  95000000, 96000000);
            // cout << snps.size() << endl;
            // snps = dbsnp->get_snps("8",  10000000, 12000000);
            // cout << snps.size() << endl;
            // exit(0);

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
        if (verbose) {
            for (int i = 0; i < (int)fasta_files.size(); i++) {
                cerr << "chromosome #" << (i + 1) << ":" << fasta_files[i]->name() << "\t" << fasta_files[i]->length() << endl;
            }
            // const chromosome_seq* chrm = fasta_files[12];
            // cout << chrm->name() << endl;
            // for (int i = 0; i < 100; i++) {
            //     int pos = i * 1534154 + 1000001;
            //     cout << pos << "\t";
            //     for (int j = 0; j < 60; j++) {
            //         cout << chrm->get_base(j + pos);
            //     }
            //     cout << endl;
            // }
            // exit(0);
        }

//        map<int,pair<int,unsigned char*> >  genome = load_genome(filename_genome);
        // if (verbose) {
        //     cerr << fasta_files.size() << " chromosomes\n";
        // }

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
        if (tkbio::check_header_consistency(num_files, headers) == false) {
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
        // 0x01   : read of next section is stored
        // 0x02   : different chromosome
        // 0x04   : touched
        // 0x80   : reach end
        bool skipping = false;
        bool chromosome_not_included = false;

        //int counter = 0;

        while (true) {
            int num_chr_finish = 0;
            vector<recfragment*> fragments;
            set<int> snp_positions;

            // load snps if dbSNP is available
            if (dbsnp != NULL && chromosome != NULL && !chromosome_not_included) {
                vector<dbsnp_locus const*> snps;
                try {
                    snps = dbsnp->get_snps(chromosome->name(), analysis_start, analysis_end);
                    for (int i = 0; i < (int)snps.size(); i++) {
                        snp_positions.insert(snps[i]->position());
                    }
                } catch (exception& e) {
                    // invalid chromosome
                    cerr << e.what() << "  \n";
                    chromosome_not_included = true;
                }
            }

            // load reads
            for (int i = 0; i < num_files; i++) {
                if ((status[i] & 0x04) != 0) {
                    bam1_t const* r = reads[i];
                    // put the stored read in range of the chromosome 
                    if (r->core.tid == current_bamchrm) {
                        if (retain_start <= r->core.pos && r->core.pos < retain_end) {
                            recfragment* frag = NULL;
                            if (dbsnp != NULL) {
                                frag = new recfragment(current_chromosome, reads[i], i, snp_positions);
                            } else if (!chromosome_not_included) {
                                frag = new recfragment(current_chromosome, reads[i], i);
                            }
                            if (frag != NULL) {
                                fragments.push_back(frag);
                            }
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
                                // if (verbose) {
                                //     cerr << "change chromosome       \n";
                                // }
                                break;
                            } else if (r->core.pos > retain_end) { // next window
                                status[i] = 0x01 | 0x04;
                                //cout << "window filled in " << i << endl;
                                break;
                            } else {
                                if (fragments.size() <= max_fragments && chromosome != NULL) {
                                    recfragment* frag = NULL;
                                    if (dbsnp != NULL) {
                                        frag = new recfragment(current_chromosome, reads[i], i, snp_positions);
                                    } else if (!chromosome_not_included) {
                                        frag = new recfragment(current_chromosome, reads[i], i);
                                    }
//                                    recfragment* frag = new recfragment(current_chromosome, reads[i], i);
                                    if (frag != NULL) {
                                        fragments.push_back(frag);
                                    }
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


            if ((int)fragments.size() > coverage && chromosome != NULL) {
                //cerr << "analyze " << __LINE__ << endl;
                // integrate pairs
                if (paired) {
                    recfragment::bundle_pairs(fragments);
                }
                // process
                if (verbose) {
                    if (chromosome == NULL) {
                        cerr << "undetermined_chrom:" << current_bamchrm;
                    } else {
                        cerr << chromosome->name();
                    }
                    cerr << ":" << analysis_start << "-" << analysis_end << " // " << fragments.size() << " // " << snp_positions.size();
                    if (dbsnp != NULL) {
                        pair<int,int> range = dbsnp->cached_range();
                        cerr << " // " << range.first << "," << range.second;
                    }
                    cerr << "       \r";

                    // if (++counter % 1000 == 0) {
                    //     cerr << "\n" << processor->to_string() << endl;
                    // }
                }
                if (!chromosome_not_included) {
                    processor->process_fragments(fragments, chromosome, analysis_start, analysis_end, *ost);
                    *ost << std::flush;
                }
            }

            if (num_chr_finish == num_files) {// || chromosome->length() <= analysis_start) { // change chromosome
                int next_chr = -1;
                chromosome_not_included = false;
                for (int i = 0; i < num_files; i++) {
                    if (status[i] == 0xff) continue;
                    bam1_t const* r = reads[i];
                    if (r->core.tid != current_bamchrm && (next_chr < 0 || r->core.tid < next_chr)) {
                        next_chr = r->core.tid;
                    }
                }
                if (next_chr < 0) {
                    break;
                }
                if (verbose) {

                    cerr << "change chrom " << headers[0]->target_name[next_chr] << "     " << endl;
                    cerr << processor->to_string();
                }
                current_bamchrm = next_chr;
                
                chromosome = NULL;
                //if (next_chr >= 0) {
                int ccode = convert_chromosome_to_code(headers[0]->target_name[next_chr]);
                for (int j = 0; j < (int)fasta_files.size(); j++) {
                    if (fasta_files[j]->code() == ccode) {
                        chromosome = fasta_files[j];
                        current_chromosome = tktools::bio::convert_chromosome_to_code(chromosome->name().c_str());
                    }
                }
                // the chromosome in bam file is not included fasta file
                //cerr << (chromosome == NULL) << endl;
                if (chromosome != NULL) {
                    skipping = false;
                    if (chromosome->name().find("_") != string::npos) {
                        chromosome_not_included = true;
                    }
                    if (chromosome->name().find("Y") != string::npos && exit_on_y) {
                        chromosome_not_included = true;
                        break;
                    }
                } else {
                    skipping = true;
                    chromosome_not_included = true;
                }
                //cerr << (chromosome == NULL) << ", " << chromosome_not_included << endl;
                //cerr << chromosome->name() << endl;
                //skipping = (chromosome == NULL);

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
                    const recfragment* f = fragments[i];
                    if (f->position3() > analysis_start) {
                        fragments[tail++] = fragments[i];
                    } else {
                        delete f;
                    }
                }
                fragments.erase(fragments.begin() + tail, fragments.end());
            }
        }

        if (verbose) {
            cerr << "finish processing\n";
        }
        *ost << processor->to_string();

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

        // processor
        delete processor;

        return 0;
    } catch (exception& e) {
        show_help();
        cerr << e.what() << endl;
        return -1;
    }
}
