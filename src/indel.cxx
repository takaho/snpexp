#include <vector>
#include <stdexcept>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

using std::string;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::exception;
using std::runtime_error;
using std::pair;
using std::make_pair;
using std::map;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::stringstream;

//using std::illegal_argument;

#include <bam.h>

#include <recrec.hxx>
#include <tktools.hxx>
#include <indel.hxx>

using namespace tkbio;
using tktools::util::get_argument_string;
using tktools::util::get_argument_integer;
using tktools::util::get_argument_float;
using tktools::util::has_option;
using tktools::split_items;
using tktools::bio::convert_chromosome_to_code;
using tktools::bio::convert_code_to_chromosome;
indel::indel() {}

indel::indel(const indel& rhs) {
    _id = rhs._id;
    _feature = rhs._feature;
    _span = rhs._span;
    _count = rhs._count;
    _coverage = rhs._coverage;
    _position = rhs._position;
}

const indel& indel::operator = (const indel& rhs) {
    if (&rhs == this) {
        return *this;
    }
    _id = rhs._id;
    _feature = rhs._feature;
    _span = rhs._span;
    _count = rhs._count;
    _coverage = rhs._coverage;
    _position = rhs._position;
    return *this;
}

indel::indel(int id, int position, int span, bool insertion) {
    _id = id;
    _position = position;
    _span = span;// > 255 ? (unsigned char)0xff : (unsigned char)_span;
    _count = 1;
    _feature = insertion ? INSERTION : 0;
    _coverage = 0;
}

double indel::heterozygosity() const {
    if (_coverage <= 0) {
        return 0.0;
    } else {
        return (double)_count / _coverage;
    }
    
}

namespace tkbio {
    bool operator == (const indel& lhs, const indel& rhs) {
        return (lhs._id == rhs._id && lhs._position == rhs._position && lhs._span == rhs._span && lhs._feature == rhs._feature);
    }
}

namespace {
    int resolve_bam_span(bam1_t const* read, int quality=0) {
        int span = 0;
        uint32_t const* cigar = bam1_cigar(read);
        int clen = read->core.n_cigar;
        //cerr << "CIGAR LENGTH " << clen << endl;
        if (clen < 1) {
            return 0;
        }
        for (int i = 0; i < clen; i++) {
            int op = bam_cigar_op(cigar[i]);
            ullong slen = (ullong)bam_cigar_oplen(cigar[i]);
            //cerr << "CIGAR SECTION SIZE " << slen << " " << __func__ << ":" << __LINE__ << endl;
            //ullong info = (slen << 32) | position;
            if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CINS) {
                span += slen;
            } else if (op == BAM_CDEL) { // D, less repeat than reference
                //
            } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) { // S or H
                break;
            } else if (op == BAM_CREF_SKIP) {
                break;
            } else if (op == BAM_CBACK) {
                break;
            }
        }
        return span;
    }
}
vector<indel> indel::parse_indel(bam1_t const* read, int quality) {
    vector <indel> indels;
    uint32_t const* cigar = bam1_cigar(read);
    int clen = read->core.n_cigar;
    if (clen == 0) {
        return indels;
    }
    ullong position = (ullong)(read->core.pos + 1);
    //bool ends_with_match = false;
    for (int i = 0; i < clen - 1; i++) {
        //cerr << i << " " << position << endl;
        //ends_with_match = false;
        int op = bam_cigar_op(cigar[i]);
        ullong slen = (ullong)bam_cigar_oplen(cigar[i]);
        //ullong info = (slen << 32) | position;
        if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) { // M
            //_sections.push_back(FEATURE_MATCH | info);
            position += slen;
            //ends_with_match = true;
        } else if (op == BAM_CINS) { // I, more repeats than reference
            if (i > 0) {
                indels.push_back(indel(0, position, slen, true));
            }
            // if (i > 0) {
            //     _sections.push_back(FEATURE_REPEAT_INSERTION | info);
            // }
            position += slen;
        } else if (op == BAM_CDEL) { // D, less repeat than reference
            if (i > 0) {
                indels.push_back(indel(0, position, slen, false));
            }
            // if (i > 0) {
            //     _sections.push_back(FEATURE_REPEAT_DELETION | info);
            // }
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) { // S or H
            break;
        } else if (op == BAM_CREF_SKIP) {
            break;
        } else if (op == BAM_CBACK) {
            break;
        }
    }
    // _stop = position + 1;
    // if (!ends_with_match) {
    //     _sections.erase(_sections.begin() + _sections.size() - 1);
    // }
    return indels;
}

string indel::to_string() const {
    stringstream ss;
    ss << id() << "\t" << position() << "\t" << span() << ":" << (is_insertion() ? "I" : "D" )<< "\t" << count() << "\t" << coverage();
    return ss.str();
}

namespace {
    // bool compare_first(const pair<int,indel>& lhs, const pair<int,indel>& rhs) {
    //     return (lhs.first < rhs.first);
    // }
}

int indel::detect_indel_polymorphism(int argc, char** argv) throw (exception) {
    try {
        const char* filename1 = get_argument_string(argc, argv, "1", "/mnt/smb/tae/stap/shira/BAM6/Sample6.bam");
        const char* filename2 = get_argument_string(argc, argv, "2", "/mnt/smb/tae/stap/shira/BAM12/Sample12.bam");
        //const char* filename_bed = get_argument_string(argc, argv, "b", NULL);
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        int coverage = get_argument_integer(argc, argv, "c", 20);
        double heterozygosity = get_argument_float(argc, argv, "z", 0.5);
        int chunk_size = get_argument_integer(argc, argv, "w", 3000);
        int margin_size = get_argument_integer(argc, argv, "m", 500);
        int size_threshold = get_argument_integer(argc, argv, "t", 6);
        int num_files = 2;
        int maximum_reads_in_window = get_argument_integer(argc, argv, "x", 0);
        bool verbose = has_option(argc, argv, "verbose");
        //bool use_preset = filename_bed != NULL;

        if (verbose) {
            cerr << "filename 1 : " << filename1 << endl;
            cerr << "filename 2 : " << filename2 << endl;
            cerr << "heterozygosity : " << heterozygosity << endl;
            cerr << "coverage   : " << coverage << endl;
            cerr << "chunk size : " << chunk_size << endl;
            cerr << "margin     : " << margin_size << endl;
            cerr << "max reads  : " << maximum_reads_in_window << endl;
            cerr << "output     : " << (filename_output == NULL ? "stdout" : filename_output) << endl;
            // if (use_preset) {
            //     cerr << "preset regions : " << filename_bed << endl;
            // }
        }

        ostream* ost = &cout;
        //str_collection** detectors = new str_collection*[num_files];
        bamFile* bamfiles = new bamFile[num_files];
        bam_header_t** headers = new bam_header_t*[num_files];
        bam1_t** reads = new bam1_t*[num_files];

        if (filename_output != NULL) {
            ofstream* fo = new ofstream(filename_output);
            if (fo->is_open() == false) {
                throw invalid_argument("cannot open output file");
            }
            ost = fo;
        }

        bamfiles[0] = bam_open(filename1, "rb");
        bamfiles[1] = bam_open(filename2, "rb");
        for (int i = 0; i < num_files; i++) {
            headers[i] = bam_header_read(bamfiles[i]);
            reads[i] = bam_init1();
            //detectors[i] = new str_collection();
        }

        if (check_header_consistency(num_files, headers) == false) {
            throw runtime_error("incompatible bam files");
        }
        int current_chromosome = -1;
        int next_chromosome = -1;
        int position = 0;
        int next_position = position + chunk_size;
        int steps = 0;
        //vector<map<int, indel> > indels;
        vector<indel> indels;
        int** covers = new int*[num_files];
        //vector<int*> covers;
        //vector<pair<int,int> > covers;
        for (int i = 0; i < num_files; i++) {
            // map<int,indel> dummy;
            // indels.push_back(dummy);
            covers[i] = new int[chunk_size];
            for (int j = 0; j < chunk_size; j++) covers[i][j] = 0;
        }
        
        for (;;) {
            bool chromosome_change = false;
            bool finished = false;
            bool active = true;
               
            for (int i = 0; i < num_files; i++) {
                for (;;) {
                    if (bam_read1(bamfiles[i], reads[i]) > 0) {
                        bam1_t const* r = reads[i];
                        //if (maximum_reads_in_window <= 0) {// || detectors[i]->size() < maximum_reads_in_window) {
                        if (active) {
                            int position5 = r->core.pos + 1;
                            int position3 = position5 + resolve_bam_span(r);
                            if (position5 < position) position5 = position;
                            if (position3 >= next_position) position3 = next_position;
                            for (int j = position5; j < position3; j++) {
                                covers[i][j - position] ++;
                            }
                            //covers[i].push_back(make_pair(position5, position3));
                            //cerr << position5 << "-" << position3 << " ==> ";
                            
                            //detectors[i]->add_read(r);
                            vector<indel> di = indel::parse_indel(r);
                            for (int j = 0; j < (int)di.size(); j++) {
                                di[j].set_id(i);
                                //int pos = di[j].position();
                                bool included = false;
                                for (int k = 0; k < (int)indels.size(); k++) {
                                    if (indels[k] == di[j]) {
                                        indels[k].set_count(indels[k].count() + 1);
                                        included = true;
                                        break;
                                    }
                                }
                                if (!included) {
                                    indels.push_back(di[j]);
                                }
                            }
                        }
                        
                        int chrm = r->core.tid;
                        int pos = r->core.pos;
                        
                        if (chrm != current_chromosome) {
                            chromosome_change = true;
                            if (next_chromosome <= 0) {
                                next_chromosome = chrm;
                            }
                            break;
                        }
                        if (pos > next_position) {
                            break;
                        }
                    } else {
                        finished = true;
                        chromosome_change = false;
                        break;
                    }
                }
            }

            // detect gaps and insertions
            
            if (current_chromosome >= 0) {
                if (verbose) {// && !use_preset) {
                    if (++steps % 1000 == 0) {
                        cerr << " " << headers[0]->target_name[current_chromosome] << ":" << position << "-" << next_position << " ";
                        cerr << "          \r";
                    }
                }
                //cerr << indels.size() << endl;
                vector<indel> accepted;
                //vector<pair<int,indel> > accepted;
                // for (int i = 0; i < (int)covers.size(); i++) {
                    
                // }
                for (int i = 0; i < (int)indels.size(); i++) {
                    indel& i1 = indels[i];
                    if (i1.is_rejected() || i1.span() < size_threshold || i1.count() <= 1) {
                        i1.reject();
                        continue;
                    }
                    int const* cv = covers[i1.id()];
                    if (position <= i1.position() && i1.position() < next_position) {
                        i1.set_coverage(cv[i1.position() - position]);
                    } else {
                        continue;
                    }
                    
                    bool rejected = false;
                    for (int j = 0; j < (int)indels.size(); j++) {
                        indel& i2 = indels[j];
                        if (i1.id() != i2.id() && i1.position() == i2.position() && i1.span() == i2.span()) {
                            i1.reject();
                            i2.reject();
                            rejected = true;
                            break;
                        }
                    }
                    //cout << i1.to_string() << " " << i1.heterozygosity() << endl;
                    if (!i1.is_rejected() && i1.coverage() >= coverage) {
                        double h = i1.heterozygosity();
                        if (heterozygosity <= h && h <= 1.0 - heterozygosity) {
                            *ost << headers[0]->target_name[current_chromosome] << "\t" << i1.to_string() << endl;
                        }
                    }
                }
                
            }
            
            if (finished) {
                break;
            }
            
            if (chromosome_change) {
                //int pos_min = numeric_limits<int>::max();
                if (verbose) {
                    //int total_reads = 0;
                    // for (int i = 0; i < num_files; i++) {
                    //     total_reads += detectors[i]->size();
                    // }
                    if (next_chromosome < 0 || next_chromosome >= headers[0]->n_targets) {
                        finished = true;
                    } else {
                        // if (use_preset == false) {
                        //     cerr << "change chromosome to " << headers[0]->target_name[next_chromosome] << endl;//<< ", sweep " << total_reads << "reads from " << position << "\r";
                        // }
                    }
                }
                if (!finished) {
                    // for (int i = 0; i < num_files; i++) {
                    //     indels[i].erase(indels[i].begin(), indels[i].end());
                    // }
                    current_chromosome = next_chromosome;
                    next_chromosome = -1;
                    position = 0;//pos_min;
                    next_position = position + chunk_size;
                } else {
                    break;
                }
            } else { // next position
                vector<indel> remnants;
                for (int i = 0; i < (int)indels.size(); i++) {
                    if (indels[i].position() >= next_position - margin_size) {
                        remnants.push_back(indels[i]);
                    }
                }
                indels = remnants;
                    
                //     //for (int i = 0; i < (int)indels.size(); i++) {
                //     map<int,indel> obj;
                //     for (map<int,indel>::iterator it = indels[i].begin(); it != indels[i].end(); it++) {
                //         if (it->first >= next_position - margin_size) {
                //             obj.insert(*it);
                //         }
                //     }
                //     indels[i] = obj;
                // }
                for (int i = 0; i < num_files; i++) {
                    //cerr << i << ":" << chunk_size << ", " << margin_size << endl;
                    memmove(covers[i], covers[i] + chunk_size - margin_size, sizeof(int) * margin_size);
                    for (int j = margin_size; j < chunk_size; j++) {
                        covers[i][j] = 0;
                    }
                }
                position = next_position;
                next_position = position + chunk_size;
            }
        }

        if (filename_output != NULL) {
            dynamic_cast<ofstream*>(ost)->close();
            delete ost;
            ost = &cout;
        }

        for (int i = 0; i < num_files; i++) {
            bam_destroy1(reads[i]);
            bam_header_destroy(headers[i]);
            bam_close(bamfiles[i]);
            delete[] covers[i];
            //delete detectors[i];
        }
        //delete[] detectors;
        delete[] covers;
        delete[] reads;
        delete[] headers;
        delete[] bamfiles;
        return 0;
    } catch (exception& e) {
        throw;
        return -1;
    }
}
