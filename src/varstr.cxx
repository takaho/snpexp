#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <set>
#include <limits>

#ifndef HAVE_BAM_H
#error You do not have bam library
#endif
#include <bam.h>

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
using std::numeric_limits;

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

bamread::bamread(bam1_t const* read) {
    _chromosome = read->core.tid;
    _start = read->core.pos;
    uint32_t const* cigar = bam1_cigar(read);
    int clen = read->core.n_cigar;
    ullong position = (ullong)(read->core.pos + 1);
    for (int i = 0; i < clen; i++) {
        int op = bam_cigar_op(cigar[i]);
        ullong slen = (ullong)bam_cigar_oplen(cigar[i]);
        ullong info = (slen << 32) | position;
        if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) { // M
            _sections.push_back(FEATURE_MATCH | info);
            position += slen;
        } else if (op == BAM_CINS) { // I
            _sections.push_back(FEATURE_REPEAT_INSERTION | info);
            position += slen;
        } else if (op == BAM_CDEL) { // D
            _sections.push_back(FEATURE_REPEAT_DELETION | info);
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) { // S or H
            break;
        } else if (op == BAM_CREF_SKIP) {
            break;
        } else if (op == BAM_CBACK) {
            break;
        }
    }
    _stop = position + 1;
}

bool bamread::shares_variation(ullong section) const {
    for (vector<ullong>::const_iterator it = _sections.begin(); it != _sections.end(); it++) {
        if (*it == section) {
            return true;
        }
    }
    return false;
}

char bamread::get_feature(ullong info) {
    char pat;
    switch (info & FEATURE_MASK) {
    case FEATURE_MATCH:
        pat = 'M'; break;
    case FEATURE_REPEAT_INSERTION:
        pat = 'I'; break;
    case FEATURE_REPEAT_DELETION:
        pat = 'D'; break;
    default:
        pat = '.'; break;
    }
    return pat;
}

string bamread::to_string() const {
    stringstream ss;
    ss << _start << "\t";
    for (vector<ullong>::const_iterator it = _sections.begin(); it != _sections.end(); it++) {
        ullong info = *it;
        ss << get_span(info) << get_feature(info);
    }
    return ss.str();
}

operator == (const str_variation& lhs, const str_variation& rhs) {
    if (lhs._position == rhs._position 
        && lhs._reference_span == rhs._reference_span 
        && lhs._read_span == rhs._read_span) {
        return true;
    } else {
        return false;
    }
}

strcollection::strcollection() {
}

strcollection::~strcollection() {
}

void strcollection::sweep(int chromosome, int start, int end) {
    if (end < 0) { // clear all
        if (chromosome < 0) {
            _reads.erase(_reads.begin(), _reads.end());
        } else {
            vector<bamread> retained;
            for (int i = 0; i < (int)_reads.size(); i++) {
                if (_reads[i].chromosome() != chromosome) {
                    retained.push_back(_reads[i]);
                }
            }
            _reads = retained;
        }
    } else {
        vector<bamread> retained;
        for (int i = 0; i < (int)_reads.size(); i++) {
            if (_reads[i].chromosome() != chromosome || _reads[i].position5() < start || _reads[i].position3() > end) {
                retained.push_back(_reads[i]);
            }
        }
        _reads = retained;
    }
}

void strcollection::add_read(bam1_t const* read) {
    _reads.push_back(bamread(read));
}

pair<int,int> strcollection::span() const {
    int pos_min = numeric_limits<int>::max();
    int pos_max = numeric_limits<int>::min();
    // if (_reads.size() == 0) {
    //     return make_pair(-1, -1);
    // }
    for (int i = 0; i < (int)_reads.size(); i++) {
        int p5 = _reads[i].position5();
        int p3 = _reads[i].position3();
        if (p5 < pos_min) pos_min = p5;
        if (p3 > pos_max) pos_max = p3;
    }
    return make_pair(pos_min, pos_max);
}

vector<str_variation> str_collection::get_variations(int coverage, double heterozygosity) const throw (exception) {
    vector<str_variation> alleles;
    for (int i = 0; i < (int)_reads.size(); i++) {
        int num = _reads[i].size();
        for (int j = 0; j < num; j++) {
            ullong section = _reads[i].get_section(j);
            ullong feature = section & bamread::FEATURE_MASK;
            if (feature != 0) {
                int span = bamread::get_span(feature);
                if ((feature & bamread::FEATURE_REPEAT_INSERTION) != 0) {
                    alleles.push_back(str_variation(position, span, 0));
                } else if ((feature & bamread::FEATURE_REPEAT_DELETION) != 0) {
                    alleles.push_back(str_variation(position, 0, span));
                }
            }
        }
    }
    for (int i = 0; i < (int)alleles.size(); i++) {
        
    }
    //for (int i = 0; i < 
}

int strcollection::detect_str(int argc, char** argv) throw (exception) {
    try {
        const char* filename1 = get_argument_string(argc, argv, "1", NULL);
        const char* filename2 = get_argument_string(argc, argv, "2", NULL);
        int coverage = get_argument_integer(argc, argv, "c", 20);
        double heterozygosity = get_argument_float(argc, argv, "z", 0.5);
        int chunk_size = get_argument_integer(argc, argv, "w", 2000);
        int margin_size = get_argument_integer(argc, argv, "m", 500);
        int num_files = 2;
        bool verbose = has_option(argc, argv, "verbose");

        if (verbose) {
            cerr << "filename 1 : " << filename1 << endl;
            cerr << "filename 2 : " << filename2 << endl;
            cerr << "heterozygosity : " << heterozygosity << endl;
            cerr << "coverage   : " << coverage << endl;
            cerr << "chunk size : " << chunk_size << endl;
            cerr << "margin     : " << margin_size << endl;
        }

        str_collection** detectors = new str_collection*[num_files];
        bamFile* bamfiles = new bamFile[num_files];
        bam_header_t** headers = new bam_header_t*[num_files];
        bam1_t** reads = new bam1_t*[num_files];

        if (verbose) {
            cerr << "opening files\n";
        }
        bamfiles[0] = bam_open(filename1, "rb");
        bamfiles[1] = bam_open(filename2, "rb");
        for (int i = 0; i < num_files; i++) {
            headers[i] = bam_header_read(bamfiles[i]);
            reads[i] = bam_init1();
            detectors[i] = new str_collection();
        }

        if (verbose) {
            cerr << "check integrity\n";
        }
        if (check_header_consistency(num_files, headers) == false) {
            throw runtime_error("incompatible bam files");
        }

        int current_chromosome = -1;
        int next_chromosome = -1;
        int position = 0;
        int next_position = position + chunk_size;
        
        if (verbose) {
            cerr << "analyzing\n";
        }
        for (;;) {
            bool chromosome_change = false;
            bool finished = false;
            for (int i = 0; i < num_files; i++) {
                for (;;) {
                    if (bam_read1(bamfiles[i], reads[i]) > 0) {
                        bam1_t const* r = reads[i];
                        detectors[i]->add_read(r);
                        if (r->core.tid != current_chromosome) {
                            chromosome_change = true;
                            next_chromosome = r->core.tid;
                            break;
                        }
                        if (r->core.pos > next_position) {
                            break;
                        }
                    } else {
                        finished = true;
                        chromosome_change = false;
                        break;
                    }
                }
            }

            if (current_chromosome >= 0) {
                cout << current_chromosome << ":" << position << "-" << next_position << "\t";
                for (int i = 0; i < num_files; i++) {
                    cout << "\t" << i << ":" << detectors[i]->size();
                }
                cout << endl;
                for (int i = 0; i < num_files; i++) {
                    vector<str_variation> gaps = detectors[i]->get_variations(coverage, heterozygosity);
                    for (int j = 0; j < (int)gaps.size(); j++) {
                        cout << gaps[j]->to_string() << endl;
                    }
                }
            }

            // for (int i = 0; i < num_files; i++) {
            //     strcollection* detector = detectors[i];
            // }
            // vector<strlocus> alleles = strcollection::get_variations(num_files, detectors, coverage, heterozygosity);
            // for (int i = 0; i < (int)alleles.size(); i++) {
            // }

            if (finished) {
                break;
            }

            if (chromosome_change) {
                int pos_min = numeric_limits<int>::max();
                for (int i = 0; i < num_files; i++) {
                    detectors[i]->sweep(current_chromosome);
                    if (!finished) {
                        if (reads[i]->core.pos < pos_min) {
                            pos_min = reads[i]->core.pos;
                        }
                        detectors[i]->add_read(reads[i]);
                    }
                }
                current_chromosome = next_chromosome;
                position = pos_min;
                next_position = position + chunk_size;
            } else {
                for (int i = 0; i < num_files; i++) {
                    detectors[i]->sweep(current_chromosome, 0, next_position - margin_size);
                }
                int pos_min = numeric_limits<int>::max();
                for (int i = 0; i < num_files; i++) {
                    pair<int,int> span = detectors[i]->span();
                    if (span.second > 0 && pos_min > span.first) {
                        pos_min = span.first;
                    }
                    detectors[i]->sweep(current_chromosome);
                }
                if (pos_min > next_position) {
                    next_position = pos_min;
                }
                position = next_position;
                next_position = position + chunk_size;
            }
        }

        for (int i = 0; i < num_files; i++) {
            bam_destroy1(reads[i]);
            bam_header_destroy(headers[i]);
            bam_close(bamfiles[i]);
        }
        delete[] reads;
        delete[] headers;
        delete[] bamfiles;
        return 0;
    } catch (exception& e) {
        throw;
        return -1;
    }
}
