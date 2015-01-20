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
//using std::sort;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::stringstream;
using std::endl;
using std::set;
using std::numeric_limits;
using std::vector;
using std::make_pair;

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
//     string resolve_cigar(const bam1_t* read) {
//         const uint32_t* cigar = bam1_cigar(read);
//         int len = read->core.n_cigar;
//         stringstream ss;
//         for (int i = 0; i < len; i++) {
//             int op = bam_cigar_op(cigar[i]);
//             int slen = bam_cigar_oplen(cigar[i]);
//             ss << slen;
//             if (op == BAM_CMATCH) {
//                 ss << 'M';
//             } else if (op == BAM_CINS) {
//                 ss << 'I';
//             } else if (op == BAM_CDEL) {
//                 ss << 'D';
//             } else if (op == BAM_CREF_SKIP) {
//                 ss << 'N';
//             } else if (op == BAM_CSOFT_CLIP) {
//                 ss << 'S';
//             } else if (op == BAM_CHARD_CLIP) {
//                 ss << 'H';
//             } else if (op == BAM_CPAD) {
//                 ss << 'P';
//             } else if (op == BAM_CEQUAL) {
//                 ss << 'M';
//             } else if (op == BAM_CDIFF) {
//                 ss << 'X';
//             } else if (op == BAM_CBACK) {
//                 break;
//             }
//         }
//         return ss.str();
//     }
// }

bamread::bamread(bam1_t const* read) {
    _chromosome = read->core.tid;
    _start = read->core.pos;
    uint32_t const* cigar = bam1_cigar(read);
    int clen = read->core.n_cigar;
    if (clen == 0) {
        return;
    }
    ullong position = (ullong)(read->core.pos + 1);
    bool ends_with_match = false;
    for (int i = 0; i < clen; i++) {
        ends_with_match = false;
        int op = bam_cigar_op(cigar[i]);
        ullong slen = (ullong)bam_cigar_oplen(cigar[i]);
        ullong info = (slen << 32) | position;
        if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) { // M
            _sections.push_back(FEATURE_MATCH | info);
            position += slen;
            ends_with_match = true;
        } else if (op == BAM_CINS) { // I, more repeats than reference
            if (i > 0) {
                _sections.push_back(FEATURE_REPEAT_INSERTION | info);
            }
            position += slen;
        } else if (op == BAM_CDEL) { // D, less repeat than reference
            if (i > 0) {
                _sections.push_back(FEATURE_REPEAT_DELETION | info);
            }
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) { // S or H
            break;
        } else if (op == BAM_CREF_SKIP) {
            break;
        } else if (op == BAM_CBACK) {
            break;
        }
    }
    _stop = position + 1;
    if (!ends_with_match) {
        _sections.erase(_sections.begin() + _sections.size() - 1);
    }
}

bool bamread::covers(int start, int end) const {
    return _start <= start && end <= _stop;
}

bool bamread::shares_variation(ullong section) const {
    for (vector<ullong>::const_iterator it = _sections.begin(); it != _sections.end(); it++) {
        if (*it == section) {
            return true;
        }
    }
    return false;
}

bool bamread::has_reference_at(int start, int end) const {
    for (int i = 0; i < (int)_sections.size(); i++) {
        ullong section = _sections[i];
        int position = get_position(section);
        if (position >= end) {
            return false;
        }
        int span = get_span(section);
        if (position <= start && end <= position + span) {
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

str_variation::str_variation(ullong code) {
    _code = code;
    _num_shares = _num_reference = _num_others = 0;
}

str_variation::str_variation(int position, int reference_span, int read_span) {
    ullong feature = 0;
    if (reference_span > 0) {
        feature = bamread::FEATURE_REPEAT_DELETION | ((ullong)reference_span << 32);
    } else if (read_span > 0) {
        feature = bamread::FEATURE_REPEAT_INSERTION | ((ullong)reference_span << 32);
    } else {
        feature = bamread::FEATURE_ERROR;
    }
    _code = feature | (ullong)position;
    _num_shares = 0;
    _num_reference = 0;
    _num_others = 0;
}

str_variation::str_variation() {
    _code = bamread::FEATURE_ERROR;
    _num_shares = 0;
    _num_reference = 0;
    _num_others = 0;
}

char str_variation::feature() const {
    switch (_code & bamread::FEATURE_MASK) {
    case bamread::FEATURE_MATCH:
        return 'M';
    case bamread::FEATURE_REPEAT_INSERTION:
        return 'I';
    case bamread::FEATURE_REPEAT_DELETION:
        return 'D';
    default:
        return '?';
    }
}

int str_variation::position() const {
    return (int)(_code & 0xffffffff);
}

int str_variation::read_span() const {
    if ((_code & bamread::FEATURE_MASK) == bamread::FEATURE_REPEAT_INSERTION) {
        return (int)(_code >> 32) & 0xffff;
    } else {
        return 0;
    }
}

int str_variation::reference_span() const {
    if ((_code & bamread::FEATURE_MASK) == bamread::FEATURE_REPEAT_DELETION) {
        return (int)(_code >> 32) & 0xffff;
    } else {
        return 0;
    }
}

void str_variation::set_counts(int share, int reference, int others) {
    _num_shares = share;
    _num_reference = reference;
    _num_others = others;
}

string str_variation::to_string() const {
    stringstream ss;
    ss << position() << ":" << feature() << ":" << (reference_span() > 0 ? reference_span() : read_span());
    ss << "\t" << _num_shares << "/" << _num_reference << "/" << _num_others;
    return ss.str();
}

namespace tkbio {
    bool operator == (const str_variation& lhs, const str_variation& rhs) {
        if (lhs._code != rhs._code) {
            cerr << lhs.to_string() << "\t" << hex << lhs._code << " ? " << rhs.to_string() << "\t" << rhs._code << dec << endl;
        }
        return lhs._code == rhs._code;
    }
}

str_collection::str_collection() {
}

str_collection::~str_collection() {
}

int str_collection::count_coverage(int start, int stop) const {
    int num = 0;
    for (vector<bamread>::const_iterator it = _reads.begin(); it != _reads.end(); it++) {
        if (it->covers(start, stop)) num++;
    }
    return num;
}

void str_collection::sweep(int chromosome, int start, int end) {
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

void str_collection::add_read(bam1_t const* read) {
    _reads.push_back(bamread(read));
}

pair<int,int> str_collection::span() const {
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
                int span = bamread::get_span(section);
                if (span < 4) continue;
                alleles.push_back(str_variation(section));
            }
        }
    }

    vector<str_variation> unique;
    for (int i = 0; i < (int)alleles.size(); i++) {
        bool overlap = false;
        for (int j = i + 1; j < (int)alleles.size(); j++) {
            if (alleles[i] == alleles[j]) {
                overlap = true;
                break;
            }
        }
        if (!overlap) {
            //cerr << alleles[i].to_string() << endl;
            unique.push_back(alleles[i]);
        }
    }

    vector<str_variation> counted;
    float lower = (float)heterozygosity;
    float upper;
    if (heterozygosity <= 0.0) {
        lower = 0.0f;
        upper = 1.0f;
    }  else {
        if (lower > 1.0f) {
            upper = lower;
            lower = 1.0f / lower;
        } else {
            upper = 1.0f / lower;
        }
    }
    for (int i = 0; i < (int)unique.size(); i++) {
        int num_identical = 0;
        int num_reference = 0;
        int num_others = 0;
        str_variation& allele = unique[i];
        int start = allele.position();
        int stop = allele.position() + allele.reference_span();
        ullong info = allele.code();
        for (int j = 0; j < (int)_reads.size(); j++) {
            const bamread& br = _reads[j];
            if (!br.covers(start, stop)) continue;
            if (br.shares_variation(info)) {
                num_identical++; // shares the same allele
            } else if (br.has_reference_at(start, stop)) {
                num_reference++; // not shares the allele
            } else {
                num_others++; // has different alleles
            }
        }
        int total = num_identical + num_reference + num_others;
        int minimum = (int)(total * lower + .5f);
        int maximum = (int)(total * upper + .5f);
        cout << allele.to_string() << "\t" << num_identical << "," << num_reference << "," << num_others << endl;
        if (total >= coverage && minimum <= num_identical && num_identical <= maximum) {
            allele.set_counts(num_identical, num_reference, num_others);
            counted.push_back(allele);
            cout << allele.to_string() << endl;
        }
    }
    return counted;
}

namespace {

    // bool compare_first_element(const pair<int,int>& lhs, const pair<int,int>& rhs) {
    //     return lhs.first < rhs.first;
    // }

    pair<int,int> find_repeat_regions(const vector<repeat_region*>& repeats, int chrom, int start, int stop) {
        int left = 0;
        int right = (int)repeats.size();
        int center = 0;
        while (left < right) {
            center = (left + right) / 2;
            const repeat_region* r = repeats[center];
            if (r->chromosome() < chrom) {
                left = center + 1;
            } else if (r->chromosome() > chrom) {
                right = center;
            } else if (r->start() > stop) {
                right = center;
            } else if (r->stop() < start) {
                left = center + 1;
            } else {
                //int index = center;
                int limit_lower = 0;
                int limit_upper = (int)repeats.size();
                for (int i = center; i >= 0; i--) {
                    r = repeats[i];
                    if (r->chromosome() != chrom || r->stop() < start) {
                        limit_lower = i + 1;
                        break;
                    }
                }
                for (int i = center + 1; i < (int)repeats.size(); i++) {
                    r = repeats[i];
                    if (r->chromosome() != chrom || r->start() > stop) {
                        limit_upper = i;
                        break;
                    }
                }
                return make_pair(limit_lower, limit_upper);
            }
        }
        return make_pair(0,0);
    }
}

int str_collection::detect_str(int argc, char** argv) throw (exception) {
    try {
        if (has_option(argc, argv, "-enum")) {
            repeat_region::enumerate_repeat_regions(argc, argv);
            return 0;
        }
        const char* filename1 = get_argument_string(argc, argv, "1", "/mnt/smb/tae/stap/shira/BAM6/Sample6.bam");
        const char* filename2 = get_argument_string(argc, argv, "2", "/mnt/smb/tae/stap/shira/BAM12/Sample12.bam");
        const char* filename_bed = get_argument_string(argc, argv, "b", NULL);
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        int coverage = get_argument_integer(argc, argv, "c", 20);
        double heterozygosity = get_argument_float(argc, argv, "z", 0.5);
        int chunk_size = get_argument_integer(argc, argv, "w", 2000);
        int margin_size = get_argument_integer(argc, argv, "m", 500);
        int num_files = 2;
        int maximum_reads_in_window = get_argument_integer(argc, argv, "x", 0);
        bool verbose = has_option(argc, argv, "verbose");
        bool use_preset = filename_bed != NULL;

        if (verbose) {
            cerr << "filename 1 : " << filename1 << endl;
            cerr << "filename 2 : " << filename2 << endl;
            cerr << "heterozygosity : " << heterozygosity << endl;
            cerr << "coverage   : " << coverage << endl;
            cerr << "chunk size : " << chunk_size << endl;
            cerr << "margin     : " << margin_size << endl;
            cerr << "max reads  : " << maximum_reads_in_window << endl;
            cerr << "output     : " << (filename_output == NULL ? "stdout" : filename_output) << endl;
            if (use_preset) {
                cerr << "preset regions : " << filename_bed << endl;
            }
        }

        ostream* ost = &cout;
        str_collection** detectors = new str_collection*[num_files];
        bamFile* bamfiles = new bamFile[num_files];
        bam_header_t** headers = new bam_header_t*[num_files];
        bam1_t** reads = new bam1_t*[num_files];
        //map<string,vector<repeat_region> > test_regions;//pair<int,int> > > test_regions;
        vector<repeat_region*> preset_regions;

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
            detectors[i] = new str_collection();
        }

        // if (verbose) {
        //     cerr << "check integrity\n";
        // }
        if (check_header_consistency(num_files, headers) == false) {
            throw runtime_error("incompatible bam files");
        }

        if (use_preset) {
            if (verbose) cerr << "loading preset regions\n";
            ifstream fi(filename_bed);
            map<string,int> c2code;
            while (!fi.eof()) {
                string line;
                getline(fi, line);
                vector<string> items = split_items(line, '\t');
                if (items.size() >= 3) {
                    const string& chrom = items[0];
                    map<string,int>::const_iterator ic = c2code.find(chrom);
                    int chromcode = -1;
                    if (ic == c2code.end()) {
                        for (int i = 0; i < headers[0]->n_targets; i++) {
                            if (chrom == headers[0]->target_name[i]) {
                                chromcode = i;
                                break;
                            }
                        }
                        c2code[chrom] = chromcode;
                    } else {
                        chromcode = ic->second;
                    }
                    if (chromcode < 0) continue;
                    int start = atoi(items[1].c_str());
                    int stop = atoi(items[2].c_str());
                    //int unitsize = 0;
                    string unitmotif;
                    if (items.size() >= 5) {
                        unitmotif = items[4];
                        //unitsize = std::atoi(items[3].c_str());
                    }
		    //cerr << line ;
                    preset_regions.push_back(new repeat_region(unitmotif.size(), unitmotif.c_str(), chromcode, start, stop));
		    //cerr << "\t***\n";
                }
            }
            fi.close();
            //cerr << "sort\n";
            std::sort(preset_regions.begin(), preset_regions.end(), repeat_region::compare_position);
            //cerr << "sorted\n";
            // for (map<string,vector<repeat_region> >::const_iterator it = test_regions.begin(); it != test_regions.end(); it++) {
            //     //for (map<string,vector<pair<int,int> > >::const_iterator it = test_regions.begin(); it != test_regions.end(); it++) {
            //     if (verbose) {
            //         cerr << it->first << "\t" << it->second.size() << endl;
            //     }
            //     std::sort(it->second.begin(), it->second.end(), repeat_region::compare_position);
            // }
            //exit(0);
        }
        
        int current_chromosome = -1;
        int next_chromosome = -1;
        int position = 0;
        int next_position = position + chunk_size;
        int steps = 0;
        //vector<repeat_region>* counting_regions = NULL;
        vector<repeat_region> target_regions;
        int preset_regions_start = 0;
        int preset_regions_stop = 0;
//        bool debug = false;
        for (;;) {
            bool chromosome_change = false;
            bool finished = false;
            bool active = true;
            if (use_preset) {
                active = false;
                pair<int,int> rindices = find_repeat_regions(preset_regions, current_chromosome, position, next_position);
                if (rindices.first < rindices.second) {
                    active = true;
                    preset_regions_start = rindices.first;
                    preset_regions_stop = rindices.second;
                } else {
                    active = false;
                }
                //cerr << position << "-" << next_position << ":" << active << "         " << endl;
            }
               
            for (int i = 0; i < num_files; i++) {
                for (;;) {
                    if (bam_read1(bamfiles[i], reads[i]) > 0) {
                        bam1_t const* r = reads[i];
                        if (maximum_reads_in_window <= 0 || detectors[i]->size() < maximum_reads_in_window) {
                            if (active) {
                                detectors[i]->add_read(r);
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
                if (verbose && !use_preset) {
                    if (++steps % 1000 == 0) {
                        cerr << " " << headers[0]->target_name[current_chromosome] << ":" << position << "-" << next_position << " ";
                        //cerr << current_chromosome << ":" << position << "-" << next_position << " ";
                        for (int i = 0; i < num_files; i++) {
                            cerr << " " << i << ":" << detectors[i]->size();
                        }
                        cerr << "     \r";
                    }
                }

                vector<vector<str_variation> > detected;
                for (int i = 0; i < num_files; i++) {
                    vector<str_variation> gaps = detectors[i]->get_variations(coverage / 2, heterozygosity);
                    detected.push_back(gaps);
                    // for (int j = 0; j < (int)gaps.size(); j++) {
                    //     cout << gaps[j].to_string() << endl;
                    // }
                }
                if (use_preset) {
                    //cerr << "check preset regions : \n";
                    for (int index = preset_regions_start; index < preset_regions_stop; index++) {
                        repeat_region* rr = preset_regions[index];
                        int cov_min = detectors[0]->count_coverage(rr->start(), rr->stop());
                        for (int i = 1; i < num_files; i++) {
                            int cov_ = detectors[i]->count_coverage(rr->start(), rr->stop());
                            if (cov_ < cov_min) {
                                cov_min = cov_;
                            }
                        }
                        if (cov_min >= coverage) {
                            rr->set_status((repeat_region::STATUS)(rr->status() | repeat_region::COVERED));
                        } else {
			  continue;
			}
                            
			//bool flag_covered = false;
			//bool flag_polymorphic = false;
                        //bool flag_accepted = false;
                        stringstream vstr;
                        //*ost << headers[0]->target_name[rr->chromosome()] << "\t" << rr->start() << "\t" << rr->stop();
                        for (int i = 0; i < (int)detected.size(); i++) {
                            //cerr << "file " << i << endl;
                            //bool has_variation = false;
                            const str_variation* sv_inside = NULL;
                            for (int j = 0; j < (int)detected[i].size(); j++) {
                                //cerr << "item " << j << endl;
                                const str_variation& sv = detected[i][j];

                                if (sv.position() >= rr->start() && sv.position() + sv.reference_span() <= rr->stop()) {
                                    int c = sv.coverage();
				    if (c >= coverage) {
				      //flag_covered = true;
				      int o = sv.occurrence();
                                    //if (c >= coverage) {
                                    //int o = sv.occurrence();
				      if (o >= c / 4 && o <= c * 4) {
                                        sv_inside = &sv;
                                        rr->set_status((repeat_region::STATUS)(repeat_region::POLYMORPHIC | rr->status()));
                                        break;                                    
				      }
				    }
                                }
                            }
                            if (sv_inside != NULL) {
			      //flag_polymorphic = true;
			      //flag_accepted = true;
			      vstr<< "\t" << sv_inside->feature() << ":" << (sv_inside->reference_span() > 0 ? sv_inside->reference_span() : sv_inside->read_span());
                            } else {
			      vstr << "\t.";
                            }
                        }
                        //*ost << "\n";
                        if (verbose && (rr->covered() || rr->polymorphic())) {//flag_accepted) {//flag_accepted) {
                            cerr << headers[0]->target_name[rr->chromosome()] << "\t" << rr->start() << "\t" << rr->stop();
                            cerr << "\t" << cov_min << "\t" << rr->covered() << "\t" << rr->polymorphic();
                            cerr << vstr.str() << flush << "\n";
                            //rr->set_accepted(true);
                        }
                    }
                } else { // without template
                    for (int i = 0; i < num_files; i++) {
                        vector<str_variation>& vi = detected[i];
                        cerr << i << ":" << vi.size() << endl;
                        for (int k = 0; k < (int)vi.size(); k++) {
                            str_variation& ak = vi[k];
                            int c = ak.coverage();
                            int o = ak.occurrence();
                            if (c < coverage || o < c / 4 || o > c * 4) {
                                continue;
                            }
                            // if (debug) {
                            cerr << "test " << i << "," << k << ":" << ak.to_string() << "\t" << endl;
                            // }
                            bool flag_shared = false;
                            cerr << num_files << endl;
                            for (int j = 0; j < num_files; j++) {
                                cerr << "i=" << i << ", j=" << j << endl;
                                if (i == j) continue;
                                vector<str_variation>& vj = detected[j];
                                for (int l = 0; l < (int)vj.size(); l++) {
                                    if (ak == vj[l]) {
                                        //if (debug) {
                                        cerr << " rejected by " << j << "," << l << ":" << vj[l].to_string() << endl;
                                        //}
                                        flag_shared = true;
                                        ak.set_counts(0,0,0);
                                        vj[l].set_counts(0,0,0);
                                        break;
                                    }
                                }
                                if (flag_shared) {
                                    cerr << ak.to_string() << " was accepted   \n";
                                    break;
                                }
                            }
                            if (!flag_shared) {
                                // test coverage in the other sequences
                                for (int j = 0; j < num_files; j++) {
                                    if (i == j) continue;
                                    if (detectors[j]->count_coverage(ak.position(), ak.position() + ak.reference_span()) >= coverage) {
                                        *ost << headers[i]->target_name[current_chromosome] 
                                             << "\t" << ak.position() 
                                             << "\t" << (ak.position() + ak.reference_span())
                                             << "\t" << ak.read_span() 
                                             << "\t" << ak.reference_span()
                                             << "\t" << ak.occurrence() << "/" << ak.opposite() 
                                             << flush << endl;
                                    }
                                }
                            }
                        }
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
                        total_reads += detectors[i]->size();
                    }
                    if (next_chromosome < 0 || next_chromosome >= headers[0]->n_targets) {
                        finished = true;
                    } else {
                        if (use_preset == false) {
                            cerr << "change chromosome to " << headers[0]->target_name[next_chromosome] << ", sweep " << total_reads << "reads from " << position << "\r";
                        }
                    }
                }
                if (!finished) {
                    for (int i = 0; i < num_files; i++) {
                        detectors[i]->sweep(current_chromosome);
                        pair<int,int> region = detectors[i]->span();
                        pos_min = pos_min < region.first ? pos_min : region.first;
                    }
                    //     if (reads[i]->core.pos < pos_min) {
                    //         pos_min = reads[i]->core.pos;
                    //     }
                    //     detectors[i]->add_read(reads[i]);
                    // }
                    if (pos_min == numeric_limits<int>::max()) {
                        pos_min = 0;
                    }
		    for (int i = 0; i < (int)preset_regions.size(); i++) {
		      repeat_region* rr = preset_regions[i];
		      if (rr->chromosome() == current_chromosome) {
			*ost << headers[0]->target_name[rr->chromosome()] << "\t" << rr->start() << "\t" << rr->stop() << "\t" << rr->covered() << "\t" << rr->polymorphic() << endl;
		      }
		    }
                    current_chromosome = next_chromosome;
                    next_chromosome = -1;
                    position = pos_min;
                    next_position = position + chunk_size;
                } else {
                    break;
                }
            } else { // next position
                for (int i = 0; i < num_files; i++) {
                    detectors[i]->sweep(current_chromosome, 0, next_position);// - margin_size);
                }
                int pos_min = numeric_limits<int>::max();
                for (int i = 0; i < num_files; i++) {
                    pair<int,int> span = detectors[i]->span();
                    if (span.second > 0 && pos_min > span.first) {
                        pos_min = span.first;
                    }
                    detectors[i]->sweep(current_chromosome);
                }
                if (pos_min == numeric_limits<int>::max()) {
                    pos_min = next_position;
                } else if (pos_min > next_position) {
                    next_position = pos_min;
                }
                position = next_position;
                next_position = position + chunk_size;
            }
        }
        if (use_preset) {
            for (int i = 0; i < (int)preset_regions.size(); i++) {
                repeat_region* rr = preset_regions[i];
                //*ost << headers[0]->target_name[rr->chromosome()] << "\t" << rr->start() << "\t" << rr->stop() << "\t" << rr->covered() << "\t" << rr->polymorphic() << endl;
                delete rr;
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
            delete detectors[i];
        }
        delete[] detectors;
        delete[] reads;
        delete[] headers;
        delete[] bamfiles;
        return 0;
    } catch (exception& e) {
        throw;
        return -1;
    }

}

repeat_region::repeat_region(const repeat_region& rhs) {
    _pattern = rhs._pattern;
    _chromosome = rhs._chromosome;
    _start = rhs._start;
    _stop = rhs._stop;
    _status = DEFAULT;
    //_accepted = rhs._accepted;
    //_coverage = rhs._coverage;
}

const repeat_region& repeat_region::operator = (const repeat_region& rhs) {
    if (this == &rhs) return *this;
    _pattern = rhs._pattern;
    _chromosome = rhs._chromosome;
    _start = rhs._start;
    _stop = rhs._stop;
    //_accepted = rhs._accepted;
    _status = rhs._status;
    //_coverage = rhs._coverage;
    return *this;
}

repeat_region::repeat_region(int size, const char* sequence, int start, int stop) {
    _pattern = string(sequence, size);
    _chromosome = -1;
    _start = start;
    _stop = stop;
    _status = DEFAULT;
    //_accepted = false;
    //_coverage = 0;
}

repeat_region::repeat_region(int size, const char* sequence, int chromosome, int start, int stop) {
    _pattern = string(sequence, size);
    _chromosome = chromosome;
    _start = start;
    _stop = stop;
    _status = DEFAULT;
    //_accepted = false;
    //_coverage = 0;
}

repeat_region::repeat_region(int start, int stop) {
    _chromosome = -1;
    _start = start;
    _stop = stop;
    _status = DEFAULT;
    //_accepted = false;
    //_coverage = 0;
}

string repeat_region::to_string() const {
    stringstream ss;
    ss << start() << "-" << stop() << "\t" << pattern() << "\t" << unitsize() << "x" << repeats() << "=" << span();//repeats();
    return ss.str();
}

repeat_region* repeat_region::analyze_str(int size, char const* buffer, int seq_pos, int unit_min, int unit_max, int required_span) throw (exception) {
    for (int u = 0; u < unit_max; u++) {
        char c = buffer[u];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            return NULL;
        }
    }
    for (int unit = unit_min; unit <= unit_max; unit++) {
        //const char* seed = buffer;
        bool remnants = true;
        for (int u = unit_min; u < unit; u++) {
            if (unit % u == 0) {
                remnants = false;
                break;
            }
        }
        if (!remnants) continue;
        bool repeated = true;
        int stop = unit;
        for (int p = unit; p < required_span; p += unit) {
            if (strncmp(buffer, buffer + p, unit) != 0) {
                repeated = false;
                break;
            }
            stop += unit;
        }
        if (repeated) {
            //int stop = required_span;
            for (;;) {
                if (strncmp(buffer, buffer + stop, unit) != 0) {
                    break;
                }
                stop += unit;
            }
            // for (int p = stop; ; p += unit) {
            //     if (strncmp(buffer, buffer + p, unit) != 0) {
            //         break;
            //     }
            //     stop += unit;
            // }
            return new repeat_region(unit, buffer, seq_pos, seq_pos + stop);
        }
    }
    return NULL;
}

// bool repeat_region::compare_position(const repeat_region& lhs, const repeat_region& rhs) {
//     if (lhs._chromosome != rhs._chromosome) {
//         return lhs._chromosome < rhs._chromosome;
//     } else {
//         return lhs._start < rhs._start;
//     }
// }

bool repeat_region::compare_position(const repeat_region* lhs, const repeat_region* rhs) {
    if (lhs->_chromosome != rhs->_chromosome) {
        return lhs->_chromosome < rhs->_chromosome;
    } else {
        return lhs->_start < rhs->_start;
    }
}

void repeat_region::enumerate_repeat_regions(int argc, char** argv) throw (exception) {
    const char* filename_genome = get_argument_string(argc, argv, "i", NULL);
    const char* filename_output = get_argument_string(argc, argv, "o", NULL);
    int str_span = get_argument_integer(argc, argv, "s", 18);
    int unit_min = get_argument_integer(argc, argv, "m", 2);
    int unit_max = get_argument_integer(argc, argv, "M", 6);
    bool verbose = has_option(argc, argv, "verbose");
    string chromosome;

    if (verbose) {
        cerr << "genome : " << filename_genome << endl;
        cerr << "span : " << str_span << endl;
        cerr << "unit : " << unit_min << ", " << unit_max << endl;
        cerr << "output : " << (filename_output != NULL ? filename_output : "stdout") << endl;
    }
    //int start = 0;
    // int str_span = 18;
    // int unit_min = 2;
    // int unit_max = 6;
    if (filename_genome == NULL || tktools::io::file_exists(filename_genome) == false) {
        throw invalid_argument(string("cannot find ") + filename_genome);
    }
    ifstream fi(filename_genome);//"/Data/mm10/mm10.fa");
    //int position = 0;

    //int margin_size = str_span * 2;
    int buffer_size = 4096;
    char* buffer = new char[buffer_size + unit_max];
    for (int i = 0; i < unit_max; i++) {
        buffer[buffer_size + i] = '\0';
    }
    string sequence;
    int bpos = 0;
    int spos = 0;
    ostream* ost = &cout;
    if (filename_output != NULL) {
        ofstream* fo = new ofstream(filename_output);
        if (fo->is_open() == false) {
            throw invalid_argument("cannot open output file");
        }
        ost = fo;
    }
    while (!fi.eof()) {
        string line;
        getline(fi, line);
        if (line.c_str()[0] == '>') {
            chromosome = line.substr(1);
            spos = 1;
            bpos = 0;
            //sequence = "";
        } else {
            sequence += line;
            for (int i = 0; i < (int)line.size(); i++) {
                char base = line.c_str()[i];
                if (base >= 'A' && base <= 'Z') {
//                if (base == 'A' || base == 'C' || base == 'G' || base == 'T' ) {
                    buffer[bpos++] = base;
                    if (bpos == buffer_size) {
                        for (int i = 0; i < buffer_size - str_span; i++) {
                            repeat_region* rep = analyze_str(buffer_size, buffer + i, spos + i, unit_min, unit_max, str_span);
                            if (rep != NULL) {
                                *ost << chromosome << "\t" << rep->start() << "\t" << rep->stop() << "\t" << rep->span() << "\t" << rep->pattern() << endl;
                                i += rep->span();
                                delete rep;
                                break;
                            }
                        }
                        memmove(buffer, buffer + buffer_size - str_span, str_span);
                        spos += buffer_size - str_span;
                        bpos = str_span;
                    }
                }
            }
        }
    }
    //         for (int i = 0; i < (int)sequence.size() - str_span; i++) {
    //             bool accepted = false;
    //             for (int j = unit_min; j < unit_max; j++) {
                    
    //             }
    //         }
    //     }
    // }
    delete[] buffer;
    if (filename_output != NULL) {
        dynamic_cast<ofstream*>(ost)->close();
        delete ost;
    }
}


