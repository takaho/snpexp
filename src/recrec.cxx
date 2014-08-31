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

using namespace std;

#include <tktools.hxx>
#include <gtf.hxx>
#include <recrec.hxx>

using namespace tktools;
using namespace tktools::io;
using namespace tktools::util;
using namespace tktools::bio;

using namespace tkbio;

#define TEST 1

// recfragment::recfragment(const string& chromosome, int position, int flag, const string& sequence, const string& cigar) {
//     _chromosome = convert_chromosome_to_code(chromosome.c_str());
//     _cigar = cigar;
//     _position = position;
//     _flag = flag;
//     _sequence = sequence;
//     _mapped = resolve_map(position, cigar, sequence);
//     _position5 = 2000000000;
//     _position3 = 0;
//     initialize(position, cigar, sequence);

//     for (vector<pair<int,char> >::const_iterator it = _mapped.begin(); it != _mapped.end(); it++) {
//         int pos = it->first;
//         if (pos < _position5) _position5 = pos;
//         if (pos > _position3) _position3 = pos;
//     }
// }

recfragment::recfragment(int chromosome, bam1_t* seq) {
    _chromosome = chromosome;
    initialize(seq);
}

// vector<pair<int,char> > recfragment::resolve_map(int position, const string& cigar, const string& sequence) {
//     recfragment f("chr1", position, 0, sequence, cigar);
//     return f.get_positions();
// }

namespace {
    char _bamnucleotide[17] = "?AC?G???T??????N";
    char _index2nucleotide[6] = "ACGT-";
}
void recfragment::initialize(bam1_t* read) {
    const uint32_t* cigar = bam1_cigar(read);
    const uint8_t* sequence = bam1_seq(read);
    int len = read->core.n_cigar;
    int spos = read->core.pos;
    int fpos = 0;
    int offset = 4;
    const char* qname = bam1_qname(read);
    int qlen = strlen(qname);
    char* buffer = new char[qlen + 1];
    memcpy(buffer, qname, sizeof(char) * qlen);
    buffer[qlen] = '\0';
    _name = buffer;
    delete[] buffer;
    _max_match_span = 0;
    _position = _position5 = spos;
    for (int i = 0; i < len; i++) {
        int op = bam_cigar_op(cigar[i]);
        int slen = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) { // M
            if (slen > _max_match_span) {
                _max_match_span = slen;
            }
            for (int j = 0; j < slen; j++) {
                //int offset = (fpos & 1) == 0 ? 4 : 0;
                uint8_t base = (sequence[fpos >> 1] >> offset) & 15;
                offset = 4 - offset;
                //cout << spos << ":" << (int)(base) << ":" << _bamnucleotide[base] << endl;
                _mapped.push_back(make_pair(spos++, _bamnucleotide[base]));
                fpos++;
            }
        } else if (op == BAM_CINS) {
            fpos += slen;
            //ss << 'I';
        } else if (op == BAM_CDEL) {
            for (int j = 0; j < slen; j++) {
                _mapped.push_back(make_pair(spos++, '.'));
            }
                                  //spos += slen;
            //ss << 'D';
        // } else if (op == BAM_CREF_SKIP) {
        //     ss << 'N';
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) { // S or H
            fpos ++;
            spos ++;
        // } else if (op == BAM_CPAD) {
        //     ss << 'P';
        // } else if (op == BAM_CEQUAL) {
        //     ss << 'M';
        // } else if (op == BAM_CDIFF) {
        //     ss << 'X';
        } else if (op == BAM_CBACK) {
            break;
        }
    }
    _position3 = spos;
    _flag = read->core.flag;
}

// void recfragment::initialize(int position, const string& sequence, const string& cigar){
//     int value = 0;
//     int spos = position;
//     int fpos = 0;
//     //vector<pair<int,char> > nucleotides;
//     const char* ptr = cigar.c_str();
//     const char* seq = sequence.c_str();
//     int max_matches = 0;
//     int left = 2000000000;
//     int right = 0;
//     for (int i = 0; i < (int)cigar.size(); i++) {
//         char c = ptr[i];
//         if (c >= '0' && c <= '9') {
//             if (value == 0) {
//                 value = atoi(ptr + i);
//             }
//         } else if (c == 'M') {
//             if (max_matches < value) max_matches = value;
//             if (left > spos) left = spos;
//             for (int j = 0; j < value; j++) {
//                 _mapped.push_back(make_pair(spos, seq[spos]));
//                 spos++;
//                 fpos++;
//             }
//             right = spos;
//         } else if (c == 'S' || c == 'H') { // soft clip & hard clip
//             fpos += value;
//             spos += value;
//         } else if (c == 'I') {
//             fpos += value;
//         } else if (c == 'D') {
//             for (int j = 0; j < value; j++) {
//                 _mapped.push_back(make_pair(spos, '-'));
//                 spos++;
//             }
//         } else {
//             throw runtime_error(string("unknown cigar character :") + c); 
//         }
//     }
//     _max_match_span = max_matches;
//     _position5 = left;
//     _position3 = right;
// }

char recfragment::get_base(int pos, int& count) const {
    int left = 0;
    int size = _mapped.size();
    int right = size;
    for (;;) {
        int center = (left + right) / 2;
        const pair<int,char>& p = _mapped[center];
        if (p.first < pos) {
            left = center + 1;
        } else if (p.first > pos) {
            right = center;
        } else {
            count = 1;
            if (center > 0 && _mapped[center - 1].first == pos) count++;
            if (center < size - 1 && _mapped[center - 1].first == pos) count++;
            return p.second;
        }
        if (left == right) {
            count = 0;
            return '\0';
        }
    }
}
void recfragment::generate_recombination_pattern(const vector<pair<int,char> >& loci, int* pattern) const {
    int index = 0;
    for (int i = 0; i < (int)loci.size(); i++) {
        int pos = loci[i].first;
        char base = loci[i].second;
        int value = 0;
        while (index < (int)_mapped.size()) {
            const pair<int,char>& locus = _mapped[index];
            if (locus.first > pos) {
                break;
            } else if (locus.first == pos) {
                value = locus.second == base ? -1 : 1;
                break;
            } else {
                index++;
            }
        }
        pattern[i] = value;
    }
}

pair<int,int> recfragment::get_recombination(const vector<pair<int,char> >& loci) const {
    const int minimum_diff = 2;
    int size = loci.size();
    int* buffer = new int[size];
    generate_recombination_pattern(loci, buffer);
    int left_score = 0;
    int right_score = 0;
    int limit_left = 0;//size;
    int limit_right = size;//0;
    // int n = 0;
    // for (int i = 0; i < size; i++) {
    //     if (buffer[i] != 0) {
    //         left_score += buffer[i];
    //         if (++n >= minimum_span) {
    //             limit_left = n;
    //             break;
    //         }
    //     }
    //     right_score += buffer[i];
    // }
    // n = 0;
    // for (int i = size - 1; i >= limit_left; i--) {
    //     if (buffer[i] != 0) {
    //         if (++n >= minimum_span) {
    //             limit_right = i;
    //             for (int j = limit_left; j < size; j++) {
    //                 right_score += buffer[j];
    //             }
    //             break;
    //         }
    //     }
    // }
    int index = limit_left;
    int max_diff = 0;
    pair<int,int> border = make_pair(-1,-1);
    while (index < limit_right) {
        //for (int i = 2; i < size - 2; i++) {
        int l1 = buffer[index];
        if (l1 != 0) {
            for (int j = index + 1; j < limit_right; j++) {
                int l2 = buffer[j];
                if (l2 != 0) {
                    if (l2 != l1) {
                        int diff = abs(left_score - right_score);
                        if (diff > max_diff && diff >= minimum_diff) {
                            max_diff = diff;
                            border = make_pair(_mapped[index].first, _mapped[j].first);
                        }
                    }
                    index = j;
                    left_score += l1;
                    right_score -= l1;
                    break;
                }
            }
        }
        //int l2 = buffer[i + 1];
    }

    delete[] buffer;
    return border;
}

string recfragment::get_recombination_pattern(const vector<pair<int,char> >& loci) const {
    int size = loci.size();
    int* buffer = new int[size];
    char* pat = new char[size + 1];
    generate_recombination_pattern(loci, buffer);
    for (int i = 0; i < size; i++) {
        if (buffer[i] == -1) {
            pat[i] = '_';
        } else if (buffer[i] == 1) {
            pat[i] = '+';
        } else {
            pat[i] = '.';
        }
    }
    pat[size] = '\0';
    string patern(pat);
    delete[] pat;
    delete[] buffer;
    return pat;
}

namespace {
    // string resolve_cigar(const bam1_t* read) {
    //     const uint32_t* cigar = bam1_cigar(read);
    //     int len = read->core.n_cigar;
    //     stringstream ss;
    //     for (int i = 0; i < len; i++) {
    //         int op = bam_cigar_op(cigar[i]);
    //         int slen = bam_cigar_oplen(cigar[i]);
    //         ss << slen;
    //         if (op == BAM_CMATCH) {
    //             ss << 'M';
    //         } else if (op == BAM_CINS) {
    //             ss << 'I';
    //         } else if (op == BAM_CDEL) {
    //             ss << 'D';
    //         } else if (op == BAM_CREF_SKIP) {
    //             ss << 'N';
    //         } else if (op == BAM_CSOFT_CLIP) {
    //             ss << 'S';
    //         } else if (op == BAM_CHARD_CLIP) {
    //             ss << 'H';
    //         } else if (op == BAM_CPAD) {
    //             ss << 'P';
    //         } else if (op == BAM_CEQUAL) {
    //             ss << 'M';
    //         } else if (op == BAM_CDIFF) {
    //             ss << 'X';
    //         } else if (op == BAM_CBACK) {
    //             break;
    //         }
    //     }
    //     return ss.str();
    // }
    
    int _get_chromosome_code(int length, const char* line) {
        const char* fpos = strstr(line, "chromosome");
        int code = -1;
        if (fpos == NULL) {
            fpos = strstr(line, "chr");
            if (fpos == NULL) {
                code = convert_chromosome_to_code(line);
            } else {
                code = convert_chromosome_to_code(fpos + 3);
            }
        } else {
            fpos += 10;
            while (true) {
                char c = *fpos;
                if (c > ' ') {
                    code = convert_chromosome_to_code(fpos);
                    break;
                } else if (c == '\0') {
                    code = -1;
                    break;
                }
                fpos++;
            }
        }
        return code;
    }

    map<int,pair<int,char*> > load_genome(const char* filename) throw (exception) {
        map<int,pair<int,char*> > chromosomes;
        ifstream fi(filename);
        if (!fi.is_open()) {
            throw invalid_argument("cannot open genome file");
        }
        size_t size_buffer = 300000000;
        size_t size_seq = 0;
        char* buffer = new char[size_buffer];
        int chrmcode = -1;
        while (!fi.eof()) {
            string line;
            getline(fi, line);
            if (line.c_str()[0] == '>') {
                if (chrmcode > 0) {
                    char* clone = new char[size_seq + 1];
                    memcpy(clone, buffer, size_seq);
                    clone[size_seq] = '\0';
                    chromosomes.insert(make_pair(chrmcode, make_pair(size_seq, clone)));
#ifdef TEST
                    break;
#endif
                }
                chrmcode = _get_chromosome_code(line.size() - 1, line.c_str() + 1);
                size_seq = 0;
            } else if (chrmcode > 0) {
                const char* ptr = line.c_str();
                for (int i = 0; i < (int)line.size(); i++) {
                    char c = ptr[i];
                    if (c >= 'a' && c <= 'z') { c -= (char)('a' - 'A'); }
                    if (c >= 'A' && c <= 'Z') {
                        buffer[size_seq++] = c;
                    }
                }
            }
        }
        if (chrmcode > 0) {
            char* clone = new char[size_seq + 1];
            memcpy(clone, buffer, size_seq);
            clone[size_seq] = '\0';
            chromosomes.insert(make_pair(chrmcode, make_pair(size_seq, clone)));
        }
        delete[] buffer;
        fi.close();
        return chromosomes;
    }
}

namespace {
    class hetero_locus {
    public:
        int position;
        int index1;
        int index2;
    public:
        hetero_locus(int position, int index1, int index2) {
            this->position = position;
            this->index1 = index1;
            this->index2 = index2;
        }
        hetero_locus(int position, const string& ref, const string& alt) {
            this->position = position;
            this->index1 = get_index(ref);
            this->index2 = get_index(alt);
        }
        bool is_available() const {
            return 0 <= index1 && index1 < 5 && 0 <= index2 && index2 < 5;
        }
    private:
        static int get_index(const string& locus) {
            if (locus == "A") {
                return 0;
            } else if (locus == "C") {
                return 1;
            } else if (locus == "G") {
                return 2;
            } else if (locus == "T") {
                return 3;
            } else if (locus == "." || locus == "-") {
                return 4;
            } else {
                return -1;
            }
        }
    };
}

namespace {
    bool compare_first_index(const pair<int,char>& lhs, const pair<int,char>& rhs) {
        return lhs.first < rhs.first;
    }
}

void recfragment::join_sequence(const recfragment* frag) {
    vector<pair<int,char> > bases = frag->_mapped;
    for (vector<pair<int,char> >::const_iterator it = frag->_mapped.begin();
         it != frag->_mapped.end(); it++) {
        _mapped.push_back(*it);
    }
    if ((_position3 > frag->_position5 && _position5 <= frag->_position3)) {
        sort(_mapped.begin(), _mapped.end(), compare_first_index);
        //vector<int> ambiguous;
        for (int i = (int)_mapped.size() - 1; i > 0; i--) {
            if (_mapped[i].first == _mapped[i-1].first) {
                if (_mapped[i].second != _mapped[i-1].second) {
                    //ambiguous.push_back(i);
                    _mapped.erase(_mapped.begin() + i, _mapped.begin() + i + 2);
                }
                i--;
            }
        }
    }
    _position3 = _position3 > frag->_position3 ? _position3 : frag->_position3;
    _position5 = _position5 < frag->_position5 ? _position5 : frag->_position5;
}

// recfragment* recfragment::join(recfragment const* f1, recfragment const* f2) {
//     if (f1->_position5 > f2->_position5) {
//         recfragment const* ft = f1;
//         f1 = f2;
//         f2 = ft;
//     }
//     recfragment* frag = new recfragment();
//     frag->_name = f1->_name;
//     frag->_flag = 0;
//     frag->_position = f1->_position;
//     int last_position = 0;
//     for (int i = 0; i < (int)f1->_mapped.size(); i++) {
//         const pair<int,char>& p = f1->_mapped[i];
//         if (p.first > last_position) {
//             last_position = p.first;
//         }
//         frag->_mapped.push_back(p);
//     }
//     for (int i = 0; i < (int)f2->_mapped.size(); i++) {
//         const pair<int,char>& p = f2->_mapped[i];
//         if (p.first < last_position) {
            
//         }
//         frag->_mapped.push_back(p);
//     }
//     frag->_position5 = f1->_position5;
//     frag->_position3 = f2->_position3;
//     frag->_max_match_span = f1->_max_match_span < f2->_max_match_span ? f2->_max_match_span : f1->_max_match_span;
//     return frag;
// }

void recfragment::bundle_pairs(vector<recfragment*>& fragments) throw (runtime_error) {
    vector<recfragment*> bundled;
    for (int i = 0; i < (int)fragments.size(); i++) {
        recfragment* p = fragments[i];
        if (p == NULL) continue;
        for (int j = i + 1; j < (int)fragments.size(); j++) {
            recfragment* q = fragments[j];
            if (p->name() == q->name()) {
                p->join_sequence(q);
                delete q;
                fragments[j] = NULL;
                break;
            }
        }
    }
}

ostream& operator << (ostream& ost, const recfragment& rhs) {
    ost << rhs.name() << " ; " << rhs.chromosome() << ":" << rhs.orientation() << ":" << rhs.position5() << "-" << rhs.position3() << " ";
    return ost;
}


namespace {
    vector<pair<int,char> > 
    scan_heterozygous_loci(vector<recfragment*>& fragments,
                           const string& chromosome_name,
                           int chromosome_length,
                           char const* chromosome_sequence,
                           int start, int end,
                           int minimum_coverage=10,
                           float hetero_threshold=0.2f) throw (exception) {
        //vector<recfragment*> bound = bundle_pairs(fragments);
        int freq[5];
        vector<pair<int,char> > hetero_loci;
        for (int pos = start; pos < end; pos++) {
            char ref = chromosome_sequence[pos];
            int refid;
            switch (ref) {
            case 'A': refid = 0; break;
            case 'C': refid = 1; break;
            case 'G': refid = 2; break;
            case 'T': refid = 3; break;
            case '-': refid = 4; break;
            default: refid = -1; break;
            }
            for (int i = 0; i < 5; i++) freq[i] = 0;
            int coverage = 0;
            for (int j = 0; j < (int)fragments.size(); j++) {
                const recfragment* fr = fragments[j];
                if (fr->position5() <= pos && pos < fr->position3()) {
                    int num = 0;
                    char base = fr->get_base(pos, num);
                    int index;
                    switch (base) {
                    case 'A': index = 0; break;
                    case 'C': index = 1; break;
                    case 'G': index = 2; break;
                    case 'T': index = 3; break;
                    case '-': index = 4; break;
                    default:
                        index = -1;
                        continue;
                    }
                    if (base < 0) continue;
                    freq[index] += num;
                    coverage++;
                }
            }
            if (coverage >= minimum_coverage) {
                int lower = (int)(coverage * hetero_threshold + 0.5f);
                int upper = coverage - lower;
                int altnum = coverage - freq[refid];
                if (lower <= altnum && altnum <= upper) {
                    int max_minor = 0;
                    int alt_index = -1;
                    for (int i = 0; i < 5; i++) {
                        if (i != refid && freq[i] > max_minor) {
                            max_minor = freq[i];
                            alt_index = i;
                        }
                    }
                    if (max_minor > 1) {
                        hetero_loci.push_back(make_pair(pos, _index2nucleotide[alt_index]));
                        cout << chromosome_name << ":" << pos << "\t" << ref << "\t" << _index2nucleotide[alt_index] << "\tA:" << freq[0] << " C:" << freq[1] << " G:" << freq[2] << " T:" << freq[3] << " -:" << freq[4] << endl;
                        
                    }
                }
            }
        }
        return hetero_loci;
    }

    void detect_recombination(vector<recfragment*>& fragments,
                              const string& chromosome_name,
                              int chromosome_length,
                              char const* chromosome_sequence,
                              int start, int end,
                              int minimum_coverage=10,
                              float hetero_threshold=0.20f) throw (exception) {
        // bind pairs
        recfragment::bundle_pairs(fragments);
        vector<pair<int,char> > hetero_loci = scan_heterozygous_loci(fragments, chromosome_name, chromosome_length, chromosome_sequence, start, end, minimum_coverage, hetero_threshold);
        for (int i = 0; i < (int)fragments.size(); i++) {
            pair<int,int> site = fragments[i]->get_recombination(hetero_loci);
            if (site.first < site.second) {
                string pattern = fragments[i]->get_recombination_pattern(hetero_loci);
                cout << pattern << "\t" << fragments[i]->name();
            }
        }
    }
}
    
int main(int argc, char** argv) {
    try {
        const char* filename = get_argument_string(argc, argv, "i", NULL);
        int coverage = get_argument_integer(argc, argv, "c", 10);
        double heterozygosity = get_argument_float(argc, argv, "H", 0.2);
        bool verbose = has_option(argc, argv, "verbose");
        int window_size = get_argument_integer(argc, argv, "w", 10000);
        int window_margin = get_argument_integer(argc, argv, "m", 2000);
        //const char* snpfile = get_argument_string(argc, argv, "V", NULL);
        const char* filename_genome = get_argument_string(argc, argv, "g", NULL);

        if (verbose) {
            cerr << "Filename         : " << filename << endl;
            cerr << "Minimum coverage : " << coverage << endl;
            cerr << "Heterozygosity   : " << heterozygosity << endl;
            //cerr << "Variation file   : " << snpfile << endl;
            cerr << "Genome file      : " << filename_genome << endl;
        }

        vector<recfragment*> fragments;
        map<string,const recfragment*> pairs;
        bam1_t* read;
        bam_header_t* header;
        bamFile bamfile;

        map<int,pair<int,char*> >  genome = load_genome(filename_genome);

        // map<int,vector<hetero_locus*> > snps;
        // if (snpfile != NULL) {
        //     ifstream fs(snpfile);
        //     if (!fs.is_open()) {
        //         throw invalid_argument("cannot open SNP file");
        //     }
        //     while (!fs.eof()) {
        //         string line;
        //         getline(fs, line);
        //         vector<string> items = split_items(line, '\t');
        //         int chrm_code = convert_chromosome_to_code(items[0].c_str());
        //         vector<hetero_locus*>& chrmsnps = snps.find(chrm_code);
        //         hetero_locus* locus = hetero_locus(atoi(items[1]), items[2], items[3]);
        //         if (locus->is_available()) {
        //             if (chrmsnps == snps.end()) {
        //                 vector<hetero_locus*> cs;
        //                 cs.push_back(locus);
        //                 snps.push_back(cs);
        //             } else {
        //                 chrmspns.push_back(locus);
        //             }
        //         } else {
        //             delete locus;
        //         }
        //     }
        //     fs.close();
        // } else {
        //     throw invalid_argument("snp file required");
        // }

        bamfile = bam_open(filename, "rb");
        header = bam_header_read(bamfile);
        read = bam_init1();

        ////int chromosome = -1;
        //int region_start = 0;
        //int region_end = region_start + window_size;

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

        //int position_last = 0;
        //int position_next = window_size;

        int retain_start = 0;
        int retain_end = window_size;
        int analysis_start = retain_start + window_margin;
        int analysis_end = retain_end - window_margin;

//        int region_start, region_end;
//        int region_start = 0;
//        int region_end = region_start + window_size;
        int current_bamchrm = -1;
        int current_chromosome = -1;
        int next_chromosome = -1;
        int chromosome_length = 0;
        string chromosome_name;
        const char* chromosome_sequence = NULL;

        while (bam_read1(bamfile, read) > 0) {
            bool processing = false;

            if (current_bamchrm != read->core.tid) {
                // finish chromosome
                //if (fragments.size() > 0) {
                processing = true;
                //    retain_end = read->core.pos + read->core.l_qseq * 2;
                //}
                const char* name = header->target_name[read->core.tid];
                next_chromosome = _get_chromosome_code(strlen(name), name);
            } else if (read->core.pos > retain_end) {
                // process current region
                if (fragments.size() > 0) {
                    processing = true;
                }
            }
            //cout << fragments.size() << endl;
            if (processing) {
                // analysis
                if (fragments.size() > coverage && chromosome_length > 0) {
#ifdef TEST
                    //cerr << "test\n";
                    cout << current_chromosome << " : " << retain_start << "," << analysis_start << "," << analysis_end << "," << retain_end << " : " << fragments.size() << endl;
                    if (fragments.size() > 1) {
                        cout << *fragments[0] << endl;
                        cout << *fragments[fragments.size() - 1] << endl;
                    }
#endif
                    detect_recombination(fragments, chromosome_name, chromosome_length, chromosome_sequence, analysis_start, analysis_end);
                }

                // change chromosome
                if (next_chromosome != current_chromosome) {
                    //cerr << "change chromosoe\n";

                    map<int, pair<int,char*> >::const_iterator _chrm = genome.find(next_chromosome);
                    if (_chrm == genome.end()) {
                        chromosome_sequence = NULL;
                        chromosome_length = 0;
                        cerr << "cannot find " << next_chromosome << endl;
                    } else {
                        chromosome_sequence = _chrm->second.second;
                        chromosome_length = _chrm->second.first;
                        chromosome_name = convert_code_to_chromosome(next_chromosome);
                        cerr << "new chromosome " << chromosome_name << " // " << chromosome_length << endl;
                    }

                    retain_start = 0;
                    analysis_start = 0;
                    analysis_end = window_size - window_margin;
                    retain_end = window_size;
                    //cerr << fragments.size() << " => ";
                    for (vector<recfragment*>::iterator it = fragments.begin(); it != fragments.end(); it++) {
                        delete *it;
                    }
                    fragments.erase(fragments.begin(), fragments.end());
                    cerr << fragments.size() << "\n";
                    current_bamchrm = read->core.tid;
                    current_chromosome = next_chromosome;
                } else {
                    //cerr << "analysis\n";

                    //int limit = analysis_end - window_margin;
                    int end = 0;
                    for (int i = 0; i < (int)fragments.size(); i++) {
                        if (fragments[i]->position3() < analysis_start) {
                            delete fragments[i];
                        } else {
                            end = i; 
                            break;
                        }
                    }
                    fragments.erase(fragments.begin(), fragments.begin() + end);
                    // for (vector<recfragment*>::iterator it = fragments.begin(); it != fragments.end(); it++) {
                    //     delete *it;
                    // }
                    analysis_start = analysis_end;
                    retain_start = analysis_end - window_margin;
                    retain_end = retain_start + window_size;
                    analysis_end = retain_end - window_margin;
                    if (fragments.size() > 0 && fragments[0]->position5() > analysis_start) {
                        analysis_start = fragments[0]->position5();
                        retain_start = analysis_start - window_margin;
                        retain_end = retain_start + window_size;
                        analysis_end = retain_end - window_margin;
                    }
                }
            }
            if (current_chromosome > 0) {
                recfragment* frag = new recfragment(current_chromosome, read);
                fragments.push_back(frag);
            }
        }

        // process remaining
        if (fragments.size() > coverage && chromosome_length > 0) {
            detect_recombination(fragments, chromosome_name, chromosome_length, chromosome_sequence, analysis_start, analysis_end);
        }


        // clean up
        for (vector<recfragment*>::iterator it = fragments.begin(); it != fragments.end(); it++) {
            delete *it;
        }
        for (map<int, pair<int,char*> >::iterator it = genome.begin(); it != genome.end(); it++) {
            delete[] it->second.second;
        }

//            if (chrm != chromosome
//            string cigar = resolve_cigar(read);
            //header;
            //fragment f(chromosome, position, flag, sequence, cigar);
            //read->core.tid; chromosome
            //}

        bam_destroy1(read);
        bam_header_destroy(header);
        bam_close(bamfile);

        return 0;
    } catch (exception& e) {
        cerr << e.what() << endl;
        return -1;
    }
}
