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

//#define TEST 1

recfragment::recfragment(int chromosome, bam1_t* seq, int num) {
    _chromosome = chromosome;
    _group = num;
    initialize(seq);
}

namespace {
    char _bamnucleotide[17] = "\?AC\?G\?\?\?T\?\?\?\?\?-N";
    unsigned char base2code(char n) {
        if (n == 'A') {
            return 0x1;
        } else if (n == 'C') {
            return 0x2;
        } else if (n == 'G') {
            return 0x4;
        } else if (n == 'T') {
            return 0x8;
        } else if (n == 'N') {
            return 0xf;
        } else if (n == '-') {
            return 0xe;
        } else {
            return 0;
        }
    }

    string resolve_cigar(const bam1_t* read) {
        const uint32_t* cigar = bam1_cigar(read);
        int len = read->core.n_cigar;
        stringstream ss;
        for (int i = 0; i < len; i++) {
            int op = bam_cigar_op(cigar[i]);
            int slen = bam_cigar_oplen(cigar[i]);
            ss << slen;
            if (op == BAM_CMATCH) {
                ss << 'M';
            } else if (op == BAM_CINS) {
                ss << 'I';
            } else if (op == BAM_CDEL) {
                ss << 'D';
            } else if (op == BAM_CREF_SKIP) {
                ss << 'N';
            } else if (op == BAM_CSOFT_CLIP) {
                ss << 'S';
            } else if (op == BAM_CHARD_CLIP) {
                ss << 'H';
            } else if (op == BAM_CPAD) {
                ss << 'P';
            } else if (op == BAM_CEQUAL) {
                ss << 'M';
            } else if (op == BAM_CDIFF) {
                ss << 'X';
            } else if (op == BAM_CBACK) {
                break;
            }
        }
        return ss.str();
    }
}

namespace {
    const unsigned char QUAL_GAP = 10;
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
    uint8_t* qptr = bam1_qual(read);
    for (int i = 0; i < len; i++) {
        int op = bam_cigar_op(cigar[i]);
        int slen = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) { // M
            if (slen > _max_match_span) {
                _max_match_span = slen;
            }
            for (int j = 0; j < slen; j++) {
                uint8_t base = (sequence[fpos >> 1] >> offset) & 15;
                unsigned char qual = qptr[fpos];
                offset = 4 - offset;
                _mapped.push_back(qseqbase(spos, _bamnucleotide[base], qual));
                //_mapped.push_back(make_pair(spos, _bamnucleotide[base]));
                spos++;
                fpos++;
            }
        } else if (op == BAM_CINS) { // I
            fpos += slen;
            if ((slen & 1) != 0) {
                offset = 4 - offset;
            }
        } else if (op == BAM_CDEL) { // D
            for (int j = 0; j < slen; j++) {
                _mapped.push_back(qseqbase(spos + j, '-', QUAL_GAP));
//make_pair(spos + j, '-'));
            }
            spos += slen;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) { // S or H
            break;
        } else if (op == BAM_CBACK) {
            break;
        } else {
            cerr << "unknown " << op << ", " << slen << " " << resolve_cigar(read) << endl;
        }
    }
    _position3 = spos;
    _flag = read->core.flag;
}

char recfragment::get_base(int pos, unsigned char qual_threshold, int& count) const {
    int left = 0;
    int size = _mapped.size();
    int right = size;
    for (;;) {
        int center = (left + right) / 2;
        const qseqbase& p = _mapped[center];
        if (p.position() < pos) {
            left = center + 1;
        } else if (p.position() > pos) {
            right = center;
        } else {
            count = 1;
            for (int i = center - 1; i <= center + 1; i++) {
                if (i >= 0 && i < size) {
                    const qseqbase& qb = _mapped[i];
                    if (qb.position() == pos && qb.quality() >= qual_threshold) {
                        count++;
                    }
                }
            }
            if (count > 0) {
                return p.nucleotide();
            } else {
                return '\0';
            }
        }
        if (left == right) {
            count = 0;
            return '\0';
        }
    }
}

char recfragment::get_base(int pos, int& count) const {
    int left = 0;
    int size = _mapped.size();
    int right = size;
    for (;;) {
        int center = (left + right) / 2;
        const qseqbase& p = _mapped[center];
        if (p.position() < pos) {
            left = center + 1;
        } else if (p.position() > pos) {
            right = center;
        } else {
            count = 1;
            if (center > 0 && _mapped[center - 1].position() == pos) count++;
            if (center < size - 1 && _mapped[center + 1].position() == pos) count++;
            return p.nucleotide();
        }
        if (left == right) {
            count = 0;
            return '\0';
        }
    }
}

// A,C,G,T,- => 0,1,2,3,4
// others => -1
int recfragment::get_base_id(int pos, int& count) const {
    char base = get_base(pos, count);
    if (base == 'A') {
        return 0;
    } else if (base == 'C') {
        return 1;
    } else if (base == 'G') {
        return 2;
    } else if (base == 'T') {
        return 3;
    } else if (base == '-') {
        count = 1;
        return 4;
    } else {
        return -1;
    }
}

int recfragment::get_base_id(int pos, unsigned char threshold, int& count) const {
    char base = get_base(pos, threshold, count);
    if (base == 'A') {
        return 0;
    } else if (base == 'C') {
        return 1;
    } else if (base == 'G') {
        return 2;
    } else if (base == 'T') {
        return 3;
    } else if (base == '-') {
        count = 1;
        return 4;
    } else {
        return -1;
    }
}

void recfragment::generate_recombination_pattern(const vector<hetero_locus*>& loci, int* pattern) const {
    int index = 0;
    for (int i = 0; i < (int)loci.size(); i++) {
        int pos = loci[i]->position();
        char ref = loci[i]->ref();//second;
        char alt = loci[i]->alt();
        int value = 0;
        while (index < (int)_mapped.size()) {
            const qseqbase& locus = _mapped[index];
            if (locus.position() > pos) break;
            if (locus.position() == pos) {
                if (locus.nucleotide() == ref) {
                    value = 1;
                } else if (locus.nucleotide() == alt) {
                    value = -1;
                } else {
                    value = 0;
                }
                break;
            } else {
                index++;
            }
        }
        pattern[i] = value;
    }
}

vector<recpattern> recfragment::get_recombination_borders(const vector<hetero_locus*>& loci, int minimum_span) const {
    vector<recpattern> borders;
    int size = loci.size();
    int* buffer = new int[size];
    generate_recombination_pattern(loci, buffer);
    int i = 0;
    while (i < size) {
        //for (int i = 0; i < size; i++) {
        int num_cont = 1;
        int gt = buffer[i];
        int last = size;
        if (gt != 0) {
            num_cont = 1;
        } else {
            continue;
        }
        for (int j = last; j < size; j++) {
            int gj = buffer[j];
            if (gj == 0) continue;
            if (gj == gt) {
                last = j;
                num_cont ++;
            } else {
                break;
            }
        }
        if (num_cont >= minimum_span) {
            int num_cont2 = 0;
            for (int j = last; j < size; j++) {
                int gj = buffer[j];
                if (gj == 0) continue;
                if (gj != gt) {
                    num_cont2 ++;
                } else {
                    break;
                }
            }
            if (num_cont2 >= minimum_span) {
                recpattern::Genotype genotype;
                if (gt == 1) {
                    genotype = recpattern::REF_ALT;
                } else {
                    genotype = recpattern::ALT_REF;
                }
                borders.push_back(recpattern(loci[i]->position(), loci[i+1]->position(), genotype, num_cont - num_cont2));
            }
        }
        i = last;
    }
    delete[] buffer;
    return borders;
}

recpattern::Genotype recfragment::get_recombination_pattern(int pos5, int pos3, int minimum_span) const {
    recpattern::Genotype genotype = recpattern::UNDETERMINED;
    //int index5 = 0;
    //int index3 = 0;
    
    // for (; index5 >= 0; index5--) {
    // }
    // for (; index3 < size; index3++) {
    // }
    return genotype;
}

pair<int,int> recfragment::get_recombination(const vector<hetero_locus*>& loci, float diff_ratio) const {
    const int minimum_diff = 2;
    int size = loci.size();
    int* buffer = new int[size];
    generate_recombination_pattern(loci, buffer);
    int left_score = 0;
    int right_score = 0;
    int limit_left = 0;
    int limit_right = size;
    int index = limit_left;
    int max_diff = 4;
    int max_index = -1;
    pair<int,int> border = make_pair(-1,-1);
    for (int i = 0; i < size; i++) {
        right_score += buffer[i];
    }
    while (index < limit_right) {
        //for (int i = 2; i < size - 2; i++) {
        int l1 = buffer[index];
        bool found = false;
        if (l1 != 0) {
            for (int j = index + 1; j < limit_right; j++) {
                int l2 = buffer[j];
                if (l2 != 0) {
                    if (l2 != l1) {
                        int diff = abs(left_score - right_score);
                        if (diff >= (abs(left_score) + abs(right_score)) * diff_ratio
                            && diff > max_diff 
                            && abs(left_score) >= minimum_diff
                            && abs(right_score) >= minimum_diff) {
                            max_diff = diff;
                            max_index = index;
                            border = make_pair(loci[index]->position(), loci[j]->position());
                        }
                    }
                    found = true;
                    break;
                }
            }
            left_score += l1;
            right_score -= l1;
        }
        index++;
    }

    delete[] buffer;
    return border;
}

//string recfragment::get_recombination_pattern(const vector<pair<int,char> >& loci) const {
string recfragment::get_recombination_pattern(const vector<hetero_locus*>& loci) const {
    int size = loci.size();
    int* buffer = new int[size];
    char* pat = new char[size + 1];
    generate_recombination_pattern(loci, buffer);
    int pos_last = 0;
    for (int i = 0; i < size; i++) {
        if (buffer[i] == -1) {
            pat[i] = 'a';
            pos_last = i;
        } else if (buffer[i] == 1) {
            pat[i] = 'A';
            pos_last = i;
        } else {
            pat[i] = '.';
        }
    }
    pat[pos_last] = '\0';
    string pattern(pat);
    delete[] pat;
    delete[] buffer;
    return pattern;
}

string hetero_locus::chromosome() const {
    return convert_code_to_chromosome(_chromosome);
}

string hetero_locus::to_string() const {
    stringstream ss;
    ss << chromosome() << ":" << position() << "\t" << ref() << "\t" << alt();
    return ss.str();
}

void chromosome_seq::set_chromosome(int num) {
    _code = num;
    _name = convert_code_to_chromosome(num);
}

unsigned char chromosome_seq::get_base_code(int pos) const {
    unsigned char code = _sequence[pos >> 1] >> ((pos & 1) == 0 ? 4 : 0) & 0x0f;
    return code;
}

namespace {
    int _base_identifiers[16] = {-1,0,1,-1,2,-1,-1,-1,3,-1,-1,-1,-1,-1,4,-1};
}

int chromosome_seq::get_base_id(int pos) const {
    unsigned char code = _sequence[pos >> 1] >> ((pos & 1) == 0 ? 4 : 0) & 0x0f;
    return _base_identifiers[code];
    //return code;
}

char chromosome_seq::get_base(int pos) const {
    //unsigned char code = _sequence[pos >> 1] >> ((pos & 1) == 0 ? 4 : 0) & 0x0f;
    return _bamnucleotide[get_base_code(pos)];
}

void chromosome_seq::mask(int start, int end) {
    if ((start & 1) != 0) {
        set_base(start, 'N');
        start++;
    }
    if ((end & 1) != 1) {
        set_base(end - 1, 'N');
        end--;
    }
    int index_start = start >> 1;
    int index_end = end >> 1;
    for (int i = index_start; i < index_end; i++) {
        _sequence[i] = 0xff;
    }
}

void chromosome_seq::set_base(int pos, char nuc) {
    int index = pos >> 1;
    bool odd = (pos & 1) != 0;
    unsigned char b;
    switch (nuc) {
    case 'N': case 'n': b = 0x0f; break;
    case 'A': case 'a': b = 0x01; break;
    case 'C': case 'c': b = 0x02; break;
    case 'G': case 'g': b = 0x04; break;
    case 'T': case 't': b = 0x08; break;
    case '-': b = 0x0e; break;
    default:
        b = 0; break;
    }
    if (b == 0) return;
    if (odd) {
        _sequence[index] = (_sequence[index] & 0xf0) | b;
    } else {
        _sequence[index] = (_sequence[index] & 0x0f) | (b << 4);
    }
}

int chromosome_seq::get_chromosome_code(int length, const char* line) {
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

void chromosome_seq::set_sequence(int length, char const* codes) {
    if (_sequence != NULL) {
        delete[] _sequence;
    }
    _length = length;
    int size = length / 2 + 1;
    _sequence = new unsigned char[size];
    for (int i = 0; i < length; i+= 2) {
        unsigned char binuc = (base2code(codes[i]) << 4) | base2code(codes[i + 1]);
        _sequence[i >> 1] = binuc;
    }
    if (length % 2 == 1) {
        _sequence[size - 1] = base2code(codes[length - 1]) << 4;
    }
}

vector<chromosome_seq*> chromosome_seq::load_genome(const char* filename) throw (exception) {
    vector<chromosome_seq*> chromosomes;
    ifstream fi(filename);
    if (!fi.is_open()) {
        throw invalid_argument("cannot open genome file");
    }
    size_t length = 0;
    int size_buffer = 300000000;
    char* buffer = new char[size_buffer];
    buffer[0] = 0;
    int chrmcode = -1;
    chromosome_seq* seq = NULL;
    while (!fi.eof()) {
        string line;
        getline(fi, line);
        if (line.c_str()[0] == '>') {
            if (seq != NULL) {
                seq->set_sequence(length, buffer);
                chromosomes.push_back(seq);
                seq = NULL;
#ifdef TEST
                break;
#endif
            }
            chrmcode = chromosome_seq::get_chromosome_code(line.size() - 1, line.c_str() + 1);
            length = 0;
            if (chrmcode > 0) {
                seq = new chromosome_seq();
                seq->_name = line.substr(1, line.size());
                seq->_code = chrmcode;
            } else {
                seq = NULL;
            }
        } else if (seq != NULL && length < size_buffer) {
            const char* ptr = line.c_str();
            for (int i = 0; i < (int)line.size(); i++) {
                char base = '\0';
                switch (ptr[i]) {
                case 'a': case 'A':
                    base = 'A'; break;
                case 'c': case 'C':
                    base = 'C'; break;
                case 'g': case 'G':
                    base = 'G'; break;
                case 't': case 'T':
                    base = 'T'; break;
                case '-': 
                    base = '-'; break;
                default: // 'N'
                    base = 'N'; break;
                }
                if (base != '\0') {
                    if (length >= size_buffer) {
                        cerr << "too much nucleotides > " << size_buffer << endl;
                        //buffer[size_buffer-1] = '\0';
                    } else {
                        buffer[length++] = base;
                    }
                }
            }
        }
    }
    if (seq != NULL) {
        seq->set_sequence(length, buffer);
        chromosomes.push_back(seq);
    }
    delete[] buffer;
    fi.close();
    return chromosomes;
}


// namespace {
//     bool compare_first_index(const pair<int,char>& lhs, const pair<int,char>& rhs) {
//         return lhs.first < rhs.first;
//     }
// }

void recfragment::join_sequence(const recfragment* frag) {
    //vector<qseqbase> bases = frag->_mapped;
    for (vector<qseqbase>::const_iterator it = frag->_mapped.begin();
         it != frag->_mapped.end(); it++) {
        _mapped.push_back(*it);
    }
    if ((_position3 > frag->_position5 && _position5 <= frag->_position3)) {
        //sort(_mapped.begin(), _mapped.end(), compare_first_index);
        sort(_mapped.begin(), _mapped.end(), qseqbase::compare_position);//compare_first_index);
        //vector<int> ambiguous;
        for (int i = (int)_mapped.size() - 1; i > 0; i--) {
            if (_mapped[i].position() == _mapped[i-1].position()) {
                if (_mapped[i].nucleotide() != _mapped[i-1].nucleotide()) {
                    //ambiguous.push_back(i);
                    _mapped.erase(_mapped.begin() + i, _mapped.begin() + i + 2);
                }
                i--;
            }
        }
    }
    _position3 = _position3 > frag->_position3 ? _position3 : frag->_position3;
    _position5 = _position5 < frag->_position5 ? _position5 : frag->_position5;
    _flag |= FLAG_JOINED;
}


bool recfragment::compare_fragment_order(const recfragment* lhs, 
                                         const recfragment* rhs) {
    if (lhs->_group != rhs->_group) {
        return lhs->_group < rhs->_group;
    } else {
        return lhs->_position5 < rhs->_position5;
    }
}

const int recfragment::FLAG_JOINED = 0x8000;

void recfragment::bundle_pairs(vector<recfragment*>& fragments) throw (runtime_error) {
    //vector<recfragment*> bundled;
    for (int i = 0; i < (int)fragments.size(); i++) {
        if (fragments[i] == NULL) {
            cerr << "null fragment at " << i << endl;
            return;
        }
    }
    sort(fragments.begin(), fragments.end(), compare_fragment_order);
    int N = fragments.size();
    for (int i = 0; i < N; i++) {
        recfragment* p = fragments[i];
        if (p == NULL || p->is_compound()) continue;
        for (int j = i + 1; j < N; j++) {
            recfragment* q = fragments[j];
            if (q == NULL) continue;
            if (p->_group != q->_group) {
                if (p->_group == q->_group && p->name() == q->name()) {
                    //cout << "JOIN:" << p->to_string() << " + " << q->to_string();
                    p->join_sequence(q);
                    //cout << " => " << p->to_string() << endl;
                    delete q;
                    fragments[j] = NULL;
                    break;
                }
            }
        }
    }

    int tail = 0;
    for (int i = 0; i < N; i++) {
        recfragment* p = fragments[i];
        if (p != NULL) {
            if (i != tail) {
                fragments[tail++] = p;
            }
        }
    }
    fragments.erase(fragments.begin() + tail, fragments.end());
}

string recfragment::to_string() const {
    stringstream ss;
    ss << _group << ":" << chromosome() << ":" << position5() << "-" << position3();
    return ss.str();
}


ostream& operator << (ostream& ost, const recfragment& rhs) {
    ost << rhs.name() << " ; " << rhs.chromosome() << ":" << rhs.orientation() << ":" << rhs.position5() << "-" << rhs.position3() << " ";
    return ost;
}


