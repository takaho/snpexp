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

recfragment::recfragment(int chromosome, bam1_t* seq) {
    _chromosome = chromosome;
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
                uint8_t base = (sequence[fpos >> 1] >> offset) & 15;
                offset = 4 - offset;
                _mapped.push_back(make_pair(spos, _bamnucleotide[base]));
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
                _mapped.push_back(make_pair(spos + j, '-'));
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
            if (center < size - 1 && _mapped[center + 1].first == pos) count++;
            return p.second;
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

void recfragment::generate_recombination_pattern(const vector<hetero_locus>& loci, int* pattern) const {
    int index = 0;
    for (int i = 0; i < (int)loci.size(); i++) {
        int pos = loci[i].position();
        char ref = loci[i].ref();//second;
        char alt = loci[i].alt();
        int value = 0;
        while (index < (int)_mapped.size()) {
            const pair<int,char>& locus = _mapped[index];
            if (locus.first > pos) break;
            if (locus.first == pos) {
                if (locus.second == ref) {
                    value = 1;
                } else if (locus.second == alt) {
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

pair<int,int> recfragment::get_recombination(const vector<hetero_locus>& loci, float diff_ratio) const {
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
                            border = make_pair(loci[index].position(), loci[j].position());
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
string recfragment::get_recombination_pattern(const vector<hetero_locus>& loci) const {
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


void recfragment::bundle_pairs(vector<recfragment*>& fragments) throw (runtime_error) {
    //vector<recfragment*> bundled;
    for (int i = 0; i < (int)fragments.size(); i++) {
        recfragment* p = fragments[i];
        if (p == NULL) continue;
        for (int j = i + 1; j < (int)fragments.size(); j++) {
            recfragment* q = fragments[j];
            if (q == NULL) continue;
            if (p->name() == q->name()) {
                p->join_sequence(q);
                delete q;
                fragments[j] = NULL;
                break;
            }
        }
    }
    int tail = 0;
    for (int i = 0; i < (int)fragments.size(); i++) {
        recfragment* p = fragments[i];
        if (p != NULL) {
            if (i != tail) {
                fragments[tail++] = p;
            }
        }
    }
    fragments.erase(fragments.begin() + tail, fragments.end());
}

ostream& operator << (ostream& ost, const recfragment& rhs) {
    ost << rhs.name() << " ; " << rhs.chromosome() << ":" << rhs.orientation() << ":" << rhs.position5() << "-" << rhs.position3() << " ";
    return ost;
}


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
