#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <sstream>

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

recfragment::recfragment(const string& chromosome, int position, int flag, const string& sequence, const string& cigar) {
    _chromosome = convert_chromosome_to_code(chromosome.c_str());
    _cigar = cigar;
    _position = position;
    _flag = flag;
    _sequence = sequence;
    _mapped = resolve_map(position, cigar, sequence);
    _position5 = 2000000000;
    _position3 = 0;
    initialize(position, cigar, sequence);

    for (vector<pair<int,char> >::const_iterator it = _mapped.begin(); it != _mapped.end(); it++) {
        int pos = it->first;
        if (pos < _position5) _position5 = pos;
        if (pos > _position3) _position3 = pos;
    }
}

vector<pair<int,char> > recfragment::resolve_map(int position, const string& cigar, const string& sequence) {
    recfragment f("chr1", position, 0, sequence, cigar);
    return f.get_positions();
}

void recfragment::initialize(const bam1_t* seq) {
    
}

void recfragment::initialize(int position, const string& sequence, const string& cigar){
    int value = 0;
    int spos = position;
    int fpos = 0;
    //vector<pair<int,char> > nucleotides;
    const char* ptr = cigar.c_str();
    const char* seq = sequence.c_str();
    int max_matches = 0;
    int left = 2000000000;
    int right = 0;
    for (int i = 0; i < (int)cigar.size(); i++) {
        char c = ptr[i];
        if (c >= '0' && c <= '9') {
            if (value == 0) {
                value = atoi(ptr + i);
            }
        } else if (c == 'M') {
            if (max_matches < value) max_matches = value;
            if (left > spos) left = spos;
            for (int j = 0; j < value; j++) {
                _mapped.push_back(make_pair(spos, seq[spos]));
                spos++;
                fpos++;
            }
            right = spos;
        } else if (c == 'S' || c == 'H') { // soft clip & hard clip
            fpos += value;
            spos += value;
        } else if (c == 'I') {
            fpos += value;
        } else if (c == 'D') {
            for (int j = 0; j < value; j++) {
                _mapped.push_back(make_pair(spos, '-'));
                spos++;
            }
        } else {
            throw runtime_error(string("unknown cigar character :") + c); 
        }
    }
    _max_match_span = max_matches;
    _position5 = left;
    _position3 = right;
}

namespace {
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
    class hetero_locus {
    private:
        int position;
        int index1;
        int index2;
    public:
        hetero_locus(int position, int index1, int index2) {
            this->position = position;
            this->index1 = index1;
            this->index2 = index2;
        }
    };
}

void recfragment::detect_recombination(const vector<recfragment*>& fragments, int start, int end) {
    int length = end - start;
    int** counts = new int*[4];
    for (int i = 0; i < 5; i++) {
        counts[i] = new int[length];
        for (int j = 0; j < length; j++) {
            counts[i][j] = 0;
        }
    }

    // reconstruct sequence
    for (int i = 0; i < (int)fragments.size(); i++) {
        const vector<pair<int,char> >& nucleotides = fragments[i]->_mapped;
        for (int j = 0; j < (int)nucleotides.size(); j++) {
            int pos = nucleotides[j].first;
            char nuc = nucleotides[j].second;
            if (start <= pos && pos < end) {
                switch (nuc) {
                case 'A':
                    index = 0; break;
                case 'C':
                    index = 1; break;
                case 'G':
                    index = 2; break;
                case 'T':
                    index = 3; break;
                case '-':
                    index = 4; break;
                default:
                    continue;
                }
                counts[i][pos - start]++;
            }
        }
    }
    // detect heterozygous position
    int coverage = 10;
    double ratio = 0.2;
    double oratio = 0.9;
    vector<hetero_locus> hetero;
    for (int i = 0; i < length; i++) {
        int total = 0;
        int types = 0;
        int max_index = -1;
        int second_index = -1;
        for (int j = 0; j < 5; j++) {
            int cnt = counts[j][i];
            if (cnt > 0) {
                total += cnt;
                types ++;
                if (cnt > max_cnt) {
                    max_cnt = cnt;
                    max_index = j;
                }
            }
        }
        if (total >= coverage) {
            int hetero_limit = coverage * ratio;
            int occupy_limit = coverage * oratio;
            int second_most = 0;
            int second_index = -1;
            for (int j = 0; j < 5; j++) {
                int cnt = counts[j][i];
                if (j != max_index && second_most < cnt) {
                    second_most = cnt;
                    second_index = j;
                }
            }
            if (second_idnex >= 0 && hetero_limit < second_most && coverage * occipy_limit < max_value + second_most) {
                hetero.push_back(hetero_locus(i, max_index, second_index));
            }
        }
    }

    for (int i = 0; i < (int)hetero.size(); i++) {
        int posi = hetero[i].position;
        for (int j = i + 1; j < (int)hetero.size(); j++) {
            int posj = hetero[j].position;
            
        }
    }

    // for (int i = 0; i < length; i++) {
    //     for (int j = 0; j 
    //              }
    //     char* sequence = new char[length];
    // }
}
    
int main(int argc, char** argv) {
    try {
        const char* filename = get_argument_string(argc, argv, "i", NULL);
        int coverage = get_argument_integer(argc, argv, "c", 10);
        double heterozygosity = get_argument_float(argc, argv, "H", 0.2);
        bool verbose = has_option(argc, argv, "verbose");
        int window_size = get_argument_integer(argc, argv, "w", 5000);
        int window_step = get_argument_integer(argc, argv, "s", 500);
        if (verbose) {
            cerr << "Filename         : " << filename << endl;
            cerr << "Minimum coverage : " << coverage << endl;
            cerr << "Heterozygosity   : " << heterozygosity << endl;
        }

        vector<recfragment*> fragments;
        map<string,const recfragment*> pairs;
        bam1_t* read;
        bam_header_t* header;
        bamFile bamfile;

        bamfile = bam_open(filename, "rb");
        header = bam_header_read(bamfile);
        read = bam_init1();

        int chromosome = -1;
        int region_start = 0;
        int region_end = region_start + window_size;
        while (bam_read1(bamfile, read) > 0) {
            recfragment* frag = new recfragment(read);
            int chrm = read->core.tid;
            bool processing = false;
            int next_position;
            if (chrm != chromosome) {
                if (chromosome >= 0) {
                    // process
                }
                chromosome = chrm;
                next_position = 0;
                processing = true;
            } else if (read->core.pos > region_end) {
                // process
                next_position = region_start + window_step;
                processing = true;
            }
            if (processing) {
                // remove outsiders
                
                chromosome = chrm;
                region_start = next_position;
                for (int i = 0; i < (int)fragments.size(); i++) {
                    recfragment* f = fragments[i];
                    if (f->chromosome() != chromosome || f->position3() < region_start) {
                        delete f;
                        fragments[i] = NULL;
                    }
                }
                int index = 0;
                for (int i = 0; i < (int)fragments.size(); i++) {
                    if (fragments[i] != NULL) {
                        fragments[index] = fragments[i];
                        fragments[i] = NULL;
                        index++;
                    }
                }
                fragments.erase(fragments.begin() + index, fragments.end());
                if (fragments.size() == 0) {
                    region_start = (frag->position5() / window_step) * window_step;
                }
                region_end = region_start + window_size;
                region_end = region_start + window_size;
            }
            fragments.push_back(frag);
        }
        if (fragments.size() > coverage) {
            // process
        }
        

//            if (chrm != chromosome
            string cigar = resolve_cigar(read);
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
