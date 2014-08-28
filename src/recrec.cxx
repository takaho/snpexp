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

fragment::fragment(const string& chromosome, int position, int flag, const string& sequence, const string& cigar) {
    _chromosome = chromosome;
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

vector<pair<int,char> > fragment::resolve_map(int position, const string& cigar, const string& sequence) {
    fragment f("chr1", position, 0, sequence, cigar);
    return f.get_positions();
}

void fragment::initialize(const bam1_t* seq) {
    
}

void fragment::initialize(int position, const string& sequence, const string& cigar){
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
            spos += value;
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


int main(int argc, char** argv) {
    try {
        const char* filename = get_argument_string(argc, argv, "i", NULL);
        int coverage = get_argument_integer(argc, argv, "c", 10);
        double heterozygosity = get_argument_float(argc, argv, "H", 0.2);
        bool verbose = has_option(argc, argv, "verbose");

        if (verbose) {
            cerr << "Filename         : " << filename << endl;
            cerr << "Minimum coverage : " << coverage << endl;
            cerr << "Heterozygosity   : " << heterozygosity << endl;
        }

        vector<fragment> fragments;
        bam1_t* read;
        bam_header_t* header;
        bamFile bamfile;

        bamfile = bam_open(filename, "rb");
        header = bam_header_read(bamfile);
        read = bam_init1();

        while (bam_read1(bamfile, read) > 0) {
            string cigar = resolve_cigar(read);
            header;
            fragment f(chromosome, position, flag, sequence, cigar);
            read->core.tid;
        }

        bam_destroy1(read);
        bam_header_destroy(header);
        bam_close(bamfile);

        return 0;
    } catch (exception& e) {
        cerr << e.what() << endl;
        return -1;
    }
}
