#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>

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
}

vector<pair<int,char> > fragment::resolve_map(int position, const string& cigar, const string& sequence) {
    int value = 0;
    int pos = 0;
    const char* ptr = cigar.c_str();
    for (int i = 0; i < (int)cigar.size(); i++) {
        char c = ptr[i];
        if (c >= '0' && c <= '9') {
            if (value == 0) {
                value = atoi(ptr + i);
            }
        } else if (c == 'M') {
            for (int j = 0; j < value; j++) {
                
                pos++;
            }
        } else if (c == 'S') {
        }
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

        

        return 0;
    } catch (exception& e) {
        cerr << e.what() << endl;
        return -1;
    }
}
