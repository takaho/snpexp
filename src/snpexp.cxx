/* 
The MIT Lincense

Copyright (c) 2014 Takaho A. Endo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */


/*
Filename 1: SRR1171556.mm10/CD45_1.bam
Filename 2: SRR1171557.mm10/CD45_2.bam
Filename 3: SRR1171560.mm10/ES_1.bam
Filename 4: SRR1171561.mm10/ES_2.bam
Filename 5: SRR1171580.mm10/STAP_1.bam
Filename 6: SRR1171581.mm10/STAP_2.bam
Filename 7: SRR1171585.mm10/STAP-SC_1.bam
Filename 8: SRR1171586.mm10/STAP-SC_2.bam
Filename 9: SRR1171565.mm10/FI-SC_1.bam
Filename 10: SRR1171566.mm10/FI-SC_2.bam
Filename 11: SRR1171590.mm10/TS_1.bam
Filename 12: SRR1171591.mm10/TS_2.bam
Filename 13: SRR1171558.mm10/EpiSC_1.bam
Filename 14: SRR1171559.mm10/EpiSC_2.bam


 */

#include <memory>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <algorithm>

using namespace std;

#ifdef HAVE_BAM_H
#include <bam.h>
#else
#error this program requires bam.h, bgzf.h and libbam
#endif

#include <tktools.hxx>
#include <gtf.hxx>
#include <snpexp.hxx>

using namespace tktools;
using namespace tktools::io;
using namespace tktools::util;
using namespace tktools::bio;

using namespace tkbio;

void base_frequency::initialize(int length) {
    _a = new unsigned short[length];
    _c = new unsigned short[length];
    _g = new unsigned short[length];
    _t = new unsigned short[length];
    for (int i = 0; i < length; i++) {
        _a[i] = _c[i] = _g[i] = _t[i] = 0;
    }
    _reference = false;
}

base_frequency::~base_frequency() {
    delete[] _a;
    delete[] _c;
    delete[] _g;
    delete[] _t;
}

void base_frequency::clear_reference() {
    if (!_reference) return;
    for (int i = 0; i < _length; i++) {
        if (_a[i] < 0) {
            _a[i] = - _a[i] - 1;
        } else if (_c[i] < 0) {
            _c[i] = - _c[i] - 1;
        } else if (_g[i] < 0) {
            _g[i] = - _g[i] - 1;
        } else if (_t[i] < 0) {
            _t[i] = - _t[i] - 1;
        }
    }
    _reference = false;
}
    
void base_frequency::set_sequence(int length, const char* sequence) {
    throw runtime_error("this method is obsolete");
    // if (_reference) {
    //     return;
    // }
    // length = length < _length ? length : _length;
    // for (int i = 0; i < length; i++) {
    //     char base = sequence[i];
    //     ushort* target;
    //     if (base == 'A') {
    //         target = _a;
    //     } else if (base == 'C') {
    //         target = _c;
    //     } else if (base == 'G') {
    //         target = _g;
    //     } else if (base == 'T') {
    //         target = _t;
    //     } else {
    //         continue;
    //     }
    //     target[i] = -1 - target[i];
    // }
    // _reference = true;
}

bool base_frequency::get_basecount(int pos, int& A, int& C, int& G, int& T) const {
    if (pos >= 0 && pos < _length) {
        A = _a[pos];
        C = _c[pos];
        G = _g[pos];
        T = _t[pos];
        return true;
    } else {
        A = C = G = T = 0;
        return false;
    }
}

//#define DEBUG

void base_frequency::add(int pos, const bam1_t* read) {
    clear_reference();
    // check cigar
    int len = read->core.n_cigar;
#ifdef DEBUG
    if (len < 4) return;
#endif
    int seqpos = 0;
    const uint32_t* cigar = bam1_cigar(read);
    const uint8_t* sequence = bam1_seq(read);
    int offset = 4;
    for (int i = 0; i < len; i++) {
        int op = bam_cigar_op(cigar[i]);
        int slen = bam_cigar_oplen(cigar[i]);
#ifdef DEBUG
        cout << pos << "\t";
        cout << op << "/" << slen << " ";
#endif
        if (slen + pos >= _length) {
            slen = _length - pos;
        }
        //char c = '\0';
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int j = 0; j < slen; j++) {
                uint8_t base = (sequence[(j+seqpos) >> 1] >> offset) & 15;
                offset = 4 - offset;
                switch (base) {
                case 15: break; // N
#ifdef DEBUG
                case 1: _a[pos + j] ++; cout << "A"; break;
                case 2: _c[pos + j] ++; cout << "C"; break;
                case 4: _g[pos + j] ++; cout << "G"; break;
                case 8: _t[pos + j] ++; cout << "T"; break;
#else
                case 1: _a[pos + j] ++; ; break;
                case 2: _c[pos + j] ++; ; break;
                case 4: _g[pos + j] ++; ; break;
                case 8: _t[pos + j] ++; ; break;
#endif
                }
            }
#ifdef DEBUG
            cout << " ";
#endif
            pos += slen;
            seqpos += slen;
        } else if (op == BAM_CINS || op == BAM_CREF_SKIP) {
#ifdef DEBUG
            cout << "=" << len << ">";
#endif
            pos += slen;
        } else if (op == BAM_CDEL) {
            break;
        } else {
            cerr << "unknown cigar character " << op << endl;
        }
    }
#ifdef DEBUG
    cout << endl;
#endif
}

void base_frequency::display_difference(base_frequency* b1, base_frequency* b2, int minimum_bases, double difference) {
    int length = b1->length() < b2->length() ? b1->length() : b2->length();
    static const char* NUCLEOTIDES = "ACGT";
    for (int i = 0; i < length; i++) {
        int num1[4], num2[4];
        num1[0] = b1->_a[i];
        num1[1] = b1->_c[i];
        num1[2] = b1->_g[i];
        num1[3] = b1->_t[i];

        num2[0] = b2->_a[i];
        num2[1] = b2->_c[i];
        num2[2] = b2->_g[i];
        num2[3] = b2->_t[i];
        int s1 = 0;
        int s2 = 0;
        int max_index = 0;
        int maximum = 0;
        for (int j = 0; j < 4; i++) {
            int n1 = num1[j];
            int n2 = num2[j];
            if (n1 + n2 > maximum) {
                max_index = j;
                maximum = n1 + n2;
            }
            s1 += n1;
            s2 += n2;
        }
        if (s1 >= minimum_bases && s2 >= minimum_bases) {
            int max1 = num1[max_index];
            int max2 = num2[max_index];
            double r1 = (double)max1 / s1;
            double r2 = (double)max2 / s2;
            double diff = fabs(r1 - r2);
            if (diff > difference) {
                double pvalue = tktools::stat::get_pvalue_exacttest(max1, s1 - max1, max2, s2 - max2);
                if (pvalue < 0.0001) {
                    cout << (i + 1);
                    for (int j = 0; j < 4; j++) {
                        cout << "\t" << NUCLEOTIDES[j] << ":" << num1[j] << "/" << num2[j];
                    }
                    cout << "\t" << max1 << "," << s1 << "," << max2 << "," << s2 << "\t" << setprecision(3) << diff << "\t" << ios::scientific << pvalue << ios::fixed << endl;
                }
            }
        }
    }
}

void base_frequency::clear_data(const vector<pair<int, int> >& mask) throw (exception) {
    int p5 = 0;
    for (vector<pair<int, int> >::const_iterator it = mask.begin(); it != mask.end(); it++) {
        int p3 = it->first;
        for (int i = p5; i < p3; i++) {
            _a[i] = _c[i] = _g[i] = _t[i] = 0;
        }
        p5 = it->second + 1;
    }
    
}

snp_allele::snp_allele(const string& name, const string& chromosome, int location, const string& reference, const string& alternative) {
    _name = name;
    _chromosome = convert_chromosome_to_code(chromosome.c_str());
    _position = location;
    _reference = reference;
    _alternative = alternative;
}

snp_allele::snp_allele(const string& name, const string& chromosome, int location) {
    _name = name;
    _chromosome = convert_chromosome_to_code(chromosome.c_str());
    _position = location;
    _reference = ".";
    _alternative = ".";
}

bool snp_allele::compare_position(const snp_allele* lhs, const snp_allele* rhs) {
    if (lhs->_chromosome != rhs->_chromosome) {
        return lhs->_chromosome < rhs->_chromosome;
    } else {
        return lhs->_position < rhs->_position;
    }
}

string snp_allele::chromosome() const {
    return convert_code_to_chromosome(_chromosome);
}

void snp_allele::set_variation(const string& variation) {
    _variation = variation;
}

namespace {
    string resolve_snp(const string& pattern, const string& reference, const string& alternative) {
        if (strncmp(pattern.c_str(), "0/0", 3) == 0) {
            return reference + reference;
        } else if (strncmp(pattern.c_str(), "1/1", 3) == 0) {
            return alternative + alternative;
        } else if (strncmp(pattern.c_str(), "0/1", 3) == 0 || strncmp(pattern.c_str(), "1/0", 3) == 0) {
            return reference + alternative;
        } else {
            return "**";
        }
    }
}

vector<snp_allele*> snp_allele::load_vcf(const char* filename, const vector<string>& strains, const char* filename_gtf, bool verbose) throw (exception) {
    ifstream fi(filename);
    if (fi.is_open() == false) {
        throw invalid_argument("no vcf file");
    }
    vector<snp_allele*> snps;
    vector<int> strain_columns;
    gtffile* gtf = NULL;
    if (filename_gtf != NULL) {
        gtf = gtffile::load_gtf(filename_gtf);
    }

    size_t num_snps = 0;
    size_t num_snps_accepted = 0;
    string prev_chrm;
    size_t filesize = 0;
    if (verbose) {
        filesize = get_file_size(filename);
    }

    if (verbose) {
        cerr << "#Chr   SNPs (accepted)   SNPs (all)\n";
    }
    while (fi.eof() == false) {
        string line;
        getline(fi, line);
        vector<string> items = split_items(line, '\t');
        if (items.size() > 4 && line.c_str()[0] != '#') {
            const string& chromosome = items[0];

            if (verbose) {
                if (prev_chrm != chromosome) {
                    if (num_snps > 0) {
                        cerr << " " << setw(3) << prev_chrm << " " << setw(8) << num_snps_accepted << " " << setw(8) << num_snps << "     " << endl;
                        num_snps = 0;
                        num_snps_accepted = 0;
                    }
                    prev_chrm = chromosome;
                }
                num_snps++;
                if (num_snps % 1000 == 0) {
                    double percentage = (double)fi.tellg() * 100.0 / filesize;
                    cerr << " " << setprecision(3) << percentage << "%  " << num_snps << "       \r";
                }
            }

            int location = atoi(items[1].c_str());
            if (gtf != NULL && gtf->contains_in_exon(chromosome, location) == false) {

                //gtf->contains_in_exon_debug(chromosome, location);
                continue;
            }
            num_snps_accepted++;

            snp_allele* snp = new snp_allele(items[2], chromosome, location, items[3], items[4]);
            if (strain_columns.size() > 0) {
                string info = "";
                for (int i = 0; i < (int)strain_columns.size(); i++) {
                    int col = strain_columns[i];
                    if (col < 0) {
                        info += "-";
                    } else {
                        info += resolve_snp(items[col], items[3], items[4]);
                    }
                }
                snp->set_variation(info);
            }
            snps.push_back(snp);
        } else if (line.find("#CHROM") == 0) {
            for (int j = 0; j < (int)strains.size(); j++) {
                int column = -1;
                for (int i = 0; i < (int)items.size(); i++) {
                    if (items[i] == strains[j]) {
                        column = i;
                        break;
                        //strain_columns.push_back(
                    }
                }
                strain_columns.push_back(column);
            }
        }
    }
    if (verbose) {
        if (num_snps > 0) {
            cerr << " " << setw(3) << prev_chrm << " " << setw(8) << num_snps_accepted << " " << setw(8) << num_snps << "     " << endl;
        }
    }
    sort(snps.begin(), snps.end(), compare_position);
    if (gtf != NULL) {
        delete gtf;
    }
    return snps;
}


namespace {

    void show_help() {
        cerr << "snpexp [options] bam1 bam2 ...\n";
        //cerr << " -1 <bam file> : 1st aligned reads\n";
        //cerr << " -2 <bam file> : 2nd aligned reads\n";
        cerr << " -V <filename> : variation file in VCF format\n";
        cerr << " -G <fileanme> : gene position in GTF format\n";
        //cerr << " -g <fasta>    : genomic sequence\n";
        cerr << " -o <filename> : output filename (default:stdout)\n";
        cerr << " -m <number>   : minimum base count to show SNPs\n";
        cerr << " -s <strain1,strain2,...> : straind such as C57BL6NJ, 129S1\n";
        cerr << " -verbose      : verbose mode\n";
        cerr << " -I            : use intergenic positions\n";
        cerr << " -h            : show this message\n";
    }

    pair<int,int> get_snps(const vector<snp_allele*>& snps, const string& chrname) {
        int left = 0;
        int right = (int)snps.size();
        int ccode = convert_chromosome_to_code(chrname.c_str());
        int index = -1;
        //cerr << "searching " << chrname << " " << ccode << endl;
        //int probe;
        int index_start, index_end;
        while (left < right) {
            int center = (left + right) >> 1;
            const snp_allele* snp = snps[center];
            if (snp->chromosome_code() < ccode) {
                left = center + 1;
            } else if (snp->chromosome_code() > ccode) {
                right = center;
            } else {
                index = center;
                break;
            }
        }
        if (index < 0) {
            return make_pair(0,0);
        }
        //cerr << "pivot " << index << endl;

        //int left_cache = left;
        //int right_cache = right;
        //probe = index;
        left = 0;
        right = index + 1;
        index_start = index_end = index;
        while (left < right) {
            int center = (left + right) >> 1;
            const snp_allele* snp = snps[center];
            if (snp->chromosome_code() < ccode) {
                left = center + 1;
            } else if (snp->chromosome_code() > ccode) {
                right = center;
            } else {
                if (center == 0) {
                    index_start = 0;
                    break;
                }
                const snp_allele* prev = snps[center - 1];
                //cerr << center << ":" << prev->chromosome_code() << endl;
                if (prev->chromosome_code() != ccode) {
                    index_start = center;
                    break;
                } else {
                    right = center;
                }
            }            
        }
        left = index;
        right = snps.size();
        for (;;) {
            int center = (left + right) >> 1;
            const snp_allele* snp = snps[center];
            if (snp->chromosome_code() < ccode) {
                left = center + 1;
            } else if (snp->chromosome_code() > ccode) {
                right = center;
            } else {
                if (center >= (int)snps.size() - 1) {
                    index_end = (int)snps.size();
                    break;
                }
                const snp_allele* post = snps[center + 1];
                //cerr << center << ":" << post->chromosome_code() << endl;
                if (post->chromosome_code() != ccode) {
                    index_end = center;
                    break;
                } else {
                    left = center + 1;
                }
            }            
        }
        //cerr << chrname << "=>" << index_start << "-" << index_end << endl;
        return make_pair(index_start, index_end);
    }
}

namespace {
    void display_snps(const string& chromosome, 
                      const vector<base_frequency*>& bfreqs,
                      const vector<snp_allele*>& snps, 
                      int minimum_basecount=10,
                      const gtffile* gtf=NULL,
                      ostream* ost=NULL) {
        //cerr << "get snps " << chromosome << endl;
        pair<int, int> chrsec = get_snps(snps, chromosome);
        int num_bams = bfreqs.size();
        if (ost == NULL) {ost = &cout; }
        for (int i = chrsec.first; i < chrsec.second; i++) {
            const snp_allele* snp = snps[i];
            //cerr << i << "   \r";
            int pos = snp->position();
            bool accepted = false;
            for (int j = 0; j < num_bams; j++) {
                const base_frequency* b = bfreqs[j];
                if (b != NULL) {
                    int A, C, G, T;
                    b->get_basecount(pos - 1, A, C, G, T);
                    if (A + C + G + T >= minimum_basecount) {
                        accepted = true;
                        break;
                    }
                }
            }
            if (accepted) {
                *ost << snp->chromosome() << "\t" << snp->position() << "\t" << snp->id() << "\t" << snp->reference() << "\t" << snp->alternative();
                const char* variation = snp->variation();
                if (strlen(variation) > 0) {
                    *ost << "\t" << variation;
                }
                for (int j = 0; j < num_bams; j++) {
                    const base_frequency* b = bfreqs[j];
                    int A, C, G, T;
                    b->get_basecount(pos - 1, A, C, G, T);
                    *ost << "\t" << A << "," << C << "," << G << "," << T;
                    //*ost << "\tA:" << A << ",C:" << C << ",G:" << G << ",T:" << T;
                }

                vector<const gtfgene*> genes;
                if (gtf != NULL) {
                    genes = gtf->find_genes(snp->chromosome(), snp->position(), snp->position());
                }
                *ost << "\t";
                if (genes.size() == 0) {
                    *ost << ".";
                } else {
                    for (int j = 0; j < (int)genes.size(); j++) {
                        if (j > 0) {
                            *ost << " /// ";
                        }
                        *ost << genes[j]->transcript_id() << " // " << genes[j]->name();
                    }
                }
                *ost << "\n";
                *ost << flush;
            }
        }
    }
}


int main(int argc, char** argv) {
    try {
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        // genomic sequence
        //const char* filename_sequence = get_argument_string(argc, argv, "g", NULL);
        //const char* filename_sequence = NULL;//get_argument_string(argc, argv, "g", NULL);
        const char* filename_gtf = get_argument_string(argc, argv, "G", NULL);
        const char* filename_vcf = get_argument_string(argc, argv, "V", NULL);
        const char* strains = get_argument_string(argc, argv, "s", NULL);
        const char* selected_chromosomes = get_argument_string(argc, argv, "C", NULL);
        int minimum_basecount = get_argument_integer(argc, argv, "m", 0);
        bool use_intergenic = has_option(argc, argv, "-I");

        vector<int> annotation_columns;
        vector<string> filenames_bam;
        for (int i = 1; i < argc; i++) {
            string fn = argv[i];
            if (fn.rfind(".bam") == fn.size() - 4 && file_exists(fn.c_str())) {
                filenames_bam.push_back(fn);
            }
        }
        int num_bams = (int)filenames_bam.size();
        bool verbose = has_option(argc, argv, "verbose");

        if (argc <= 1 || has_option(argc, argv, "-h")) {
            show_help();
            return 0;
        }
        if (num_bams <= 0) {
            throw invalid_argument("at least 1 bam file required");
        }

        if (verbose) {
            for (int i = 0; i < (int)filenames_bam.size(); i++) {
                cerr << "Filename " << (i + 1) << ": " << filenames_bam[i] << endl;
            }
            //cerr << "Filename 2: " << filename_bam2 << endl;
            //cerr << "Genome    : " << (filename_sequence == NULL ? "not_used" : filename_sequence) << endl;
            cerr << "Output    : " << (filename_output == NULL ? "stdout" : filename_output) << endl;
            if (filename_gtf != NULL) {
                cerr << "GTF file  : " << filename_gtf << endl;
            }
            if (filename_vcf != NULL) {
                cerr << "VCF        : " << (filename_vcf == NULL ? "not_used" : filename_vcf) << endl;
                cerr << "Minimum    : " << minimum_basecount << endl;
                cerr << "Strains    : " << (strains == NULL ? "not_used" : strains) << endl;
            }
            if (selected_chromosomes != NULL) {
                cerr << "Chromosomes: " << selected_chromosomes << endl;
            }
        }

        bamFile* bamfiles = new bamFile[num_bams];
        bam1_t** reads = new bam1_t*[num_bams];
        bam_header_t** headers = new bam_header_t*[num_bams];
        ostream* ost = &cout;
        //fasta_sequence* chromosome = NULL;
        gtffile* gtf = NULL;
        vector<snp_allele*> snps;

        // setup bam structs
        for (int i = 0; i < num_bams; i++) {
            bamfiles[i] = bam_open(filenames_bam[i].c_str(), "rb");
            headers[i] = bam_header_read(bamfiles[i]);
            reads[i] = bam_init1();
        }

        if (filename_output != NULL) {
            ofstream* fout = new ofstream(filename_output);
            if (fout->is_open() == false) {
                throw runtime_error(string("cannot open ") + filename_output);
            }
            ost = fout;
        }

        // Gene file
        if (filename_gtf != NULL) {
            if (verbose) {
                cerr << "loading GTF\n";
            }
            gtf = gtffile::load_gtf(filename_gtf);
        }

        // VCF file
        vector<string> strain_labels;
        if (filename_vcf != NULL) {
            if (verbose) { 
                cerr << "loading " << filename_vcf << endl;
            }
            if (strains != NULL) {
                strain_labels = split_items(strains, ',');
            }
            snps = snp_allele::load_vcf(filename_vcf, strain_labels, use_intergenic ? NULL : filename_gtf, verbose);
            if (verbose) {
                cerr << snps.size() << " alleles loaded\n";
            }
        }

        for (int i = 0; i < num_bams; i++) {
            bamfiles[i] = bam_open(filenames_bam[i].c_str(), "rb");
            reads[i] = bam_init1();
            headers[i] = bam_header_read(bamfiles[i]);
        }


        int current_chromosome = -1;
        vector<base_frequency*> bfreqs;
        vector<string> bam_labels;
        set<string> accepted_chromosomes;
        if (selected_chromosomes != NULL) {
            vector<string> chrms = split_items(string(selected_chromosomes), ',');
            for (int i = 0; i < (int)chrms.size(); i++) {
                accepted_chromosomes.insert(chrms[i]);
            }
        }
        //accepted_chromosomes.insert("chrX");
        //accepted_chromosomes.insert("chr9");

        int* chromosomes = new int[num_bams];
        for (int i = 0; i < num_bams; i++) {
            //mode |= (1 << i);
            bfreqs.push_back(NULL);
            chromosomes[i] = -1;
        }

        // print header
        *ost << "#Chr\tPosition\tID_REF\tRef\tAlt";
        if (strains != NULL) *ost << "\t" << strains;

        for (int i = 0; i < num_bams; i++) {
            string fn = filenames_bam[i];
            size_t pos = fn.rfind(tktools::io::file_separator());
            if (pos == string::npos) {
                pos = 0;
            } else {
                pos ++;
            }
            size_t dot = fn.rfind('.');
            if (dot == string::npos || dot < pos) {
                dot = fn.size();
            }
            bam_labels.push_back(fn.substr(pos, dot - pos));
            *ost << "\t" << fn.substr(pos, dot - pos);
        }

        if (filename_gtf != NULL) {
            *ost << "\tGenes";
        }
        *ost << "\n";
        // enuemrate SNPs
        for (;;) {
            bool flag_terminate = false;
            int max_chromosome = chromosomes[0];
            if (current_chromosome >= 0) {
                for (int i = 0; i < num_bams; i++) {
                    if (reads[i] != NULL) {
                        chromosomes[i] = reads[i]->core.tid;
                    } else {
                        chromosomes[i] = -1;
                        flag_terminate = true;
                    }
                }
            } else { // first chromosome
                for (int i = 0; i < num_bams; i++) {
                    if (bam_read1(bamfiles[i], reads[i]) > 0) {
                        chromosomes[i] = reads[i]->core.tid;
                    } else {
                        flag_terminate = true;
                        chromosomes[i] = -1;
                    }
                }
            }
            max_chromosome = chromosomes[0];
            for (int i = 0; i < num_bams; i++) {
                if (max_chromosome < chromosomes[i]) {
                    max_chromosome = chromosomes[i];
                }
            }
            current_chromosome = max_chromosome;
            // skip unaligned chromosome
            for (int i = 0; i < num_bams; i++) {
                if (chromosomes[i] != current_chromosome) {
                    for (;;) {
                        if (bam_read1(bamfiles[i], reads[i]) <= 0) {
                            flag_terminate = true;
                            break;
                        } else if (reads[i]->core.tid == current_chromosome) {
                            chromosomes[i] = current_chromosome;
                            break;
                        }
                    }
                }
            }

            // read
            bool flag_skipping = false;
            for (int i = 0; i < num_bams; i++) {
                if (chromosomes[i] != current_chromosome) {
                    bfreqs[i] = NULL;
                    continue;
                }
                string chromosome_name = headers[i]->target_name[chromosomes[i]];
                if (accepted_chromosomes.size() == 0 || accepted_chromosomes.find(chromosome_name) != accepted_chromosomes.end()) {
                    bfreqs[i] = new base_frequency(chromosomes[i], headers[i]->target_name[chromosomes[i]], headers[i]->target_len[chromosomes[i]]);
                    bfreqs[i]->add(reads[i]->core.pos, reads[i]);
                    if (verbose) {
                        cerr << " reading " << i << ":" << bam_labels[i] << ":"
                             << chromosome_name << "         \r";
                        //<< headers[i]->target_name[chromosomes[i]] << "            \r";
                    }
                } else {
                    if (!flag_skipping) {
                        for (int j = 0; j < i; j++) {
                            delete bfreqs[i];
                            bfreqs[i] = NULL;
                        }
                        if (verbose) {
                            cerr << "skipping chromosome " << chromosome_name << endl;
                        }
                    }
                    flag_skipping = true;
                    bfreqs[i] = NULL;
                }
                if (bfreqs[i] == NULL) { // skip
                    for (;;) {
                        if (bam_read1(bamfiles[i], reads[i]) <= 0) {
                            flag_terminate = true;
                            break;
                        } else if (reads[i]->core.tid != chromosomes[i]) {
                            break;
                        }
                    }
                } else {
                    for (;;) {
                        if (bam_read1(bamfiles[i], reads[i]) <= 0) {
                            flag_terminate = true;
                            break;
                        } else if (reads[i]->core.tid != chromosomes[i]) {
                            break;
                        } else {
                            bfreqs[i]->add(reads[i]->core.pos, reads[i]);
                        }
                    }
                }
            }

            // Clear integenic regions
            if (current_chromosome >= 0 && gtf != NULL && !flag_skipping) {
                vector<gtfgene*> genes = gtf->get_genes();
                vector<pair<int,int> > active_range = gtf->get_exonregion(headers[0]->target_name[current_chromosome]);
                for (int i = 0; i < num_bams; i++) {
                    bfreqs[i]->clear_data(active_range);
                }
            }

            // display SNPs
            if (snps.size() > 0 && !flag_skipping) {
                string chrname = headers[0]->target_name[current_chromosome];
                if (verbose) {
                    cerr << "display variance " << chrname << "    \r";
                }
                display_snps(chrname, bfreqs, snps, minimum_basecount, gtf, ost);
            }

            for (int i = 0; i < num_bams; i++) {
                delete bfreqs[i];
                bfreqs[i] = NULL;
            }

            if (flag_terminate) {
                break;
            }
        }

        // cleanup
        if (gtf != NULL) {
            delete gtf;
        }

        for (int i = 0; i < num_bams; i++) {
            delete bfreqs[i];
        }
        delete[] chromosomes;
        for (vector<snp_allele*>::iterator it = snps.begin(); it != snps.end(); it++) {
            delete *it;
        }
        for (int i = 0; i < num_bams; i++) {
            bam_destroy1(reads[i]);
            bam_header_destroy(headers[i]);
            bam_close(bamfiles[i]);
        }
        delete[] reads;
        delete[] bamfiles;
        delete[] headers;

        if (filename_output != NULL) {
            dynamic_cast<ofstream*>(ost)->close();
            delete ost;
        }

        return 0;
    } catch (exception& e) {
        cerr << e.what() << endl;
        show_help();
        return -1;
    }
}
