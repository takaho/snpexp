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

#include <bam.h>

#include <tktools.hxx>
#include <seq_gene.hxx>
#include <gtf.hxx>
#include <snpexp.hxx>

using namespace tktools;
using namespace tktools::io;
using namespace tktools::util;
using namespace tktools::bio;

using namespace tkbio;

void BaseFrequency::initialize(int length) {
    _a = new short[length];
    _c = new short[length];
    _g = new short[length];
    _t = new short[length];
    for (int i = 0; i < length; i++) {
        _a[i] = _c[i] = _g[i] = _t[i] = 0;
    }
    _reference = false;
}

BaseFrequency::~BaseFrequency() {
    delete[] _a;
    delete[] _c;
    delete[] _g;
    delete[] _t;
}

void BaseFrequency::clear_reference() {
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
    
void BaseFrequency::set_sequence(int length, const char* sequence) {
    if (_reference) {
        return;
    }
    length = length < _length ? length : _length;
    for (int i = 0; i < length; i++) {
        char base = sequence[i];
        short* target;
        if (base == 'A') {
            target = _a;
        } else if (base == 'C') {
            target = _c;
        } else if (base == 'G') {
            target = _g;
        } else if (base == 'T') {
            target = _t;
        } else {
            continue;
        }
        target[i] = -1 - target[i];
    }
    _reference = true;
}

bool BaseFrequency::get_basecount(int pos, int& A, int& C, int& G, int& T) const {
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

void BaseFrequency::add(int pos, const bam1_t* read) {
    clear_reference();
    // check cigar
    int len = read->core.n_cigar;
#ifdef DEBUG
    if (len < 4) return;
#endif
    int seqpos = 0;
    //int relpos = 0;
    const uint32_t* cigar = bam1_cigar(read);
    const uint8_t* sequence = bam1_seq(read);
    //int start = pos;
    //int end = pos;
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
            // pos += len;
            // seqpos += len;
            // if ((seqpos & 1) != 0) {
            //     offset = 4 - offset;
            // }                
            //pos -= len;
        } else {
            cerr << "unknown cigar character " << op << endl;
        }
    }
#ifdef DEBUG
    cout << endl;
#endif
}

void BaseFrequency::display_difference(BaseFrequency* b1, BaseFrequency* b2, int minimum_bases, double difference) {
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
        //int max1 = 0;
        //int max2 = 0;
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
                //double pvalue = 0;
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
        // int s1 = a1 + c1 + g1 + t1;
        // int s1 = a1 + c1 + g1 + t1;
        
        // int a2 = b2->_a[i];
        // int c2 = b2->_c[i];
        // int g2 = b2->_g[i];
        // int t2 = b2->_t[i];
        // int N2 = a2 + c2 + g2 + t2;

        // if (N1 >= minimum_bases && N2 >= minimum_bases) {
            
        // }
        
    }
}

void BaseFrequency::clear_data(const vector<pair<int, int> >& mask) throw (exception) {
    //int index = 0;
    // if (range.size() == 0) {
    //     return;
    // }

    int p5 = 0;
    for (vector<pair<int, int> >::const_iterator it = mask.begin(); it != mask.end(); it++) {
        int p3 = it->first;
        for (int i = p5; i < p3; i++) {
            _a[i] = _c[i] = _g[i] = _t[i] = 0;
        }
        p5 = it->second + 1;
    }
    
    // int p5 = range[0].first - 1;
    // int p3 = range[0].second - 1;
    // for (int i = 0; i < _length; i++) {
    //     if (i < p5) {
    //         cerr << "clear " << i << " " << p3 << "-" << p5 << endl;
    //         for (; i < p5; i++) {
    //             _a[i] = _c[i] = _g[i] = _t[i] = 0;
    //         }
    //         i = p3;
    //         index ++;
    //     } else {
    //         i = p3;
    //         index++;
    //     }
    //     if (index >= (int)range.size()) {
    //         for (int j = i; j < _length; j++) {
    //             _a[i] = _c[i] = _g[i] = _t[i] = 0;
    //         }
    //         break;
    //     } else {
    //         p5 = range[index].first - 1;
    //         p3 = range[index].second - 1;
    //     }
    // }
}

// namespace {
//     bool paircomp(const pair<int, int>& lhs, const pair<int, int>& rhs) {
//         return lhs.first < rhs.first;
//     }

//     // vector<pair<int,int> > merge_range(vector<pair<int,int> >& range) {
//     //     vector<pair<int,int> > available;
//     //     for (int i = (int)range.size() - 1; i > 0; i--) {
//     //         pair<int, int>& r = range[i];
//     //         pair<int, int>& q = range[i - 1];
//     //         if (r.first >= q.second) {
//     //             q.second = r.second < q.second ? q.second : r.second;
//     //             //r.second = r.first;
//     //         } else {
//     //             available.push_back(r);
//     //         }
//     //     }
//     //     available.push_back(range[0]);
//     //     return available;
//     // }

//     vector<pair<int,int> > determine_exon_regions(const gtffile* gtf, string chromosome) {
//         if (chromosome.find("chr") == 0) {
//             chromosome = chromosome.substr(3);
//         }
//         const vector<gtfgene*>& genes = gtf->get_genes();
//         vector<pair<int,int> > regions;
//         for (int i = 0; i < (int)genes.size(); i++) {
//             const gtfgene* gene = genes[i];
//             if (gene->chromosome() == chromosome) {
//                 vector<gtfexon> exons = gene->exons();
//                 for (int j = 0; j < (int)exons.size(); j++) {
//                     regions.push_back(make_pair(exons[j].position5() - 1, exons[j].position3()));
//                 }
//             }
//         }
//         if (regions.size() <= 1) {
//             return regions;
//         }
//         sort(regions.begin(), regions.end(), paircomp);
//         vector<pair<int,int> > merged;
//         int start = regions[0].first;
//         int end = regions[1].second;
//         for (int i = 1; i < (int)regions.size(); i++) {
//             int s = regions[i].first;
//             int e = regions[i].second;
//             if (start < e && s <= end) {
//                 end = e < end ? end : e;
//             } else {
//                 merged.push_back(make_pair(start, end));
//                 start = s;
//                 end = e;
//             }
//         }
//         merged.push_back(make_pair(start, end));
//         return merged;
//     }
// }

SNPAllele::SNPAllele(const string& name, const string& chromosome, int location, const string& reference, const string& alternative) {
    _name = name;
    _chromosome = convert_chromosome_to_code(chromosome.find("chr") == 0 ? chromosome.c_str() + 3 : chromosome.c_str());
    _position = location;
    _reference = reference;
    _alternative = alternative;
}

SNPAllele::SNPAllele(const string& name, const string& chromosome, int location) {
    _name = name;
    _chromosome = convert_chromosome_to_code(chromosome.find("chr") == 0 ? chromosome.c_str() + 3 : chromosome.c_str());
    _position = location;
    _reference = ".";
    _alternative = ".";
}

bool SNPAllele::compare_position(const SNPAllele* lhs, const SNPAllele* rhs) {
    if (lhs->_chromosome != rhs->_chromosome) {
        return lhs->_chromosome < rhs->_chromosome;
    } else {
        return lhs->_position < rhs->_position;
    }
}

string SNPAllele::chromosome() const {
    return convert_code_to_chromosome(_chromosome);
}

void SNPAllele::set_variation(const string& variation) {
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

vector<SNPAllele*> SNPAllele::load_vcf(const char* filename, const vector<string>& strains) throw (exception) {
    ifstream fi(filename);
    if (fi.is_open() == false) {
        throw invalid_argument("no vcf file");
    }
    vector<SNPAllele*> snps;
    vector<int> strain_columns;
    while (fi.eof() == false) {
        string line;
        getline(fi, line);
        vector<string> items = split_items(line, '\t');
        if (items.size() > 4 && line.c_str()[0] != '#') {
            SNPAllele* snp = new SNPAllele(items[2], items[0], atoi(items[1].c_str()), items[3], items[4]);
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
    sort(snps.begin(), snps.end(), compare_position);
    return snps;
}


namespace {

    void show_help() {
        cerr << "snpexp [options] bam1 bam2 ...\n";
        //cerr << " -1 <bam file> : 1st aligned reads\n";
        //cerr << " -2 <bam file> : 2nd aligned reads\n";
        cerr << " -V <vcf>      : VCF file\n";
        cerr << " -g <fasta>    : genomic sequence\n";
        cerr << " -o <filename> : output filename (default:stdout)\n";
        cerr << " -m <number>   : minimum base count to show SNPs\n";
        cerr << " -s <strain1,strain2,...> : straind such as C57BL6NJ, 129S1\n";
    }

    pair<int,int> get_snps(const vector<SNPAllele*>& snps, const string& chrname) {
        int left = 0;
        int right = (int)snps.size();
        int ccode = convert_chromosome_to_code(chrname.find("chr") == 0 ? chrname.c_str() + 3 : chrname.c_str());
        int index = -1;
        //cerr << "searching " << chrname << " " << ccode << endl;
        //int probe;
        int index_start, index_end;
        while (left < right) {
            int center = (left + right) >> 1;
            const SNPAllele* snp = snps[center];
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
            const SNPAllele* snp = snps[center];
            if (snp->chromosome_code() < ccode) {
                left = center + 1;
            } else if (snp->chromosome_code() > ccode) {
                right = center;
            } else {
                if (center == 0) {
                    index_start = 0;
                    break;
                }
                const SNPAllele* prev = snps[center - 1];
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
            const SNPAllele* snp = snps[center];
            if (snp->chromosome_code() < ccode) {
                left = center + 1;
            } else if (snp->chromosome_code() > ccode) {
                right = center;
            } else {
                if (center >= (int)snps.size() - 1) {
                    index_end = (int)snps.size();
                    break;
                }
                const SNPAllele* post = snps[center + 1];
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
                      const vector<BaseFrequency*>& bfreqs,
                      const vector<SNPAllele*>& snps, 
                      int minimum_basecount=10,
                      ostream* ost=NULL) {
        //cerr << "get snps " << chromosome << endl;
        pair<int, int> chrsec = get_snps(snps, chromosome);
        int num_bams = bfreqs.size();
        if (ost == NULL) {ost = &cout; }
        for (int i = chrsec.first; i < chrsec.second; i++) {
            const SNPAllele* snp = snps[i];
            //cerr << i << "   \r";
            int pos = snp->position();
            bool accepted = false;
            for (int j = 0; j < num_bams; j++) {
                const BaseFrequency* b = bfreqs[j];
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
                    const BaseFrequency* b = bfreqs[j];
                    int A, C, G, T;
                    b->get_basecount(pos - 1, A, C, G, T);
                    *ost << "\t" << A << "," << C << "," << G << "," << T;
                    //*ost << "\tA:" << A << ",C:" << C << ",G:" << G << ",T:" << T;
                }
                *ost << "\n";
            }
        }
    }
}


int main(int argc, char** argv) {
    try {
        //const char* filename_bam1 = get_argument_string(argc, argv, "1", NULL);
        //const char* filename_bam2 = get_argument_string(argc, argv, "2", NULL);
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        // genomic sequence
        const char* filename_sequence = get_argument_string(argc, argv, "g", NULL);
        const char* filename_gtf = get_argument_string(argc, argv, "G", NULL);
        const char* filename_vcf = get_argument_string(argc, argv, "V", NULL);
        const char* strains = get_argument_string(argc, argv, "s", NULL);
        int minimum_basecount = get_argument_integer(argc, argv, "m", 0);

        vector<int> annotation_columns;
        vector<string> filenames_bam;
        for (int i = 1; i < argc; i++) {
            string fn = argv[i];
            if (fn.rfind(".bam") == fn.size() - 4 && file_exists(fn.c_str())) {
                filenames_bam.push_back(fn);
            }
        }
        int num_bams = (int)filenames_bam.size();
        if (strains != NULL) {
        }

        //fasta_sequence* seq = load_fasta_sequence(filename_sequence);
        //cout << seq->name() << " " << seq->length() << endl;
        //delete seq;
        //exit(0);

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
            cerr << "Genome    : " << (filename_sequence == NULL ? "not_used" : filename_sequence) << endl;
            cerr << "Output    : " << (filename_output == NULL ? "stdout" : filename_output) << endl;
            if (filename_gtf != NULL) {
                cerr << "GTF file  : " << filename_gtf << endl;
            }
            if (filename_vcf != NULL) {
                cerr << "VCF        : " << (filename_vcf == NULL ? "not_used" : filename_vcf) << endl;
                cerr << "Minimum    : " << minimum_basecount << endl;
                cerr << "Strains    : " << (strains == NULL ? "not_used" : strains) << endl;
            }
        }

        bamFile* bamfiles = new bamFile[num_bams];
        bam1_t** reads = new bam1_t*[num_bams];
        bam_header_t** headers = new bam_header_t*[num_bams];
        ostream* ost = &cout;
        //fasta_sequence* chromosome = NULL;
        gtffile* gtf = NULL;
        vector<SNPAllele*> snps;

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
            gtf = gtffile::load(filename_gtf);
        }

        // VCF file
        if (filename_vcf != NULL) {
            if (verbose) { 
                cerr << "loading " << filename_vcf << endl;
            }
            vector<string> strain_labels;
            if (strains != NULL) {
                strain_labels = split_items(strains, ',');
            }
            snps = SNPAllele::load_vcf(filename_vcf, strain_labels);
            if (verbose) {
                cerr << snps.size() << " alleles loaded\n";
            }
        }

        for (int i = 0; i < num_bams; i++) {
            bamfiles[i] = bam_open(filenames_bam[i].c_str(), "rb");
            reads[i] = bam_init1();
            headers[i] = bam_header_read(bamfiles[i]);
            // for (int j = 0; j < headers[i]->n_targets; j++) {
            //     cerr << i << ":" << j << " " << headers[i]->target_name[j] << endl;
            // }
        }
        //exit(0);

        int mode = 0; // bit 1 : read , 0: use current
        int current_chromosome = -1;
        vector<BaseFrequency*> bfreqs;
        //BaseFrequency** bfreqs = new BaseFrequency*[num_bams];
        int* chromosomes = new int[num_bams];
        for (int i = 0; i < num_bams; i++) {
            mode |= (1 << i);
            bfreqs.push_back(NULL);
            chromosomes[i] = -1;
        }

        for (;;) {
            bool flag_terminate = false;
            //bool flag_break = false;
            //bool flag_different = false;
            int max_chromosome = chromosomes[0];
            if (current_chromosome >= 0) {
                for (int i = 0; i < num_bams; i++) {
                    //bfreqs[i]->add(reads[i]->core.pos, reads[i]);
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
                //cerr << i << " " << chromosomes[i] << " : " << headers[i]->target_name[chromosomes[i]] << endl;
                if (max_chromosome < chromosomes[i]) {
                    max_chromosome = chromosomes[i];
                }
            }
            current_chromosome = max_chromosome;
            // skip unaligned chromosome
            for (int i = 0; i < num_bams; i++) {
                if (chromosomes[i] != current_chromosome) {
                    // if (verbose) {
                    //     cerr << "skip " << i << " : " << headers[i]->target_name[chromosomes[i]] << endl;
                    // }
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
            for (int i = 0; i < num_bams; i++) {
                if (chromosomes[i] != current_chromosome) {
                    bfreqs[i] = NULL;
                    continue;
                }
                bfreqs[i] = new BaseFrequency(chromosomes[i], headers[i]->target_name[chromosomes[i]], headers[i]->target_len[chromosomes[i]]);
                bfreqs[i]->add(reads[i]->core.pos, reads[i]);
                if (verbose) {
                    cerr << " " << i << " reading " << headers[i]->target_name[chromosomes[i]] << "            \r";
                }
                for (;;) {
                    if (bam_read1(bamfiles[i], reads[i]) <= 0) {
                        flag_terminate = true;
                        break;
                    } else if (reads[i]->core.tid != chromosomes[i]) {
                        //cerr << i << " finished " << headers[i]->target_name[chromosomes[i]] << ", " << chromosomes[i] << " / " << reads[i]->core.tid << endl;
                        break;
                    } else {
                        bfreqs[i]->add(reads[i]->core.pos, reads[i]);
                    }
                }
            }
            //cerr << "setup finished \n";

            // Clear integenic regions
            if (current_chromosome >= 0 && gtf != NULL) {
                vector<gtfgene*> genes = gtf->get_genes();
                vector<pair<int,int> > active_range = gtf->get_exonregion(headers[0]->target_name[current_chromosome]);
                for (int i = 0; i < num_bams; i++) {
                    bfreqs[i]->clear_data(active_range);
                }
            }

            // display SNPs
            if (snps.size() > 0) {
                string chrname = headers[0]->target_name[current_chromosome];
                if (verbose) {
                    cerr << "display variance " << chrname << "    \r";
                }
                display_snps(chrname, bfreqs, snps, minimum_basecount, ost);
            }

            // for (int i = 0; i < num_bams; i++) {
            //     delete bfreqs[i];
            //     bfreqs[i] = NULL;
            // }
            // if (chromosome != NULL) {
            //     delete chromosome;
            // }
            
            // if (verbose) {
            //     cerr << "prepare chromosome " << headers[0]->target_name[max_chromosome] << endl;
            // }
            for (int i = 0; i < num_bams; i++) {
                delete bfreqs[i];
                bfreqs[i] = NULL;
                //bfreqs[i] = new BaseFrequency(max_chromosome, headers[i]->target_name[max_chromosome], headers[i]->target_len[max_chromosome]);
            }
            // if (filename_sequence != NULL) {
            //     chromosome = load_fasta_sequence(filename_sequence, headers[0]->target_name[max_chromosome]);
            // }
            //current_chromosome = max_chromosome;
            if (flag_terminate) {
                break;
            }
            // mode = 0;
            // for (int i = 0; i < num_bams; i++) {
            //     if (chromosomes[i] != max_chromosome) {
            //         mode |= (1 << i);
            //     }
            // }
        }
        if (gtf != NULL) {
            delete gtf;
        }

        for (int i = 0; i < num_bams; i++) {
            delete bfreqs[i];
        }
        delete[] chromosomes;
        //delete chromosome;

        //delete bfreq1;
        //delete bfreq2;
        //delete chromosome;

        for (vector<SNPAllele*>::iterator it = snps.begin(); it != snps.end(); it++) {
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

        // bam_destroy1(read1);
        // bam_destroy1(read2);
        // bam_header_destroy(header1);
        // bam_header_destroy(header2);
        // bam_close(bamfile1);
        // bam_close(bamfile2);

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
