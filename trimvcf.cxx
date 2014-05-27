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


namespace {
    bool check_region(int pos, const vector<pair<int,int> >& regions) {
        int left = 0;
        int right = regions.size();
        while (left != right) {
            int center = (left + right) >> 1;
            const pair<int,int>& p = regions[center];
            if (pos < p.first) {
                right = center;
            } else if (p.second < pos) {
                left = center + 1;
            } else {
                return true;
            }
        }
        return false;
    }

    void annotate_rsid(int argc, char** argv) throw (exception) {
        try {
            const char* gtf = get_argument_string(argc, argv, "G", "mm10_genes.gtf");
            const char* diff = get_argument_string(argc, argv, "e", NULL);

            const char* filename = get_argument_string(argc, argv, "i", argv[1]);

            map<string,string> expression;

            if (diff != NULL) {
                ifstream fe(diff);
                while (!fe.eof()) {
                    string line;
                    getline(fe, line);
                    vector<string> items = split_items(line, '\t');
                    if (items.size() > 13 && items[6] == "OK") {
                        if (items[13] == "yes") {
                            expression[items[0]] = items[9];
                        } else {
                            expression[items[0]] = "(" + items[9] + ")";//"0";
                        }
                    }
                }
            }

            // const char* strains = get_argument_string(argc, argv, "c", NULL);
            // const char* vcf = get_argument_string(argc, argv, "V", NULL);

            // vector<int> columns;
            // if (strains != NULL) {
            //     vector<string> cstr = split_items(string(strains), ',');
                
                
            // }

            ifstream fi(filename);
            if (fi.is_open() == false) {
                throw invalid_argument("cannot open input file");
            }
            gtffile* genes = gtffile::load(gtf);
            while (!fi.eof()) {
                string line;
                getline(fi, line);
                vector<string> items = split_items(line, '\t');
                if (items.size() <= 2) continue;
                vector<const gtfgene*> selected = genes->find_genes(items[0], atoi(items[1].c_str()));
                if (selected.size() > 0) {
                    string gstr = "";
                    set<string> touched;
                    for (int i = 0; i < (int)selected.size(); i++) {
                        const gtfgene* g = selected[i];
                        if (touched.find(g->name()) == touched.end()) {
                            if (gstr != "") {
                                gstr += " /// ";
                            }
                            gstr += g->name() + " // " + g->transcript_id();
                            if (expression.size() > 0) {
                                if (expression.find(g->name()) == expression.end()) {
                                    gstr += " // 0";
                                } else {
                                    gstr += " // " + expression[g->name()];
                                }
                            }
                            touched.insert(g->name());
                        }
                    }
                    cout << line << "\t" << gstr << endl;
                } else {
                    cout << line << "\t" << "." << endl;
                }
            }
            fi.close();
            delete genes;
        } catch (exception& e) {
            throw;
        }
    }
}


int main(int argc, char** argv) {
    try {

        if (true) {
            annotate_rsid(argc, argv);
            return 0;
        }

        const char* vcf = get_argument_string(argc, argv, "i", NULL);
        const char* gtf = get_argument_string(argc, argv, "G", NULL);
        const char* output = get_argument_string(argc, argv, "o", NULL);

        bool verbose = has_option(argc, argv, "verbose");
        bool rsonly = has_option(argc, argv, "r");
        if (vcf == NULL || gtf == NULL) {
            throw invalid_argument("VCF and GTF files are required");
        }
        ostream* ost = &cout;
        if (output != NULL) {
            ofstream* fo = new ofstream(output);
            if (fo->is_open() == false) {
                throw invalid_argument("cannot open output file");
            }
            ost = fo;
        }
        size_t filesize = get_file_size(vcf);
        ifstream fi(vcf);
        if (fi.is_open() == false) {
            throw runtime_error("cannot open VCF file");
        }
        gtffile* genes = gtffile::load(gtf);
        size_t num_lines = 0;
        string chromosome;
        vector<pair<int,int> > mask;
        while (!fi.eof()) {
            string line;
            getline(fi, line);
            if (verbose) {
                num_lines++;
                if (++num_lines % 1000 == 0) {
                    size_t cur = fi.tellg();
                    double percent = (double)cur * 100.0 / filesize;
                    cerr << (num_lines / 1000) << " " << setprecision(2) << percent << "         \r";
                }
            }
            bool accepted = false;
            const char* ptr = line.c_str();
            if (ptr[0] == '#') {
                accepted = true;
            } else {
                size_t span = line.size();
                int column = 0;
                string chrm;
                int start = 0;
                for (int i = 0; i < (int)span; i++) {
                    char c = ptr[i];
                    if (c == '\t') {
                        if (column == 0) {
                            chrm = line.substr(0, i);
                            if (chrm != chromosome) {
                                mask = genes->get_exonregion(chrm);
                                chromosome = chrm;
                            }
                        } else if (column == 1) {
                            int pos = atoi(ptr + start);
                            accepted = check_region(pos, mask);
                            if (!accepted) break;
                            if (!rsonly) break;
                        } else if (column == 2) {
                            accepted = ptr[start] != '.';
                            break;
                        }
                        start = i + 1;
                        column ++;
                    }
                }
            }
            if (accepted) {
                *ost << line << "\n";
            }
        }
        if (output != NULL) {
            dynamic_cast<ofstream*>(ost)->close();
            delete ost;
            ost = &cout;
        }
        delete genes;
        return 0;
    } catch (exception& e) {
        cerr << e.what() << endl;
        return -1;
    }
}
