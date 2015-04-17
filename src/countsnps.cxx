#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

#include <tktools.hxx>
#include <distsnp.hxx>
#include <gtf.hxx>

using namespace tktools;
using namespace tktools::io;
using namespace tktools::util;
using namespace tkbio;

namespace {
    void display(dbsnp_file const* dbsnp, int const* counts, ostream& ost, const char* target=NULL) {
        if (target != NULL) {
            for (int i = 0; i < dbsnp->strain_size(); i++) {
                if (dbsnp->get_strain(i) == target) {
                    ost << target << ":" << counts[i];
                }
            }
        } else {
            for (int i = 0; i < dbsnp->strain_size(); i++) {
                if (i == 0) {
                    //
                } else {
                    ost << "\t";
                }
                ost << dbsnp->get_strain(i) << ":" << counts[i];
            }
        }
    }
}
int main(int argc, char** argv) {
    try {
        const char* filename_dbsnp = get_argument_string(argc, argv, "d", "/mnt/smb/tae/snps/mgp.v3.snps.rsIDdbSNPv137.vcf");
        //const char* strain = get_argument_string(argc, argv, 
        dbsnp_file* dbsnp = dbsnp_file::load_dbsnp(filename_dbsnp);
        int step = 100000;
        int num_strains = dbsnp->strain_size();
        int* snpcounts = new int[num_strains];
        int length_x = 166650297;
        for (int i = 0; i < num_strains; i++) {
            snpcounts[i] = 0;
        }
        const char* chromosome = "X";
        for (int pos = 0; pos < length_x; pos += step) {
            int tail = pos + step < length_x ? pos + step: length_x;
            
            vector<dbsnp_locus const*> snps = dbsnp->get_snps(chromosome, pos, tail);
            for (int i = 0; i < (int)snps.size(); i++) {
                for (int j = 0; j < num_strains; j++) {
                    unsigned char genotype = snps[i]->get_genotype(j);
                    if (genotype == (unsigned char)0x22) {
                        snpcounts[j]++;
                    }
                }
            }
            cerr << chromosome << ":1-" << (pos + step) << " ";
            display(dbsnp, snpcounts, cerr, "CASTEiJ");
            cerr << "        \r";
        }
        cerr << endl;
        display(dbsnp, snpcounts, cout);
        return 0;
    } catch (exception& e) {
        cerr << e.what() << endl;
        return -1;
    }
}
    
