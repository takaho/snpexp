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
        const char* filename_output = get_argument_string(argc, argv, "o", NULL);
        const char* target = get_argument_string(argc, argv, "t", "CASTEiJ");
        dbsnp_file* dbsnp = dbsnp_file::load_dbsnp(filename_dbsnp);
        int step = 100000;
        int num_strains = dbsnp->strain_size();
        int* snpcounts = new int[num_strains];
        int length_x = 166650297;
        int target_index = -1;
        for (int i = 0; i < num_strains; i++) {
            snpcounts[i] = 0;
            if (dbsnp->get_strain(i) == target) {
                target_index = i;
            }
        }
        const char* chromosome = "X";
        ofstream fo;
        if (filename_output != NULL) {
            fo.open(filename_output);
            if (fo.is_open() == false) {
                throw exception(string("cannot open ") + filename_output);
            }
        }
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
                if (target_index >= 0 && filename_output != NULL) {
                    unsigned char genotype = snps[i]->get_genotype(target_index);
                    if (genotype == (unsigned char)0x22) {
                        fo << chromosome << "\t" << snps[i]->position();
                    }
                }
            }
            cerr << chromosome << ":1-" << (pos + step) << " ";
            display(dbsnp, snpcounts, cerr, target);//"CASTEiJ");
            cerr << "        \r";
        }
        cerr << endl;
        display(dbsnp, snpcounts, cout);
        if (filename_output != NULL) {
            fo.close();
        }
        return 0;
    } catch (exception& e) {
        cerr << e.what() << endl;
        return -1;
    }
}
    
