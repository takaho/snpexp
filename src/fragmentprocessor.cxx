#include <iostream>
#include <fstream>
#include <vector>

using std::cout;
using std::endl;
using std::vector;
using std::ostream;

#include <fragmentprocessor.hxx>
#include <recrec.hxx>
#include <distsnp.hxx>

using namespace tkbio;

fragment_processor::fragment_processor() {
    _variation_db = NULL;
}
fragment_processor::~fragment_processor() {
}

vector<hetero_locus*> recombination_detector
::scan_heterozygous_loci(const vector<recfragment*>& fragments, 
                         chromosome_seq const* chromosome, int start, int end) const
    throw (exception) {
    set<int> pos;
    return scan_heterozygous_loci(fragments, chromosome, start, end, pos);
}

vector<hetero_locus*> recombination_detector
::scan_heterozygous_loci(const vector<recfragment*>& fragments, 
                         chromosome_seq const* chromosome, int start, int end, 
                         set<int> accepted) const throw (exception) {
    int freq[5];
    vector<hetero_locus*> candidates;
    bool gapped = false;
    for (int pos = start; pos < end; pos++) {
        if (accepted.size() > 0 && accepted.find(pos) == accepted.end()) {
            continue;
        }
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
        if (coverage >= _coverage && refnum > 1 && altnum > 1) {
            int threshold = (int)(coverage * _minimum_minor_ratio + 0.5);
            if (threshold <= altnum && threshold <= refnum) {
                candidates.push_back(new hetero_locus(chromosome->code(), pos, refid, refnum, altid, altnum));
            }
        }
    }
    return candidates;
}

recombination_detector::recombination_detector(int coverage, float hetero_threshold){
    _coverage = coverage;
    _minimum_minor_ratio = hetero_threshold;
}

void recombination_detector::process_fragments(const vector<recfragment*>& fragments,
                                               chromosome_seq const* chromosome,
                                               int start, int end, ostream& ost) throw (exception) {
    if (chromosome != NULL) {
        std::cout << chromosome->name() << ";" << chromosome->length() << "\t" << start << "\t" << end << "\t" << fragments.size() << std::endl;
        return;
    }

    vector<hetero_locus*> loci;
    if (_variation_db == NULL) {
        loci = scan_heterozygous_loci(fragments, chromosome, start, end);
    } else {
        vector<dbsnp_locus const*> snps = _variation_db->get_snps(chromosome->name(), start, end);
        set<int> pos;
        for (int i = 0; i < (int)snps.size(); i++) {
            pos.insert(snps[i]->position());
        }
        if (pos.size() > 1) {
            loci = scan_heterozygous_loci(fragments, chromosome, start, end, pos);
        }
    }
}


