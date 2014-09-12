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
    _quality_threshold = 10;
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
                         const set<int>& accepted) const throw (exception) {
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
                int index = fr->get_base_id(pos, _quality_threshold, num);
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

recombination_detector::recombination_detector(int coverage, float hetero_threshold,
                                               float recom_threshold){
    _coverage = coverage;
    _minimum_minor_ratio = hetero_threshold;
    _minimum_recombination_reads = coverage;
    _minimum_recombination_composition = recom_threshold;
}

bool recombination_detector::check_acceptable_recombination(int counts[4]) const {
    int total = counts[0] + counts[1] + counts[2] + counts[3];
    int thr = (int)(total * _minimum_recombination_composition + 0.5);
    for (int i = 0; i < 4; i++) {
        if (counts[i] < thr) {
            return false;
        }
    }
    return true;
}

void recombination_detector::process_fragments(const vector<recfragment*>& fragments,
                                               chromosome_seq const* chromosome,
                                               int start, int end, ostream& ost) throw (exception) {
    // if (chromosome == NULL) {
    //     std::cout << chromosome->name() << ";" << chromosome->length() << "\t" << start << "\t" << end << "\t" << fragments.size() << std::endl;
    //     return;
    // }
    //cout << chromosome->name() << ":" << start << "-" << end << endl;
    vector<hetero_locus*> loci;
    set<int> pos;
    if (_variation_db != NULL) {
        vector<dbsnp_locus const*> snps = _variation_db->get_snps(chromosome->name(), start, end);
        for (int i = 0; i < (int)snps.size(); i++) {
            pos.insert(snps[i]->position());
        }
        // cout << chromosome->name() << " " << start << "-" << end << " ";
        // cout << pos.size() << "SNPs\n";
    // } else {
    //     cout << "vcf null\n";
    }
    loci = scan_heterozygous_loci(fragments, chromosome, start, end);
    int** patterns = new int*[fragments.size()];
    for (int i = 0; i < (int)fragments.size(); i++) {
        patterns[i] = new int[loci.size()];
        fragments[i]->generate_recombination_pattern(loci, patterns[i]);
    }

    int tolerance = 2;
    int counts[4];
    for (int i = 0; i < (int)loci.size(); i++) {
        int loops = i + tolerance < (int)loci.size() ? tolerance : (int)loci.size() - i;
        for (int l = 1; l < loops; l++) {
            counts[0] = counts[1] = counts[2] = counts[3] = 0;
            //for (int j = 0; j < 4; j++) counts[j] = 0;
            for (int j = 0; j < (int)fragments.size(); j++) {
                int p0 = patterns[j][i];
                int p1 = patterns[j][i + l];
                if (p0 == 0 || p1 == 0) continue;
                int gn = 0;
                if (p0 == 1) {
                    gn = (p1 == 1) ? 3 : 2;
                } else {
                    gn = (p1 == 1) ? 1 : 0;
                }
                counts[gn]++;
            }
            if (counts[0] == 0 || counts[1] == 0 
                || counts[2] == 0 || counts[3] == 0
                || (counts[0] + counts[1] + counts[2] + counts[3] < _minimum_recombination_reads) ) {
                break;
            }
            if (check_acceptable_recombination(counts)) {
                ost << chromosome->name() << ":" << loci[i]->position() << "-" << loci[i+l]->position() << "\t" << loci[i]->ref() << ";" << loci[i]->alt() << "\t" << loci[i+l]->ref() << ";" << loci[i+l]->alt() << "\t";
                for (int i = 0; i < 4; i++) {
                    ost << "\t" << counts[i];
                }
                ost << "\n";
                i += (l - 1);
                break;
            }
        }
    }
    for (int i = 0; i < (int)loci.size(); i++) {
        delete loci[i];
    }
    for (int i = 0; i < (int)fragments.size(); i++) {
        delete[] patterns[i];
    }
    delete[] patterns;
}


