#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

using std::cout;
using std::endl;
using std::vector;
using std::ostream;
using std::stringstream;
using std::make_pair;

#include <fragmentprocessor.hxx>
#include <recrec.hxx>
#include <distsnp.hxx>
#include <tktools.hxx>

using tktools::split_items;
using namespace tkbio;

namespace {
    // template<typename T> T min(T a, T b) {
    //     return a < b ? a : b;
    // }
    // template<typename T> T max(T a, T b) {
    //     return a > b ? a : b;
    // }
}

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
    _snp_stretches = 1;
    _gap_tolerance = 0;
}

bool recombination_detector::check_acceptable_recombination(int counts[4]) const {
    int total = counts[0] + counts[1] + counts[2] + counts[3];
    if (total < _minimum_recombination_reads) {
        return false;
    }
    int thr = (int)(total * _minimum_recombination_composition + 0.5);
    for (int i = 0; i < 4; i++) {
        if (counts[i] < thr) {
            return false;
        }
    }
    return true;
}

namespace {
    int get_base_id(const string& pattern) {
        if (pattern == "A") {
            return 0;
        } else if (pattern == "C") {
            return 1;
        } else if (pattern == "G") {
            return 2;
        } else if (pattern == "T") {
            return 3;
        } else if (pattern == "-") {
            return 4;
        } else {
            return -1;
        }
    }
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
    if (_variation_db != NULL) {
        set<int> pos;
        vector<dbsnp_locus const*> snps = _variation_db->get_snps(chromosome->name(), start, end);
        for (int i = 0; i < (int)snps.size(); i++) {
            //pos.insert(snps[i]->position());
            int refid = get_base_id(snps[i]->reference());
            int altid = get_base_id(snps[i]->alternative());
            if (refid >= 0 && altid >= 0) {//snps[i]->reference().size() == 1 && snps[i]->alternative().size() == 1) {
                //int chromosome_code; // chromosome->name()
                loci.push_back(new hetero_locus(chromosome->code(), snps[i]->position(), refid, 0, altid, 0));//snps[i]->reference().c_str()[0], 0, snps[i]->alternative().c_str()[0], 0));
            }
        }
        if (loci.size() < 2) {
            return;
        }
        // cout << chromosome->name() << " " << start << "-" << end << " ";
        // cout << pos.size() << "SNPs\n";
    // } else {
    //     cout << "vcf null\n";
                //loci = scan_heterozygous_loci(fragments, chromosome, start, end, pos);
    } else {
        loci = scan_heterozygous_loci(fragments, chromosome, start, end);
    }
    int** patterns = new int*[fragments.size()];
    for (int i = 0; i < (int)fragments.size(); i++) {
        patterns[i] = new int[loci.size()];
        fragments[i]->generate_recombination_pattern(loci, patterns[i]);
    }

    //int tolerance = 2;
    int counts[4];


    for (int i = _snp_stretches; i < (int)loci.size() - _snp_stretches; i++) {
        for (int j = 0; j < 4; j++) counts[j] = 0;
        int head = i - _snp_stretches - _gap_tolerance;
        if (head < 0) {
            head = 0;
        }
        int tail = i + _snp_stretches + _gap_tolerance;
        if (tail > (int)loci.size()) {
            tail = (int)loci.size();
        }
        for (int j = 0; j < (int)fragments.size(); j++) {
            int const* pattern = patterns[j];
            int backward = 0;
            int forward = 0;
            int span_forward = 0;
            int span_backward = 0;
            // for (int k = i - _snp_stretches; k < i + _snp_stretches; k++) {
            //     if (pattern[k] == 1) {
            //         cout << "A";
            //     } else if (pattern[k] == -1) {
            //         cout << "a";
            //     } else {
            //         cout << ".";
            //     }
            // }
            // cout << "\t";
            for (int k = i - 1; k >= head; k--) {
                int p = pattern[k];
                if (p != 0) {
                    if (backward == 0) {
                        backward = p;
                    } else if (p + backward == 0) {//backward != p) {
                        //if (backward != p) {
                        backward = 0;
                        break;
                    }
                    span_backward ++;
                    if (span_backward >= _snp_stretches) {
                        break;
                    }
                }
            }
            for (int k = i; k < tail; k++) {
                int p = pattern[k];
                if (p != 0) {
                    if (forward == 0) {
                        forward = p;
                    } else if (forward + p == 0) {
                        forward = 0;
                        break;
                    }
                    span_forward ++;
                    if (span_forward >= _snp_stretches) {
                        break;
                    }
                }
            }
            //cout << i << ":" << j << "\t" << backward << "/" << forward << endl;
            if (forward != 0 && backward != 0) {
                counts[(forward == 1 ? 1 : 0) | (backward == 1 ? 2 : 0)] ++;
            }
        }
        if (check_acceptable_recombination(counts)) {
            ost << chromosome->name() << ":" << loci[std::max(0, i - _snp_stretches)]->position() << "-" << loci[std::min((int)loci.size() - 1, i + _snp_stretches - 1)]->position();
// << "\t" << loci[i]->ref() << ";" << loci[i]->alt() << "\t" << loci[i+l]->ref() << ";" << loci[i+l]->alt() << "\t";
            for (int j = 0; j < 4; j++) {
                ost << "\t";
                for (int k = i - _snp_stretches; k < i + _snp_stretches; k++) {
                    hetero_locus const* locus = loci[k];//- _snp_stretches + k];
                    if (j == 0) {
                        ost << locus->ref();
                    } else if (j == 3) {
                        ost << locus->alt();
                    } else {
                        if (k < i) {
                            ost << locus->ref();
                        } else {
                            ost << locus->alt();
                        }
                    }
                }
                ost << ";" << counts[j];
                //ost << loci[i - _snp_stretches + k - 1]
                //ost << "\t" << counts[j];
            }
            ost << "\n";
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


snp_enumerator::snp_enumerator(int coverage, float hetero_threshold) {
    _coverage = coverage;
    _minimum_minor_ratio = hetero_threshold;
    _display_mode = 1;
}

namespace {
    inline int get_index(char c) {
        switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case '-': return 4;
        }
        return -1;
    }
}

string snp_enumerator::get_genotype_symbol(dbsnp_locus const* snp, int counts[5]) const {//throw (exception) {
    int total = counts[0] + counts[1] + counts[2] + counts[3] + counts[4];
    if (total < _coverage) {
        return ".";
    }
    const string& reference = snp->reference();
    const string& alternative = snp->alternative();
    if (reference.size() == 1 && alternative.size() == 1) {
        int refind = get_index(reference.c_str()[0]);
        int altind = get_index(alternative.c_str()[0]);
        if (refind < 0 || altind < 0) {
            return "-";//undetermined";
        }
        int n0 = counts[refind];
        int n1 = counts[altind];
        if (n0 + n1 < _coverage) {
            return ".";
        }
        int thr = (int)(_minimum_minor_ratio * (n0 + n1) + .5);
        if (n0 < thr) {
            return "1/1";
        } else if (n1 < thr) {
            return "0/0";
        } else {
            return "0/1";
        }
    } else {
        if (reference.size() != 1) {
            return "-";//undetermined";
        }
        int refind = get_index(snp->reference().c_str()[0]);
        if (refind < 0) {
            return "-";//undetermined";
        }
        vector<string> items = split_items(alternative, ',');
        vector<pair<int,int> > nums;
        nums.push_back(make_pair(refind, counts[refind]));
        int total_alleles = counts[refind];
        for (int i = 0; i < (int)items.size(); i++) {
            if (items[i].size() == 1) {
                int index = get_index(items[i].c_str()[0]);
                if (index >= 0) {
                    nums.push_back(make_pair(index, counts[index]));
                    total_alleles += counts[index];
                } else {
		  nums.push_back(make_pair(-1, 0));
		}
            }
        }
	if (total_alleles < _coverage) {
	  return "-";
	}
        stringstream ss;
        int threshold = (int)(_minimum_minor_ratio * total_alleles + 0.5);
        //bool flag = false;
	int num_alleles = 0;

// 	cout << endl;
//         for (int i = 0; i < (int)nums.size(); i++) {
// 	  cout << i << ":" << nums[i].first << ", " << nums[i].second << endl;
// 	}
// 	cout << endl;

        for (int i = 0; i < (int)nums.size(); i++) {
            int max_index = -1;
            int max_counts = threshold;
	    int array_index = -1;
            for (int j = 0; j < (int)nums.size(); j++) {
                pair<int,int>& value = nums[j];
                if (value.second >= max_counts) {
		  array_index = j;
                    max_index = value.first;
                    max_counts = value.second;
                }
            }
            if (max_index >= 0) {//max_counts >= threshold) {
	      //cout << i << "::" << max_index << ", " << max_counts << endl;
	      if (num_alleles > 0) {
		ss << "/";
	      }
	      num_alleles ++;
	      ss << array_index;
	      nums[array_index].second = 0;
            }
        }
	if (num_alleles == 0) {
	  return "-";
	} else if (num_alleles == 1) {
	  ss << "/" + ss.str();
	}
        return ss.str();
    }
}

namespace {
    char nucleotides[6] = "ACGT-";
} 

void fragment_processor::set_display_mode(int mode) throw (std::invalid_argument) {
  if (mode < 0 || mode > 3) {
    throw std::invalid_argument("display mode is 0,1,2 or 3");
  }
  _display_mode = mode;
}

void snp_enumerator::process_fragments(const vector<recfragment*>& fragments,
                                       chromosome_seq const* chromosome,
                                       int start, int end, ostream& ost) 
    throw (exception){
    if (_variation_db != NULL) {
        //set<int> pos;
        vector<dbsnp_locus const*> snps = _variation_db->get_snps(chromosome->name(), start, end);
        for (int i = 0; i < (int)snps.size(); i++) {
            int position = snps[i]->position();
            if (position < start || position > end) continue;
            //pos.insert(snps[i]->position());
            int freq[5];
            int total = 0;
            freq[0] = freq[1] = freq[2] = freq[3] = freq[4] = 0;
            for (int j = 0; j < (int)fragments.size(); j++) {
                const recfragment* frag = fragments[j];
                int num;
                int index = frag->get_base_id(position, _quality_threshold, num);
                if (index >= 0) {
                    freq[index] += num;
                    total += num;
                }
            }
	    if (_display_mode > 0 && total < _coverage) continue;
	    string symbol = get_genotype_symbol(snps[i], freq);
	    if (_display_mode > 1 && (symbol == "." || symbol == "-" || symbol == "0/0")) {
                continue;
	    }
            ost << chromosome->name() << "\t" << position << "\t" << snps[i]->reference() << "\t" << snps[i]->alternative();
            ost << "\t" << symbol;
            for (int j = 0; j < 5; j++) {
                ost << "\t" << freq[j];
            }
            ost << "\n";
        }
    } else {
        for (int pos = start; pos < end; pos++) {
            char ref = chromosome->get_base(pos);
            int refind = get_index(ref);
            if (refind < 0) continue;
            int freq[5];
            freq[0] = freq[1] = freq[2] = freq[3] = freq[4] = 0;
            int total = 0;
            for (int j = 0; j < (int)fragments.size(); j++) {
                const recfragment* frag = fragments[j];
                int num;
                int index = frag->get_base_id(pos, _quality_threshold, num);
                if (index >= 0) {
                    freq[index] += num;
                    total += num;
                }
            }
            if (total < _coverage) continue;
            int max_minor = 0;
            int max_ind = -1;
            for (int j = 0; j < 5; j++) {
                if (j != refind && max_minor < freq[j]) {
                    max_minor = freq[j];
                    max_ind = j;
                }
            }
            if (max_ind >= 0 && max_minor >= (int)((freq[refind] + max_minor) * _minimum_minor_ratio + 0.5)) {
                ost << chromosome->name() << "\t" << pos << "\t" << nucleotides[refind] << "\t" << nucleotides[max_ind];
                for (int j = 0; j < 5; j++) {
                    ost << "\t" << freq[j];
                }
                ost << "\n";
            }
        }
    }
}
