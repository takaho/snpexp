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
    _allele_balance = 0.0;
    _recombination_mode = DSBR;
}

void recombination_detector::set_detection_mode(Mode mode) {
    _recombination_mode = mode;
}

bool recombination_detector::check_acceptable_recombination(int counts[4]) const {
    int total = counts[0] + counts[1] + counts[2] + counts[3];
    if (total < _minimum_recombination_reads) {
        return false;
    }
    int a1 = counts[0] + counts[1];
    int A1 = counts[2] + counts[3];
    int a2 = counts[0] + counts[2];
    int A2 = counts[1] + counts[3];
    if (a1 * _allele_balance > A1 || a1 < _allele_balance * A1 
        || a2 * _allele_balance > A2 || a2 < _allele_balance * A2) {
        return false;
    }
    int thr = (int)(total * _minimum_recombination_composition + 0.5);
    int num_zeros = 0;
    int num_detected = 0;
    for (int i = 0; i < 4; i++) {
        if (counts[i] < thr) {
            if (counts[i] == 0) {
                num_zeros ++;
            }
        } else {
            num_detected ++;
        }
    }
    if (_recombination_mode == DSBR) {
        return (num_detected == 3 && num_zeros == 1);
    } else if (_recombination_mode == Meiotic) {
        return num_detected == 4;
    } else {
        return (num_detected == 3 && num_zeros == 1) || (num_detected == 4);
    }
}

void recombination_detector::set_allele_balance(double ratio) {
    if (ratio <= 0.0) {
        _allele_balance = 0.0;
    } else if (ratio > 1.0) {
        _allele_balance = 1.0 / ratio;
    } else {
        _allele_balance = ratio;
    }
}

void recombination_detector::set_haplotype_parameters(int stretches, int gaps) {
    _snp_stretches = stretches;
    if (gaps >= 0) {
        _gap_tolerance = gaps;
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
            int refid = hetero_locus::get_base_id(snps[i]->reference());
            int altid = hetero_locus::get_base_id(snps[i]->alternative());
            if (refid >= 0 && altid >= 0) {//snps[i]->reference().size() == 1 && snps[i]->alternative().size() == 1) {
                //int chromosome_code; // chromosome->name()
                hetero_locus* hl = new hetero_locus(chromosome->code(), snps[i]);
                if (hl->is_available()) {
                    loci.push_back(hl);
                } else {
                    delete hl;
                }
                //                loci.push_back(new hetero_locus(chromosome->code(), snps[i]->position(), refid, 0, altid, 0));//snps[i]->reference().c_str()[0], 0, snps[i]->alternative().c_str()[0], 0));
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
    int haplo_counts[4];
    int* genotype_cache = new int[fragments.size()];
    int tail = (int)loci.size() - _snp_stretches * 2;
    for (int site_start = 0; site_start < tail; site_start++) {
        int gap_limit = std::min(_gap_tolerance, (int)loci.size() - site_start - _snp_stretches * 2);
        for (int gap = 0; gap <= gap_limit; gap++) {
            // clear counts
            for (int i = 0; i < 4; i++) haplo_counts[i] = 0;

            // set genotype [site_start:site_start + _snp_stretches] when gap is 0
            if (gap == 0) {
                for (int i = 0; i < (int)fragments.size(); i++) {
                    int const* pattern = patterns[i];
                    int genotype = pattern[site_start];
                    if (genotype != 0) {
                        for (int j = 1; j < _snp_stretches; j++) {
                            int p = pattern[site_start + j];
                            if (p != genotype) {
                                genotype = 0;
                                break;
                            }
                        }
                    }
                    genotype_cache[i] = genotype;
                }
            }

            // set forward genotype
            int num_available = 0;
            int forward_start = site_start + _snp_stretches + gap;

            int steps = std::min(_snp_stretches, (int)loci.size() - forward_start);
            for (int i = 0; i < (int)fragments.size(); i++) {
                if ((int)fragments.size() - i + num_available < _minimum_recombination_reads) {
                    break;
                }
                int genotype_backward = genotype_cache[i];
                if (genotype_backward == 0) {
                    continue;
                }
                int const* pattern = patterns[i];
                int genotype = pattern[forward_start];
                if (genotype != 0) {
                    for (int j = 1; j < steps; j++) {
                        if (genotype != pattern[forward_start + j]) {
                            genotype = 0;
                            break;
                        }
                    }
                }
                // if genotype is determined, count up the haplotype
                if (genotype != 0) {
                    num_available++;
                    int genotype_forward = genotype;
                    haplo_counts[(genotype_backward > 0 ? 2 : 0) | (genotype_forward > 0 ? 1 : 0)]++;
                }
            }
            if (check_acceptable_recombination(haplo_counts)) {
                hetero_locus const* snp_5 = loci[site_start];// + _snp_stretches - 1];
                hetero_locus const* snp_3 = loci[forward_start + _snp_stretches - 1];
                ost << chromosome->name() << ":" << snp_5->position() << "-" << snp_3->position() << "\t" << snp_5->id() << "-" << snp_3->id();
                for (int j = 0; j < 4; j++) {
                    ost << "\t";
                    for (int k = 0; k < _snp_stretches; k++) {
                        hetero_locus const* locus = loci[site_start + k];
                        ost << (((j & 2) != 0) ? locus->ref() : locus->alt());
                    }
                    for (int k = 0; k < _snp_stretches; k++) {
                        hetero_locus const* locus = loci[forward_start + k];
                        ost << (((j & 1) != 0) ? locus->ref() : locus->alt());
                    }
                    ost << ";" << haplo_counts[j];
                }
                ost << "\n";
                site_start += _snp_stretches + gap;
                break;
            }
        }
    }
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
	int num_alleles = 0;

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
