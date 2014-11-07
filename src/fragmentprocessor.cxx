#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <limits>

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
#include <gtf.hxx>

using tktools::split_items;
using namespace tkbio;

fragment_processor::fragment_processor() {
    _variation_db = NULL;
    _quality_threshold = 10;
    _gtffile = NULL;
}

fragment_processor::~fragment_processor() {
}

string fragment_processor::to_string() const {
    return string("");
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

void recombination_detector::set_gtf(gtffile const* gtf) throw (exception) {
    throw runtime_error("recombination detection does not support GTF file");
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

    // detect heterozygous loci
    vector<hetero_locus*> loci;
    if (_variation_db != NULL) {
        set<int> pos;
        vector<dbsnp_locus const*> snps = _variation_db->get_snps(chromosome->name(), start, end);
        for (int i = 0; i < (int)snps.size(); i++) {
            int refid = hetero_locus::get_base_id(snps[i]->reference());
            int altid = hetero_locus::get_base_id(snps[i]->alternative());
            if (refid >= 0 && altid >= 0) {//snps[i]->reference().size() == 1 && snps[i]->alternative().size() == 1) {
                hetero_locus* hl = new hetero_locus(chromosome->code(), snps[i]);
                if (hl->is_available()) {
                    loci.push_back(hl);
                } else {
                    delete hl;
                }
            }
        }
        if (loci.size() < 2) {
            return;
        }
    } else {
        // de novo detection of SNPs
        loci = scan_heterozygous_loci(fragments, chromosome, start, end);
    }

    // get genotypes as integer arrays
    int** patterns = new int*[fragments.size()];
    for (int i = 0; i < (int)fragments.size(); i++) {
        patterns[i] = new int[loci.size()];
        fragments[i]->generate_recombination_pattern(loci, patterns[i]);
    }
    int* genotype_cache = new int[fragments.size()];

    int haplo_counts[4];
    int tail = (int)loci.size() - _snp_stretches * 2;
    for (int site_start = 0; site_start < tail; site_start++) {
        int gap_limit = std::min(_gap_tolerance, (int)loci.size() - site_start - _snp_stretches * 2);
        for (int gap = 0; gap <= gap_limit; gap++) {
            // clear counts
            for (int i = 0; i < 4; i++) haplo_counts[i] = 0;

            if (gap == 0) { // set backward haplotypes
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

            // set forward haplotypes
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

            // output results if recombination is detected
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

    // clean up
    delete[] genotype_cache;
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

namespace {
    bool compare_by_second(const pair<int,int>& lhs, const pair<int,int>& rhs) {
        return lhs.second < rhs.second;
    }
    int get_next_position(const vector<pair<int,int> >& exons, int position) {
        if (exons.size() == 0) {
            return std::numeric_limits<int>::max();
        }
        int left = 0;
        int right = (int)exons.size();
        for (;;) {
            int center = (left + right) >> 1;
            const pair<int,int>& exon = exons[center];
            if (exon.first > position) {
                left = center + 1;
            } else if (exon.second < position) {
                right = center;
            } else {
                return position;
            }
            if (left == right) {
                return std::numeric_limits<int>::max();
            }
        }
    }
}

namespace {
    bool check_heterozygosity(int counts[5], double hvalue, int reference) {
        if (reference < 0 || reference > 5) {
            return false;
        }
        int num = counts[reference];
        //int max_index = -1;
        int max_counts = 0;
        for (int i = 0; i < 4; i++) {
	  if (i != reference && counts[i] > (int)max_counts) {
                max_counts = counts[i];
                //max_index = i;
            }
        }
        int total = max_counts + num;
        double ratio = (double)num / total;
        if (ratio >= hvalue && hvalue <= 1.0 - hvalue) {
            return true;
        } else {
            return false;
        }
    }
}

// namespace {
//     vector<pair<int,int> > get_exon_regions(
// }

void snp_enumerator::process_fragments(const vector<recfragment*>& fragments,
                                       chromosome_seq const* chromosome,
                                       int start, int end, ostream& ost) 
    throw (exception){

    // use GTF
    vector<pair<int,int> > exons;
    if (_gtffile != NULL) {
        vector<const gtfgene*> genes = _gtffile->find_genes(chromosome->name(), start, end);
        for (int i = 0; i < (int)genes.size(); i++) {
            const vector<gtfexon>& geneexons = genes[i]->exons();
            for (int j = 0; j < (int)geneexons.size(); j++) {
                pair<int,int> epos = geneexons[j].position();
                if (epos.second < start || epos.first > end) continue;
                exons.push_back(epos);
            }
        }
        if (exons.size() == 0) {
            return;
        }
        sort(exons.begin(), exons.end(), compare_by_second);
        for (int i = (int)exons.size() - 1; i >= 0; i--) {
            if (exons[i].first == exons[i].second) continue;
            for (int j = i - 1; j >= 0; j--) {
                if (exons[j].second < exons[i].first) {
                    break;
                } else {
                    exons[i].first = std::min(exons[i].first, exons[j].first);
                    exons[i].second = std::max(exons[i].second, exons[j].second);
                    exons[j].first = exons[j].second;
                }
            }
        }
        int index = 0;
        for (int i = 0; i < (int)exons.size(); i++) {
            pair<int,int>& e = exons[i];
            if (e.first != e.second) {
                exons[index] = e;
                index++;
            }
        }
        exons.erase(exons.begin() + index, exons.end());
    }

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
            if (_display_mode == 0 || (_display_mode > 0 && check_heterozygosity(freq, _minimum_minor_ratio, get_index(snps[i]->reference().c_str()[0])))) {
                ost << chromosome->name() << "\t" << position << "\t" << snps[i]->reference() << "\t" << snps[i]->alternative();
                ost << "\t" << symbol;
                for (int j = 0; j < 5; j++) {
                    ost << "\t" << freq[j];
                }
                ost << "\n";
            }
        }
    } else {
        for (int pos = start; pos < end; pos++) {
            if (_gtffile != NULL) {
                pos = get_next_position(exons, pos);
                if (pos >= end) {
                    break;
                }
            }

            char ref = chromosome->get_base(pos);
            int refind = get_index(ref);
            if (refind < 0) continue;
            int freq[5];
            for (int i = 0; i < 5; i++) { freq[i] = 0; }
            //freq[0] = freq[1] = freq[2] = freq[3] = freq[4] = 0;
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
                ost << chromosome->name() << "\t" << pos << "\t" << ref << "\t" << nucleotides[max_ind];
                for (int j = 0; j < 5; j++) {
                    ost << "\t" << freq[j];
                }
                ost << "\n";
            }
        }
    }
}

strain_estimator::strain_estimator() throw (exception) {
    _coverage = 50;
    _hetero_thr = 0.05f;
    _num_strains = 0;
    _matrix = NULL;
    _variation_db = NULL;
}

strain_estimator::strain_estimator(dbsnp_file const* vardb) throw (exception) {
    _coverage = 50;
    _hetero_thr = 0.05f;
    initialize_matrix(vardb);
}

/***
    strain_1/strain_2  0/0 0/1 1/0 1/1
    genotype           0/0 0/1 1/1
 => 4 * 3 = 12
 */
namespace {
    const int STRAIN_GROUPS = 12;
}
void strain_estimator::initialize_matrix(dbsnp_file const* vardb) {
    _variation_db = vardb;
    _num_strains = _variation_db->strain_size();
    int matsize = _num_strains * (_num_strains - 1) / 2;
    _matrix = new int*[matsize];
    int index = 0;
    for (int i = 0; i < _num_strains; i++) {
        for (int j = 0; j < i; j++) {
            _matrix[index] = new int[STRAIN_GROUPS];
            for (int k = 0; k < STRAIN_GROUPS; k++) {
                _matrix[index][k] = 0;
            }
        }
    }
}

strain_estimator::~strain_estimator() {
    int index = 0;
    for (int i = 0; i < _num_strains; i++) {
        for (int j = 0; j < i; j++) {
            delete[] _matrix[index++];
        }
    }
    delete[] _matrix;
}

void strain_estimator::process_fragments(const vector<recfragment*>& fragments,
                                         chromosome_seq const* chromosome,
                                         int start, int end, ostream& ost) throw (exception) {


    vector<dbsnp_locus const*> snps = _variation_db->get_snps(chromosome->name(), start, end);
    //int num_slots = _num_strains * (_num_strains - 1) / 2;
    //int* indices = new int[num_slots];
    //for (int i = 0; i < num_slots; i++) indices[i] = -1;
    float thr_lower = 0.05f;
    float thr_upper = 1.0f - thr_lower;
    for (int i = 0; i < (int)snps.size(); i++) {
        // get nucleotide and reject non-standards
        int refid = hetero_locus::get_base_id(snps[i]->reference());
        int altid = hetero_locus::get_base_id(snps[i]->alternative());
        if (refid < 0 || altid < 0) continue;
        
        // count bases and reject few coverage
        int position = snps[i]->position();
        if (position < start || position >= end) continue;
        int freq[5];
        //int total = 0;
        freq[0] = freq[1] = freq[2] = freq[3] = freq[4] = 0;
        for (int j = 0; j < (int)fragments.size(); j++) {
            const recfragment* frag = fragments[j];
            int num;
            int index = frag->get_base_id(position, _quality_threshold, num);
            if (index >= 0) {
                freq[index] += num;
                //total += num;
            }
        }
        if (freq[refid] + freq[altid] < _coverage) continue;
        float ratio = (float)freq[altid] / (freq[refid] + freq[altid]);

        int group = 0;
        if (ratio < thr_lower) { // ref/ref
            group = 0;
        } else if (ratio >= thr_upper) { // alt/alt
            group = 2; 
        } else {
            group = 1;
        }
        //if (total < _coverage) continue;

        //
        int num = 0;
        for (int j = 0; j < _num_strains; j++) {
            unsigned char g1 = snps[i]->get_genotype(j);
            int index = 0;
            if (g1 == (unsigned char)0x00) {
                index = 0;
            } else if (g1 == (unsigned char)0x11) {
                index = 6;
            } else{
                index = -1;
                continue;
            }
            for (int k = 0; k < j; k++) {
                unsigned char g2 = snps[i]->get_genotype(k);
                if (g2 == (unsigned char)0x00) {
                    //
                } else if (g2 == (unsigned char)0x11) {
                    index += 3;
                } else{
                    index = -1;
                    continue;
                }
                _matrix[num++][index + group] ++;
                //indices[num++] = index;
            }
        }
    }
//    delete[] indices;
}

string strain_estimator::to_string() const {
    stringstream ss;
    ss << "#Strain1\tStrain2";
    for (int i = 0; i < STRAIN_GROUPS; i++) {
        ss << "\t";
        ss << (i / 6) << "/" << ((i / 3) % 2) << ":";
        switch (i % 3) {
        case 0:
            ss << "0/0"; break;
        case 1:
            ss << "0/1"; break;
        case 2:
            ss << "1/1"; break;
        }
    }
    ss << endl;
    int index = 0;
    for (int i = 0; i < _num_strains; i++) {
        ss << _variation_db->get_strain(i);
        for (int j = 0; j < i; j++) {
            ss << "\t" << _variation_db->get_strain(j);
            for (int k = 0; k < STRAIN_GROUPS; k++) {
                ss << "\t" << _matrix[index][k];
            }
            ss << endl;
            index++;
        }
    }
    return ss.str();
}
