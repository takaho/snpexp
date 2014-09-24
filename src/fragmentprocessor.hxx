#ifndef TKBIO_BAMPROCESSOR_H
#define TKBIO_BAMPROCESSOR_H

#include <stdexcept>
#include <iosfwd>
#include <vector>
#include <set>
#include <string>

using std::exception;
using std::ostream;
using std::ostream;
using std::vector;
using std::set;
using std::string;

//#include <recrec.hxx>

//using namespace tkbio;

namespace tkbio {
    class recfragment;
    class chromosome_seq;
    class recpattern;
    class dbsnp_file;
    class hetero_locus;
    class dbsnp_locus;

    class fragment_processor {
    protected:
        dbsnp_file const* _variation_db;
        unsigned char _quality_threshold;
        int _display_mode; // 0: all, 1: covered, 2:non-standard, 3: hetero
    public:
        fragment_processor();
        virtual ~fragment_processor();
        virtual void process_fragments(const vector<recfragment*>& fragments,
                                       chromosome_seq const* chromosome,
                                       int start, int end, ostream& ost) throw (exception) = 0;
        void set_vcf(dbsnp_file const* vardb) { _variation_db = vardb; }
        void set_quality_threshold(unsigned char thr) {
            _quality_threshold = thr;
        }
        virtual void set_display_mode(int mode) throw (std::invalid_argument);
    };
    
    class recombination_detector : public fragment_processor {
        typedef enum Mode {DSBR=1, Meiotic=2, BOTH=3} Mode;
        int _coverage;
        float _minimum_minor_ratio;
        int _minimum_recombination_reads;
        float _minimum_recombination_composition;
        int _snp_stretches; // how many snps are required to define genotype (1)
        int _gap_tolerance; // how many
        double _allele_balance;
        Mode _recombination_mode; // 3:all 1:double-strand break 2:meitotic recombination
    public:
        recombination_detector(int coverage=25, float hetero_threshold=0.25f,
                               float recomb_threshold=0.05f);

        void set_allele_balance(double ratio);
        Mode get_detection_mode() const { return _recombination_mode; }
        void set_detection_mode(Mode mode);
        double allele_balance() const { return _allele_balance; }
        int coverage_threshold() const { return _coverage; }
        float heterozygosity_threshold() const { return _minimum_minor_ratio; }
        void set_haplotype_parameters(int stretches, int gaps=-1);
        int get_haplotype_streches() const { return _snp_stretches; }
        int get_gap_tolerance() const { return _gap_tolerance; }
        virtual void process_fragments(const vector<recfragment*>& fragments,
                                       chromosome_seq const* chromosome,
                                       int start, int end, ostream& ost) throw (exception);
        virtual bool check_acceptable_recombination(int counts[4]) const;
        vector<hetero_locus*> scan_heterozygous_loci(const vector<recfragment*>& fragments, chromosome_seq const* chromosome, int start, int end, const set<int>& accepted) const throw (exception);
        vector<hetero_locus*> scan_heterozygous_loci(const vector<recfragment*>& fragments, chromosome_seq const* chromosome, int start, int end) const throw (exception);
    };

    class snp_enumerator : public fragment_processor {
        int _coverage;
        float _minimum_minor_ratio;
        int _display_mode; // 0: all, 1: covered, 2:non-standard, 3: hetero
    private:
        string get_genotype_symbol(dbsnp_locus const* snp, int counts[5]) const;
    public:
        int mode() const { return _display_mode; }
        snp_enumerator(int coverage, float hetero_threshold=0.0f);
        virtual void process_fragments(const vector<recfragment*>& fragments,
                                       chromosome_seq const* chromosome,
                                       int start, int end, ostream& ost) throw (exception);

    };
}


#endif
