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
    };
    
    class recombination_detector : public fragment_processor {
        int _coverage;
        float _minimum_minor_ratio;
        int _minimum_recombination_reads;
        float _minimum_recombination_composition;
    public:
        recombination_detector(int coverage=25, float hetero_threshold=0.25f,
                               float recomb_threshold=0.05f);

        int coverage_threshold() const { return _coverage; }
        float heterozygosity_threshold() const { return _minimum_minor_ratio; }
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
    private:
        string get_genotype_symbol(dbsnp_locus const* snp, int counts[5]) const;
    public:
        snp_enumerator(int coverage, float hetero_threshold=0.0f);
        virtual void process_fragments(const vector<recfragment*>& fragments,
                                       chromosome_seq const* chromosome,
                                       int start, int end, ostream& ost) throw (exception);

    };
}


#endif
