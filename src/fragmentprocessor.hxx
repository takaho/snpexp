#ifndef TKBIO_BAMPROCESSOR_H
#define TKBIO_BAMPROCESSOR_H

#include <stdexcept>
#include <iosfwd>
#include <vector>
#include <set>

using std::exception;
using std::ostream;
using std::ostream;
using std::vector;
using std::set;

//#include <recrec.hxx>

//using namespace tkbio;

namespace tkbio {
    class recfragment;
    class chromosome_seq;
    class recpattern;
    class dbsnp_file;
    class hetero_locus;

    class fragment_processor {
    protected:
        dbsnp_file const* _variation_db;
    public:
        fragment_processor();
        virtual ~fragment_processor();
        virtual void process_fragments(const vector<recfragment*>& fragments,
                                       chromosome_seq const* chromosome,
                                       int start, int end, ostream& ost) throw (exception) = 0;
        void set_vcf(dbsnp_file const* vardb) { _variation_db = vardb; }
    };
    
    class recombination_detector : public fragment_processor {
        int _coverage;
        float _minimum_minor_ratio;

    public:
        recombination_detector(int coverage=20, float hetero_threshold=0.25f);

        int coverage_threshold() const { return _coverage; }
        float heterozygosity_threshold() const { return _minimum_minor_ratio; }
        virtual void process_fragments(const vector<recfragment*>& fragments,
                                       chromosome_seq const* chromosome,
                                       int start, int end, ostream& ost) throw (exception);
        vector<hetero_locus*> scan_heterozygous_loci(const vector<recfragment*>& fragments, chromosome_seq const* chromosome, int start, int end, set<int> accepted) const throw (exception);
        vector<hetero_locus*> scan_heterozygous_loci(const vector<recfragment*>& fragments, chromosome_seq const* chromosome, int start, int end) const throw (exception);
    };

    class snp_enumerator : public fragment_processor {
        int _coverage;
    public:
        snp_enumerator(int coverage);
    };
}


#endif
