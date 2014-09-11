#ifndef TKBIO_BAMPROCESSOR_H
#define TKBIO_BAMPROCESSOR_H

#include <stdexcept>
#include <iosfwd>
#include <vector>

using std::exception;
using std::ostream;
using std::ostream;
using std::vector;

//#include <recrec.hxx>

//using namespace tkbio;

namespace tkbio {
    class recfragment;
    class chromosome_seq;
    class recpattern;
    class dbsnp_file;

    class fragment_processor {
    protected:
        dbsnp_file* _variation_db;
    public:
        virtual ~fragment_processor() {};
        virtual void process_fragments(const vector<recfragment*>& fragments,
                               chromosome_seq const* chromosome,
                               int start, int end, ostream& ost) throw (exception) = 0;
        void set_vcf(dbsnp_file* vardb) { _variation_db = vardb; }
    };
    
    class recombination_detector : public fragment_processor {
        int _coverage;
        float _threshold;
    public:
        recombination_detector(int coverage=20, float hetero_threshold=0.25f);

        int coverage_threshold() const { return _coverage; }
        float heterozygosity_threshold() const { return _threshold; }
        virtual void process_fragments(const vector<recfragment*>& fragments,
                               chromosome_seq const* chromosome,
                                       int start, int end, ostream& ost) throw (exception);
    };

    class snp_enumerator : public fragment_processor {
        int _coverage;
    public:
        
    };
}


#endif
