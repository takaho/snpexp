#ifndef TKBIO_VARSTR_H
#define TKBIO_VARSTR_H
#include <string>
#include <stdexcept>
#include <vector>
#include <map>

using std::string;
using std::vector;
using std::exception;

#include <bam.h>

#include <recrec.hxx>

using namespace tkbio;

namespace tkbio {

    typedef unsigned long long ullong;

    class bamread {

        const static ullong FEATURE_MASK             = 0x3000000000000LL;
        const static ullong FEATURE_MATCH            = 0x0000000000000LL;
        const static ullong FEATURE_REPEAT_INSERTION = 0x1000000000000LL; 
        const static ullong FEATURE_REPEAT_DELETION  = 0x2000000000000LL;

        vector<unsigned long long> _sections;

        int _chromosome;
        int _start;
        int _stop;
    public:
        bamread(const bam1_t* read);
        int chromosome() const { 
            return _chromosome;
        }
        int position5() const { return _start; }
        int position3() const { return _stop; }
        int size() const { return _sections.size(); }
        ullong get_section(int index) const { return _sections[index]; }
        static int get_position(ullong info) { return (int)info & 0xffffffff; }
        static int get_span(ullong info) { return (int)((info >> 32) & 0xffff); }
        static char get_feature(ullong info);// { return info & FEATURE_MASK; }
        bool covers(int start, int end) const { return start <= _start && end <= _stop; }
        bool shares_variation(ullong section) const;

        string to_string() const;
        friend class strcollection;
    };

    class str_variation {
        int _position;
        int _reference_span;
        int _read_span;
        int _num_shares;
        int _num_reference;
        int _num_others;
    public:
        str_variation(int position, int reference_span, int read_span);
        void set_counts(int share, int reference, int others);

        string to_string() const;
        friend class strcollection;
        friend ::operator == (const str_variation& lhs, const str_variation& rhs);
    };

    class str_collection {
    private:
        vector<bamread> _reads;
    public:
        str_collection();
        ~str_collection();
        void sweep(int chromosome, int start=0, int end=-1);
        void add_read(bam1_t const* read);
        int size() const { return _reads.size(); }
        vector<strvariation> get_variations(int coverage, double heterozygosity) const throw (exception);
        pair<int,int> span() const;
        static int detect_str(int argc, char** argv) throw (exception);
    };
}
#endif
