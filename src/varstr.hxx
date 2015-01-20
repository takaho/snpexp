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

    class bamread;
    class str_collection;
    class str_variation;

    class bamread {
    public:
        const static ullong FEATURE_MASK             = 0x3000000000000LL;
        const static ullong FEATURE_MATCH            = 0x0000000000000LL;
        const static ullong FEATURE_REPEAT_INSERTION = 0x1000000000000LL; 
        const static ullong FEATURE_REPEAT_DELETION  = 0x2000000000000LL;
        const static ullong FEATURE_ERROR            = 0x3000000000000LL;

    private:
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
        bool covers(int start, int end) const;
        bool shares_variation(ullong section) const;
        bool has_reference_at(int start, int end) const;

        string to_string() const;
        friend class strcollection;
    };

    class str_variation {
        ullong _code;
        int _num_shares;
        int _num_reference;
        int _num_others;
    public:
        str_variation();
        str_variation(ullong section_code);
        str_variation(int position, int reference_span, int read_span);
        void set_counts(int share, int reference, int others);
        int occurrence() const { return _num_shares; }
        int opposite() const { return _num_reference; }
        int coverage() const { return _num_shares + _num_reference + _num_others; }
        string to_string() const;
        ullong code() const { return _code; }
        int reference_span() const;
        int read_span() const;
        int position() const;
        char feature() const;
        friend bool operator == (const str_variation& lhs, const str_variation& rhs);
        friend class str_collection;
    };

    class repeat_region {
    public:
        typedef enum STATUS {DEFAULT=0, COVERED=1, POLYMORPHIC=2} STATUS;
    private:
        int _chromosome;
        string _pattern;
        int _start;
        int _stop;
        //bool _accepted;
        STATUS _status;
        //int _coverage;
    public:
        repeat_region(const repeat_region& rhs);
        const repeat_region& operator = (const repeat_region& rhs);
        repeat_region(int size, const char* sequence, int start, int stop);
        repeat_region(int size, const char* sequence, int chromosome, int start, int stop);
        repeat_region(int start=0, int stop=0);
        //string chromosome() const;
        int chromosome() const { return _chromosome; }
        int unitsize() const { return _pattern.size(); }
        int start() const { return _start; }
        int stop() const { return _stop; }
        int span() const { return _stop - _start; }
        int repeats() const { return span() / unitsize(); }
        bool accepted() const { return (_status & COVERED & POLYMORPHIC) != 0; }
        bool covered() const { return (_status & COVERED) != 0; }
        bool polymorphic() const { return (_status & POLYMORPHIC) != 0; }
        STATUS status() const { return _status; }
        void set_status(STATUS status) { _status = status; }
        //void set_accepted(bool flag) { _accepted = flag; }
        //int coverage() const { return _coverage; }
        //void set_coverage(int coverage) { _coverage = coverage; }
        const string& pattern() const { return _pattern; }
        string to_string() const;
        static repeat_region* analyze_str(int size, char const* buffer, int seq_pos, int unit_min, int unit_max, int required_span) throw (exception);
        static void enumerate_repeat_regions(int argc, char** argv) throw (exception);
        //static bool compare_position(const repeat_result& lhs, const repeat_result& rhs);
        static bool compare_position(const repeat_region* lhs, const repeat_region* rhs);
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
        vector<str_variation> get_variations(int coverage, double heterozygosity) const throw (exception);
        pair<int,int> span() const;
        int count_coverage(int start, int stop) const;
        static int detect_str(int argc, char** argv) throw (exception);
    };

}
#endif
