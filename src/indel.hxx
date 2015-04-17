#include <vector>
#include <stdexcept>
#include <string>

using std::vector;
using std::exception;
using std::string;

typedef unsigned char uchar;
typedef unsigned long long ullong;
typedef unsigned short ushort;

namespace tkbio {
    class indel {
    public:
        const static uchar INSERTION = (uchar)0x01;
        const static uchar DELETION  = (uchar)0x00;
        const static uchar REJECTED = (uchar)0x02;
    private:
        uchar _feature;
        uchar _id;
        ushort _span;
        ushort _count;
        ushort _coverage;
        int _position;
    public:
        indel();
        indel(const indel& rhs);
        const indel& operator = (const indel& rhs);
        indel(int id, int position, int span, bool insertion=true);
        //indel(bam1_t const* read);
        int id() const { return _id; }
        int count() const { return _count; }
        int coverage() const { return _coverage; }
        int position() const { return _position; }
        int span() const { return _span; }
        void reject() { _feature |= REJECTED; }
        bool is_insertion() const { return (_feature & INSERTION) != 0; }
        bool is_rejected() const { return (_feature & REJECTED) != 0; }

        void set_id(int id) { _id = id; }
        void set_coverage(int c) { _coverage = c; }
        void set_count(int c) { _count = c; }

        double heterozygosity() const;
        string to_string() const;

        static vector<indel> parse_indel(bam1_t const* read, int quality=0);
        static int detect_indel_polymorphism(int argc, char** argv) throw (exception);
        static int detect_multiple_polymorphism(int argc, char** argv) throw (exception);
        friend bool operator == (const indel& lhs, const indel& rhs);
    };
}    
