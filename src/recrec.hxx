#ifndef TKBIO_H
#define TKBIO_H

#include <string>
#include <bam.h>
#include <iosfwd>

//using namespace std;
using std::string;
using std::ostream;

typedef unsigned int uint;
namespace tkbio {
    class recfragment;

    class hetero_locus {
    public:
        int _chromosome;
        int _position;
        unsigned char _reference;
        unsigned char _alt;
        int _refcount;
        int _altcount;
        // int position;
        // int index1;
        // int index2;
        // int count1;
        // int count2;
    public:
        // hetero_locus(int position, int index1, int index2) {
        //     this->position = position;
        //     this->index1 = index1;
        //     this->index2 = index2;
        //     this->count1 = this->count2 = 0;
        // }
        hetero_locus(int chromosome, int position, unsigned char reference, int reference_count, unsigned char alt, int alt_count) {
            this->_chromosome = chromosome;
            this->_position = position;
            this->_reference = reference;
            this->_alt = alt;
            this->_refcount = reference_count;
            this->_altcount = alt_count;
        }
        int position() const { return _position; }
        int chromosome_code() const { return _chromosome; }
        string chromosome() const;
        int count_alt() const { return _altcount; }
        int count_ref() const { return _refcount; }
        unsigned char get_ref() const { return _reference; }
        unsigned char get_alt() const { return _alt; }
        char ref() const { return "ACGT-"[_reference]; }
        char alt() const { return "ACGT-"[_reference]; }
        string to_string() const;
        // hetero_locus(int position, const string& ref, const string& alt) {
        //     this->position = position;
        //     this->index1 = get_index(ref);
        //     this->index2 = get_index(alt);
        // }
        bool is_available() const {
            // return 0 <= index1 && index1 < 5 && 0 <= index2 && index2 < 5;
            return _refcount > 0 && _altcount > 0 && _reference < 5 && _altcount < 5;
        }
        friend class recfragment;
    // public:
    //     static int get_index(const string& locus) {
    //         if (locus == "A") {
    //             return 0;
    //         } else if (locus == "C") {
    //             return 1;
    //         } else if (locus == "G") {
    //             return 2;
    //         } else if (locus == "T") {
    //             return 3;
    //         } else if (locus == "." || locus == "-") {
    //             return 4;
    //         } else {
    //             return -1;
    //         }
    //     }
    };

    class chromosome_seq {
    public:
        static const unsigned char A;
        static const unsigned char C;
        static const unsigned char G;
        static const unsigned char T;
        static const unsigned char N;
        static const unsigned char INS;
    private:
        int _length;
        unsigned char* _sequence;
        string _name;
        int _code;
        int _bam_code;
    private:
        const chromosome_seq& operator = (const chromosome_seq& rhs);
        chromosome_seq(const chromosome_seq& rhs);
    public:
        chromosome_seq() {
            _sequence = NULL;
            _code = _bam_code = -1;
            _length = 0;
        }
        ~chromosome_seq() {
            delete[] _sequence;
        }
        int length() const { return _length; }
        int code() const { return _code; }
        int bam_id() const { return _bam_code; }
        const string& name() const { return _name; }
        void set_chromosome(int num);
        void set_sequence(int length, unsigned char* codes) {
            if (_sequence != NULL) {
                delete[] _sequence;
            }
            _length = length;
            int size = length / 2 + 1;
            _sequence = new unsigned char[size];
            memcpy(_sequence, codes, sizeof(unsigned char) * size);
        }
        char get_base(int pos) const;
        unsigned char get_base_code(int pos) const {
            unsigned char code = _sequence[pos >> 1] >> ((pos & 1) == 0 ? 4 : 0) & 0x0f;
            return code;
        }
        void set_bam_id(int num) {
            _bam_code = num;
        }
        static int get_chromosome_code(int length, const char* line);
        static vector<chromosome_seq*> load_genome(const char* filename) throw (exception);
    };

    class recfragment {
        string _name;
        string _sequence;
        string _cigar;
        int _chromosome;
        int _position;
        int _flag;
        vector<pair<int,char> > _mapped;
        int _position5;
        int _position3;
        int _max_match_span;
    private:
        // inhibit default and copy constructors
        recfragment();
        const recfragment& operator = (const recfragment&);
        recfragment(const recfragment& rhs);

        // initializers
        void initialize(int position, const string& sequence, const string& cigar);
        void initialize(bam1_t* sequence);

        //
        vector<pair<int,char> > get_positions() const;
        void join_sequence(const recfragment* frag);
        //void generate_recombination_pattern(const vector<pair<int,char> >& loci, int* pattern) const;
        void generate_recombination_pattern(const vector<hetero_locus>& loci, int* pattern) const;

    public:
        //recfragment(const string& chromosome, int position, int flag, const string& sequence, const string& cigar);
        // constructor with a BAM read
        recfragment(int chromosome, bam1_t* sequence);

        // properties
        const string& name() const { return _name; }
        int chromosome() const { return _chromosome; }
        int position() const { return _position; }
        int span() const { return _position3 - _position5 + 1; }
        int position5() const { return _position5; }
        int position3() const { return _position3; }
        bool is_temporary() const {return _cigar == "";}
        char orientation() const { return (_flag & 16) == 0 ? '+' : '-'; }

        // getters
        char get_base(int position, int& count) const;
        void test_heterozygosity(int num, char const* reference, char* const alt, int* result) const;

        // recombination test
        string get_recombination_pattern(const vector<hetero_locus>& loci) const;
        pair<int,int> get_recombination(const vector<hetero_locus>& loci, float diff_ratio=0.0) const;
        
        //string get_recombination_pattern(const vector<pair<int,char> >& loci) const;
        //pair<int,int> get_recombination(const vector<pair<int,char> >& loci, float diff_ratio=0.0) const;

        // join pairs
        static void bundle_pairs(vector<recfragment*>& fragments) throw (runtime_error);

        static vector<pair<int,char> > resolve_map(int position, const string& cigar, const string& sequence);
        static void detect_recombination(const vector<recfragment*>& frags,
                                         int clength, char const* cseq,
                                         int start, int end) throw (exception);
    };

}

#endif
