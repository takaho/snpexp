#ifndef TKBIO_H
#define TKBIO_H

#include <string>
#include <bam.h>
#include <iosfwd>

//using namespace std;
using std::string;
using std::ostream;
using std::vector;
using std::pair;
using std::exception;
using std::invalid_argument;
using std::runtime_error;

typedef unsigned int uint;
namespace tkbio {
    class recfragment;

    //// Chromosome
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
        void set_sequence(int length, char const* codes);
        void set_base(int position, char base);
        void mask(int start, int end);
        char get_base(int pos) const;
        unsigned char get_base_code(int pos) const;
        int get_base_id(int pos) const;
        void set_bam_id(int num) {
            _bam_code = num;
        }
        static int get_chromosome_code(int length, const char* line);
        static vector<chromosome_seq*> load_genome(const char* filename) throw (exception);
    };

    /// Heterozygous locus
    /// containing chromosome information and reference/alt SNPs
    class hetero_locus {
    public:
        int _chromosome;
        int _position;
        int _reference;
        int _alt;
        int _refcount;
        int _altcount;
    public:
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
        int get_ref() const { return _reference; }
        int get_alt() const { return _alt; }
        char ref() const { return "ACGT-"[_reference]; }
        char alt() const { return "ACGT-"[_alt];}
        string to_string() const;
        bool is_available() const {
            // return 0 <= index1 && index1 < 5 && 0 <= index2 && index2 < 5;
            return _refcount > 0 && _altcount > 0 && _reference < 5 && _altcount < 5;
        }
        friend class recfragment;
    };

    class recpattern {
    public:
        typedef enum Genotype {REF_REF=0, REF_ALT=1, ALT_REF=2, ALT_ALT=3, UNDETERMINED=-1} Genotype;
    private:
        int _position5;
        int _position3;
        int _score;
        Genotype _genotype;
    public:
        recpattern(int pos5, int pos3, Genotype gt, int score=0) {
            _position5 = pos5;
            _position3 = pos3;
            _genotype = gt;
            _score = score;
        }
        Genotype genotype() const { return _genotype; }
        int position5() const { return _position5; }
        int position3() const { return _position3; }
        int score() const { return _score; }
    };

    ///// sequence read to detect recombination
    class recfragment {
        string _name;     // name of the read
        string _sequence; // read sequence
        string _cigar;    // CIGAR
        int _group;       // number to distinguish multiple bams
        int _chromosome;  // chromosome coded by tktools
        int _position;    // chromosomal position
        int _flag;        // SAM file flag
        int _position5;   // upper margin of alignment
        int _position3;   // lower margin of alignment
        int _max_match_span;  // maximum match length
        vector<pair<int,char> > _mapped;  // nucleotides
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
        void generate_recombination_pattern(const vector<hetero_locus>& loci, int* pattern) const;

    public:
        //recfragment(const string& chromosome, int position, int flag, const string& sequence, const string& cigar);
        // constructor with a BAM read
        recfragment(int chromosome, bam1_t* sequence, int group=0);

        // properties
        const string& name() const { return _name; }
        int chromosome() const { return _chromosome; }
        int group() const { return _group; }
        int position() const { return _position; }
        int span() const { return _position3 - _position5 + 1; }
        int position5() const { return _position5; }
        int position3() const { return _position3; }
        bool is_temporary() const {return _cigar == "";}
        char orientation() const { return (_flag & 16) == 0 ? '+' : '-'; }

        // getters
        char get_base(int position, int& count) const;
        int get_base_id(int position, int& count) const;
        void test_heterozygosity(int num, char const* reference, char* const alt, int* result) const;
        string to_string() const;

        // recombination test
        string get_recombination_pattern(const vector<hetero_locus>& loci) const;
        vector<recpattern> get_recombination_borders(const vector<hetero_locus>& loci, int minimum_span=2) const;
        recpattern::Genotype get_recombination_pattern(int pos5, int pos3, int minimum_span=2) const;
        pair<int,int> get_recombination(const vector<hetero_locus>& loci, float diff_ratio=0.0) const;

        // join pairs
        static void bundle_pairs(vector<recfragment*>& fragments) throw (runtime_error);

        // execute recombination detection
        static void detect_recombination(const vector<recfragment*>& frags,
                                         int clength, char const* cseq,
                                         int start, int end) throw (exception);
    };

}

#endif
