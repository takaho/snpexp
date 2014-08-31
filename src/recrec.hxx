#ifndef TKBIO_H
#define TKBIO_H

#include <string>
#include <bam.h>
#include <iosfwd>

using namespace std;
typedef unsigned int uint;
namespace tkbio {
    class recfragment {
        string _name;
        string _sequence;
        string _cigar;
        //string _chromosome;
        int _chromosome;
        int _position;
        int _flag;
        vector<pair<int,char> > _mapped;
        int _position5;
        int _position3;
        int _max_match_span;
    private:
        void initialize(int position, const string& sequence, const string& cigar);
        void initialize(bam1_t* sequence);
        vector<pair<int,char> > get_positions() const;
        //static recfragment* join_pair(recfragment* f1, recfragment* f2);
        void join_sequence(const recfragment* frag);
        void generate_recombination_pattern(const vector<pair<int,char> >& loci, int* pattern) const;
        recfragment();
    public:
        recfragment(const string& chromosome, int position, int flag, const string& sequence, const string& cigar);
        recfragment(int chromosome, bam1_t* sequence);
        const string& name() const { return _name; }
        int chromosome() const { return _chromosome; }
        int position() const { return _position; }
        int span() const { return _position3 - _position5 + 1; }
        int position5() const { return _position5; }
        int position3() const { return _position3; }
        bool is_temporary() const {return _cigar == "";}
        
        //int length() const { return _length; }
        char orientation() const { return (_flag & 16) == 0 ? '+' : '-'; }
        char get_base(int position, int& count) const;
        void test_heterozygosity(int num, char const* reference, char* const alt, int* result) const;

        string get_recombination_pattern(const vector<pair<int,char> >& loci) const;
        pair<int,int> get_recombination(const vector<pair<int,char> >& loci) const;

        // join pairs
        static void bundle_pairs(vector<recfragment*>& fragments) throw (runtime_error);

        static vector<pair<int,char> > resolve_map(int position, const string& cigar, const string& sequence);
        static void detect_recombination(const vector<recfragment*>& frags,
                                         int clength, char const* cseq,
                                         int start, int end) throw (exception);
//                                         int start, int end) throw (exception);
        //friend ostream& ::operator << (ostream& ost, const recfragment& rhs);
    };

}

#endif
