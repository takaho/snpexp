#ifndef TKBIO_H
#define TKBIO_H

#include <string>
#include <bam.h>

using namespace std;
typedef unsigned int uint;
namespace tkbio {
    class recfragment {
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
        void initialize(bam1_t const* sequence);
        vector<pair<int,char> > get_positions() const;

    public:
        recfragment(const string& chromosome, int position, int flag, const string& sequence, const string& cigar);
        recfragment(bam1_t const* sequence);
        int chromosome() const { return _chromosome; }
        int position() const { return _position; }
        int span() const { return _position3 - _position5 + 1; }
        int position5() const { return _position5; }
        int position3() const { return _position3; }
        //int length() const { return _length; }
        char orientation() const { return (_flag & 16) == 0 ? '+' : '-'; }
        char get_base(int position) const;

        static vector<pair<int,char> > resolve_map(int position, const string& cigar, const string& sequence);
        static void detect recombination(const vector<recfragment*>& frags,
                                         int start, int end) throw (exception);
    };

}

#endif
