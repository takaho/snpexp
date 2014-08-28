#ifndef TKBIO_H
#define TKBIO_H

#include <string>

using namespace std;

namespace tkbio {

    class fragment {
        string _sequence;
        string _cigar;
        string _chromosome;
        int _position;
        int _flag;
        vector<pair<int,char> > _mapped;
        int _position5;
        int _position3;
        int _max_match_span;
    private:
        void initialize(int position, const string& sequence, const string& cigar);
        vector<pair<int,char> > get_positions() const;

    public:
        fragment(const string& chromosome, int position, int flag, const string& sequence, const string& cigar);

        int position() const { return _position; }
        int span() const { return _position3 - _position5 + 1; }
        //int length() const { return _length; }
        char orientation() const { return (_flag & 16) == 0 ? '+' : '-'; }
        char get_base(int position) const;

        static vector<pair<int,char> > resolve_map(int position, const string& cigar, const string& sequence);

    };

}

#endif
