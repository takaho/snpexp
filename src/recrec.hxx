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
        

        
    public:
        fragment(const string& chromosome, int position, int flag, const string& sequence, const string& cigar);

        int position() const { return _position; }
        int length() const { return _length; }
        char orientation() const { return (flag & 16) == 0 : '+' : '-'; }
        char get_base(int position) const;

        static vector<pair<int,char> > resolve_map(int position, const string& cigar, const string& sequence);

    };

}

#endif
