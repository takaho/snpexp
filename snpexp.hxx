/* 
The MIT Lincense

Copyright (c) 2014 Takaho A. Endo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */

#ifndef TKBIO_SNPEXP_HXX
#define TKBIO_SNPEXP_HXX

namespace tkbio {
    class BaseFrequency {
    private:
        short* _a;
        short* _c;
        short* _g;
        short* _t;
        //int _size;
        //int _position;
        int _chid;
        string _chname;
        int _length;
        bool _reference;
    private:
        void initialize(int length);
        void clear_reference();
    public:
        BaseFrequency(int chromosome, const string& name, int length) {
            _length = length;
            _chid = chromosome;
            _chname = name;
            initialize(length);
        }
        ~BaseFrequency();
        bool has_reference() const { return _reference; }
        int length() const { return _length; }
        int chromosome_id() const { return _chid; }
        const string& chromosome_name() const { return _chname; }
        void add(int position, const bam1_t* read);
        void set_sequence(int length, const char* sequence);
        void add(int position, char base) throw (out_of_range) {
            if (position < 0 || position >= _length) {
                throw out_of_range("invalid base position");
            }
            switch (base) {
            case 'A': _a[position]++; break;
            case 'C': _c[position]++; break;
            case 'G': _g[position]++; break;
            case 'T': _t[position]++; break;
            }
        }
        bool get_basecount(int pos, int& A, int& C, int& G, int& T) const;
        void clear_data(const vector<pair<int, int> >& sorted_range) throw (exception);
        static void display_difference(BaseFrequency* b1, BaseFrequency* b2, int minimum_bases=100, double ratio_diff=0.20);
    };

    class SNPAllele {
        string _name;
        int _chromosome;
        int _position;
        string _reference;
        string _alternative;
        string _variation;
        SNPAllele();
        const SNPAllele& operator = (const SNPAllele& );
        SNPAllele(const SNPAllele&);
    public:
        SNPAllele(const string& name, const string& chromosome, int location, const string& reference, const string& alternative);
        SNPAllele(const string& name, const string& chromosome, int location);
        void set_variation(const string& variation);
        const char* variation() const { return _variation.c_str(); }
        const char* reference() const { return _reference.c_str(); }
        const char* alternative() const { return _alternative.c_str(); }
        const string& id() const { return _name; }
        int position() const { return _position; }
        int chromosome_code() const { return _chromosome; }
        string chromosome() const;
        static bool compare_position(const SNPAllele* lhs, const SNPAllele* rhs);
        static vector<SNPAllele*> load_vcf(const char* filename, const vector<string>& strains) throw (exception);
    };
}

#endif
