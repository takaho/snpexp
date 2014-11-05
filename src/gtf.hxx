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

#ifndef GTF_H
#define GTF_H

#include <string>
#include <vector>
#include <stdexcept>
#include <iosfwd>
#include <map>

using namespace std;
using std::string;
using std::map;

namespace tkbio {
    class gtfexon;
    class gtfgene;
    class geffile;

    /// Exon 
    ///
    class gtfexon {
    public:
        enum Feature {UNKNOWN=-1, OTHERS=0, EXON=1, CDS=2, START_CODON=3, STOP_CODON=4, INTRON=5};
    private:
        static string* _features;
        static const int _num_features;
        int _pos5;
        int _pos3;
        Feature _feature;
    public:
        gtfexon(int pos5, int pos3, const char* feature) {
            _pos5 = pos5;
            _pos3 = pos3;
            _feature = encode_feature(feature);
        }
        gtfexon(int pos5, int pos3, Feature feat) {
            _pos5 = pos5;
            _pos3 = pos3;
            _feature = feat;
        }
        static const string& get_feature(Feature _feature); // feature of the exon
        static Feature encode_feature(const char* feature); // convert label to feature in integer
        int position5() const { return _pos5; } // minimum limit of position
        int position3() const { return _pos3; } // maximum limit of position
        pair<int,int> position() const { return make_pair(_pos5, _pos3); } // range
        Feature feature() const { return _feature; } // feature in integer
        const string& feature_name() const { return get_feature(_feature); }
        friend class gtfgene;
        friend class gtffile;
    };

    /// Gene
    class gtfgene {
    private:
        string _name; // name of the gene
        string _transcript_id; // identifier of each transcripts, typically GenBank ID
        string _tss_id; // identifier of transcription start site
        string _chromosome; // chromosome 
        int _chromosome_code;
        char _orientation; // orientation
        int _position5; // minimum limit of position
        int _position3; // maximum limit of position
        vector<gtfexon> _exons; // exons
        
        gtfgene(const string& transcript_id, const string& name, const string& tss_id);
    private:
        void determine_span();
    public:
        const string& chromosome() const { return _chromosome; }
        char orientation() const { return _orientation; }
        void insert(const string& chromosome, char orientation, int pos5, int pos3, const char* feature);
        void insert(int pos5, int pos3, const char* feature);
        const vector<gtfexon>& exons() const { return _exons; }
        const string& transcript_id() const { return _transcript_id; }
        const string& name() const { return _name; }
        const string& tss_id() const { return _tss_id; }
        bool contains(const string& chromosome, int location) const;
        bool contains_in_exon(const string& chromosome, int location) const;
        bool contains_in_exon(int chromosome_code, int location) const;
        bool contains_in_exon_debug(int chromosome_code, int location) const;
        bool contains_in_exon(int chromosome_code, int start, int end) const;
        bool contains_in_exon_debug(int chromosome_code, int start, int end) const;
        pair<int,int> get_span() const;

        string to_string() const;

        static bool compare_position(const gtfgene* lhs, const gtfgene* rhs);

        friend class gtffile;
    };

    // Genes
    class gtffile {
        map<string, gtfgene*> _genes; // genbank to gene, unique
        map<string, vector<gtfgene*> > _tssid; // tss_id to genbank
        map<string, vector<gtfgene*> > _symbol; // symbol to genbank
        vector<gtfgene*> _sorted_genes; // sorted genes by their positions
        gtffile(const gtffile& rhs);
        const gtffile& operator = (const gtffile& rhs);
        bool _sorted;
        int _maximum_gene_span;
    private:
        void sort_genes();
    public:
        gtffile();
        ~gtffile();
        void add(const char* transcript_id, const char* name, const char* tss_id, const char* feature, const char* chromosome, char orientation, int position5, int position3) throw (exception);
        int size() const { return _genes.size(); }
        vector<gtfgene*> get_tss() const;
        vector<gtfgene*> get_genes() const;
        vector<string> get_symbols() const;
        vector<const gtfgene*> find_genes(const string& chromosome, int start, int end=-1) const;
        vector<const gtfgene*> find_genes(int chromosome, int start, int end) const;
        bool contains_in_exon(const string& chromosome, int location) const;
        bool contains_in_exon_debug(const string& chromosome, int location) const;
        vector<pair<int,int> > get_exonregion(const string& chromosome, bool utr=true) const;
        string to_string() const;

        static gtffile* load_gtf(const char* fileanme) throw (exception);
    };
}

#endif
