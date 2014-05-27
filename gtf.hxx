#ifndef GTF_H
#define GTF_H

#include <string>
#include <vector>
#include <stdexcept>
#include <iosfwd>
#include <map>

using namespace std;

namespace tkbio {
    class gtfexon;
    class gtfgene;
    class geffile;

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
        static const string& get_feature(Feature _feature);
        static Feature encode_feature(const char* feature);
        int position5() const { return _pos5; }
        int position3() const { return _pos3; }
        pair<int,int> position() const { return make_pair(_pos5, _pos3); }
        Feature feature() const { return _feature; }
        const string& feature_name() const { return get_feature(_feature); }
        friend class gtfgene;
        friend class gtffile;
    };


    class gtfgene {
    private:
        string _name;
        string _transcript_id;
        string _tss_id;
        string _chromosome;
        char _orientation;
        int _position5;
        int _position3;
        vector<gtfexon> _exons;
        
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
        pair<int,int> get_span() const;
        friend class gtffile;
        static bool compare_position(const gtfgene* lhs, const gtfgene* rhs);
    };

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
        // list of representative genes having tss
        vector<gtfgene*> get_tss() const;
        vector<gtfgene*> get_genes() const;
        vector<string> get_symbols() const;
        vector<const gtfgene*> find_genes(string chromosome, int start, int end=-1) const;
        vector<pair<int,int> > get_exonregion(string chromosome, bool utr=true) const;
        
        static gtffile* load(const char* fileanme) throw (exception);
    };
}

#endif
