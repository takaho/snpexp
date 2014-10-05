#include <iosfwd>
#include <string>
#include <vector>
#include <stdexcept>
#include <cfloat>

using namespace std;

namespace tkbio {
    class seq_gene;

    class exon {
    public:
        enum Feature {CDS=1, UTR=2, OTHERS=-1};
    private:
        Feature _feature;
        int _start;
        int _end;
    public:
        exon(Feature feature, int start, int end) {
            _feature = feature;
            _start = start;
            _end = end;
        }
        exon(const exon& rhs) {
            _feature = rhs._feature;
            _start = rhs._start;
            _end = rhs._end;
        }
        const exon& operator = (const exon& rhs) {
            if (this == &rhs) return *this;
            _feature = rhs._feature;
            _start = rhs._start;
            _end = rhs._end;
            return *this;
        }
        friend class seq_gene;
        char feature() const { 
            switch (_feature) {
            case CDS: return 'C';
            case UTR: return 'U';
            case OTHERS: return 'x';
            }
            return '?';
        }
        int start() const { return _start; }
        int end() const { return _end; }
    };

    class seq_gene {
    private:
        string _id;
        string _symbol;
        int _geneid;
        int _chrm;
        int _start;
        int _end;
        char _orientation;
        //vector<const target_probe*> probes;
        vector<exon> _exons;
        bool _pseudo;
    public:
        seq_gene(const string& name, const string& symbol, int geneid);
        virtual ~seq_gene();

        void set_location(const char* chrm, char orientation, int start, int end);
        void set_location(int chrm, char orientation, int start, int end);
        void set_symbol(const string& sym);

        int get_distance(int chromosome, int position) const;
        int get_size() const { return _end - _start + 1; };
        const string& symbol() const { return _symbol; }
        int geneid() const { return _geneid; }
        string chromosome() const;
        const string& genbank() const { return _id; }
        int chromosome_code() const { return _chrm; }
        int start() const { return _start; }
        int end() const { return _end; }
        char orientation() const { return _orientation; }
        bool is_pseudo() const { return _pseudo; }
        void set_pseudo(bool flag=true) { _pseudo = flag; }
        void add_exon(exon::Feature feature, int start, int end);
        const vector<exon>& exons() const { return _exons; }
        int count_exons() const;
        int count_bases() const;
        bool contains(int chromosome_code, int start, int end) const;
        bool contains(const string& chromosome, int start, int end) const;
        ///const exon& get_exon(int position) const;
        int get_exon_index(int pos) const;//, bool utr_sensitive=false) const;

        virtual void output(ostream& ost) const;
        virtual void serialize(ostream& ost) const;
        static seq_gene* deserialize(const string& line);

        static vector<seq_gene*> load(const char* filename, const char* assembly, bool verbose=false) throw (logic_error);

        static bool compare_position(const seq_gene* lhs, const seq_gene* rhs);
        static bool compare_tss_position(const seq_gene* lhs, const seq_gene* rhs);
        static vector<const seq_gene*> find_genes(const vector<seq_gene*>& genes, const string& chromosome, int start, int end);
        static vector<const seq_gene*> find_genes(const vector<seq_gene*>& genes, int chromosome_code, int start, int end);
            
        static int get_maximum_size() { return maximum_size; }
    private:
        static int maximum_size;
    };


    class spliced_gene;
    class splice_section {
    public:
        typedef enum section_feature {unknown=-1, exon=1, intron=2, utr=3, rna=4} section_feature;
        //typedef enum feature {exon:1, intron:2, utr:3} feature;
    private:
        int _start;
        int _end;
        int _value;
        section_feature _feature;
        splice_section() {_start = _end = _value = -1; _feature = unknown; }
    public:
        splice_section(int start, int end, const char* feature) throw (invalid_argument);
        splice_section(int start, int end, section_feature feature) throw (invalid_argument);
        int start() const { return _start; }
        int end() const { return _end; }
        int value() const { return _value; }
        const char* feature() const;
        section_feature feature_code() const { return _feature; }
        void set_value(int value) { _value = value; }
        void add_value(int value) { _value += value; }
    private:
        static section_feature get_feature_code(const char* feature);
        static const char* get_feature_label(section_feature code);

        friend class spliced_gene;
        friend bool operator == (const splice_section& lhs, const splice_section& rhs);
    public:
        static splice_section out_of_range;
    };

    bool operator == (const splice_section& lhs, const splice_section& rhs);

    class spliced_gene : public seq_gene {
    private:
        vector<splice_section> _sections;
    public:
        spliced_gene(const string& name, const string& symbol, int geneid);// throw (invalid_argument) ;
        virtual ~spliced_gene() {}
        void set_section(int start, int end, const char* feature);
        void set_section(int start, int end, splice_section::section_feature feature);
        void set_introns();
        int count() const { return (int)_sections.size(); }
        const splice_section& get_feature(int index) const { return _sections[index]; }
        splice_section& get_section_at_index(int index);
        splice_section& get_section_at_position(int position);

        static vector<spliced_gene*> load(const char* filename, const char* assembly, bool verbose=false) throw (logic_error);

        virtual void output(ostream& ost) const {
            seq_gene::output(ost);
            ost << "\t" << _sections.size();
        }
        virtual void serialize(ostream& ost) const;
        static seq_gene* deserialize(const string& line);

        static int maximum_size;
    private:
//            static vector<spliced_gene*> load_cache(const char* filename_cache, const char* filename, const char* assembly) throw (logic_error);
    };
}

ostream& operator << (ostream& ost, const tkbio::seq_gene& gene);
