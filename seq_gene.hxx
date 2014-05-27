#include <iosfwd>
#include <string>
#include <vector>
#include <stdexcept>
#include <cfloat>

using namespace std;

namespace tkbio {
    class seq_gene;

    class target_probe {
        int _chromosome;
        int _position;
    public:
        target_probe(const char* chromosome, int position);
        target_probe(int chromosome, int position);
        //target_probe(const string& name, const string& chromosome, int position);
        virtual const string& id() const = 0;
        virtual ~target_probe() {}
        int chromosome_code() const { return _chromosome; }
        string chromosome() const;
        int position() const { return _position; }
        virtual void output(ostream& ost) const = 0;
    };

    class exon {
    public:
        enum Feature {CDS=1, UTR=2, OTHERS=-1};
    public:
        virtual ~exon() {};
        virtual int start() const = 0;
        virtual int end() const = 0;
        virtual int length() const = 0;
        virtual string feature() const = 0;
        friend class seq_gene;
        virtual void output(ostream& ost) const {
            ost << feature() << "[" << start() << ":" << end() << "]";
        }
        static bool compare_exon_position(const exon* lhs, const exon* rhs) {
            return lhs->start() < rhs->start();
        }
    };

        
    class singletonexon : public exon {
    private:
        Feature _feature;
        int _start;
        int _end;

        singletonexon();
    public:
        singletonexon(char feature, int start, int end);
        singletonexon(Feature feature, int start, int end);
        singletonexon(const singletonexon& rhs);

        const singletonexon& operator = (const singletonexon& rhs);

        virtual string feature() const;
        virtual int start() const { return _start; }
        virtual int end() const { return _end; }
        virtual int length() const { return _end - _start + 1; }
    };

    class featuredexon : public exon {
    private:
        int _num_exons;
        int* _borders;
        string _features;
        //exon::Feature* _features;
    private:
        featuredexon();
    public:
        featuredexon(const exon& rhs);
        featuredexon(const featuredexon& rhs);
        featuredexon& operator = (const featuredexon& rhs);
        ~featuredexon();
        virtual int start() const {return _borders[0];}
        virtual int end() const { return _borders[_num_exons]; }
        virtual int length() const { return end() - start() + 1; }
        virtual string feature() const {return _features; }
        void add(const exon& e) throw (invalid_argument);
        int count() const { return _num_exons; }
        int border(int index) const { return _borders[index]; }
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
        vector<const target_probe*> probes;
        vector<exon*> _exons;
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
        void add_exon(exon* e);
        const vector<exon*>& exons() const { return _exons; }
        int count_exons() const { return _exons.size(); }
        int count_bases() const;
        ///const exon& get_exon(int position) const;
        int get_exon_index(int pos) const;//, bool utr_sensitive=false) const;

        virtual void output(ostream& ost) const;
        virtual void serialize(ostream& ost) const;
        static seq_gene* deserialize(const string& line);

        static vector<seq_gene*> load(const char* filename, const char* assembly, bool verbose=false) throw (logic_error);
        static void insert_probe(vector<seq_gene*>& genes, const target_probe* probe, int distance_up=-10000, int distance_down=10000);

        static bool compare_position(const seq_gene* lhs, const seq_gene* rhs);
        static bool compare_tss_position(const seq_gene* lhs, const seq_gene* rhs);
            
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

    class twocolor_probe : public target_probe {
        string _id;
        float _values[4];
        static double _inv_log_2;
        static const float NO_DATA;// = -FLT_MAX;
    public:
        twocolor_probe(const string& name, const char* chromosome, int position);
        twocolor_probe(const string& name, int chromosome, int position);
        virtual ~twocolor_probe() {}
        const string& id() const { return _id; }
        float cy3value() const { return _values[0]; }
        float cy5value() const { return _values[1]; }
        float cy3call() const { return _values[2]; }
        float cy5call() const { return _values[3]; }

        float cy3value(float val) { _values[0] = val; return val; }
        float cy5value(float val) { _values[1] = val; return val; }
        float cy3call(float val)  { _values[2] = val; return val; }
        float cy5call(float val)  { _values[3] = val; return val; }
        bool has_stat() const { return _values[2] != NO_DATA && _values[3] != NO_DATA; }

        void cy3(float value, float call) { _values[0] = value; _values[1] = call; }
        void cy5(float value, float call) { _values[2] = value; _values[3] = call; }

        void output(ostream& ost) const;
        static vector<twocolor_probe*> load(const char* filename, bool verbose=false);
        static void normalize(vector<twocolor_probe*> probes);
    };

}
ostream& operator << (ostream& ost, const tkbio::seq_gene& gene);
