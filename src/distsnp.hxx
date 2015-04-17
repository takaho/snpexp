#ifndef TKBIO_DISTSNP_HXX
#define TKBIO_DISTSNP_HXX
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_BAM_H
#include <bam.h>
#endif

#include <map>
#include <string>

using std::map;
using std::string;
using std::exception;
using std::invalid_argument;
using std::runtime_error;
using std::out_of_range;
using std::logic_error;

#include <recrec.hxx>


namespace tkbio {
    class dbsnp_file;
    class dbsnp_locus;

    class cache_position {
    public:
        string chromosome;
        size_t position;
        size_t file_position;
    public:
        cache_position(const string& chromosome, size_t location, size_t file_pos) {
            this->chromosome = chromosome;
            this->position = location;
            this->file_position = file_pos;
            
        }

        int compare_position(const string& chromosome, size_t location, size_t file_pos) {
            if (this->chromosome < chromosome) {
                return -1;
            } else if (this->chromosome > chromosome) {
                return 1;
            }
            if (this->position < location) {
                return -1;
            } else if (this->position > location) {
                return 1;
            } else {
                return 0;
            }
        }
    };

    class dbsnp_locus {
        size_t _position;
        string _reference;
        string _alternative;
        int _num_strains;
        string _rsid;
        unsigned char* _strains;
    private:
        dbsnp_locus();
        const dbsnp_locus& operator = (const dbsnp_locus& rhs);
        dbsnp_locus(const dbsnp_locus& rhs);
    public:
        static const unsigned char NO_INORMATION;
        static const unsigned char REFERENCE_ALLELE;
        static const unsigned char REFERENCE_ALLELE_HOMO;
        dbsnp_locus(size_t position, const string& reference, const string& alternative);
        dbsnp_locus(size_t position, const string& rsid, const string& reference, const string& alternative, int num_strains=0);
        ~dbsnp_locus();
        void set_genotype(int strain_index, char const* info);
        unsigned char get_genotype(int strain_index) const;
        friend class dbsnp_file;
        const string& reference() const { return _reference; }
        const string& alternative() const { return _alternative; }
        const string& rsid() const { return _rsid; }
        size_t position() const { return _position; }
        string to_string() const;
        string to_string(const string& chromosome) const;
    };

    class dbsnp_file {
        string _filename;
        vector<string> _strains;
        vector<cache_position*>  _indicators;//dbsnp_locus*> > _indicators;
        vector<dbsnp_locus*> _cache;
        map<string,pair<int,int> > _cover_range;

        string _cache_chromosome;
        int _cache_range_start;
        int _cache_range_end;
    private:
        static const int CACHE_RETAIN_TOLERANCE;
        static const size_t CACHE_LINE_INTERVAL;
    private:
        static string get_cache_filename(const char* filename);
        void save_cache(const char* filename) const throw (exception);
        static dbsnp_file* load_cache(const char* filename) throw (exception);
        void add_cache(const string& chromosome, size_t pos, size_t fpos);
        void load_snps(const string& chromosome, int start, int end);
    public:
        dbsnp_file(const char* filename, const vector<string>& strains);
        ~dbsnp_file();
        void remove_strains() { _strains.erase(_strains.begin(), _strains.end()); }
        int strain_size() const { return _strains.size(); }
        const string& get_strain(int index) const { return _strains[index]; }
        const string& cached_chromosome() const { return _cache_chromosome; }
        std::pair<int,int> cached_range() const { return std::make_pair(_cache_range_start, _cache_range_end); }
        vector<dbsnp_locus const*> get_snps(string chromosome, int start, int end) const throw (exception);
        int get_strain_index(const string& strain) const;
        vector<string> get_strain() const { return _strains; }
        dbsnp_locus const* get_snp(const string& chromosome, int position) const;
        void add_file_position(const string& chrom, size_t pos, size_t fpos) {
            _indicators.push_back(new cache_position(chrom, pos, fpos));
        }
        static dbsnp_file* load_dbsnp(const char* filename, bool force_uncached=false) throw (exception);
    };

    class snp_distance {
    public:
        static void main(int argc, char** argv) throw (exception);
    };

    class gtffile;
    class denovo_snp;

    class polymorphic_allele {
        int _chromosome;
        int _position;
        int _bases[8];
        char _reference;
        char _alternative;
        static polymorphic_allele NO_DATA;
    private:
        void resolve_variation();
    public:
        polymorphic_allele(int chromosome, int position);
        polymorphic_allele(int chromosome, int position, int const* bases1, int const* bases2);
        void set_variation(char ref, char alt);
        char reference() const;
        char alternative() const;
        void set_bases(int slot, int a, int c, int g, int t) throw (out_of_range);
        int position() const { return _position; }
        void get_bases(int slot, int*& bases) const throw (out_of_range);
        int chromosome() const { return _chromosome; }
        int coverage(int slot) const throw (out_of_range);
        string to_string() const;
        float frequency(int slot) const throw (out_of_range);
        friend class denovo_snp;
        friend bool operator == (const polymorphic_allele& lhs, const polymorphic_allele& rhs);
    };

    class allele_count {
    public:
        int homo_homo;
        int homo_hetero;
        int hetero_homo;
        int hetero_hetero;
        int hetero_hetero_diff;
        int uncountable;
        allele_count() {
            homo_homo = 0;
            homo_hetero = 0;
            hetero_homo = 0;
            hetero_hetero = 0;
            hetero_hetero_diff = 0;
            uncountable = 0;
        }
        double homozygosity(bool whole=true) const;
        string to_string() const;
    };

    class denovo_snp {
        const gtffile* _genes;
        int _chromosome;
        int _start;
        int _stop;
        int _size;
        int* _position;
        int** _count1;
        int** _count2;
        int _mapped1;
        int _mapped2;
        int _reserved;
        int _quality;
    private:
        denovo_snp(const denovo_snp& rhs);
        void initialize_buffer_without_gtf(int chromosome, int start, int stop) throw (std::logic_error);
        void initialize_buffer(int chromosome, int start, int stop) throw (std::logic_error);
        //void initialize_buffer(const gtffile* gtf, int chromosome, int start, int stop) throw (std::logic_error);
        //void set_read(int slot, bam1_t const* read);
        int get_index(int position) const;
        void release_buffer();
        const denovo_snp& operator = (const denovo_snp& rhs);
    public:
        denovo_snp(const gtffile* gtf, int chromosome, int start, int stop) throw (std::logic_error);
        denovo_snp(int chromosome, int start, int stop) throw (std::logic_error);
        ~denovo_snp();
        int start() const { return _start; }
        int stop() const { return _stop; }
        int has_buffer() const { return _size > 0; }
        int mapped(int slot) const { return slot == 0 ? _mapped1 : _mapped2; }
        int chromosome() const { return _chromosome; }
        int quality() const { return _quality; }
        void set_quality(int qual);
        void add_read(int slot, bam1_t const* read) throw (out_of_range);
        void set_scope(int chromosome, int start, int stop);
        void set_scope_without_gtf(int chromosome, int start, int stop);
        
        polymorphic_allele get_allele(int position) const;
        vector<polymorphic_allele> get_polymorphism(int coverage, double hetero_threshold) const;///, double tolerance=0.0) const;
        vector<polymorphic_allele> get_polymorphism(int coverage, double hetero_threshold, int start, int stop) const;///, double tolerance=0.0) const;
        allele_count count_allele_types(int coverage, double heterozygosity, int start, int stop, chromosome_seq const* chromosome) const;
        //void add(in position, );
        //void clear();
        std::string to_string() const;
        static void enumerate_hetero(int argc, char** argv) throw (exception);
        static void detect_heterozygous(int argc, char** argv) throw (exception);
    };
}

#endif
