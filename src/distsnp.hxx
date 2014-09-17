#ifndef TKBIO_DISTSNP_HXX
#define TKBIO_DISTSNP_HXX

#include <map>
#include <string>

using std::map;
using std::string;
using std::exception;
using std::invalid_argument;
using std::runtime_error;

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
        unsigned char* _strains;
    public:
        dbsnp_locus(size_t position, string reference, string alternative, int num_strains=1);
        ~dbsnp_locus();
        void set_genotype(int strain_index, char const* info);
        friend class dbsnp_file;
        const string& reference() const { return _reference; }
        const string& alternative() const { return _alternative; }
        size_t position() const { return _position; }
        string to_string() const;
        string to_string(const string& chromosome) const;
    };

    class dbsnp_file {
        string _filename;
        vector<string> _strains;
        vector<cache_position*>  _indicators;//dbsnp_locus*> > _indicators;
        vector<dbsnp_locus*> _cache;
        string _current_chromosome;
    private:
        static string get_cache_filename(const char* filename);
        void save_cache(const char* filename) const throw (exception);
        static dbsnp_file* load_cache(const char* filename) throw (exception);
        void add_cache(const string& chromosome, size_t pos, size_t fpos);
        void load_snps(const string& chromosome, int start, int end);
    public:
        dbsnp_file(const char* filename, const vector<string>& strains);
        ~dbsnp_file();
        int strain_number() const { return _strains.size(); }
        const string& get_strain(int index) const { return _strains[index]; }
        vector<dbsnp_locus const*> get_snps(string chromosome, int start, int end) const throw (exception);
        void add_file_position(const string& chrom, size_t pos, size_t fpos) {
            _indicators.push_back(new cache_position(chrom, pos, fpos));
        }
        static dbsnp_file* load_dbsnp(const char* filename, bool force_uncached=false) throw (exception);
    };

    // class snp_result {
    //     string _chromosome;
    //     int _length;
    //     vector<string> _strains;
    //     int ** _
    // };
}

#endif
