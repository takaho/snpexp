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

#ifndef TKTOOLS_HXX
#define TKTOOLS_HXX
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdio>
#include <set>
#include <iosfwd>

using std::vector;
using std::string;
using std::exception;
using std::runtime_error;
using std::invalid_argument;
using std::ostream;
using std::set;


namespace tktools {
    vector<string> split_items( const string& line, char delimiter = '\t' );
    string strip( const string& line, const char* characters = "\n\r" );

    namespace util {
        bool has_argument( int argc, char** argv, const char* option ) throw ( std::invalid_argument );
        bool has_option( int argc, char** argv, const char* option ) throw ( std::invalid_argument );
        int get_argument_integer( int argc, char** argv, const char* option, int defval = 0 ) throw ( std::invalid_argument );
        double get_argument_float( int argc, char** argv, const char* option, double defval = 0.0 ) throw ( std::invalid_argument );

        const char* get_argument_string( int argc, char** argv, const char* option, const char* defval = "" ) throw ( std::invalid_argument );

        string get_log_string( const char* __filename, int __linenumber );

        long get_usec();
    }

    namespace io {
        bool file_exists( const string& filename );
        vector<string> get_filenames( const string& directory ) throw ( exception );
        vector<string> get_filenames( const string& directory, const string& extension) throw ( exception );
        vector<string> get_filenames( const string& directory, const set<string>& extensions ) throw ( exception );
        string get_file_extension( const string& filename );
        size_t get_file_size( const string& filename );

        bool is_directory( const char* filepath );
        const char* file_separator();
        //bool is_file( const char* filepath );
    }

#ifdef HAVE_ZLIB_H
    namespace zip {
        // load compressed data from the current cursor of the file
        unsigned char* load_and_decompress( FILE* input, size_t size = 0 ) throw ( std::logic_error );
        //unsigned char* load_and_decompress( std::istream& input, size_t size = 0 ) throw ( std::logic_error );
        // save content in ZIP compression
        size_t compress_and_save( size_t size, unsigned char* input, FILE* output ) throw ( std::logic_error );
        //size_t compress_and_save( size_t size, unsigned char* input, std::ostream& output ) throw ( std::logic_error );


        size_t decompress_buffer( size_t size_org, unsigned char* data_cmp,
                                  unsigned char*& dst ) throw ( std::logic_error );
        size_t compress_buffer( size_t size_org, unsigned char* data_org,
                                unsigned char*& dst ) throw ( std::logic_error );
        
        //unsigned char* compress_buffer
    }
#endif

    namespace graphics {
        void save_png_graphics( char const* filename, int width, int height, char const* const* image ) throw ( exception );
    }

    namespace stat {
        typedef enum KS_SIGN {
            DPLUS = 1,
            DMINUS = -1,
            DMAX = 0
        } KS_DIRECTION;
        double get_pvalue_of_gaussian(double Z);
        double get_pvalue_of_wilcoxontest(int num, const double* sample1, const double* sample2);
        double get_pvalue_ttest( int n1, double avg1, double var1,
                                 int n2, double avg2, double var2 );
        double get_pvalue_ttest( int n1, double avg1, double var1, double mu );
        // f-test
        double get_pvalue_ftest(int n1, double var1, int n2, double var2);
        // t-test
        double get_pvalue_ttest( int num_samples, const double* data1, const double* data2 );
        // hypergeometric distribution
        double get_pvalue_exacttest( int n00, int n01, int n10, int n11 );
        // Kolmogoroc-Smirnov test
        double get_pvalue_ks_uniform(int n, double const* values, double lower, double upper);
        double get_pvalue_ks_norm(int n, double const* values, KS_SIGN = DMAX);
        double get_pvalue_ks_norm(int n, double const* values, double mean, double sigma, KS_SIGN=DMAX);
        double get_pvalue_ks_pair(int n, double const* x, int m, double const* y);
        double get_log_factorial( int n );
    }

    namespace bio {
        class fasta_sequence {
        private:
            string _name;
            int _length;
            char* _sequence;
        private:
            fasta_sequence();
            const fasta_sequence& operator = (const fasta_sequence& rhs);
            fasta_sequence(const fasta_sequence& rhs);
        public:
            //fasta_sequence(const char* filename, const char* name=NULL) throw (exception);
            ~fasta_sequence();
            int length() const { return _length; }
            const string& name() const { return _name; }
            const char* sequence() const { return _sequence; }
            string get_subsequence( int start, int end, bool complementary=false ) const;
            void set_sequence(int length, const char* sequence);
            void reverse();

            //static vector<fasta_sequence*> load_fasta_sequences(const char* filename) throw (std::exception);
            static fasta_sequence* load_file(const char* filename, const char* name=NULL) throw (std::exception);
        };

        fasta_sequence* load_fasta_sequence(const char* filename, const char* name=NULL) throw ( std::exception );

        int convert_chromosome_to_code(const char* name);
        string convert_code_to_chromosome(int code);

        // returned value is the length of converted sequence
        int convert_solid_color_to_sequence(int length, const char* codes, char* sequence, char top='\0');
        // returned value is the length of converted sequence
        int convert_sequence_to_solid_color(int length, const char* sequence, char* codes, char top='\0');
    }

//     template <typename T> pair<T,T> detect_percentile( int width, int heihgt, T const* const* values, float minimum_threshold = 0.01, float maximum_threshold = 0.99 );
//     template <typename T> pair<T,T> detect_percentile( int size, T const* values, float minimum_threshold = 0.01, float maximum_threshold = 0.99 );
}

#endif
