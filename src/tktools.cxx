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

#include <fstream>
#include <iostream>
#include <string>
#include <dirent.h>
#include <cmath>
#include <cstring>
#include <memory>
#include <cstdlib>
#include <algorithm>

#include "tktools.hxx"

#ifdef HAVE_ZLIB_H
#include <zlib.h>
#endif

using namespace std;
using namespace tktools;

#ifdef _WIN32
const char* tktools::io::file_separator() { return "\\"; } 
#else
const char* tktools::io::file_separator() { return "/"; } 
#endif

string tktools::strip( const string& line, const char* characters ) {
    int len = line.size();
    int start = len;
    int end = 0;
    for ( int i = 0; i < len; i++ ) {
        char c = line.c_str()[ i ];
        if ( strchr( characters, c ) == NULL ) {
            start = i;
            break;
        }
    }
    for ( int i = len - 1; i >= start; i-- ) {
        char c = line.c_str()[ i ];
        if ( strchr( characters, c ) == NULL ) {
            end = i + 1;
            break;
        }
    }
    return line.substr( start, end - start );
}

vector<string> tktools::split_items( const string& line, char delimiter ) {
    vector <string> items;
    const char* ptr = line.c_str();
    int size = line.size();
    int pos_prev = 0;
    for ( int i = 0; i <= size; i++ ) {
        char c = ptr[ i ];
        if ( c == delimiter || c == '\0' ) {
            items.push_back( line.substr( pos_prev, i - pos_prev ) );
            pos_prev = i + 1;
        }
    }
    return items;
}

size_t tktools::io::get_file_size( const string& filename ) {
    ifstream file_in( filename.c_str() );
    if ( file_in.is_open() ) {
        file_in.seekg( 0, ios::end );
        size_t end = file_in.tellg();
        file_in.close();
        return end;
    }
    return 0;
}

bool tktools::io::file_exists( const string& filename ) {
    ifstream file_in( filename.c_str() );
    if ( file_in.is_open() ) {
        file_in.close();
        return true;
    }
    return false;
}

string tktools::io::get_file_extension( const string& filename ) {
    int slen = strlen( filename.c_str() );
    for ( int i = slen - 1; i >= 0; i-- ) {
        if ( filename[ i ] == '.' ) {
            return string( filename.c_str() + i + 1 );
        }
    }
    return string( "" );
}

bool tktools::io::is_directory( const char* filepath ) {
    DIR* dir = opendir( filepath );
    if ( dir == NULL ) return false;
    closedir( dir );
    return true;
}
vector<string> tktools::io::get_filenames( const string& directory, const string& extension) throw ( exception ) {
    set<string> exts;
    exts.insert(extension);
    return tktools::io::get_filenames(directory, exts);
}

vector<string> tktools::io::get_filenames( const string& directory, const set<string>& extensions ) throw ( exception ) {
    vector<string> filenames;
    DIR* dir = opendir( directory.c_str() );
    dirent* dent;
    if ( dir == NULL ) {
        throw logic_error( string( "cannot open " ) + directory );
    }
    while ( ( dent = readdir( dir ) ) != NULL ) {
        const char* name = dent->d_name;
        if ( name[ 0 ] == '.' ) continue;
        string ext = tktools::io::get_file_extension( name );
        if ( extensions.size() == 0 || extensions.find( ext ) != extensions.end() ){
            //if ( ext != "fa" && ext != "txt" && ext != "fasta" ) continue;
            string filepath = string( directory ) + string( "/" ) + string( dent->d_name );
            filenames.push_back( filepath );
        }
    }
    closedir( dir );
    return filenames;
}

vector<string> tktools::io::get_filenames( const string& directory ) throw ( exception ) {
    vector<string> filenames;
    DIR* dir = opendir( directory.c_str() );
    dirent* dent;
    if ( dir == NULL ) {
        throw logic_error( string( "cannot open " ) + directory );
    }
    while ( ( dent = readdir( dir ) ) != NULL ) {
        const char* name = dent->d_name;
        if ( name[ 0 ] == '.' ) continue;
        //string ext = tktools::io::get_file_extension( name );
        //if ( extensions.size() == 0 || extensions.find( ext ) != extensions.end() ){
            //if ( ext != "fa" && ext != "txt" && ext != "fasta" ) continue;
        string filepath = string( directory ) + string( "/" ) + string( dent->d_name );
        filenames.push_back( filepath );
        //}
    }
    closedir( dir );
    return filenames;
}

//#define HAVE_PNG_H
#ifdef HAVE_PNG_H
#define USE_PNG
#endif
#ifdef HAVE__USR_X11_INCLUDE_LIBPNG12_PNG_H
//#define USE_PNG
#endif
#ifdef USE_PNG
#include <png.h>

void tktools::graphics::save_png_graphics( char const* filename, int width, int height, char const* const* image ) throw ( exception ) {
    FILE* fp;
    png_structp png_ptr;
    png_infop info_ptr;

    fp = fopen( filename, "wb" );
    if ( !fp ) {
        cerr << __FILE__ << " " << __LINE__ << endl;
        throw invalid_argument( string( "cannot open " ) + string( filename ) );
    }
    png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
    if ( !png_ptr ) {
        cerr << __FILE__ << " " << __LINE__ << endl;
        throw logic_error( "cannt instanciate png_ptr" );
    }
    info_ptr = png_create_info_struct( png_ptr );
    if ( !info_ptr ) {
        cerr << __FILE__ << " " << __LINE__ << endl;
        throw logic_error( "cannt instanciate png_infop" );
    }
    png_init_io( png_ptr, fp );
    png_set_IHDR( png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT );
    png_write_info( png_ptr, info_ptr );
    png_write_image( png_ptr, (png_byte**)image );
    png_write_end( png_ptr, info_ptr );
    png_destroy_write_struct( &png_ptr, &info_ptr );
    fclose( fp );
    return;
}
#endif


namespace {
    double* _exacttest_cache_log          = NULL;
    unsigned int _exacttest_cache_log_reserved = 0;
    unsigned int _exacttest_cache_log_size     = 0;
    unsigned int _exacttest_cache_log_max      = 512;
    
    static double _coeff_log_gamma[ 6 ] = { 
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5
    };
    
    
    double log_gamma( double x ) {
        int i;
        double y, tmp, ser;
        y = x;
        tmp = x + 5.5;
        tmp -= ( x + .5 ) * log( tmp );
        ser = 1.000000000190015;
        for ( i = 0; i <= 5; i++ ) {
            ser += _coeff_log_gamma[ i ] /(++y);
        }
        return -tmp + log( 2.5066282746310005 * ser / x );
    }

    double log_factorial( unsigned int n ) {
        unsigned int i;
        if ( _exacttest_cache_log == 0 ) {
            _exacttest_cache_log_reserved = 16;
            _exacttest_cache_log = (double*)malloc( sizeof( double ) * _exacttest_cache_log_reserved );
            _exacttest_cache_log_size = 2;
            _exacttest_cache_log[ 0 ] = 0.0;
            _exacttest_cache_log[ 1 ] = 0.0;
        }
        /* return cached value if stored */
        if ( n < _exacttest_cache_log_size ) {
            return _exacttest_cache_log[ n ];
        }
        /* calculate approximation if the number is large */
        if ( n > _exacttest_cache_log_max ) {
            return log_gamma( n + 1 );
        }
        /* expand buffer size */
        if ( n >= _exacttest_cache_log_reserved ) {
            unsigned int size_exp = _exacttest_cache_log_size << 1;
            double* buf;
            while ( size_exp <= n ) {
                size_exp <<= 1;
            }
            buf = (double*)malloc( sizeof( double ) * size_exp );
            memcpy( buf, _exacttest_cache_log, 
                    sizeof( double ) * _exacttest_cache_log_reserved );
            free( _exacttest_cache_log );
            _exacttest_cache_log = buf;
            _exacttest_cache_log_reserved = size_exp;
        }
        /* Fill cache */
        {
            double di = (double) ( _exacttest_cache_log_size + 1 );
            for ( i = _exacttest_cache_log_size; i <= n; i++ ) {
                double lg = log_gamma( di ++ );
                _exacttest_cache_log[ i ] = lg;
            }
            _exacttest_cache_log_size = n + 1;
        }
        return _exacttest_cache_log[ n ];
    }
    
    double SQR( double a ) {
        return a == 0.0 ? 0.0 : a * a;
    }


    double betacf( double a, double b, double x ) {
        double qap, qam, qab, em, tem, d;
        double bz, bm=1.0, bp, bpp;
        double az=1.0,am=1.0, ap, app, aold;
        int m;
        double md;
        
        qab = a + b;
        qap = a + 1.0;
        qam = a - 1.0;
        bz  = 1.0 - qab * x / qap;
        md = 1.0;
        for ( m = 1; m <= 100; m++, md += 1.0 ) {/* 100 is the maximum count of iteration */
            em = md;
            tem = em + em;
            d = em * (b-em) * x / ( ( qam + tem ) * ( a + tem ) );
            ap = az + d * am;
            bp = bz + d * bm;
            d = - ( a + em ) * ( qab + em ) * x / ( ( qap + tem ) * ( a + tem ) );
            app = ap + d * az;
            bpp = bp + d * bz;
            aold = az;
            am = ap / bpp;
            bm = bp / bpp;
            az = app /bpp;
            bz = 1.0;
            /* 1e-7 is the cutoff threshold to finish approximation */
            if ( fabs( az - aold ) < ( 1e-7 * fabs( az ) ) ) return az;
        }
        fprintf( stderr, "a or b too big, or the number of iteration is too small in %s:%d", __FILE__, __LINE__ );
        return az;
        
    }

    double betai( double a, double b, double x ) {
        double bt;
        if ( x < 0.0 || x > 1.0 ) {
            fprintf( stderr, "invalid value %f in %s:%d\n", x, __FILE__, __LINE__ );
            return 0.0;
        }
        if ( x == 0.0 || x == 1.0 ) {
            bt = 0.0;
        } else {
            bt = exp( log_gamma( a + b ) - log_gamma( a ) - log_gamma( b ) + a * log( x ) + b * log( 1.0 - x ) );
        }
        if ( x < ( a + 1.0 ) / ( a + b + 2.0 ) ) {
            return bt * betacf( a, b, x ) /a;
        } else {
            return 1.0 - bt * betacf( b, a, 1.0 - x ) / b;
        }
    }

}

double tktools::stat::get_log_factorial( int n ) {
    if ( n <= 0 ) return 0.0;
    return log_factorial( (unsigned int)n );
}


double tktools::stat::get_pvalue_of_gaussian(double Z) {
    double pvalue = 1.0 - 0.5 * ( 1.0 + erf( Z * .70710678118654752440 ) );
    return pvalue;
}

namespace {
/** Quick sort for U-test
 */
void qsort_pair( int left, int right, double* values, int* flags ) {
    int scan_left, scan_right;
    double pivot;
    if ( right - left <= 1 ) return;

    scan_left  = left;
    scan_right = right - 1;
    pivot = values[ ( left + right ) >> 1 ];

    while ( scan_left <= scan_right ) {
        while ( scan_left <= right && values[ scan_left ] < pivot ) scan_left ++;
        while ( scan_right >= left && values[ scan_right ] > pivot ) scan_right --;
        if ( scan_left > scan_right ) break;
        { /* swap */
            double v = values[ scan_left ];
            int    b = flags[ scan_left ];
            values[ scan_left ] = values[ scan_right ];
            values[ scan_right ] = v;
            flags[ scan_left ] = flags[ scan_right ];
            flags[ scan_right ] = b;
        }
        scan_left ++;
        scan_right --;
    }
    if ( scan_left >= right ) return;
    qsort_pair( left, scan_left, values, flags );
    qsort_pair( scan_left, right, values, flags );
}
}

double tktools::stat::get_pvalue_of_wilcoxontest(int num, const double* sample1, const double* sample2) {
    int i;
    int num_effectives = 0;
    double pvalue;
    double  tvalue = 0.0;
    double* diffs = (double*)malloc( num * sizeof( double ) );
    double* ranks = (double*)malloc( num * sizeof( double ) );
    int*    flags = (int*)malloc( num * sizeof( int ) );
    for ( i = 0; i < num; i++ ) {
        double dif;
        dif = sample1[ i ] - sample2[ i ];
        if ( dif == 0.0 ) continue;
        diffs[ num_effectives ] = fabs( dif );
        flags[ num_effectives ] = dif > 0 ? 1 : 0;
        num_effectives ++;
    }
    qsort_pair( 0, num_effectives, diffs, flags );
    i = 0;
    while ( i < num_effectives ) {
        int j;
        int end = num_effectives;
        double rank = 0.0;
        double rank_accum = (double)( i + 1 );
        
        for ( j = i + 1; j < num_effectives; j++ ) {
            if ( diffs[ i ] != diffs[ j ] ) {
                end = j;
                break;
            }
            rank_accum += (double)( j + 1 );
        }
        rank = rank_accum / ( end - i );
        for ( j = i; j < end; j++ ) ranks[ j ] = rank;
        i = end;
    }
    for ( i = 0; i < num_effectives; i++ ) {
/*         fprintf( stderr, "%d\t%.2f\t%s%.2f\n", i, ranks[ i ], flags[ i ] ? "+" : "-", diffs[ i ] ); */
        if ( flags[ i ] != 0 ) {
            tvalue += ranks[ i ];
        }
    }
    {
        double avg = ( num_effectives * ( num_effectives + 1 ) ) * 0.25;
        double sd;

        if ( num_effectives < 1000 ) {
            sd = sqrt( (double)( num_effectives * ( num_effectives + 1 ) * ( 2 * num_effectives + 1 ) ) / 24.0 );
        } else {
            double ne = (double)num_effectives;
            sd = sqrt( ne * ( ne + 1 ) * ( 2 * ne + 1 ) / 24.0 );
        }
        //fprintf( stderr, "N = %d\n", num_effectives );

        double z = ( tvalue - avg ) / sd;
        //fprintf( stderr, "t:%.2f\tavg:%.2f\tsd:%.2f\tz=%.4f\n", tvalue, avg, sd, z );
        pvalue = get_pvalue_of_gaussian( z );
        if ( z < 0.0 ) {
            pvalue = ( 1.0 - pvalue ) * 2.0;
        } else {
            pvalue *= 2.0;
        }
    }

    free( diffs );
    free( ranks );
    free( flags );
    return pvalue;
}

double tktools::stat::get_pvalue_ftest(int n1, double var1, int n2, double var2) {
    double pvalue;
    int df1, df2;
    double f;
    if (var1 > var2) {
        f = var1 / var2;
        df1 = n1 - 1;
        df2 = n2 - 1;
    } else {
        f = var2 / var1;
        df1 = n2 - 1;
        df2 = n1 - 1;
    }
    pvalue = 2.0 * betai(0.5 * df2, 0.5 * df1, df2 / (df2 + df1 * f));
    if (pvalue > 1.0) pvalue = 2.0 - pvalue;
    return pvalue;
}

double tktools::stat::get_pvalue_ttest( int n1, double avg1, double var1,
                                         int n2, double avg2, double var2 ) {
    double v1   = var1 / n1;
    double v2   = var2 / n2;
    double svar = sqrt( v1 + v2 );
    double t    = fabs( avg1 - avg2 ) / svar;
    /* degree of freedom */
    double df = SQR( v1 + v2 ) 
        / ( SQR( var1 / n1 ) / ( n1 - 1 ) + SQR( var2 / n2 ) / ( n2 - 1 ) );
    double prob = betai( 0.5 * df, 0.5, df / ( df + SQR( t ) ) );
#ifdef DEBUG
    fprintf( stderr, "Sample1 ( %d, %f, %f ) vs Sample2 ( %d, %f, %f )\n", N1, avg1, var1, N2, avg2, var2 ); 
    fprintf( stderr, "t0 = %.4f, df = %.4f,  p => %.3e: %d\n", t, df, prob, __LINE__ );
#endif
    return prob;
}

double tktools::stat::get_pvalue_ttest( int n, double avg, double var, double mean ) {
    double t = ( avg - mean ) / ( sqrt( var / ( n - 1 ) ) );
    double df = n - 1;
    double p = betai( 0.5 * df, 0.5, df / ( df + t * t ) );
    return p;
}

double tktools::stat::get_pvalue_exacttest( int n00, int n01, int n10, int n11 ) {
    int div = n01;
    double prob_base = log_factorial( n00 + n01 ) + log_factorial( n10 + n11 ) 
        + log_factorial( n00 + n10 ) + log_factorial( n01 + n11 )
        - log_factorial( n00 + n01 + n10 + n11 );
    double prob_start = prob_base 
        - log_factorial( n00 ) - log_factorial( n01 ) 
        - log_factorial( n10 ) - log_factorial( n11 );
    double pvalue = exp( prob_start );

#ifdef DEBUG
    printf( "%3d  %3d : %3d \n%3d  %3d : %3d\n--------------\n%3d  %3d : %3d\n\n", 
            n00, n01, ( n00 + n01 ), n10, n11, ( n10 + n11 ), ( n00 + n10 ), ( n01 + n11 ),
            ( n00 + n01 + n10 + n11 ) );
#endif    
    n00 += div;
    n01 -= div;
    n10 -= div;
    n11 += div;

    if ( n01 < 0 || n10 < 0 ) {
        div = n01 < n10 ? n01 : n10;
        n00 += div;
        n10 -= div;
        n01 -= div;
        n11 += div;
    }

    while ( n00 >= 0 && n11 >= 0 ) {
        double p = prob_base
            - log_factorial( n00-- ) - log_factorial( n01++ ) 
            - log_factorial( n10++ ) - log_factorial( n11-- );
        if ( p < prob_start ) {
            pvalue += exp( p );
        }
    }
#ifdef DEBUG
    printf( "P value = %.5e\n", pvalue );
#endif
    return pvalue;
}


double tktools::stat::get_pvalue_ttest( int N, const double* values1, const double* values2 ) {
    double sum_d = 0.0;
    double sum_d2 = 0.0;
    int i;
    double p;
    for ( i = 0; i < N; i++ ) {
        double d = values1[ i ] - values2[ i ];
        sum_d  += d;
        sum_d2 += d * d;
    }
    {
        double avg_d = sum_d / N;
        double vd = ( sum_d2 - sum_d * avg_d );
        double t0 = fabs( avg_d ) / sqrt( vd / ( N * ( N - 1 )  ) );
        double df = N - 1;
        p = betai( 0.5 * df, 0.5, df / ( df + t0 * t0 ) );
    }
    return p;
}

namespace {
    const int MAX_LOOP = 4096;
    double get_pvalue_ks( double K ) {
//        double pvalue = 1.0 - 0.5 * ( 1.0 + erf( Z * .70710678118654752440 ) );
//        return pvalue;
        if (K <= 0.0 ) return 1.0;
        int j = 1;
        double accum_p = 0.0;
        double e2_ = - 2 * K * K;
        while (true) {
            double p = ((j & 1) == 0) ? -exp(e2_ * j * j) : exp(e2_ * j * j);
            accum_p += p;
            if (j > 100 && (fabs(p) < 1e-200 || j > MAX_LOOP)) break;
            j ++;
        }
        double pvalue = 2.0 * accum_p;
        return pvalue;
    }
}

double tktools::stat::get_pvalue_ks_uniform(int n, double const* values, double lower, double upper) {
    double* x = new double[n];
    memcpy(x, values, sizeof(double) * n);
    sort(x, x + n);
    double dmax = 0;
    for ( int i = 0; i < n; i++ ) {
        double v = x[i];
        if (i < n - 1 && x[i + 1] == v) continue;
        double d = fabs( (double)(i + 1) / n - (v - lower) / (upper - lower) );
        //cout << "i=" << i << ", F=" << (double)(i + 1) / n << ", f=" << (v-lower)/(upper-lower) << ", d=" << d << endl;
        if ( d > dmax ) dmax = d;
    }
    delete[] x;
    //cerr << "D=" << dmax << endl;
    return get_pvalue_ks( dmax * sqrt((double) n ) );
}

double tktools::stat::get_pvalue_ks_norm(int n, double const* values, KS_SIGN sign) {
    double s1 = 0.0;
    double s2 = 0.0;
    for ( int i = 0; i < n; i++ ) {
        double v = values[i];
        s1 += v;
        s2 += v * v;
    }
    double avg = s1 / n;
    return get_pvalue_ks_norm(n, values, avg, sqrt((s2 - s1 * avg) / (n - 1)), sign);
}

double tktools::stat::get_pvalue_ks_norm(int n, double const* values, double mean, double sigma, KS_SIGN sign) {
    double* x = new double[n];
    memcpy(x, values, sizeof(double) * n);
    sort(x, x + n);
    double dmax = 0.0;
    double dmin = 0.0;
    double accum = 0.0;
    double div = 1.0 / n;
    double coeff = .70710678118654752440 / sigma;
    double K;
    //int num_samples = 0;
    for ( int i = 0; i < n; i++ ) {
        accum += div;
        double d = accum - 0.5 * (1.0 + erf( (x[i] - mean) * coeff ) );
        if (d > dmax) dmax = d;
        if (d < dmin) dmin = d;
    }
    delete[] x;
    if (sign == DPLUS) {
        K = dmax * sqrt(n);
        return exp( -2 * K * K);
    } else if (sign == DMINUS) {
        K = -dmin * sqrt(n);
        return exp( -2 * K * K);
    } else {
        K = (dmax > -dmin ? dmax : -dmin) * sqrt(n);
        cerr << "n=" << n << ", u=" << mean << ", sd=" << sigma << ", D+=" << dmax << ", D-=" << (-dmin) << ", K=" << K << ", " << endl;
    return get_pvalue_ks( K );
    }
}

double tktools::stat::get_pvalue_ks_pair(int n, double const* x, int m, double const* y) {
    double* X = new double[n];
    double* Y = new double[m];
    memcpy(X, x, sizeof(double) * n);
    memcpy(Y, y, sizeof(double) * m);
    sort(X, X + n);
    sort(Y, Y + m);
    double dmax = 0.0;
    double dn = 1.0 / (double)n;
    double dm = 1.0 / (double)m;
    {
        int i = 0;
        int j = 0;
        double xi = X[0];
        double yj = Y[0];
        for (;;) {
            if (xi > yj) {
                double t0 = (double)(i + 1) * dn;
                while (j < m) {
                    yj = Y[j];
                    if (xi <= yj) break;
                    double d = fabs(t0 - (double)(j + 1) * dm);
                    if (d > dmax) dmax = d;
                    j++;
                }
                if (j >= m) break;
            } else {
                double t1 = (double)(j + 1) * dm;
                while (i < n) {
                    xi = X[i];
                    if (xi > yj) break;
                    double d = fabs((double)(i + 1) * dn - t1);
                    if (d > dmax) dmax = d;
                    i++;
                }
                if (i >= n) break;
            }
        }
    }
    delete[] X;
    delete[] Y;
    double K = dmax * sqrt( m * n / (m + n) );
    return get_pvalue_ks( K );
}

namespace {
    bool option_duplicated( int argc, char** argv, const char* option = NULL ) {
        if ( option == NULL || strlen( option ) == 0) {//<= 1 ) {
            for ( int i = 1; i < argc; i++ ) {
                if ( argv[ i ][ 0 ] == '-' ) {
                    for ( int j = i + 1; j < argc; j++ ) {
                        if ( strcmp( argv[ i ], argv[ j ] ) == 0 ) {
                            return true;
                        }
                    }
                }
            }
            return false;
        } else {
            for ( int i = 1; i < argc; i++ ) {
                if ( argv[ i ][ 0 ] == '-' && strncmp( argv[ i ] + 1, option, strlen( option ) ) == 0  ) {
                    for ( int j = i + 1; j < argc; j++ ) {
                        if ( argv[ j ][ 0 ] == '-' && strncmp( argv[ j ] + 1, option, strlen( option ) ) == 0  ) {
                            return true;
                        }
                    }
                }
            }
            return false;
        }
    }

    const char* get_argument_pointer( int argc, char** argv, const char* option ) throw ( invalid_argument ) {
        if ( option_duplicated( argc, argv, option ) == true ) {
            throw invalid_argument( string( "option " ) + string( option ) + string( " duplicated" ) );
        }
        if ( option == NULL ) {
            return argv[ 1 ];
        }
        if ( strlen( option ) == 1 ) {
            for ( int i = 1; i < argc; i++ ) {
                //cout << option << " " << i << " : " << argv[ i ] << endl;
                if ( argv[ i ][ 0 ] == '-' ) {
                    if ( argv[ i ][ 1 ] == option[ 0 ] ) {
                        if ( strlen( argv[ i ] ) > 2 ) {
                            return argv[ i ] + 2;
                        } else if ( i < argc - 1 ) {
                            //cout << "return " << argv[ i + 1 ] << " " << __FILE__ << " " << __LINE__ << endl;
                            return argv[ i + 1 ];
                        } else {
                            throw invalid_argument( string( option ) + " does not have value" );
                        }
                    }
                }
            }
            //cout << "not found " << endl;
            return NULL;
        } else {
            for ( int i = 1; i < argc; i++ ) {
                if ( strcmp( argv[ i ] + 1, option ) == 0 ) {
                    if ( i < argc - 1 ) {
                        return argv[ i + 1 ];
                    } else {
                        throw invalid_argument( string( option ) + " does not have value" );
                    }
                }
            }
            return NULL;
        }
    }
}

bool tktools::util::has_argument( int argc, char** argv, const char* option ) throw ( invalid_argument ) {
    if ( option_duplicated( argc, argv, option ) ) {
        throw invalid_argument( string( "option " ) + string( option ) + string( " duplicated" ) );
    }
    if ( strlen( option ) == 1 ) {
        for ( int i = 1; i < argc; i++ ) {
            if ( argv[ i ][ 0 ] == '-' && argv[ i ][ 1 ] == option[ 0 ] ) {
                return true;
            }
        }
    } else {
        for ( int i = 1; i < argc; i++ ) {
            if ( strcmp( option, argv[ i ] + 1 ) == 0 ) {
                return true;
            }
        }
    }
    return false;
}

bool tktools::util::has_option( int argc, char** argv, const char* option ) throw ( invalid_argument ) {
    if ( option_duplicated( argc, argv, option ) ) {
        throw invalid_argument( string( "option " ) + string( option ) + string( " duplicated" ) );
    }
    for ( int i = 1; i < argc; i++ ) {
        if ( strcmp( argv[ i ] + 1, option ) == 0 ) {
            return true;
        }
    }
    return false;
}

int tktools::util::get_argument_integer( int argc, char** argv, const char* option, int defval )  throw ( invalid_argument ) {
    try {
        const char* ptr = get_argument_pointer( argc, argv, option );
        //cout << "PTR = " << ( ptr == NULL ) << " : " << (void*)ptr << " , " << defval << endl;
        if ( ptr == NULL ) {
            //cout << "returning " << defval << endl;
            return defval;
        }
        return atoi( ptr );
    } catch ( invalid_argument& e ) {
        throw;
    }
}

double tktools::util::get_argument_float( int argc, char** argv, const char* option, double defval )  throw ( invalid_argument ) {
    try {
        const char* ptr = get_argument_pointer( argc, argv, option );
        if ( ptr == NULL ) return defval;
        return atof( ptr );
    } catch ( invalid_argument& e ) {
        throw;
    }
}

const char* tktools::util::get_argument_string( int argc, char** argv, const char* option, const char* defval ) throw ( invalid_argument ) {
    try {
        const char* ptr = get_argument_pointer( argc, argv, option );
        if ( ptr == NULL ) return defval;
        //cout << "option = " << option << ", PTR = " << ptr << endl;
        return ptr;
    } catch ( invalid_argument& e ) {
        throw;
    } catch (...) {
        throw;
    }
}

string tktools::util::get_log_string( const char* filename, int linenumber ) {
    char buffer[ 1024 ];
    snprintf( buffer, 1024, "at %d of %s", linenumber, filename );
    return string( buffer );
}

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
long tktools::util::get_usec() {
    struct timeval tv;
    struct timezone tz;
    gettimeofday( &tv, &tz );
    time_t sec = tv.tv_sec;
    suseconds_t usec = tv.tv_usec;
    return (long)usec + (long)sec * 1000;
}
#else
#include <ctime>
long tktools::util::get_usec() {
    return (long)( time( NULL ) * 1000 );
}
#endif

#ifdef HAVE_ZLIB_H
namespace {
    size_t expand_content_buffer( size_t current_size, unsigned char*& buffer ) {
        size_t next_size;
        if ( current_size < 100000000 ) {
            next_size = current_size * 2;
        } else {
            next_size = current_size + 100000000;
        }
        //cerr << current_size << " => " << next_size << " " << __FILE__ << " " << __LINE__ << endl;
        if ( next_size < 16 ) next_size = 16;
        unsigned char* heap = new unsigned char[ next_size ];
        memcpy( heap, buffer, current_size );
        delete[] buffer;
        buffer = heap;
        return next_size;
    }
}

size_t tktools::zip::compress_and_save( unsigned long size, unsigned char* input, FILE* output ) throw ( logic_error ) {
    unsigned char* cmpbuf;
    unsigned long cmpsize = 65536;
    z_stream zcmp;
    int status;
    size_t size_stream = 0;

    //cmpsize = size * 1.1;

    zcmp.zalloc = Z_NULL;
    zcmp.zfree  = Z_NULL;
    zcmp.opaque = Z_NULL;

    if ( deflateInit( &zcmp, Z_DEFAULT_COMPRESSION ) != Z_OK ) {
        cerr << ( zcmp.msg != NULL ? zcmp.msg : "???  " ) << __FILE__ << " " << __LINE__ << endl;
        throw logic_error( "Z-compression initialization failed" );
    }

    cmpbuf = new unsigned char[ cmpsize ];

    zcmp.next_in   = input;
    zcmp.avail_in = size;
    zcmp.next_out  = cmpbuf;
    zcmp.avail_out = cmpsize;

    for ( ;; ) {
        if ( zcmp.avail_in == 0 ) {
            status = deflate( &zcmp, Z_FINISH );
        } else {
            status = deflate( &zcmp, Z_NO_FLUSH );
        }
        if ( status == Z_STREAM_END ) break;
        if ( status != Z_OK ) {
            delete[] cmpbuf;
            cerr << ( zcmp.msg != NULL ? zcmp.msg : "???  " ) << __FILE__ << " " << __LINE__ << endl;
            throw logic_error( "compression failed" );
        }
        if ( zcmp.avail_out == 0 ) {
            size_t out_size = cmpsize;
            //cerr << "writing " << out_size << " bytes" << endl;
            if ( fwrite( cmpbuf, 1, out_size, output ) != out_size ) {
                delete[] cmpbuf;
                throw logic_error( "output failed" );
            }
            size_stream += out_size;
            zcmp.next_out = cmpbuf;
            zcmp.avail_out = cmpsize;
        }
    }    

    size_t size_remain = cmpsize - zcmp.avail_out;
    if ( size_remain > 0 ) {
        //cerr << "writing " << size_remain << " bytes" << endl;
        if ( fwrite( cmpbuf, 1, size_remain, output ) != size_remain ) {
            delete[] cmpbuf;
            throw logic_error( "output failed" );
        }
        size_stream += size_remain;
    }

    delete[] cmpbuf;
    return size_stream;
}

size_t tktools::zip::compress_buffer( size_t size, unsigned char* data_org, 
                                      unsigned char*& dst ) throw ( logic_error ) {
    unsigned char* cmpbuf;
    unsigned long cmpsize = 65536;
    z_stream zcmp;
    int status;
    size_t size_stream = 0;
    size_t dstsize = cmpsize;
    size_t offset = 0;
    dst = new unsigned char[ dstsize ];

    zcmp.zalloc = Z_NULL;
    zcmp.zfree  = Z_NULL;
    zcmp.opaque = Z_NULL;

    if ( deflateInit( &zcmp, Z_DEFAULT_COMPRESSION ) != Z_OK ) {
        cerr << ( zcmp.msg != NULL ? zcmp.msg : "???  " ) << __FILE__ << " " << __LINE__ << endl;
        throw logic_error( "Z-compression initialization failed" );
    }

    cmpbuf = new unsigned char[ cmpsize ];

    zcmp.next_in   = data_org;
    zcmp.avail_in = size;
    zcmp.next_out  = cmpbuf;
    zcmp.avail_out = cmpsize;

    for ( ;; ) {
        if ( zcmp.avail_in == 0 ) {
            status = deflate( &zcmp, Z_FINISH );
        } else {
            status = deflate( &zcmp, Z_NO_FLUSH );
        }
        if ( status == Z_STREAM_END ) break;
        if ( status != Z_OK ) {
            delete[] cmpbuf;
            cerr << ( zcmp.msg != NULL ? zcmp.msg : "???  " ) << __FILE__ << " " << __LINE__ << endl;
            throw logic_error( "compression failed" );
        }
        if ( zcmp.avail_out == 0 ) {
            size_t out_size = cmpsize;
            if ( out_size + offset >= dstsize ) {
                dstsize = expand_content_buffer( dstsize, dst );
            }
            memcpy( dst + offset, cmpbuf, out_size );
            offset += out_size;
            //            cerr << "writing " << out_size << " bytes " << offset << " " << __FILE__ << " " << __LINE__ << endl;
//             if ( !output.write( (const char*)cmpbuf, out_size ) ) {//== NULL ) {
//                 //if ( fwrite( cmpbuf, 1, out_size, output ) != out_size ) {
//                 delete[] cmpbuf;
//                 throw logic_error( "output failed" );
//             }
            size_stream += out_size;
            zcmp.next_out = cmpbuf;
            zcmp.avail_out = cmpsize;
        }
    }    

//     cout << __FILE__ << " " << __func__ << " " << __LINE__ << endl;
//     cout << (void*)cmpbuf << endl;
    size_t size_remain = cmpsize - zcmp.avail_out;
//    cout << size_remain << endl;
    if ( size_remain > 0 ) {
//         cout << __FILE__ << " " << __func__ << " " << __LINE__ << endl;
        //cerr << "writing " << size_remain << " bytes" << endl;
        if ( size_remain + offset >= dstsize ) {
            dstsize = expand_content_buffer( dstsize, dst );
        }
        memcpy( dst + offset, cmpbuf, size_remain );
        offset += size_remain;
//         if ( !output.write( (const char*)cmpbuf, size_remain ) ) {//fwrite( cmpbuf, 1, size_remain, output ) ) {//!= size_remain ) {
//             delete[] cmpbuf;
//             throw logic_error( "output failed" );
//         }
        size_stream += size_remain;
    }
    //cout << (void*)cmpbuf << endl;
    //cout << __FILE__ << " " << __func__ << " " << __LINE__ << endl;

    delete[] cmpbuf;
    return size_stream;
}


// namespace {
//     unsigned char* decompress_buffer( size_t size, const unsigned char* compressed ) {
        
//     }
// }

size_t tktools::zip::decompress_buffer( size_t buffer_size, unsigned char* buffer, unsigned char*& dst ) throw ( logic_error ) {
    z_stream zcmp;
    int status;
    size_t content_size = (buffer_size * 2 + 3) & 0xfffffffc;
    unsigned char* contents = new unsigned char[ content_size ];

    zcmp.zalloc = Z_NULL;
    zcmp.zfree  = Z_NULL;
    zcmp.opaque = Z_NULL;

    if ( inflateInit( &zcmp ) != Z_OK ) {
        cerr << ( zcmp.msg != NULL ? zcmp.msg : "???  " ) << __FILE__ << " " << __LINE__ << endl;
        throw logic_error( "Z-decompression initialization failed" );
    }


    zcmp.next_out = contents;
    zcmp.avail_out = content_size;
    status = Z_OK;

    zcmp.next_in = buffer;
    zcmp.avail_in = buffer_size;

    try {
        for ( ;; ) {

            //cerr << "IN " << zcmp.avail_in << " / OUT " << zcmp.avail_out << " / ftell " << ftell( input ) << "  " << __func__ << " " << __LINE__ << endl;

    //cout << "avail_in " << zcmp.avail_in << ", avail_out " << zcmp.avail_out << endl;
            status = inflate( &zcmp, Z_NO_FLUSH );
            if ( status == Z_STREAM_END ) break;
            if ( status != Z_OK ) {
                throw logic_error( "Z-decompression failed" );
            }
            if ( zcmp.avail_out == 0 ) {
                size_t expanded_size = expand_content_buffer( content_size, contents );
                zcmp.next_out = contents + content_size;
                zcmp.avail_out = expanded_size - content_size;
                content_size = expanded_size;
            }
        }

        //cout << "avail_in " << zcmp.avail_in << ", avail_out " << zcmp.avail_out << endl;
        if ( inflateEnd( &zcmp ) != Z_OK ) {
            cerr << ( zcmp.msg != NULL ? zcmp.msg : "???  " ) << __FILE__ << " " << __LINE__ << endl;
            throw logic_error( "decompression end failed" );
        }

    } catch ( exception & e ) {
        //delete[] buffer;
        delete[] contents;
        throw;
    }
    //if ( deflateInit( &zcmp, Z_DEFAULT_COMPRESSION ) != Z_OK )

    //delete[] buffer;
    dst = contents;
    //return contents;
//    cerr << content_size << " ;;; " << ((long long)zcmp.next_out - (long long)contents) << endl;
    return (size_t)(zcmp.next_out) - (size_t)(contents);
//    return content_size;
}

unsigned char* tktools::zip::load_and_decompress( FILE* input, size_t buffer_size ) throw ( logic_error ) {
    if ( buffer_size == 0 ) { // until end
        size_t current_pos = ftell( input );
        fseek( input, 0, SEEK_END );
        size_t filesize = ftell( input );
        fseek( input, current_pos, SEEK_SET );
        buffer_size = filesize - current_pos;
    }

    unsigned char* buffer;

    buffer = new unsigned char[ buffer_size ];
    size_t size_read = fread( buffer, 1, buffer_size, input );
    if (size_read == 0) {
      return NULL;
    }
    unsigned char* contents;
    decompress_buffer( buffer_size, buffer, contents );
    delete[] buffer;
    return contents;
}

#endif

tktools::bio::fasta_sequence::fasta_sequence() {
    _name = "";
    _length = 0;
    _sequence = NULL;
}

//tktools::bio::fasta_sequence::fasta_sequence(const char* filename, const char* name) throw ( exception ) {
//    load_file(filename, name);
//}

namespace {
    class fasta_info;
    class fasta_cache;

    class fasta_info {
    private:
        size_t _position;
        size_t _size;
        size_t _bases;
        string _name_org;
        string _chromosome;
    public:
        fasta_info();
        fasta_info(size_t pos, size_t size, size_t bases, const string& name, const string& chrm);
        size_t position() const { return _position; }
        size_t size() const { return _size; }
        size_t bases() const { return _bases; }
        const string& name() const { return _name_org; }
        const string& chromosome() const { return _chromosome; }
        friend class fasta_cache;
    };

    fasta_info::fasta_info(size_t pos, size_t size, size_t bases, const string& name, const string& chrm) {
        _position = pos;
        _size = size;
        _bases = bases;
        _name_org = name;
        _chromosome = chrm;
    }

    class fasta_cache {
        int _num;
        vector<fasta_info> _info;
        string _original_file;
        size_t _original_size;
        string _cache_file;
    private:
        void generate_cache();
        void load_from_cache();
        void save_cache();
        void initialize_with_file(const char* filename) throw (runtime_error);
        
    public:
        fasta_cache(const char* filename_fasta) throw (runtime_error);
        const fasta_info& get_info(int num) const throw (runtime_error);
        const fasta_info& get_info(const char* name) const;
        int get_size() const { return _info.size(); }
        const string& filename() const { return _original_file; }
    private:
        static const fasta_info not_found;
    };

    const fasta_info fasta_cache::not_found = fasta_info();

    fasta_info::fasta_info() {
        _position = 0;
        _size = 0;
        _bases = 0;
    }

    fasta_cache::fasta_cache(const char* filename) throw (runtime_error) {
        initialize_with_file(filename);
    }


    vector<string> __chromosome_keywords;

    void fasta_cache::initialize_with_file(const char* filename) throw (runtime_error){
        if (__chromosome_keywords.size() == 0) {
            __chromosome_keywords.push_back("chromosome");
            __chromosome_keywords.push_back("Chromosome");
            __chromosome_keywords.push_back("chr");
            __chromosome_keywords.push_back("Chr");
        }

        //char chromosome_keywords*[4] = {};
        _original_file = filename;
        size_t pos = _original_file.rfind(tktools::io::file_separator());
        if (pos == string::npos) {
            _cache_file = string(".");
            pos = 0;
        } else {
            _cache_file = _original_file.substr(0, pos + 1) + ".";
            pos ++;
        }
        _cache_file += _original_file.substr(pos, _original_file.size() - pos) + ".fcache";
        //cerr << pos << " " << _cache_file << endl;
        //exit(0);
        if (tktools::io::file_exists(_cache_file.c_str())) {
            ifstream fi(_cache_file.c_str());
            while (fi.eof() == false) {
                string line;
                getline(fi, line);
                vector<string> items = split_items(line, '\t');
                if (items.size() > 4) {
                    _info.push_back(fasta_info(atoi(items[0].c_str()), atoi(items[1].c_str()), atoi(items[2].c_str()), items[3], items[4]));
                }
            }
            fi.close();
        } else {
            ifstream fi(filename);
            if (fi.is_open() == false) {
                throw runtime_error(string("cannot open fasta file " ) + filename);
            }

            fasta_info info;
            size_t position = 0;
            //size_t size = 0;
            size_t num_bases = 0;
            string name;
            string chromosome;
            size_t filesize = tktools::io::get_file_size(filename);
            while (!fi.eof()) {
                string line;
                size_t p1 = fi.tellg();
                getline(fi, line);
                if (line.c_str()[0] == '>') {
                    if (num_bases > 0) {
                        _info.push_back(fasta_info(position, p1 - position, num_bases, name, chromosome));
                    }
                    name = line.substr(1, line.size() - 1);
                    num_bases = 0;
                    position = p1;
                    //cout << "revise pos " << position << endl;
                    bool assigned = false;
                    const char* ptr = line.c_str();
                    for (int i = 0; i < (int)__chromosome_keywords.size(); i++) {
                        const char* cptr = strstr(ptr, __chromosome_keywords[i].c_str());
                        if (cptr != NULL) {
                            ptr = cptr + __chromosome_keywords[i].size();
                            for (;;ptr++) {
                                char c = *ptr;
                                if (c == '\0') break;
                                if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z')){
                                    int span = 1;
                                    for (;; span++) {
                                        c = ptr[span];
                                        if (c == '\0') break;
                                        if ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z')){
                                            continue;
                                        } else {
                                            break;
                                        }
                                    }
                                    chromosome = string(ptr, span);
                                    assigned = true;
                                    break;
                                }
                            }
                                    
                        }
                        if (assigned) break;
                    }
                    if (!assigned) {
                        chromosome = name;
                    }
                } else {
                    const char* ptr = line.c_str();
                    for (int i = 0; i < (int)line.size(); i++) {
                        char c = ptr[i];
                        if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) {
                            num_bases ++;
                        }
                    }
                }
            }
            if (num_bases > 0) {
                _info.push_back(fasta_info(position, filesize - position, num_bases, name, chromosome));
            }
            fi.close();

            // dump cache
            ofstream fo(_cache_file.c_str());
            for (vector<fasta_info>::const_iterator it = _info.begin(); it != _info.end(); it++) {
                fo << it->position() << "\t" << it->size() << "\t" << it->bases() << "\t" << it->name() << "\t" << it->chromosome() << "\n";
            }
            fo.close();
        }
    }

    const fasta_info& fasta_cache::get_info(int num) const throw (runtime_error){
      if (num < 0 || num > (int)_info.size()) {
            throw out_of_range("");
        }
        return _info[num];
    }

    const fasta_info& fasta_cache::get_info(const char* name) const {
        if (name == NULL) {
            return get_info(0);
        }
        for (int i = 0; i < (int)_info.size(); i++) {
            if (_info[i].name() == name || _info[i].chromosome() == name) {
                return _info[i];
            }
        }
        return not_found;
    }

}

tktools::bio::fasta_sequence* tktools::bio::fasta_sequence::load_file(const char* filename, const char* name) throw (exception) {
    try {
        fasta_cache fdata(filename);
        fasta_info info = fdata.get_info(name);
        ifstream fi(filename);
        fi.seekg(info.position());
        size_t size = info.size();
        char* buffer = new char[size];
        fi.read(buffer, info.size());
        fi.close();
        if (buffer[0] == '>') {
            //const char* ptr = const_cast<const char*>(buffer);
            fasta_sequence* fasta = NULL;
            size_t pos;
            for (pos = 0; pos < size; pos++) {
                char c = buffer[pos];
                if (c == '\n') {
                    buffer += pos + 1;
                    fasta = new fasta_sequence();
                    fasta->_name = info.chromosome();
                    fasta->_length = info.size();
                    fasta->_sequence = new char[fasta->_length + 1];
                    break;
                }
            }
            if (fasta == NULL) {
                delete[] buffer;
                throw runtime_error(string("not found ") + name);
            }
            char* ptr_write = fasta->_sequence;
            //const char* buffer;
            for (; pos < size; pos++) {
                char c = *buffer;
                if (c >= 'a' && c <= 'z') {
                    c -= (char)('a' - 'A');
                }
                if (c >= 'A' && c <= 'Z' ) {
                    *ptr_write = c;
                    ptr_write++;
                }
                buffer++;
            }
            *ptr_write = '\0';
            fasta->_length = info.bases();
            return fasta;
        } else {
            throw runtime_error("missing header");
        }
    } catch (exception& e) {
        cerr << e.what() << endl;
        throw;
    }
}

tktools::bio::fasta_sequence* tktools::bio::load_fasta_sequence(const char* filename, const char* name) throw ( exception ) {
    try {
        return tktools::bio::fasta_sequence::load_file(filename, name);
    } catch ( exception& e ) {
        cerr << __FILE__ << " " << __LINE__ << endl;
        throw;
    }
}

void tktools::bio::fasta_sequence::set_sequence(int length, const char* sequence) {
    if (length > _length) {
        delete[] _sequence;
        _sequence = new char[length + 1];
    }
    _length = length;
    char* ptr_write = _sequence;
    for ( int i = 0; i < length; i++ ) {
        char base = sequence[i];
        if (base >= 'a') base -= (char)('a' - 'A');
        if (base >= 'A' && base <= 'Z') {
            *ptr_write++ = base;
        }
    }
    *ptr_write = '\0';
    _length = strlen(_sequence);
}

int tktools::bio::convert_chromosome_to_code(const char* name) {
    int code;
    if (strncmp(name, "chr", 3) == 0) {
        name += 3;
    }
    if ( name[ 0 ] == 'X' ) {//strcmp( name, "X" ) == 0 ) {
        code = 64;
    } else if ( name[ 0 ] == 'Y' ) {// strcmp( name, "Y" ) == 0 ) {
        code = 65;
    } else if ( name[ 0 ] == 'W' ) {
        code = 66;
    } else if ( name[ 0 ] == 'Z' ) {
        code = 67;
    } else if ( name[ 0 ] == 'M' || strcmp( name, "mitochondria" ) == 0 ) {
        code = 128;
    } else if ( name[ 0 ] == 'C' || strcmp( name, "chroloplast" ) == 0 ) {
        code = 129;
    } else if ( name[ 0 ] == 'U' ) {
        code = 192;
    } else if ( strcmp(name, "I") == 0) {
        code = 1;
    } else if ( strcmp(name, "II") == 0) {
        code = 2;
    } else if ( strcmp(name, "III") == 0) {
        code = 3;
    } else if ( strcmp(name, "IV") == 0) {
        code = 4;
    } else if ( strcmp(name, "V") == 0) {
        code = 5;
    } else {
        code = atoi(name);
    }

    if ( code > 0 ) {
        if ( strstr( name, "andom" ) != NULL ) {
            code = -code;
        }
        return code;
    } else {
        cerr << "undetermined chromosome " << name << endl;
        return 0;
    }
}

string tktools::bio::convert_code_to_chromosome(int code) {
    if ( code > 0 && code < 64 ) {
        char buffer[256];
        snprintf( buffer, 255, "%d", code );
        string cstr(buffer);
        return cstr;
        //return string(buffer);
    } else if ( code == 64 ) {
        return string( "X" );
    } else if ( code == 65 ) {
        return string( "Y" );
    } else if ( code == 66 ) {
        return string( "W" );
    } else if ( code == 67 ) {
        return string( "Z" );
    } else if ( code == 128 ) {
        return string( "M" );
    } else if ( code == 129 ) {
        return string( "C" );
    } else if ( code == 192 ) {
        return string( "Un" );
    } else if ( code < 0 ) {
        return convert_code_to_chromosome(-code ) + "_random";
    } else {
        return string( "undetermined" );
    }
}

tktools::bio::fasta_sequence::~fasta_sequence() {
    delete[] _sequence;
}

namespace {
    inline char reverse_base(char base) {
        switch ( base ) {
        case 'a': base = 't'; break;
        case 'c': base = 'g'; break;
        case 'g': base = 'c'; break;
        case 't': base = 'a'; break;
        case 'A': base = 'T'; break;
        case 'C': base = 'G'; break;
        case 'G': base = 'C'; break;
        case 'T': base = 'A'; break;
        }
        return base;
    }
}

void tktools::bio::fasta_sequence::reverse() {
    int center = _length / 2;
    for ( int i = 0; i < center; i++ ) {
        char bt = reverse_base(_sequence[i]);
        _sequence[i] = reverse_base(_sequence[_length - 1 - i ]);
        _sequence[_length - i - 1] = bt;
    }
    if ( (_length & 1) != 0 ) {
        _sequence[center] = reverse_base(_sequence[center]);
    }
}
namespace {
    char complementary_base(char base) {
        switch (base) {
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:
            break;
        }
        return base;
    }
}
string tktools::bio::fasta_sequence::get_subsequence( int start, int end, bool complementary ) const {
    if ( start > end ) {
        int s = start;
        start = end;
        end = s;
    }
    char* buffer = new char[ end - start + 1 ];
    if (complementary == false) {
        memcpy( buffer, _sequence + start, end - start );
    } else {
        for (int i = start; i < end; i++) {
            buffer[end - i - 1] = complementary_base(*(_sequence + i));
        }
    }
    buffer[ end - start ] = '\0';
    string subseq = buffer;
    delete[] buffer;
    return subseq;
}

int tktools::bio::convert_solid_color_to_sequence(int length, const char* codes, char* sequence, char top) {
    if (top == '\0') {
        top = codes[0];
        codes++;
        length--;
    }
    static char solid_converter[18] = "ACGTCATGGTACTGCA\0";
    char base_current = top;
    for (int i = 0; i < length; i++) {
        int num = codes[i] - '0';
        switch (base_current) {
        case 'A':
            break;
        case 'C':
            num += 4;
            break;
        case 'G':
            num += 8;
            break;
        case 'T':
            num += 12;
            break;
        default:
            num = 16;
            length = i;
            break;
        }
        base_current = solid_converter[num];
        sequence[i] = base_current;
    }
    sequence[length] = '\0';
    return length;
}

int tktools::bio::convert_sequence_to_solid_color(int length, const char* sequence, char* codes, char top) {
    static char pattern[18] = "0123103223013210\0";
    int shift;
    if (top == '\0') {
        top = sequence[0];
        length --;
    }
    codes[0] = top;
    switch (top) {
    case 'A':
        shift = 0; break;
    case 'C':
        shift = 4; break;
    case 'G': 
        shift = 8; break;
    case 'T':
        shift = 12; break;
    default:
        codes[1] = '\0';
        return 0;
    }
    for (int i = 0; i < length; i++) {
        char base = sequence[i];
        int position = shift;
        switch (base) {
        case 'A': 
            position = shift; shift = 0; break;
        case 'C': 
            position = shift + 1; shift = 4; break;
        case 'G':
            position = shift + 2; shift = 8; break;
        case 'T': 
            position = shift + 3; shift = 12; break;
        default:
            length = i;
            break;
        }
        codes[i + 1] = pattern[position];
    }
    codes[length + 1] = '\0';
    return length + 1;
}
