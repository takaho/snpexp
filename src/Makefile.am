bin_PROGRAMS = snpexp trimvcf recrec
#countsnps
recrec_SOURCES = recrec.cxx recrec.hxx tktools.hxx tktools.cxx gtf.hxx gtf.cxx main.cxx fragmentprocessor.hxx fragmentprocessor.cxx distsnp.cxx distsnp.hxx varstr.cxx varstr.hxx indel.hxx indel.cxx
snpexp_SOURCES = snpexp.cxx snpexp.hxx tktools.hxx tktools.cxx gtf.hxx gtf.cxx 
trimvcf_SOURCES = trimvcf.hxx trimvcf.cxx gtf.hxx gtf.cxx tktools.hxx tktools.cxx seq_gene.hxx seq_gene.cxx
#countsnps_SOURCES = gtf.hxx gtf.cxx distsnp.hxx distsnp.cxx countsnps.cxx recrec.hxx recrec.cxx tktools.hxx tktools.cxx
CXX = c++
AM_CXXFLAGS = -O3 -pedantic -Wall -Wno-long-long -m64 -I/usr/local/include
AM_LDFLAGS = -L/usr/local/lib
LIBS = -lz -lbam -lpthread
