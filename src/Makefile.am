bin_PROGRAMS = snpexp trimvcf
snpexp_SOURCES = snpexp.cxx snpexp.hxx tktools.hxx tktools.cxx gtf.hxx gtf.cxx bam.h bgzf.h
trimvcf_SOURCES = trimvcf.hxx trimvcf.cxx gtf.hxx gtf.cxx tktools.hxx tktools.cxx
CXX = c++
AM_CXXFLAGS = -O3 -pedantic -Wall -Wno-long-long -m64 -I/usr/local/include
AM_LDFLAGS = -L/usr/local/lib
LIBS = -lz -lbam
