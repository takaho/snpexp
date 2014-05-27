README for SnpExp
======

# Introduction

SnpExp is an allele frequency counter with BAM file and VCF file.

# Install

Before building this project, you have to install libbam provided by Li et al. (http://samtools.sourceforge.net/). 
If you build SAMtools bam.h, bgzf.h and libbam.o files will appear in the directory.
These files should be copied to paths for header files and library files typically at /usr/local/include/ and /usr/local/lib/.

If these files are correctly installed, you can setup and build the program.

    ./configure
    make

Built binaries can be installed your executable paths.

# Usage of SnpExp

    snpexp [-h] [--version] 
           <options> bam1 bam2 ...

##Options
###-V FILENAME
Filename of variants provided in VCF format

###-o FILENAME
Filename of output. If no filename is given, the result will be outputted to standard output.

###-G FILENAME
GTF file of gene annotation which includes genic positions.

###-I
Retain intergenic SNPs.

###-verbose
Verbose mode.

###-s strain1,strain2
Straind such as C57BL6NJ, 129S1 for mice

###-m number (default 0)
Minimum number of detected bases. If all bam files do not include the given number of bases, the SNPs will not be outputted.

### bam files
BAM files should be sorted before processing. They should also share the same set of chromosomes because the program analyze data chromosome by chromosome.

### -h
Show command help.

