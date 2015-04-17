AC_INIT([snpexp], [0.4.1], [takaho.endo@riken.jp])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([snpexp.cxx])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
#AC_PROG_CC
AC_CHECK_HEADERS([zlib.h bam.h pthread.h])
AC_PROG_CXX
#AC_CHECK_SIZEOF(void*, 8)
#AC_CHECK_SIZEOF(int, 4)
#AC_CHECK_LIB(pthread, pthread_create)
#AC_CHECK_LIB(bam, bam_read)
#if test -f /usr/local/include/bam.h
#then
#	AC_DEFINE(HAVE_BAM_H)
#fi
AC_SEARCH_LIBS([deflate], [z], [], [])
AC_OUTPUT
#AC_OUTPUT(Makefile)
