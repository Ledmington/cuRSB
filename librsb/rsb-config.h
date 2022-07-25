/*

Copyright (C) 2008-2022 Michele Martone

This file is part of librsb.

librsb is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

librsb is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with librsb; see the file COPYING.
If not, see <http://www.gnu.org/licenses/>.

*/
/* This header file is not intended to be included librsb programs: it is only for inspection. */
#ifndef RSB_CONFIG_H_INCLUDED
#define RSB_CONFIG_H_INCLUDED
/* rsb-config.h.  Generated from rsb-config.h.in by configure.  */
/* rsb-config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* C compiler. */
#define RSB_CC "gcc"

/* Compilation flags. */
#define RSB_CFLAGS "-g -O2 -std=c99"

/* */
#define RSB_COPYRIGHT_STRING "Copyright (c) 2008-2022 Michele Martone"

/* Compilation flags. */
#define RSB_CXXFLAGS "-g -O2 -fopenmp"

/* Define to 1 if you have the <assert.h> header file. */
#define RSB_HAVE_ASSERT_H 1

/* Define to 1 if you have the <complex.h> header file. */
#define RSB_HAVE_COMPLEX_H 1

/* Define to 1 if you have the <ctype.h> header file. */
#define RSB_HAVE_CTYPE_H 1

/* Define to 1 if you have the <dirent.h> header file. */
#define RSB_HAVE_DIRENT_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define RSB_HAVE_DLFCN_H 1

/* Define to 1 if you have the <dmalloc.h> header file. */
/* #undef HAVE_DMALLOC_H */

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* Define to 1 if you have the `dup' function. */
#define RSB_HAVE_DUP 1

/* Define to 1 if you have the <execinfo.h> header file. */
#define RSB_HAVE_EXECINFO_H 1

/* fileno(): C FILE to posix file descriptor. */
#define RSB_HAVE_FILENO 1

/* Define to 1 if you have the `fread' function. */
#define RSB_HAVE_FREAD 1

/* Define to 1 if you have the `fwrite' function. */
#define RSB_HAVE_FWRITE 1

/* Get an environment variable. */
#define RSB_HAVE_GETENV 1

/* If present, will give us host name. */
#define RSB_HAVE_GETHOSTNAME 1

/* getopt */
#define RSB_HAVE_GETOPT 1

/* Define to 1 if you have the <getopt.h> header file. */
#define RSB_HAVE_GETOPT_H 1

/* getopt_long is GNU candy */
#define RSB_HAVE_GETOPT_LONG 1

/* gettimeofday */
#define RSB_HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <gmock/gmock.h> header file. */
/* #undef HAVE_GMOCK_GMOCK_H */

/* Define to 1 if you have the <gsl/gsl_sort.h> header file. */
/* #undef HAVE_GSL_GSL_SORT_H */

/* Define to 1 if you have the <gtest/gtest.h> header file. */
/* #undef HAVE_GTEST_GTEST_H */

/* used by RSB_WANT_EXPERIMENTAL_BINARY_COO, and not present on older zlib */
#define RSB_HAVE_GZFREAD 1

/* Define to 1 if you have the <hwloc.h> header file. */
#define RSB_HAVE_HWLOC_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define RSB_HAVE_INTTYPES_H 1

/* Define to 1 if you have the `isatty' function. */
#define RSB_HAVE_ISATTY 1

/* Define to 1 if you have the <limits.h> header file. */
#define RSB_HAVE_LIMITS_H 1

/* Define to 1 if you have the <malloc.h> header file. */
#define RSB_HAVE_MALLOC_H 1

/* Define to 1 if you have the <math.h> header file. */
#define RSB_HAVE_MATH_H 1

/* This function is obsolete. */
#define RSB_HAVE_MEMALIGN 1

/* Define to 1 if you have the <memory.h> header file. */
#define RSB_HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define RSB_HAVE_MEMSET 1

/* Define to 1 if you have the <mkl/mkl.h> header file. */
/* #undef HAVE_MKL_MKL_H */

/* If present, the mlockall function makes all allocations memory resident. */
#define RSB_HAVE_MLOCKALL 1

/* Define to 1 if you have the <omp.h> header file. */
#define RSB_HAVE_OMP_H 1

/* Define to 1 if you have the <oski/oski.h> header file. */
/* #undef HAVE_OSKI_OSKI_H */

/* Define to 1 if you have the <papi.h> header file. */
/* #undef HAVE_PAPI_H */

/* The POSIX aligned memory allocator.(The function posix_memalign() is
   available since glibc 2.1.91) */
#define RSB_HAVE_POSIX_MEMALIGN 1

/* Define to 1 if you have the <pthread.h> header file. */
#define RSB_HAVE_PTHREAD_H 1

/* Define to 1 if you have the `rand' function. */
#define RSB_HAVE_RAND 1

/* Define to 1 if you have the <regex.h> header file. */
#define RSB_HAVE_REGEX_H 1

/* rindex */
#define RSB_HAVE_RINDEX 1

/* Define to 1 if you have the <rpc/xdr.h> header file. */
#define RSB_HAVE_RPC_XDR_H 1

/* Define to 1 if you have the `sched_getaffinity' function. */
#define RSB_HAVE_SCHED_GETAFFINITY 1

/* Define to 1 if you have the <sched.h> header file. */
#define RSB_HAVE_SCHED_H 1

/* setenv */
#define RSB_HAVE_SETENV 1

/* Define to 1 if you have the <signal.h> header file. */
#define RSB_HAVE_SIGNAL_H 1

/* Define to 1 if you have the <stdarg.h> header file. */
#define RSB_HAVE_STDARG_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define RSB_HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define RSB_HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define RSB_HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define RSB_HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define RSB_HAVE_STRING_H 1

/* Define to 1 if you have the `strncmp' function. */
#define RSB_HAVE_STRNCMP 1

/* strrchr */
#define RSB_HAVE_STRRCHR 1

/* If present, the sysconf function gives lots of system info. */
#define RSB_HAVE_SYSCONF 1

/* Define to 1 if you have the <sys/mman.h> header file. */
#define RSB_HAVE_SYS_MMAN_H 1

/* Define to 1 if you have the <sys/resource.h> header file. */
#define RSB_HAVE_SYS_RESOURCE_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define RSB_HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/systemcfg.h> header file. */
/* #undef HAVE_SYS_SYSTEMCFG_H */

/* Define to 1 if you have the <sys/time.h> header file. */
#define RSB_HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define RSB_HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <sys/utsname.h> header file. */
#define RSB_HAVE_SYS_UTSNAME_H 1

/* times */
#define RSB_HAVE_TIMES 1

/* Define to 1 if you have the <times.h> header file. */
/* #undef HAVE_TIMES_H */

/* Define to 1 if you have the <time.h> header file. */
#define RSB_HAVE_TIME_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define RSB_HAVE_UNISTD_H 1

/* Define to 1 if you have the `vprintf' function. */
#define RSB_HAVE_VPRINTF 1

/* Define to 1 if you have the <zlib.h> header file. */
#define RSB_HAVE_ZLIB_H 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define RSB_LT_OBJDIR ".libs/"

/* OSKI path to installed lua modules. User set OSKI_LUA_PATH environment
   variable at runtime will override this one, however. */
/* #undef OSKI_LUA_PATH */

/* Name of package */
#define RSB_PACKAGE "librsb"

/* Define to the address where bug reports for this package should be sent. */
#define RSB_PACKAGE_BUGREPORT "michelemartone_AT_users_DOT_sourceforge_DOT_net"

/* Define to the full name of this package. */
#define RSB_PACKAGE_NAME "librsb"

/* Define to the full name and version of this package. */
#define RSB_PACKAGE_STRING "librsb ec50843"

/* Define to the one symbol short name of this package. */
#define RSB_PACKAGE_TARNAME "librsb"

/* Define to the home page for this package. */
#define RSB_PACKAGE_URL ""

/* Define to the version of this package. */
#define RSB_PACKAGE_VERSION "ec50843"

/* Extra (undocumented) developer oriented control switches. */
/* #undef RSB_ALLOW_INTERNAL_GETENVS */

/* If set, the library will use smaller indices in blocks. */
#define RSB_BLOCK_SMALL_INDICES 1

/* Maximal number of supported threads (default 128). */
#define RSB_CONST_MAX_SUPPORTED_THREADS 128

/* If not null, the library will rely on this for memory hierarchy info,
   unless RSB_USER_SET_MEM_HIERARCHY_INFO is set. */
#define RSB_DETECTED_MEM_HIERARCHY_INFO "L3:11/64/19712K,L2:16/64/1024K,L1:8/64/32K"

/* If defined, will not account for internally used memory. */
#define RSB_DISABLE_ALLOCATOR_WRAPPER 1

/* Performance Application Programming Interface. */
/* #undef RSB_HAVE_PAPI */

/* Will include mkl/mkl.h */
/* #undef RSB_INCLUDE_MKL_MKL_H */

/* Inner error verbosity (internal debug level). */
#define RSB_INT_ERR_VERBOSITY 0

/* Extra internal memory checks (for debugging). */
/* #undef RSB_MEM_DBG */

/* Error verbosity (often known as debug level). */
#define RSB_OUT_ERR_VERBOSITY 0

/* If set, sort operations will happen in place. */
#define RSB_SORT_IN_PLACE 0

/* If not null, the library will rely on this for memory hierarchy info. */
#define RSB_USER_SET_MEM_HIERARCHY_INFO ""

/* If undefined, NDEBUG will be defined. */
/* #undef RSB_USE_ASSERT */

/* Usable rusage and getrusage? */
#define RSB_USE_GETRUSAGE 1

/* Use librsbpp. */
#define RSB_USE_LIBRSBPP 1

/* Enable calling MKL from RSB (internal, deprecated). */
/* #undef RSB_USE_MKL */

/* experimental. */
#define RSB_WANT_ACTION_SIGNAL 1

/* If 1, will allow the user to set hard limits to the memory allocated by
   librsb. Trespass attempts will fail. */
#define RSB_WANT_ALLOCATOR_LIMITS 0

/* No ARMPL support wanted in the benchmarking program. */
#define RSB_WANT_ARMPL 0

/* */
#define RSB_WANT_DMALLOC 0

/* On some architectures (notably modern Intel), floating point computations
   on non double aligned data make loose some clock cycle. */
#define RSB_WANT_DOUBLE_ALIGNED 1

/* "Will use Google Test" */
/* #undef RSB_WANT_GTEST */

/* Supported input/output functionality. */
#define RSB_WANT_IO_LEVEL 7

/* If set, RSB_WANT_KERNELS_DEBUG will enable comparative consistency checking
   of the multiplying kernels against a naive, trusted implementation. */
#define RSB_WANT_KERNELS_DEBUG 0

/* Enabling collection of time statistics in librsb operations (this
   introduces an overhead). */
/* #undef RSB_WANT_LIBRSB_STATS */

/* long types for rsb_coo_idx_t and rsb_nnz_idx_t */
/* #undef RSB_WANT_LONG_IDX */

/* Looping kernels. */
/* #undef RSB_WANT_LOOPING_KERNELS */

/* No MKL support wanted in the benchmarking program. */
#define RSB_WANT_MKL 0

/* Support for reading matrices in parallel (Experimental, untested). */
#define RSB_WANT_OMPIO_SUPPORT 0

/* Recursive kernels parallelized with OpenMP. */
#define RSB_WANT_OMP_RECURSIVE_KERNELS 1

/* OSKI comparative benchmarking. */
/* #undef RSB_WANT_OSKI_BENCHMARKING */

/* Performance Counters. */
/* #undef RSB_WANT_PERFORMANCE_COUNTERS */

/* Enabling experimental RSB_NUM_THREADS environment variable. */
#define RSB_WANT_RSB_NUM_THREADS 1

/* If set, a reference, unoptimized Sparse BLAS Level 1 interface will be
   functional. */
#define RSB_WANT_SPARSE_BLAS_LEVEL_1 1

/* If set, the library will be much more verbose. Should be enabled for
   debugging purposes only. */
#define RSB_WANT_VERBOSE_MESSAGES 0

/* experimental. */
#define RSB_WANT_XDR_SUPPORT 1

/* Support for reading gzipped matrices. */
#define RSB_WANT_ZLIB_SUPPORT 1

/* HWLOC API support. */
#define RSB_WITH_HWLOC 0

/* LIKWID marker API support. */
#define RSB_WITH_LIKWID 0

/* Sparse BLAS interface compilation. */
#define RSB_WITH_SPARSE_BLAS_INTERFACE 1

/* The size of `char', as computed by sizeof. */
#define RSB_SIZEOF_CHAR 1

/* The size of `complex', as computed by sizeof. */
#define RSB_SIZEOF_COMPLEX 16

/* The size of `double', as computed by sizeof. */
#define RSB_SIZEOF_DOUBLE 8

/* The size of `double complex', as computed by sizeof. */
#define RSB_SIZEOF_DOUBLE_COMPLEX 16

/* The size of `float', as computed by sizeof. */
#define RSB_SIZEOF_FLOAT 4

/* The size of `float complex', as computed by sizeof. */
#define RSB_SIZEOF_FLOAT_COMPLEX 8

/* The size of `int', as computed by sizeof. */
#define RSB_SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define RSB_SIZEOF_LONG 8

/* The size of `long double', as computed by sizeof. */
#define RSB_SIZEOF_LONG_DOUBLE 16

/* The size of `long double complex', as computed by sizeof. */
#define RSB_SIZEOF_LONG_DOUBLE_COMPLEX 32

/* The size of `long int', as computed by sizeof. */
#define RSB_SIZEOF_LONG_INT 8

/* The size of `long long int', as computed by sizeof. */
#define RSB_SIZEOF_LONG_LONG_INT 8

/* The size of `short int', as computed by sizeof. */
#define RSB_SIZEOF_SHORT_INT 2

/* The size of `size_t', as computed by sizeof. */
#define RSB_SIZEOF_SIZE_T 8

/* The size of `void *', as computed by sizeof. */
#define RSB_SIZEOF_VOID_P 8

/* Define to 1 if you have the ANSI C header files. */
#define RSB_STDC_HEADERS 1

/* VCS REVISION */
#define RSB_VCS_REVISION "ec50843"

/* Version number of package */
#define RSB_VERSION "ec50843"

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#define RSB_restrict __restrict
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
#endif /* RSB_CONFIG_H_INCLUDED */
