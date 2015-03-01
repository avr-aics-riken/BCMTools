##############################
#
# BCMTools
#
# Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
##############################

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.57])
AC_INIT([BCMTools], [0.9.5], [jonishi@iis.u-tokyo.ac.jp], [BCMTools])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/BlockManager.cpp])
AC_CONFIG_HEADERS([config.h])


#
# revision No.
#
BCM_REVISION="20140323_1600"
AC_SUBST(BCM_REVISION)


#
# Checks for programs.
#
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB


#
# Checks for libraries.
#


#
# Checks for header files.
#
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stdint.h stdlib.h string.h sys/time.h unistd.h])


#
# Checks for typedefs, structures, and compiler characteristics.
#
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_MODE_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T


#
# Checks for library functions.
#
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday memset mkdir sqrt strdup])



#
# default prefix setting
#
AC_PREFIX_DEFAULT([/usr/local/BCMTools])
if test "$prefix" = "NONE" ; then
  BCM_INST_DIR=/usr/local/BCMTools
else
  BCM_INST_DIR="$prefix"
fi


#
# MPI environment
#
enable_mpi="none"


BCM_CXX="$CXX"
AC_SUBST(BCM_CXX)

case "$BCM_CXX" in
  mpicxx)   enable_mpi="wrappermpi" ;;
  mpicc)    enable_mpi="wrappermpi" ;;
  mpiCC)    enable_mpi="wrappermpi" ;;
  mpic++)   enable_mpi="wrappermpi" ;;
  mpiFCCpx) enable_mpi="wrappermpi" ;;
esac


AC_SUBST(MPI_DIR)
AC_SUBST(MPI_CFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)



# 
# MPI library
#
AC_ARG_WITH(MPI, [AC_HELP_STRING([--without-MPI],[do not use MPI library])], , with_MPI=yes)
if test "$with_MPI" = "no" ; then
  CXXFLAGS=-D_BCM_WITHOUT_MPI_
else

  # openmpi
  AC_ARG_WITH( ompi, [AC_HELP_STRING([--with-ompi=dir],[Specify OpenMPI install directory])], , with_ompi=none)
  if test "$with_ompi" != "none" ; then
    MPI_DIR=$with_ompi;
    MPI_CFLAGS=-I$MPI_DIR/include
    MPI_LDFLAGS=-L$MPI_DIR/lib
    MPI_LIBS=-lmpi
    enable_mpi="ompi"
  fi
fi


#
# TextParser
#
AC_ARG_WITH(parser, [AC_HELP_STRING([--with-parser=dir],[Specify TextParser install directory])], , with_parser=none)
AC_SUBST(TP_DIR)
AC_SUBST(TP_CFLAGS)
AC_SUBST(TP_LDFLAGS)

if test "$with_parser" != "none" ; then
  TP_DIR=$with_parser;
  TP_CFLAGS=`$TP_DIR/bin/tp-config --cflags`
  TP_LDFLAGS=`$TP_DIR/bin/tp-config --libs`
else
  echo "  Error: TextParser library is not specified."
  echo "         Specify --with-parser option."
  exit 1
fi


#
# Polylib
#
AC_ARG_WITH(polylib, [AC_HELP_STRING([--with-polylib=dir],[Specify Polylib install directory])], , with_polylib=none)
AC_SUBST(PL_DIR)
AC_SUBST(PL_CFLAGS)
AC_SUBST(PL_LDFLAGS)

if test "$with_polylib" != "none" ; then
  PL_DIR=$with_polylib;
  PL_CFLAGS=`$PL_DIR/bin/polylib-config --cflags`
  PL_LDFLAGS=`$PL_DIR/bin/polylib-config --libs`
else
  echo "  Error: Polylib is not specified."
  echo "         Specify --with-polylib option."
  exit 1
fi


#
# Warning
#

if test "$with_MPI" = "yes" ; then
  if test "$enable_mpi" != "wrappermpi" ; then
    if test "$with_ompi" = "none" ; then
      echo "  Warning: OpenMPI library is not specified."
      echo "           If compilation error occurs about MPI, specify --with-ompi option."
    fi
  fi
fi



#
# BCM special flags
#
BCM_CFLAGS="-I$BCM_INST_DIR/include"
AC_SUBST(BCM_CFLAGS)

BCM_LDFLAGS="-L$BCM_INST_DIR/lib"
AC_SUBST(BCM_LDFLAGS)

BCM_LIBS="-lBCM"
AC_SUBST(BCM_LIBS)


#
# MAKE_SUB_DIRS
#
AC_SUBST(MAKE_SUB_DIRS)
MAKE_SUB_DIRS="src doc doxygen"





AC_CONFIG_FILES([Makefile \
                 src/Makefile \
                 bcm-config \
                 include/bcm_Version.h \
                 doc/Makefile \
                 doxygen/Makefile])
AC_OUTPUT

chmod +x ./bcm-config
