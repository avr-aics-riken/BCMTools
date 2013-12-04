###################################################################
#
# FFV : Frontflow / violet BCM
#
# Copyright (c) 2013 All right reserved.
#
# Institute of Industrial Science, University of Tokyo, Japan. 
#
###################################################################

########################
# FX10
########################

AR          = ar cr
RANLIB      = ranlib
RM          = \rm -f
MPI_DIR     =
BCM_DIR     = $(HOME)/BCMTools
BCM_FILEIO_DIR     = $(HOME)/BCMTools/examples/FileIO
CPM_DIR     = $(HOME)/MyLib/CPMlib
TP_DIR      = $(HOME)/MyLib/TPlib
POLYLIB_DIR = $(HOME)/MyLib/Polylib_2_0_3_Rel
CUTLIB_DIR  = $(HOME)/MyLib/Cutlib-2.0.5
PMLIB_DIR   = $(HOME)/MyLib/PMlib-1.5
XML2_DIR    = $(HOME)/MyLib/libxml2-2.7.8_lib
OMP_FLAGS   = -Kopenmp
DEFINES     = -DDEBUG
UDEF_OPT    = -D__K_FPCOLL
UDEF_OPT    = -D__K_FAPP -D__K_FPCOLL
XG          = -Xg
CC          = mpifccpx
CFLAGS      = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,optmsg=1 $(OMP_FLAGS) -V -Nsrc $(XG)
CXX         = mpiFCCpx
CXXFLAGS    = $(CFLAGS)
FC          = mpifrtpx
FCFLAGS     = -Kfast,ocl,preex,simd=2,uxsimd,array_private,auto,parallel $(OMP_FLAGS) -Cpp -Qt -x-
F90         = mpifrtpx
F90FLAGS    = $(FCFLAGS)
LDFLAGS     = --linkfortran $(OMP_FLAGS)
LIBS        =
UDEF_LIB_PATH_SPEC =
UDEF_LIBS_SPEC     = -lz

## iff double
CFLAGS     += -D_REAL_IS_DOUBLE_
CXXFLAGS   += -D_REAL_IS_DOUBLE_
FCFLAGS    += -D_REAL_IS_DOUBLE_ -Ad
F90FLAGS   += -D_REAL_IS_DOUBLE_ -Ad

## iff large block
CFLAGS     += -D_LARGE_BLOCK
CXXFLAGS   += -D_LARGE_BLOCK
FCFLAGS    += -D_LARGE_BLOCK
F90FLAGS   += -D_LARGE_BLOCK

## iff different restart with staging
#CFLAGS     += -D_STAGING_
#CXXFLAGS   += -D_STAGING_
#FCFLAGS    += 
#F90FLAGS   += 
