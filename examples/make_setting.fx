##############################
#
# BCMTools
#
# Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
##############################

########################
# FX10
########################

MPI_DIR			= $(HOME)/MyLib/MPIlib
POLYLIB_DIR = $(HOME)/MyLib/Polylib-2.2
CUTLIB_DIR  = $(HOME)/MyLib/Cutlib-3.0.0_beta_mod
CUTLIB_DIR  = $(HOME)/MyLib/Cutlib-3.0.0_beta
TP_DIR      = $(HOME)/MyLib/TextParser-1.2
PMLIB_DIR   = $(HOME)/MyLib/PMlib-1.5
XML2_DIR    = $(HOME)/MyLib/libxml2-2.7.8_lib
BCM_DIR     = $(HOME)/BCMTools
BCM_FILEIO  = $(HOME)/BCMTools/FileIO
BCM_UTILS   = $(HOME)/BCMTools/Utils

CC          = mpifccpx
CXX         = mpiFCCpx
FC          = mpifrtpx
F90         = mpifrtpx

OPT_FLAGS   = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,optmsg=1 
OPT_FLAGS_F = -Kfast,ocl,preex,simd=2,uxsimd,array_private,parallel,auto,optmsg=1 
OMP_FLAGS   = -Kopenmp
DEFINES     = #-DDEBUG
UDEF_OPT    = #-D__K_FAPP_L0 -D__K_FAPP_L3
XG          = -Xg

CFLAGS      = $(OPT_FLAGS)   $(OMP_FLAGS) $(XG) -V -Nsrc
FCFLAGS     = $(OPT_FLAGS_F) $(OMP_FLAGS) $(XG) -Cpp -Qt -x-
CXXFLAGS    = $(CFLAGS)
F90FLAGS    = $(FCFLAGS)
LDFLAGS     = --linkfortran $(OMP_FLAGS)
LIBS        =

UDEF_LIB_PATH_SPEC =
UDEF_LIBS_SPEC     = -lz

AR          = ar cr
RANLIB      = ranlib
RM          = \rm -f

## iff double
CFLAGS     += -D_REAL_IS_DOUBLE_
CXXFLAGS   += -D_REAL_IS_DOUBLE_
FCFLAGS    += -D_REAL_IS_DOUBLE_ -CcdRR8
F90FLAGS   += -D_REAL_IS_DOUBLE_ -CcdRR8

## iff large block
CFLAGS     += -D_BLOCK_IS_LARGE
CXXFLAGS   += -D_BLOCK_IS_LARGE
FCFLAGS    += -D_BLOCK_IS_LARGE
F90FLAGS   += -D_BLOCK_IS_LARGE

## iff different restart with staging
#CFLAGS     += -D_STAGING_
#CXXFLAGS   += -D_STAGING_
#FCFLAGS    += 
#F90FLAGS   += 

