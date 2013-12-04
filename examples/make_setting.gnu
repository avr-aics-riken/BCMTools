##############################
#
# BCMTools
#
# Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
##############################

########################
# INTEL (gcc, build openmpi with gcc)
########################

MPI_DIR			= $(HOME)/MyLib/MPIlib
POLYLIB_DIR = $(HOME)/MyLib/Polylib-2.2
CUTLIB_DIR  = $(HOME)/MyLib/Cutlib-3.0.0_beta
CUTLIB_DIR  = $(HOME)/MyLib/Cutlib-3.0.0_beta_mod
TP_DIR      = $(HOME)/MyLib/TextParser-1.2
PMLIB_DIR   = $(HOME)/MyLib/PMlib-1.6
XML2_DIR    = $(HOME)/MyLib/libxml2-2.7.8
BCM_DIR     = $(HOME)/BCMTools
BCM_FILEIO  = $(HOME)/BCMTools/FileIO
BCM_UTILS   = $(HOME)/BCMTools/Utils

CC          = gcc
CXX         = g++
FC          = gfortran
F90         = gfortran

OPT_FLAGS   = -O3
OMP_FLAGS   = -fopenmp
DEFINES     = #-DDEBUG
UDEF_OPT    =

CFLAGS      = $(OPT_FLAGS) $(OMP_FLAGS) $(DEFINES) -g -Wall -Wno-sign-compare -Wunknown-pragmas 
CXXFLAGS    = $(CFLAGS)
FCFLAGS     = $(OPT_FLAGS) $(OMP_FLAGS) $(DEFINES) -cpp -ffree-form 
F90FLAGS    = $(FCFLAGS)
LDFLAGS     = $(OMP_FLAGS)
LIBS        = -lgfortran -lmpi -lmpi_cxx -lmpi_mpifh
LIBS        = -lgfortran -lmpi -lmpi_cxx -lmpi_f77 -lmpi_f90
LIBS        = -lgfortran -lmpi -lmpi_cxx -lmpi_f77 -lmpi_f90 -L/opt/intel/lib/intel64 -lifport -lifcore -lirc -limf

UDEF_LIB_PATH_SPEC = 
UDEF_LIBS_SPEC     =

AR          = ar cru
RANLIB      = ranlib
RM          = \rm -f

## iff double
CFLAGS     += -D_REAL_IS_DOUBLE_
CXXFLAGS   += -D_REAL_IS_DOUBLE_
FCFLAGS    += -D_REAL_IS_DOUBLE_ -fdefault-real-8
F90FLAGS   += -D_REAL_IS_DOUBLE_ -fdefault-real-8

## iff large block
CFLAGS     += -D_BLOCK_IS_LARGE
CXXFLAGS   += -D_BLOCK_IS_LARGE
FCFLAGS    += -D_BLOCK_IS_LARGE
F90FLAGS   += -D_BLOCK_IS_LARGE

