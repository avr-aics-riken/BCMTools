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
# INTEL (gcc, build openmpi with gcc)
########################

AR          = ar cru
RANLIB      = ranlib
RM          = \rm -f
MPI_DIR	    = /usr/lib/openmpi
#MPI_DIR     = /opt/mpi/openmpi-1.6.2-gcc-v4.4.6
#MPI_DIR	    = /usr/local
BCM_DIR     = $(HOME)/BCMTools
BCM_FILEIO_DIR     = $(HOME)/BCMTools/examples/FileIO
CPM_DIR     = $(HOME)/MyLib/CPMlib
TP_DIR      = $(HOME)/MyLib/TPlib
POLYLIB_DIR = $(HOME)/MyLib/Polylib_2_0_3_Rel
CUTLIB_DIR  = $(HOME)/MyLib/Cutlib-2.0.5
PMLIB_DIR   = $(HOME)/MyLib/PMlib-1.5
XML2_DIR    = $(HOME)/MyLib/libxml2-2.7.8_lib
DEFINES     = -DDEBUG
OMP_FLAGS   = -fopenmp 
UDEF_OPT    =
CC          = mpicc
CFLAGS      = $(DEFINES) -DTIMING -O3 -ftree-vectorize -g -Wall -Wno-sign-compare -Wunknown-pragmas $(OMP_FLAGS)
CXX         = mpicxx
CXXFLAGS    = $(CFLAGS)
FC          = mpif90
FCFLAGS     = -O3 -ffree-form $(OMP_FLAGS) -cpp
F90         = mpif90
F90FLAGS    = $(FCFLAGS)
LDFLAGS     = $(OMP_FLAGS)
LIBS        = -lgfortran -L/opt/intel/lib/intel64 -lifport -lifcore -lirc -limf
#LIBS        = -lgfortran 
UDEF_LIB_PATH_SPEC = 
UDEF_LIBS_SPEC     =

## iff double
CFLAGS     += -D_REAL_IS_DOUBLE_
CXXFLAGS   += -D_REAL_IS_DOUBLE_
FCFLAGS    += -D_REAL_IS_DOUBLE_ -fdefault-real-8
F90FLAGS   += -D_REAL_IS_DOUBLE_ -fdefault-real-8

