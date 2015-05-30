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

# Default environment
FFVC=intel

### compiler detection

ifeq ($(FFVC),gnu)
include ../make_setting.gnu
endif

ifeq ($(FFVC),fx)
include ../make_setting.fx
endif

ifeq ($(FFVC),ibm)
include ../make_setting.ibm
endif

ifeq ($(FFVC),intel)
include ../make_setting.intel
endif

#####

#UTIL    = ../util

UDEF_INC_PATH=-I$(POLYLIB_DIR)/include \
              -I$(CUTLIB_DIR)/include \
              -I$(PMLIB_DIR)/include \
              -I$(TP_DIR)/include \
							-I$(XML2_DIR)/include/libxml2 \
              -I$(MPI_DIR)/include \
 							-I$(BCM_DIR)/include \
							-I$(BCM_FILEIO)/include \
							-I$(BCM_UTILS)/include 
             
UDEF_LIB_PATH=-L$(POLYLIB_DIR)/lib \
              -L$(CUTLIB_DIR)/lib \
              -L$(PMLIB_DIR)/lib \
              -L$(TP_DIR)/lib \
							-L$(BCM_DIR)/lib \
							-L$(XML2_DIR)/lib \
              -L$(MPI_DIR)/lib
              
UDEF_LIBS = -lbcm -lMPIPolylib -lcutlib -lpmlib -lTextParser_mpi -lxml2 -lmpi 

#####

PROGRAM = a.out

GIT_VERSION = $(shell sh -c 'git describe --always --dirty')
#GIT_DATE    = $(shell sh -c 'date --iso=seconds')
BUILD_DATE    = $(shell sh -c 'date +%s')

VERSION_FLAGS = -D__GIT_VERSION=\"$(GIT_VERSION)\" -D__BUILD_DATE=\"$(BUILD_DATE)\"
CXXFLAGS += $(VERSION_FLAGS)
CFLAGS   += $(VERSION_FLAGS)

vpath %.cpp $(BCM_FILEIO)/src
vpath %.cpp $(BCM_UTILS)/src

CSRC := $(wildcard *.cpp)
CHDR := $(wildcard *.h  )
FSRC := $(wildcard *.f90)
OBJS = $(CSRC:.cpp=.o) $(FSRC:.f90=.o)
BCM_UTILS_OBJS = ConfigFile.o ConfigBase.o BCMOctree.o BCMPolylib.o
BCM_FILEIO_OBJS = BitVoxel.o ErrorUtil.o BCMFileSaver.o LeafBlockSaver.o FileSystemUtil.o

PROGRAM: $(OBJS) $(BCM_UTILS_OBJS) $(BCM_FILEIO_OBJS)
	$(CXX) \
	$(LDFLAGS) \
	$(UDEF_OPT) \
	-o $(PROGRAM) $(OBJS) $(BCM_UTILS_OBJS) $(BCM_FILEIO_OBJS) \
	$(UDEF_LIB_PATH) $(UDEF_LIBS) \
	$(UDEF_LIB_PATH_SPEC) $(UDEF_LIBS_SPEC) \
	$(LIBS) 

.cpp.o:
	$(CXX) $(CXXFLAGS) $(UDEF_INC_PATH) $(UDEF_OPT) -o $@ -c $<

.SUFFIXES: .f90
.f90.o:
	$(F90) $(F90FLAGS) $(UDEF_INC_PATH) $(UDEF_OPT) -o $@ -c $<

clean:
	rm -rf *.o $(PROGRAM) depend.inc

cleandata:
	rm -rf BIN STL VTK PLOT3D BCM_OUT data-*

depend:
	$(CXX) -MM $(CSRC) $(CXXFLAGS) $(UDEF_INC_PATH) $(UDEF_OPT) > depend.inc

-include depend.inc

