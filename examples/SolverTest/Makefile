###################################################################
#
# FFV : Frontflow / violet BCM
#
# Copyright (c) 2013 All right reserved.
#
# Institute of Industrial Science, University of Tokyo, Japan. 
#
###################################################################

# Default environment
FFVC=intel

### compiler detection

ifeq ($(FFVC),gnu)
include ../../make_setting.gnu
endif

ifeq ($(FFVC),fx)
include ../../make_setting.fx
endif

ifeq ($(FFVC),ibm)
include ../../make_setting.ibm
endif

ifeq ($(FFVC),intel)
include ../../make_setting.intel
endif

#####

UTIL    = ../util

UDEF_INC_PATH=-I$(CUTLIB_DIR)/include \
              -I$(POLYLIB_DIR)/include \
              -I$(PMLIB_DIR)/include \
              -I$(TP_DIR)/include \
							-I$(BCM_DIR)/include \
							-I$(BCM_FILEIO_DIR)/include \
							-I$(XML2_DIR)/include/libxml2 \
							-I$(UTIL)/include \
              -I$(MPI_DIR)/include
              
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

vpath %.cpp $(UTIL)/src
vpath %.cpp $(BCM_FILEIO_DIR)/src

CSRC := $(wildcard *.cpp)
CHDR := $(wildcard *.h  )
FSRC := $(wildcard *.f90)
OBJS = $(CSRC:.cpp=.o) $(FSRC:.f90=.o)
UTIL_OBJS = ConfigFile.o ConfigBase.o BCMOctree.o BCMPolylib.o
BCM_FILEIO_OBJS = BitVoxel.o ErrorUtil.o BCMFileSaver.o LeafBlockSaver.o FileSystemUtil.o

PROGRAM: $(OBJS) $(UTIL_OBJS) $(BCM_FILEIO_OBJS)
	$(CXX) \
	$(LDFLAGS) \
	$(UDEF_OPT) \
	-o $(PROGRAM) $(OBJS) $(UTIL_OBJS) $(BCM_FILEIO_OBJS) \
	$(UDEF_LIB_PATH) $(UDEF_LIBS) \
	$(UDEF_LIB_PATH_SPEC) $(UDEF_LIBS_SPEC) \
	$(LIBS) 

.cpp.o:
	$(CXX) $(CXXFLAGS) $(UDEF_INC_PATH) $(UDEF_OPT) $(XG) -o $@ -c $<

.SUFFIXES: .f90
.f90.o:
	$(F90) $(F90FLAGS) $(UDEF_OPT) $(XG) -o $@ -c $<

clean:
	rm -rf *.o $(PROGRAM) depend.inc

cleandata:
	rm -rf BIN STL VTI data-*

depend:
	$(CXX) $(XG) -MM $(CSRC) $(CXXFLAGS) $(UDEF_INC_PATH) $(UDEF_OPT) > depend.inc

-include depend.inc

