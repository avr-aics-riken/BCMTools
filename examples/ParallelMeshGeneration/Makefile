include ../../Makefile.inc

POLYLIB := $(HOME)/FFV-exp/Polylib
TP      := $(HOME)/FFV-exp/TextParser
BCMLIB  := $(HOME)/FFV-exp/BCMTools
CUTLIB  := $(HOME)/FFV-exp/Cutlib
MPI_DIR := /opt/openmpi
SILO    := /usr
UTIL    := ../../Utils

PROGRAM = a.out

POLYLIB_INC = -I$(POLYLIB)/include -I$(POLYLIB)/include/common -D_REAL_IS_DOUBLE_
POLYLIB_LIB = -L$(POLYLIB)/lib -lPOLYmpi
TP_INC = -I$(TP)/include
TP_LIB = -L$(TP)/lib -lTPmpi
BCM_INC = -I$(BCMLIB)/include
BCM_LIB = -L$(BCMLIB)/lib -lBCM 
CUTLIB_INC = -I$(CUTLIB)/include -DCUTINFO_NORMAL_DOUBLE
CUTLIB_LIB = -L$(CUTLIB)/lib -lCUT
SILO_INC = -I$(SILO)/include
SILO_LIB = -L$(SILO)/lib -lsiloh5 -lhdf5 -lz

vpath %.cpp $(UTIL)/src

OBJS = main.o 
UTIL_OBJS = ConfigFile.o ConfigBase.o

CFLAGS += -I../../include -I$(UTIL)/include $(SILO_INC) $(POLYLIB_INC) $(TP_INC) $(BCM_INC) $(CUTLIB_INC) -I$(MPI_DIR)/include $(OPTIONS)
LIBS = -lifport -lifcore                    $(SILO_LIB) $(POLYLIB_LIB) $(TP_LIB) $(BCM_LIB) $(CUTLIB_LIB) -L$(MPI_DIR)/lib -lmpi -lmpi_cxx -lmpi_f77 -lmpi_f90

PROGRAM: $(OBJS) $(UTIL_OBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(UTIL_OBJS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(PROGRAM) depend.inc

cleandata:
	rm -rf *.silo *.stlb *.id polylib_config_*

depend: $(OBJS:.o=.cpp) $(UTIL_OBJS:.o=.cpp)
	@ rm -rf depend.inc
	@ for i in $^; do\
		$(CXX) $(CFLAGS) -MM $$i >> depend.inc;\
	done

-include depend.inc



