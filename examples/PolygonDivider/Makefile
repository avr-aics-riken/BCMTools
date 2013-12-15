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

PROGRAM = a.out

BCMLIB = ../../lib/libbcm.a
UTIL = ../../Utils

POLYLIB = $(HOME)/opt_intel/Polylib
TP = $(HOME)/opt_intel/TextParser
SILO = $(HOME)/opt_intel

POLYLIB_INC = -I$(POLYLIB)/include
POLYLIB_LIB = -L$(POLYLIB)/lib -lPOLY
TP_INC = -I$(TP)/include
TP_LIB = -L$(TP)/lib -lTP
SILO_INC = -I$(SILO)/include
SILO_LIB = -L$(SILO)/lib -lsiloh5 -lhdf5 -lz

vpath %.cpp $(UTIL)/src

OBJS = main.o Solver.o 
UTIL_OBJS = ConfigFile.o ConfigBase.o BCMOctree.o

include ../../Makefile.inc

CFLAGS += -I../../include -I$(UTIL)/include $(SILO_INC) $(POLYLIB_INC) \
           $(TP_INC) $(OPTIONS)

LDFLAGS = -L$(POLYLIB)/lib
LIBS = $(BCMLIB) $(SILO_LIB) $(POLYLIB_LIB) $(TP_LIB) -lmpi -lmpi_cxx

PROGRAM: $(OBJS) $(UTIL_OBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(UTIL_OBJS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(PROGRAM) depend.inc

depend: $(OBJS:.o=.cpp) $(UTIL_OBJS:.o=.cpp)
	@ rm -rf depend.inc
	@ for i in $^; do\
		$(CXX) $(CFLAGS) -MM $$i >> depend.inc;\
	done

-include depend.inc
