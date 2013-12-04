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
POLYLIB = ../Polylib_2_0_3_Rel_search_nearest
XML2LIB = /home/jonishi/MyLib/libxml2-2.7.8

vpath %.cpp $(UTIL)/src

OBJS = main.o Solver.o 
UTIL_OBJS = ConfigFile.o ConfigBase.o BCMOctree.o

include ../../Makefile.inc

CFLAGS += -I../../include -I$(UTIL)/include $(SILO_INC) -I$(XML2LIB)/include \
         -I$(POLYLIB)/include -I/usr/include/libxml2 $(OPTIONS)

LDFLAGS = -L$(POLYLIB)/lib
LIBS = $(BCMLIB) $(SILO_LIB) -L$(POLYLIB)/lib -L$(XML2LIB)/lib -lPolylib -lxml2 -lmpi -lmpi_cxx -lsiloh5

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
