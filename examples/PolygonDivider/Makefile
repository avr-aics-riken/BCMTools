PROGRAM = a.out

BCMLIB = ../../lib/libbcm.a
UTIL = ../util
POLYLIB = ../Polylib_2_0_3_Rel_search_nearest

vpath %.cpp $(UTIL)/src

OBJS = main.o Solver.o 
UTIL_OBJS = ConfigFile.o ConfigBase.o BCMOctree.o

include ../../Makefile.inc

CFLAGS = -I../../include -I$(UTIL)/include $(SILO_INC) \
         -I$(POLYLIB)/include -I/usr/include/libxml2 $(OPTIONS)

LDFLAGS = -L$(POLYLIB)/lib
LIBS = $(BCMLIB) $(SILO_LIB) -lPolylib -lxml2

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
