POLYLIB = $(HOME)/FFV/Polylib
TP      = $(HOME)/FFV/TextParser
SILO    = /usr

PROGRAM = a.out

BCMLIB = ../../lib/libbcm.a
UTIL = ../../Utils

POLYLIB_INC = -I$(POLYLIB)/include
POLYLIB_LIB = -L$(POLYLIB)/lib -lMPIPOLY
TP_INC = -I$(TP)/include
TP_LIB = -L$(TP)/lib -lMPITP
SILO_INC = -I$(SILO)/include
SILO_LIB = -L$(SILO)/lib -lsiloh5 -lhdf5 -lz

vpath %.cpp $(UTIL)/src

OBJS = main.o 
UTIL_OBJS = ConfigFile.o ConfigBase.o BCMOctree.o BCMPolylib.o

include ../../Makefile.inc

CFLAGS += -I../../include -I$(UTIL)/include $(SILO_INC) $(POLYLIB_INC) $(TP_INC) -I$(MPI_DIR)/include $(OPTIONS)
LIBS = $(BCMLIB) -lifport -lifcore          $(SILO_LIB) $(POLYLIB_LIB) $(TP_LIB) -L$(MPI_DIR)/lib -lmpi -lmpi_cxx -lmpi_f77 -lmpi_f90

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
