PROGRAM = a.out

BCMLIB = ../../lib/libbcm.a


OBJS = main.o Config.o ConfigFile.o Octree.o MakeTree.o \
  BlockFactory.o Solver.o


include ../../Makefile.inc

CFLAGS = -I../../include $(SILO_INC) $(OPTIONS)

LDFLAGS =
LIBS = $(BCMLIB) $(SILO_LIB)


PROGRAM: $(OBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(PROGRAM) depend.inc

depend: $(OBJS:.o=.cpp)
	@ rm -rf depend.inc
	@ for i in $^; do\
		$(CXX) $(CFLAGS) -MM $$i >> depend.inc;\
	done

-include depend.inc
