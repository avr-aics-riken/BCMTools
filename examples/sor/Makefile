PROGRAM = a.out

BCMLIB = ../../lib/libbcm.a

OBJS = main.o Config.o ConfigFile.o Octree.o MakeTree.o \
  BlockFactory.o Solver.o Timing.o


include ../../Makefile.inc

CFLAGS = -I../../include $(OPTIONS)


PROGRAM: $(OBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(BCMLIB) $(LIBS)

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
