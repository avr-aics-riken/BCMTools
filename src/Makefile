BCMLIB = ../lib/libbcm.a

DATACLASS_LIST = ./DataClass.list

all: $(BCMLIB)

OBJS = CommBuffer/CommBuffer.o Block/BlockManager.o \
  DataClass/Scalar3D.o DataClass/Vector3D.o DataClass/Scalar3DUpdater.o DataClass/Vector3DUpdater.o

include ../Makefile.inc


#AR = ar cr
AR = ar cru
RANLIB = ranlib


DATACLASS_OPT = $(shell awk '/^ *[^\#]/{ printf " -DUSE_%s", $$1}' $(DATACLASS_LIST))

CFLAGS = -I. -I../include $(OPTIONS) $(DATACLASS_OPT)

DataClasss/Scalar3D.o: $(DATACLASS_LIST)
DataClasss/Vector3D.o: $(DATACLASS_LIST)
DataClasss/Scalar3DUpdater.o: $(DATACLASS_LIST)
DataClasss/Vector3DUpdater.o: $(DATACLASS_LIST)


$(BCMLIB): $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

.cpp.o:
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf $(OBJS) $(BCMLIB) depend.inc

depend: $(OBJS:.o=.cpp)
	@ rm -rf depend.inc
	@ for i in $^; do\
		$(CXX) $(CFLAGS) -MM $$i >> depend.inc;\
	done

-include depend.inc
