##############################
#
# BCMViewer - BCM mesh viewer
#
# Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
##############################

include ../Makefile.inc

BCMLIB = ../lib/libbcm.a
CFLAGS += -I. -I../include 

all: $(BCMLIB)

OBJS = CommBuffer/CommBuffer.o Block/BlockManager.o \
  DataClass/Scalar3D.o DataClass/Vector3D.o DataClass/Scalar3DUpdater.o DataClass/Vector3DUpdater.o

#include ../Makefile.inc
#DATACLASS_OPT = $(shell awk '/^ *[^\#]/{ printf " -DUSE_%s", $$1}' $(DATACLASS_LIST))
#CFLAGS = -I. -I../include $(OPTIONS) $(DATACLASS_OPT)
#
#DATACLASS_LIST = ./DataClass.list
#DataClasss/Scalar3D.o: $(DATACLASS_LIST)
#DataClasss/Vector3D.o: $(DATACLASS_LIST)
#DataClasss/Scalar3DUpdater.o: $(DATACLASS_LIST)
#DataClasss/Vector3DUpdater.o: $(DATACLASS_LIST)

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

