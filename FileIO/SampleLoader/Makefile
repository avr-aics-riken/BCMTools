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

#
# Makefile for BCMTools FileIO/SampleLoader
#


include ../../../Makefile.inc

RM         = rm -rf

PROGRAM    = main

BCMTOOLS   = ../../..
FILEIO     = ..
UTIL       = ../../util
TEXTPARSER = ../../../../ext/TextParser

INC_DIR   += -I$(BCMTOOLS)/include
INC_DIR   += -I$(UTIL)/include
INC_DIR   += -I$(TEXTPARSER)/include
INC_DIR   += -I$(FILEIO)/include

LIB_DIR   += -L$(BCMTOOLS)/lib
LIB_DIR   += -L$(LIBXML2)/lib64
LIB_FLAGS += -lbcm
LIB_FLAGS += $(TEXTPARSER)/lib/libTextParser_mpi.a

CFLAGS    += $(DEFINES)
CFLAGS    += $(INC_DIR)
CFLAGS    += $(OPTIONS)

LFLAGS    += $(LIB_DIR)
LFLAGS    += $(LIB_FLAGS)
LFLAGS    += $(OPTIONS)

# UTIL
OBJECTS   += $(UTIL)/src/BCMOctree.o

# FILEIO
OBJECTS   += $(FILEIO)/src/BitVoxel.o
OBJECTS   += $(FILEIO)/src/ErrorUtil.o
OBJECTS   += $(FILEIO)/src/BCMFileLoader.o
OBJECTS   += $(FILEIO)/src/LeafBlockLoader.o
OBJECTS   += $(FILEIO)/src/BCMFileSaver.o
OBJECTS   += $(FILEIO)/src/LeafBlockSaver.o
OBJECTS   += $(FILEIO)/src/FileSystemUtil.o

# 
OBJECTS   += main.o


all : $(PROGRAM)

$(PROGRAM) : $(OBJECTS)
	$(CXX) $(OBJECTS) $(LFLAGS) -o $(PROGRAM)

$(OBJECTS) : %.o : %.cpp
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	$(RM) $(OBJECTS) $(PROGRAM)
