##############################
#
# BCMViewer - BCM mesh viewer
#
# Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
##############################

UNAME = $(shell uname -s)
PTYPE = $(shell uname -p) 

#CXX        = llvm-g++
CXX        = g++

RM         = rm -rf

TARGET     = main

UTIL       = ./UTIL
BCMTools   = ./BCMTools
TEXTPARSER = ../../ext/TextParser

DEFINES   += -D_DEBUG
OMP_FLAG   = -fopenmp

INC_DIR   += -I./
INC_DIR   += -I$(UTIL)
INC_DIR   += -I$(BCMTools)
INC_DIR   += -I$(TEXTPARSER)/include

LIB_DIR   += 

ifeq ($(UNAME), Linux)
LIB_FLAGS += -lglut
LIB_FLAGS += -lGL
LIB_FLAGS += -lGLU
else
LIB_FLAGS += -framework OpenGL
LIB_FLAGS += -framework GLUT
endif
LIB_FLAGS += $(TEXTPARSER)/lib/libTextParser.a

LIB_FLAGS += 

CFLAGS    += $(OMP_FLAG)
CFLAGS    += $(DEFINES)
CFLAGS    += $(INC_DIR)
CFLAGS    += -g

LFLAGS    += $(OMP_FLAG)
LFLAGS    += $(LIB_DIR)
LFLAGS    += $(LIB_FLAGS)

OBJECTS   += $(BCMTools)/BCMOctree.o
OBJECTS   += $(BCMTools)/LeafBlockLoader.o
OBJECTS   += $(BCMTools)/BitVoxel.o
OBJECTS   += $(BCMTools)/ErrorUtil.o
OBJECTS   += $(UTIL)/ShaderProgramObject.o
OBJECTS   += $(UTIL)/TextureObject.o
OBJECTS   += main.o
OBJECTS   += GridBCM.o
OBJECTS   += Render.o

all : $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CXX) $(OBJECTS) $(LFLAGS) -o $(TARGET)

$(OBJECTS) : %.o : %.cpp
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	$(RM) $(TARGET) $(OBJECTS)

