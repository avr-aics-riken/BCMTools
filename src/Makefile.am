###
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
###

lib_LIBRARIES = libBCM.a

libBCM_a_CXXFLAGS = @MPI_CFLAGS@ -I$(top_builddir)/include

libBCM_a_SOURCES = \
    BlockManager.cpp \
    CommBuffer.cpp \
    Scalar3D.cpp \
    Scalar3DUpdater.cpp \
    Vector3D.cpp \
    Vector3DUpdater.cpp

BCM_includedir = $(includedir)
BCM_include_HEADERS = \
    $(top_builddir)/include/BCMTools.h \
    $(top_builddir)/include/BlockBase.h \
    $(top_builddir)/include/BlockManager.h \
    $(top_builddir)/include/CommBuffer.h \
    $(top_builddir)/include/DataClass.h \
    $(top_builddir)/include/Index3D.h \
    $(top_builddir)/include/NeighborInfo.h \
    $(top_builddir)/include/PointerSetter.h \
    $(top_builddir)/include/Scalar3D.h \
    $(top_builddir)/include/Scalar3DUpdater.h \
    $(top_builddir)/include/VCUpdater.h \
    $(top_builddir)/include/Vec3.h \
    $(top_builddir)/include/Vector3D.h \
    $(top_builddir)/include/Vector3DUpdater.h \
    $(top_builddir)/include/bcm_Version.h
    
EXTRA_DIST = Makefile_hand
