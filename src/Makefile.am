###
##############################
#
# BCMTools
#
# Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
##############################
###

lib_LIBRARIES = libBCM.a

libBCM_a_CXXFLAGS = @MPI_CFLAGS@ @PL_CFLAGS@ @TP_CFLAGS@ @BLOCK_OPT@ -I$(top_builddir)/include -I../Utils/include

libBCM_a_SOURCES = \
    BlockManager.cpp \
    CommBuffer.cpp \
    Scalar3D.cpp \
    Scalar3DUpdater.cpp \
    Vector3D.cpp \
    Vector3DUpdater.cpp \
		../Utils/src/BCMOctree.cpp \
		../Utils/src/BCMPolylib.cpp

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
    $(top_builddir)/include/bcm_Version.h \
    $(top_builddir)/Utils/include/BBDivider.h \
    $(top_builddir)/Utils/include/BCMOctree.h \
    $(top_builddir)/Utils/include/BCMPolylib.h \
    $(top_builddir)/Utils/include/Block.h \
    $(top_builddir)/Utils/include/BlockFactory.h \
    $(top_builddir)/Utils/include/BoundaryConditionSetterBase.h \
    $(top_builddir)/Utils/include/BoundaryInfo.h \
    $(top_builddir)/Utils/include/BoundingBox.h \
    $(top_builddir)/Utils/include/ConfigBase.h \
    $(top_builddir)/Utils/include/ConfigFile.h \
    $(top_builddir)/Utils/include/CylinderDivider.h \
    $(top_builddir)/Utils/include/Divider.h \
    $(top_builddir)/Utils/include/MultiRootDivider.h \
    $(top_builddir)/Utils/include/Node.h \
    $(top_builddir)/Utils/include/Partition.h \
    $(top_builddir)/Utils/include/Pedigree.h \
    $(top_builddir)/Utils/include/PolygonBBoxDivider.h \
    $(top_builddir)/Utils/include/PolygonDivider.h \
    $(top_builddir)/Utils/include/RectangleDivider.h \
    $(top_builddir)/Utils/include/RootGrid.h \
    $(top_builddir)/Utils/include/SiloWriter.h \
    $(top_builddir)/Utils/include/SimpleDivider.h \
    $(top_builddir)/Utils/include/SphereDivider.h
    
EXTRA_DIST = Makefile_hand
