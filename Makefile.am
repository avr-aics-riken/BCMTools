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

SUBDIRS = @MAKE_SUB_DIRS@

bin_SCRIPTS= bcm-config
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = bcm-config.in LICENSE Utils/src examples BCMViewer Makefile.inc README.md 
dist_data_DATA = AUTHORS COPYING ChangeLog README LICENSE NEWS

clean-local:
	find ./ -name *.obj | xargs rm -f
	find ./ -name *.exe | xargs rm -f
	find ./ -name *.lib | xargs rm -f

distclean-local: clean-local
	rm -f ./bin/*
	rm -f ./lib/*

