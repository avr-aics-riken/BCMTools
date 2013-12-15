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

SUBDIRS = @MAKE_SUB_DIRS@

bin_SCRIPTS= bcm-config
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = bcm-config.in LICENSE Utils examples FileIO BCMViewer Makefile.inc
dist_data_DATA = AUTHORS COPYING ChangeLog README README.md LICENSE NEWS

clean-local:
	find ./ -name *.obj | xargs rm -f
	find ./ -name *.exe | xargs rm -f
	find ./ -name *.lib | xargs rm -f

distclean-local: clean-local
	rm -f ./bin/*
	rm -f ./lib/*

