###################################################################################
#
# BCMTools
#
# Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2016 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################

# - Try to find BCMTools
# Once done, this will define
#
#  BCM_FOUND - system has BCMTools
#  BCM_INCLUDE_DIRS - BCMTools include directories
#  BCM_LIBRARIES - link these to use BCMTools

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(BCM_PKGCONF BCM)

if(CMAKE_PREFIX_PATH)
  set(BCM_CANDIDATE_PATH ${CMAKE_PREFIX_PATH})
  file(GLOB tmp "${CMAKE_PREFIX_PATH}/[Jj][Hh][Pp][Cc][Nn][Dd][Ff]*/")
  list(APPEND BCM_CANDIDATE_PATH ${tmp})
endif()

# Include dir
find_path(BCM_INCLUDE_DIR
  NAMES BCMTools.h
  PATHS ${BCM_ROOT} ${BCM_PKGCONF_INCLUDE_DIRS} ${BCM_CANDIDATE_PATH}
  PATH_SUFFIXES include
)

# Finally the library itself
find_library(BCM_LIBRARY
  NAMES BCM
  PATHS ${BCM_ROOT} ${BCM_PKGCONF_LIBRARY_DIRS} ${BCM_CANDIDATE_PATH}
  PATH_SUFFIXES lib
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(BCM_PROCESS_INCLUDES BCM_INCLUDE_DIR)
set(BCM_PROCESS_LIBS BCM_LIBRARY)
libfind_process(BCM)
