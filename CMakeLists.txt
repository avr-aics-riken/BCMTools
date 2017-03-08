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
#
#-*- CMAKE -*-#
# Process this file with cmake to produce a configure script.

# Following parameters must be given at cmake command line

# Build options, see Readme.md for details.
#
# -D INSTALL_DIR=install_directory
#
# -D with_MPI={yes|no}
#
# -D with_TP=installed_directory
#
# -D with_PL=installed_directory
#
# -D with_example={no|yes}
#
# -D enable_OPENMP={no|yes}
#
# -D real_type={float|double}
#
# -D enable_LARGE_BLOCK={no|yes}
#

cmake_minimum_required(VERSION 2.6)


#######
# Build options
#######

option (with_MPI "Enable MPI" "ON")
option (enable_OPENMP "Enable OpenMP" "OFF")
option (with_example "Compiling examples" "OFF")
option (with_PL "Enable Polylib" "OFF")
option (with_TP "Enable TextParser" "OFF")
option (enable_LARGE_BLOCK "Optimizing block size" "OFF")
option (real_type "Type of floating point" "float")


#######
# Project setting
#######

project(BCMTools C CXX)

if(INSTALL_DIR STREQUAL "")
  set(INSTALL_DIR "/usr/local/BCMTools")
endif()

# project()でCMAKE_INSTALL_PREFIX=/usr/localとなるので、再定義
set (CMAKE_INSTALL_PREFIX "${INSTALL_DIR}")

#set(CMAKE_BUILD_TYPE Release)



#######
set(PROJECT_VERSION "1.0.1")
set(LIB_REVISION "20170308_1400")
#######



#######
# Cmake module
#######

set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)
include(${PROJECT_SOURCE_DIR}/cmake/CompileOptionSelector.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/bcm_package.cmake)




#######
# Options
#######

AddOptimizeOption()

checkOpenMP()

if(enable_LARGE_BLOCK)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_BLOCK_IS_LARGE_")
endif()

precision()


#######
# Display options
#######

message(" ")
message( STATUS "Destination PATH: "               ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_C_COMPILER: "               ${CMAKE_C_COMPILER})
message( STATUS "CMAKE_C_FLAGS: "                  ${CMAKE_C_FLAGS})
message( STATUS "CMAKE_CXX_COMPILER: "             ${CMAKE_CXX_COMPILER})
message( STATUS "CMAKE_CXX_FLAGS: "                ${CMAKE_CXX_FLAGS})

if(enable_Fortran)
  message( STATUS "CMAKE_Fortran_COMPILER: "         ${CMAKE_Fortran_COMPILER})
  message( STATUS "CMAKE_Fortran_FLAGS: "            ${CMAKE_Fortran_FLAGS})
endif()

message(" ")
message( STATUS "Type of floating point: "      ${real_type})
message( STATUS "OpenMP support        : "      ${enable_OPENMP})
message( STATUS "MPI support           : "      ${with_MPI})
message( STATUS "TextParser support    : "      ${with_TP})
message( STATUS "Polylib support       : "      ${with_PL})
message( STATUS "Example               : "      ${with_example})
message( STATUS "LARGE_BLOCK support   : "      ${enable_LARGE_BLOCK})
message(" ")

if(CMAKE_CXX_COMPILER MATCHES ".*FCCpx$")
else()
  if(with_MPI)
    include(FindMPI)

    if(MPI_CXX_FOUND)
      message("MPI-CXX-CMP  : ${MPI_CXX_COMPILER}")
      message("MPI-CXX-FLAGS: ${MPI_CXX_COMPILE_FLAGS}")
      message("MPI-CXX-INC  : ${MPI_CXX_INCLUDE_PATH}")
      message("MPI-LNK-FLAGS: ${MPI_CXX_LINK_FLAGS}")
      message("MPI-LIBS     : ${MPI_CXX_LIBRARIES}")
    endif(MPI_CXX_FOUND)
  endif()
endif()


message(" ")
message("PROJECT : BCMTools " ${PROJECT_VERSION} " : " ${LIB_REVISION})
message(" ")



#######
# Find Libraries to depend on
#######

# TextParser
IF(with_TP)
  SET(TP_DIR "${with_TP}")
  SET(TP_INC "${TP_DIR}/include")
  SET(TP_LIB "${TP_DIR}/lib")
ELSE()
  MESSAGE("Error: can not find TextParser.")
ENDIF()


# Polylib
IF(with_PL STREQUAL "no")
ELSE()
  SET(PL_DIR "${with_PL}")
  SET(PL_INC "${PL_DIR}/include")
  SET(PL_LIB "${PL_DIR}/lib")
ENDIF()


add_definitions(-DHAVE_CONFIG_H)



#######
# Example
#######

if (with_example)
  enable_testing()
  message("Code test will be performed")
else()
  message("Code test is skipped")
endif()



#######
# Check header files
#######

include(CheckIncludeFiles)

CHECK_INCLUDE_FILES(inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(strings.h HAVE_STRINGS_H)
CHECK_INCLUDE_FILES(string.h HAVE_STRING_H)
CHECK_INCLUDE_FILES(sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD_H)



#######
# Sub Directories
#######

add_subdirectory(src)

add_subdirectory(doc)

if(with_example)
  add_subdirectory(examples)
endif()


#######
# configure files
#######

configure_file( config.h.in.cmake config.h )
configure_file( include/bcm_Version.h.in ${PROJECT_BINARY_DIR}/include/bcm_Version.h @ONLY)


#######
# Install files
#######

install(FILES  ChangeLog.md
               Readme.md
               License.txt
        DESTINATION share)


#######
# Installer
#######

bcm_package()
