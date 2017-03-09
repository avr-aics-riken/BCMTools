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

include_directories(
       ${PROJECT_BINARY_DIR}/include          # bcm_Version.h
       ${PROJECT_SOURCE_DIR}/include
       ${PROJECT_SOURCE_DIR}/Utils/include
)


set(conf_files
    src/ConfigBase.cpp
    src/ConfigFile.cpp
)

add_library(BCMconfig STATIC ${conf_files})
install(TARGETS BCMconfig DESTINATION lib)
