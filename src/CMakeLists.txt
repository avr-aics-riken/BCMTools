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
       ${TP_INC}
       ${PL_INC}
)


link_directories(
      ${TP_LIB}
      ${PL_LIB}
)


set(bcm_files
    BlockManager.cpp
    CommBuffer.cpp
    Scalar3D.cpp
    Scalar3DUpdater.cpp
    Vector3D.cpp
    Vector3DUpdater.cpp
    ../Utils/src/BCMOctree.cpp
    ../Utils/src/BCMPolylib.cpp
)


if(with_MPI)
  add_library(BCMmpi STATIC ${bcm_files})
  target_link_libraries(BCMmpi -lPOLYmpi -lTPmpi)
  install(TARGETS BCMmpi DESTINATION lib)
else()
  add_library(BCM STATIC ${bcm_files})
  target_link_libraries(BCM -lPOLY -lTP)
  install(TARGETS BCM DESTINATION lib)
endif()


install(FILES
        ${PROJECT_SOURCE_DIR}/include/BCMTools.h
        ${PROJECT_SOURCE_DIR}/include/BlockBase.h
        ${PROJECT_SOURCE_DIR}/include/BlockManager.h
        ${PROJECT_SOURCE_DIR}/include/CommBuffer.h
        ${PROJECT_SOURCE_DIR}/include/DataClass.h
        ${PROJECT_SOURCE_DIR}/include/Index3D.h
        ${PROJECT_SOURCE_DIR}/include/NeighborInfo.h
        ${PROJECT_SOURCE_DIR}/include/PointerSetter.h
        ${PROJECT_SOURCE_DIR}/include/Scalar3D.h
        ${PROJECT_SOURCE_DIR}/include/Scalar3DUpdater.h
        ${PROJECT_SOURCE_DIR}/include/VCUpdater.h
        ${PROJECT_SOURCE_DIR}/include/Vec3.h
        ${PROJECT_SOURCE_DIR}/include/Vector3D.h
        ${PROJECT_SOURCE_DIR}/include/Vector3DUpdater.h
        ${PROJECT_BINARY_DIR}/include/bcm_Version.h
        DESTINATION include
)
