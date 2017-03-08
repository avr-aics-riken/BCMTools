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
      ${PROJECT_BINARY_DIR}/include
      ${PROJECT_SOURCE_DIR}/include
      ${PROJECT_SOURCE_DIR}/Utils/include
      ${TP_INC}
      ${PL_INC}
      ${CT_INC}
      ${SILO_INC}
      ${H5_INC}
)

link_directories(
      ${PROJECT_BINARY_DIR}/src
      ${PROJECT_BINARY_DIR}/Utils
      ${CT_LIB}
      ${PL_LIB}
      ${TP_LIB}
      ${SILO_LIB}
      ${H5_LIB}
)



#### Sample1 : Mesh Generation

add_executable(mesh-gen
               ParallelMeshGeneration/main.cpp
               ${PROJECT_SOURCE_DIR}/Utils/src/ConfigBase.cpp
               ${PROJECT_SOURCE_DIR}/Utils/src/ConfigFile.cpp
)

if(with_MPI)
  target_link_libraries(mesh-gen -lBCMmpi -lCUT -lPOLYmpi -lTPmpi -lsiloh5 -lhdf5)
  set (test_parameters -np 4
                      "mesh-gen"
                      "${PROJECT_SOURCE_DIR}/examples/ParallelMeshGeneration/test.conf"
  )
  add_test(NAME TEST_1 COMMAND "mpirun" ${test_parameters}
  )
else()
#  target_link_libraries(mesh-gen -lPOLY -lTP)
endif()
