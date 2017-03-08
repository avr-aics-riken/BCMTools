/*
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
*/

#ifndef SOLVER_H
#define SOLVER_H

#include "mpi.h"
#include <string>
#include "BCMTools.h"
#include "Config.h"
#include "Block.h"
#include "BlockManager.h"
#include "Scalar3D.h"
#include "Vector3D.h"
#include "Vec3.h"

using namespace Vec3class;

class BlockManager;

class Solver {

  BlockManager& blockManager;

  const MPI::Intracomm& comm;

  const static int vc = 1;

  std::string outputFile;

  Vec3i size;
  int nx, ny, nz;

  int id_s;
  int id_v;


public:

  Solver(const Config& conf);

  ~Solver();

  void initialize();

  void run();

private:

#if 0
  void setInitialCondition();

  void setBoundaryCondition();

  void calcSorInBlock(int nx, int ny, int nz,
                      REAL_TYPE* fData, Index3DS fIndex,
                      REAL_TYPE c0, REAL_TYPE cx, REAL_TYPE cy, REAL_TYPE cz);

  void setInitialConditionInBlock(int nx, int ny, int nz,
                                  REAL_TYPE* fData, Index3DS fIndex,
                                  REAL_TYPE value);

  void setBoundaryConditionInBlock(int nx, int ny, int nz,
                                   REAL_TYPE* fData, Index3DS fIndex,
                                   const BoundaryInfo* boundaryInfo);

  void setDirichletBoundaryInBlock(int nx, int ny, int nz,
                                  REAL_TYPE* fData, Index3DS fIndex,
                                  Face face, REAL_TYPE value);

  void dumpDataClass(int dataClassID);

#endif
};

#endif // SOLVER_H
