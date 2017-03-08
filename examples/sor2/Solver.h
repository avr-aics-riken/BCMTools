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

  const int vc;

  const std::vector<REAL_TYPE>& boundaryValue;

  int nLoopInner;
  int nLoopOuter;

  REAL_TYPE omega;

  bool separateVCUpdate;

  const std::string& outputFileName;

  Vec3i size;
  int nx, ny, nz;

  int id_f;
  static const int tag_f = 100;

  Scalar3D<REAL_TYPE>* work;

public:

  Solver(const Config& conf, const std::vector<REAL_TYPE>& boundaryValue);

  ~Solver();

  void initialize();

  void run();

  void checkResult(bool verbose);

private:

  void setInitialCondition();

  void setBoundaryCondition();

  void calcSorInBlock(int nx, int ny, int nz,
                      REAL_TYPE* fData, const REAL_TYPE* sData, Index3DS fIndex,
                      REAL_TYPE omega,
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

  void setSource(int nx, int ny, int nz,
                       REAL_TYPE* sData, Index3DS sIndex,
                       const Vec3r& orig, const Vec3r& delta);

  void dumpDataClass(int dataClassID);


};

#endif // SOLVER_H
