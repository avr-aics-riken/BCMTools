/*
 * BCMTools
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef SOLVER_H
#define SOLVER_H

#include <string>
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
                      double* fData, Index3DS fIndex,
                      double c0, double cx, double cy, double cz);

  void setInitialConditionInBlock(int nx, int ny, int nz,
                                  double* fData, Index3DS fIndex,
                                  double value);

  void setBoundaryConditionInBlock(int nx, int ny, int nz,
                                   double* fData, Index3DS fIndex,
                                   const BoundaryInfo* boundaryInfo);

  void setDirichletBoundaryInBlock(int nx, int ny, int nz,
                                  double* fData, Index3DS fIndex,
                                  Face face, double value);

  void dumpDataClass(int dataClassID);

#endif
};

#endif // SOLVER_H
