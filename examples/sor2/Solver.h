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

  const std::vector<double>& boundaryValue;

  int nLoopInner;
  int nLoopOuter;

  double omega;

  bool separateVCUpdate;

  const std::string& outputFileName;

  Vec3i size;
  int nx, ny, nz;

  int id_f;
  static const int tag_f = 100;

  Scalar3D<double>* work;

public:

  Solver(const Config& conf, const std::vector<double>& boundaryValue);

  ~Solver();

  void initialize();

  void run();

  void checkResult(bool verbose);

private:

  void setInitialCondition();

  void setBoundaryCondition();

  void calcSorInBlock(int nx, int ny, int nz,
                      double* fData, const double* sData, Index3DS fIndex,
                      double omega,
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

  void setSource(int nx, int ny, int nz,
                       double* sData, Index3DS sIndex,
                       const Vec3d& orig, const Vec3d& delta);

  void dumpDataClass(int dataClassID);


};

#endif // SOLVER_H
