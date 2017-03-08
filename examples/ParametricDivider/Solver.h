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

};

#endif // SOLVER_H
