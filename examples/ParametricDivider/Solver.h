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
