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

#include <iostream>
#include "mpi.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "BCMTools.h"
#include "Config.h"

#include "RootGrid.h"
#include "BCMOctree.h"
#include "SimpleDivider.h"
#include "Partition.h"

#include "BoundaryConditionSetter.h"

#include "BlockFactory.h"

#include "Block.h"
#include "BlockManager.h"

#include "Solver.h"

#include "Timing.h"

int main(int argc, char** argv)
{
  MPI::Init(argc, argv);
  MPI::Comm& comm = MPI::COMM_WORLD;
  int myrank = comm.Get_rank();

  Config conf;
  if (myrank == 0) {
    if (argc != 2) {
      std::cout << "usage: " << argv[0] << " configfile" << std::endl;
      comm.Abort(EX_USAGE);
    }

    std::cout << "# of MPI processes = " << comm.Get_size() << std::endl;
#ifdef _OPENMP
    std::cout << "# of OpenMP threads = " << omp_get_max_threads() << std::endl;
#endif

    std::cout <<  std::endl << "Configuration file: " << argv[1] << std::endl;
  }

  conf.load(argv[1]);
  if (myrank == 0) conf.print();

  RootGrid* rootGrid = new RootGrid(1, 1, 1);
  if (conf.type == 'x') {
    rootGrid->setPeriodicY();
    rootGrid->setPeriodicZ();
  }
  if (conf.type == 'y') {
    rootGrid->setPeriodicZ();
    rootGrid->setPeriodicX();
  }
  if (conf.type == 'z') {
    rootGrid->setPeriodicX();
    rootGrid->setPeriodicY();
  }

  Divider* divider;
  if (conf.treeType == "flat") {
    divider = new FlatDivider(rootGrid, conf.maxLevel);
  }
  else if (conf.treeType == "simple") {
    divider = new SimpleDivider(rootGrid, conf.minLevel, conf.maxLevel);
  }
  else {
    exit(EX_READ_CONFIG);
  }

  BCMOctree::Ordering ordering;
  if (conf.ordering == "Z") {
    ordering = BCMOctree::Z;
  }
  else if (conf.ordering == "Hilbert") {
    ordering = BCMOctree::HILBERT;
  }
  else if (conf.ordering == "random") {
    ordering = BCMOctree::RANDOM;
  }
  else {
    exit(EX_READ_CONFIG);
  }

  BCMOctree* tree = new BCMOctree(rootGrid, divider, ordering);

  int numLeafNode = tree->getNumLeafNode();

  std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();

  Partition* partition = new Partition(comm.Get_size(), numLeafNode);

  if (myrank == 0) {
    std::cout << std::endl << "Partitioning" << std::endl;
    partition->print();
  }

  // ブロック内のセル数
  Vec3i size(conf.size, conf.size, conf.size);

  BoundaryConditionSetter* boundaryConditionSetter = new BoundaryConditionSetter(&conf);

  BlockFactory* blockFactory = new BlockFactory(tree, partition, boundaryConditionSetter, size);

  BlockManager& blockManager = BlockManager::getInstance();

  for (int id = partition->getStart(myrank); id < partition->getEnd(myrank); id++) {
    Node* node = leafNodeArray[id];
    Block* block = blockFactory->makeBlock(node);
    blockManager.registerBlock(block);
  }

  blockManager.endRegisterBlock();

  blockManager.printBlockLayoutInfo();

  delete tree;
  delete partition;
  delete boundaryConditionSetter;

  std::vector<double> boundaryValue(1);
  boundaryValue[0] = conf.b0;
  boundaryValue[1] = conf.b1;

  Solver solver(conf, boundaryValue);

  Timing::start(INIT);
  solver.initialize();
  Timing::stop(INIT);

  Timing::start(RUN);
  solver.run();
  Timing::stop(RUN);

  solver.checkResult(conf.type, conf.verbose);

  if (myrank == 0) {
    std::cout << std::endl << "Timings" << std::endl;
    Timing::print(INIT,     "  Solver::initialize           ");
    Timing::print(RUN,      "  Solver::run                  ");
#ifdef TIMING
    Timing::print(SOR,      "  SOR in block                 ");
    Timing::print(VCUPDATE, "  VCUpdate + BC                ");
    Timing::print(BC,       "  BC                           ");
#endif
  }

  MPI::Finalize();

  return EX_SUCCESS;
}
