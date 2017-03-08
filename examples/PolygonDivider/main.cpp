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

#include "mpi.h"
#include <iostream>
#ifdef _OPENMP
#include "omp.h"
#endif

#include "Polylib.h"

#include "Solver.h"
#include "BCMTools.h"
#include "Config.h"

#include "RootGrid.h"
#include "BCMOctree.h"
#include "PolygonDivider.h"
#include "Partition.h"

#include "BoundaryConditionSetter.h"

#include "BlockFactory.h"

#include "Block.h"
#include "BlockManager.h"

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

  RootGrid* rootGrid = new RootGrid(conf.rootN);

  PolylibNS::Polylib* polylib = PolylibNS::Polylib::get_instance();
  polylib->load(conf.polylibConf);

  Divider* divider = new PolygonDivider(conf.origin, conf.rootLength, rootGrid,
                                        conf.minLevel, conf.maxLevel,
                                        polylib, conf.polygonGroup,
                                        conf.polygonInsideOut);


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

  BlockFactory* blockFactory = new BlockFactory(tree, partition, boundaryConditionSetter,
                                                conf.origin, conf.rootLength, size);

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

  Solver solver(conf);

  solver.initialize();

  solver.run();

  MPI::Finalize();

  return EX_SUCCESS;
}
