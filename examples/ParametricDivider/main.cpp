#include <iostream>
#include "mpi.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "Solver.h"
#include "BCMTools.h"
#include "Config.h"

#include "RootGrid.h"
#include "BCMOctree.h"
#include "SimpleDivider.h"
#include "SphereDivider.h"
#include "CylinderDivider.h"
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
  
  Divider* divider;
  if (conf.treeType == "flat") {
    divider = new FlatDivider(rootGrid, conf.maxLevel);
  }
  else if (conf.treeType == "simple") {
    divider = new SimpleDivider(rootGrid, conf.minLevel, conf.maxLevel);
  }
  else if (conf.treeType == "sphere") {
    divider = new SphereDivider(rootGrid, conf.minLevel, conf.maxLevel,
                                conf.cx, conf.cy, conf.cz, conf.r);
//                              conf.cx, conf.cy, conf.cz, conf.r, 1.0e-6);
  }
  else if (conf.treeType == "cylinder") {
    divider = new CylinderDivider(rootGrid, conf.minLevel, conf.maxLevel, 
                                  conf.cx, conf.cy, conf.r);
//                                conf.cx, conf.cy, conf.r, 1.0e-6);
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


  Solver solver(conf);

  solver.initialize();

  solver.run();

  MPI::Finalize();

  return EX_SUCCESS;
}
