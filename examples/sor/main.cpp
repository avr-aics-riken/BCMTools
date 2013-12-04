#include <iostream>
#include "mpi.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "Solver.h"
#include "BCMTools.h"
#include "Config.h"
#include "Octree.h"
#include "MakeTree.h"
#include "Partition.h"
#include "Block.h"
#include "BlockFactory.h"
#include "BlockManager.h"
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
    conf.load(argv[1]);
    conf.print();
  }
  conf.bcast();
  
  Octree* tree = makeTree(conf.treeType, conf.level);

  int numLeafNode = tree->getNumLeafNode();

  if (conf.randomShuffle) tree->randomShuffle();

  std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();

  Partition partition(comm.Get_size(), numLeafNode);

  if (myrank == 0) {
    std::cout << std::endl << "Partitioning" << std::endl;
    partition.print();
  }

  // ブロック内のセル数
  Vec3i size(conf.size, conf.size, conf.size);

  BlockManager& blockManager = BlockManager::getInstance();

  for (int id = partition.getStart(myrank); id < partition.getEnd(myrank); id++) {
    Node* node = leafNodeArray[id];
    int level = node->pedigree.level;
    Vec3r origin = BlockFactory::getOrigin(node);
    Vec3r blockSize = BlockFactory::getBlockSize(node);
    BoundaryInfo* boundaryInfo = BlockFactory::makeBoundaryInfo(node, tree, &conf);
    NeighborInfo* neighborInfo = BlockFactory::makeNeighborInfo(node, tree,
                                                                boundaryInfo,
                                                                &partition);
    Block* block = new Block(size, origin, blockSize, level,
                             neighborInfo, boundaryInfo);
    blockManager.registerBlock(block);

//  BoundaryInfo::print(boundaryInfo);
  }

  blockManager.endRegisterBlock();

  blockManager.printBlockLayoutInfo();

  delete tree;

  std::vector<double> boundaryValue(2);
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
