#include <iostream>
#include "mpi.h"
#include "BCMTools.h"
#include "Config.h"
#include "RootGrid.h"
#include "BCMPolylib.h"
#include "BCMOctree.h"
#include "PolygonBBoxDivider.h"
#include "BoundaryConditionSetter.h"
#include "BlockFactory.h"
#include "Block.h"
#include "BlockManager.h"
#include "BlockBoundingBox.h"
#include "SiloWriter.h"

int main(int argc, char** argv)
{
  MPI::Init(argc, argv);
  MPI::Comm& comm = MPI::COMM_WORLD;
  int myRank = comm.Get_rank();

  Config conf;
  BCMOctree* tree;
  PolylibNS::BCMPolylib* pl = new PolylibNS::BCMPolylib;

  if (myRank == 0) {
    if (argc != 2) {
      std::cout << "usage: " << argv[0] << " configfile" << std::endl;
      comm.Abort(EX_USAGE);
    }

    std::cout << "# of MPI processes = " << comm.Get_size() << std::endl;

    std::cout <<  std::endl << "Configuration file: " << argv[1] << std::endl;

  }

  conf.load(argv[1]);

  if (myRank == 0) {
    conf.print();

    RootGrid* rootGrid = new RootGrid(conf.rootN);

    if (conf.polygonGroupList.size() > 0) {
      if (pl->load(conf.polylibConfig) != PLSTAT_OK) {
        comm.Abort(EX_OPEN_FILE);
      }
    }

  //// 仮想セル領域も考慮せずにポリゴン検索
  //Divider* divider = new PolygonBBoxDivider(conf.origin, conf.rootLength,
  //                                          rootGrid, conf.minLevel, pl,
  //                                          conf.polygonGroupList,
  //                                          conf.boundingBoxList);
    // 仮想セル領域も考慮してポリゴン検索
    Divider* divider = new PolygonBBoxDivider(conf.origin, conf.rootLength,
                                              rootGrid, conf.minLevel, pl,
                                              conf.polygonGroupList,
                                              conf.boundingBoxList,
                                              (double)conf.vc/conf.size);

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

    tree = new BCMOctree(rootGrid, divider, ordering);
    tree->broadcast();
  }
  else {
    tree = BCMOctree::ReceiveFromMaster();
  }

  int numLeafNode = tree->getNumLeafNode();
  std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();

  Partition* partition = new Partition(comm.Get_size(), numLeafNode);
  if (myRank == 0) {
    std::cout << std::endl << "Partitioning" << std::endl;
    partition->print();
  }

  // ブロック内のセル数
  Vec3i size(conf.size, conf.size, conf.size);

  BoundaryConditionSetter* boundaryConditionSetter = new BoundaryConditionSetter(&conf);

  BlockFactory* blockFactory = new BlockFactory(tree, partition, boundaryConditionSetter,
                                                conf.origin, conf.rootLength, size);

  BlockManager& blockManager = BlockManager::getInstance();

  for (int id = partition->getStart(myRank); id < partition->getEnd(myRank); id++) {
    Node* node = leafNodeArray[id];
    Block* block = blockFactory->makeBlock(node);
    blockManager.registerBlock(block);
  }

  blockManager.endRegisterBlock();

  blockManager.printBlockLayoutInfo();

  SiloWriter writer(conf.output, "mesh");
  writer.writeDomain("block_mesh", "domain");


  // ポリゴンの分配
  if (myRank == 0) {
    BlockBoundingBox bbb(tree, conf.origin, conf.rootLength, size, conf.vc);
    for (int iRank = 0; iRank < comm.Get_size(); iRank++) {

      // ランクiRankに属するブロック群のバインディングボックスを決定
      BoundingBox box;
      for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
        Node* node = leafNodeArray[id];
        box.addBox(bbb.getBoundingBox(node));
      }
      // std::cout << iRank << ": " << box << std::endl;

      // Polylibにバウンディングボックス情報をと登録
      pl->set_bounding_box(iRank, box);
    }

    // ポリゴンデータを各ランクに送信
    pl->send_to_all();
  }
  else {
    // ポリゴンデータを受信して、Polylibを初期化
    pl->load_from_rank0();
  }

  // 各ランクで独立にPolylibが使えることのテスト
  std::string file;
  pl->save_parallel(&file, "stl_b");

  delete pl;
  delete tree;
  delete partition;
  delete boundaryConditionSetter;


  MPI::Finalize();

  return EX_SUCCESS;
}
