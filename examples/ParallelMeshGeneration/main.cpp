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

#include "Cutlib.h"
#include "CutInfo/CutInfo.h"
#include "CutInfo/CutNormalArray.h"
#include "GridAccessor/Cell.h"

int main(int argc, char** argv)
{
  MPI::Init(argc, argv);
  MPI::Comm& comm = MPI::COMM_WORLD;
  int myRank = comm.Get_rank();

  Config conf;
  BCMOctree* tree;

  if (myRank == 0) {
    if (argc != 2) {
      std::cout << "usage: " << argv[0] << " configfile" << std::endl;
      comm.Abort(EX_USAGE);
    }

    std::cout << "# of MPI processes = " << comm.Get_size() << std::endl;

    std::cout <<  std::endl << "Configuration file: " << argv[1] << std::endl;

  }

  conf.load(argv[1]);

  PolylibNS::BCMPolylib* pl = new PolylibNS::BCMPolylib;
  if (conf.polygonGroupList.size() > 0) {
    if (pl->load(conf.polylibConfig) != PLSTAT_OK) {
      comm.Abort(EX_OPEN_FILE);
    }
  }

  if (myRank == 0) {
    conf.print();

    RootGrid* rootGrid = new RootGrid(conf.rootN);


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
  } else {
    // ポリゴンデータを受信して、Polylibを初期化
    pl->load_from_rank0();
  }

  // 各ランクで独立にPolylibが使えることのテスト
  std::string file;
  pl->save_parallel(&file, "stl_b");


  // Cutlibのテスト
	for(int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		Vec3i size      = block->getSize();
		Vec3r origin    = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize  = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {conf.vc};

		double bpos[3] = {origin.x, origin.y, origin.z};
		unsigned int bbsize[3] = {size.x, size.y, size.z};
		unsigned int gcsize[3] = {conf.vc, conf.vc, conf.vc};
		double dx[3] = {cellSize.x, cellSize.x, cellSize.x};
		size_t ncell[3];
		double org[3];
		for(int i=0; i<3; i++) {
			ncell[i] = bbsize[i] + 2*gcsize[i];
			org[i] = bpos[i] - gcsize[i]*dx[i];
		}

		cutlib::GridAccessor*   grid   = new cutlib::Cell(org, dx);
		cutlib::CutPosArray*    cutPos = new cutlib::CutPos32Array(ncell);
		cutlib::CutBidArray*    cutBid = new cutlib::CutBid5Array(ncell);
		cutlib::CutNormalArray* cutNormal = new cutlib::CutNormalArray(ncell);

		int ret = CalcCutInfo(grid, pl, cutPos, cutBid, cutNormal);

		int nNormal = cutNormal->getNumNormal();

		std::cout << "CalcCutInfo: " << ret << " " << nNormal << " " << myRank << std::endl;
	}


  delete pl;
  delete tree;
  delete partition;
  delete boundaryConditionSetter;


  MPI::Finalize();

  return EX_SUCCESS;
}
