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

#include <mpi.h>
#include "Solver.h"

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "real.h"
#include "comm.h"

#include "RootGrid.h"
#include "BCMOctree.h"
#include "SimpleDivider.h"
#include "Partition.h"
#include "PolygonBBoxDivider.h"
#include "BBDivider.h"

#include "SphereDivider.h"
#include "SphereDivider2.h"
#include "SphereDivider3.h"

#include "RectangleDivider.h"
//#include "BoundaryConditionSetter.h"
#include "BlockFactory.h"
#include "Block.h"
#include "BlockManager.h"
#include "BlockBoundingBox.h"
#include "Polylib.h"
#include "Cutlib.h"
#include "CutInfo/CutInfo.h"
#include "GridAccessor/Cell.h"

#include "bcut.h"
#include "bstl.h"
#include "PM.h"

#include "outputVtk.h"

#define GLOBAL_VALUE_DEFINE
#include "gv.h"
#undef GLOBAL_VALUE_DEFINE

void PrintInt32t(int32_t i) {
	int m[32];
	for(int n=0; n<32; n++) {
		m[n] = (i >> n)%2;
	}
	for(int n=0; n<32; n++) {
		std::cout << m[31-n];
	}
}

Solver::Solver()
	: blockManager(BlockManager::getInstance()) {
}

Solver::~Solver() {
	MPI::Finalize();
}

int Solver::Init(int argc, char** argv){
/////////////////////////////////////////////
// Init MPI
/////////////////////////////////////////////
	MPI::Init(argc, argv);
	MPI::Comm& comm = MPI::COMM_WORLD;
	this->myrank = comm.Get_rank();
	if( argc != 2 ) {
		PrintLog(0, "usage: %s configfile", argv[0]);
		comm.Abort(EX_USAGE);
	}

	PrintLog(1, GetSolverName().c_str());
	PrintLog(2, "%-20s : %s", "Git hash", __GIT_VERSION);
	PrintLog(2, "%-20s : %s", "Build date", __BUILD_DATE);
	PrintLog(2, "%-20s : %s", "Configuration file", argv[1]);
	PrintLog(2, "%-20s : %d", "MPI processes", comm.Get_size());
#ifdef _OPENMP
	PrintLog(2, "%-20s : %d", "OpenMP threads", omp_get_max_threads());
#endif
/////////////////////////////////////////////


/////////////////////////////////////////////
// Init Config
/////////////////////////////////////////////
	PrintLog(1, "Loading configuration file");
	conf.load(argv[1]);
	g_pconf = &conf;
/*
	if(myrank == 0) {
		conf.print();
	}
*/

	PrintLog(2, "%-20s : %d", "min level", conf.minLevel);
	PrintLog(2, "%-20s : %d", "max level", conf.maxLevel);
	PrintLog(2, "%-20s : %s", "tree type", conf.treeType.c_str());
	PrintLog(2, "%-20s : %s", "ordering",  conf.ordering.c_str());
	PrintLog(2, "%-20s : %d", "block size", conf.size);
	PrintLog(2, "%-20s : %d", "vc width", conf.vc);

	g_pPM = new pm_lib::PerfMonitor;
	g_pPM->initialize(tm_END);
	g_pPM->setRankInfo(this->myrank);
	g_pPM->setProperties(tm_Init_LoadSTL,      "LoadSTL",      pm_lib::PerfMonitor::COMM, true);
	g_pPM->setProperties(tm_Init_DivideDomain, "DivideDomain", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Init_OrderBlock,   "OrderBlock",   pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_Update,   "Update",   pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateT,  "UpdateT",  pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateUX, "UpdateUX", pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateUY, "UpdateUY", pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateUZ, "UpdateUZ", pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateP,  "UpdateP",  pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_UpdateU,  "UpdateU",  pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_DOT,      "DOT",      pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_DOT_Calc, "DOT_Calc", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_DOT_Comm, "DOT_Comm", pm_lib::PerfMonitor::COMM, true);
	g_pPM->setProperties(tm_JacobiSmoother,      "JacobiSmoother",      pm_lib::PerfMonitor::CALC, false);
	g_pPM->setProperties(tm_JacobiSmoother_Calc, "JacobiSmoother_Calc", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_JacobiSmoother_Swap, "JacobiSmoother_Swap", pm_lib::PerfMonitor::CALC, true);
	g_pPM->setProperties(tm_JacobiSmoother_Comm, "JacobiSmoother_Comm", pm_lib::PerfMonitor::COMM, true);
	int nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
#endif
	g_pPM->setParallelMode("Hybrid", nThreads, comm.Get_size());
/////////////////////////////////////////////


/////////////////////////////////////////////
// Init BCM
/////////////////////////////////////////////
	PrintLog(1, "Loading STL file(s)");

	PM_Start(tm_Init_LoadSTL, 0, 0, false);
	PolylibNS::BCMPolylib* pl = new PolylibNS::BCMPolylib;
	if (conf.polygonGroupList.size() > 0) {
		if (pl->load(conf.polylibConfig) != PLSTAT_OK) {
			comm.Abort(EX_OPEN_FILE);
		}
	}
	PM_Stop(tm_Init_LoadSTL);

	if(myrank == 0) {
		PrintLog(1, "Building BCMOctree");

		PM_Start(tm_Init_DivideDomain, 0, 0, false);
		rootGrid = new RootGrid(conf.rootN);
		if( conf.periodicX ) {
			rootGrid->setPeriodicX();
		}
		if( conf.periodicY ) {
			rootGrid->setPeriodicY();
		}
		if( conf.periodicZ ) {
			rootGrid->setPeriodicZ();
		}

		Divider* divider;
		if (conf.treeType == "flat") {
			divider = new FlatDivider(
												rootGrid,
												conf.maxLevel);
		} else if (conf.treeType == "simple") {
			divider = new SimpleDivider(
												rootGrid,
												conf.minLevel,
												conf.maxLevel);
		} else if (conf.treeType == "polygon") {
			divider = new PolygonBBoxDivider(
												conf.origin,
												conf.rootLength,
												rootGrid,
												conf.minLevel,
												pl,
												conf.polygonGroupList,
												conf.boundingBoxList,
												(REAL_TYPE)conf.vc/conf.size);
		} else if (conf.treeType == "bb") {
			divider = new BBDivider(
												conf.origin,
												conf.rootLength,
												rootGrid,
												conf.minLevel,
												pl,
												conf.polygonGroupList,
												conf.boundingBoxList,
												conf.sphericalBoxList,
												(REAL_TYPE)conf.vc/conf.size);
		} else if (conf.treeType == "sphere3") {
			divider = new SphereDivider3(
												conf.origin,
												conf.rootLength,
												rootGrid,
												conf.minLevel,
												pl,
												conf.polygonGroupList,
												conf.boundingBoxList,
												conf.sphericalBoxList,
												(REAL_TYPE)conf.vc/conf.size);
		} else if(conf.treeType == "sphere") {
			divider = new SphereDivider2(
												rootGrid,
												conf.minLevel,
												conf.maxLevel,
												conf.dividerOx,
												conf.dividerOy,
												conf.dividerOz,
												conf.dividerR,
												conf.dividerDR,
												conf.dividerHollow);
		} else if(conf.treeType == "cylinder") {
		} else if(conf.treeType == "rectangle") {
			divider = new RectangleDivider(
												rootGrid,
												conf.minLevel,
												conf.maxLevel,
												conf.dividerX0,
												conf.dividerY0,
												conf.dividerZ0,
												conf.dividerX1,
												conf.dividerY1,
												conf.dividerZ1);
		} else {
			exit(EX_READ_CONFIG);
		}
		PM_Stop(tm_Init_DivideDomain);

		PM_Start(tm_Init_OrderBlock, 0, 0, false);
		BCMOctree::Ordering ordering;
		if (conf.ordering == "Z") {
			ordering = BCMOctree::Z;
		} else if (conf.ordering == "Hilbert") {
			ordering = BCMOctree::HILBERT;
		} else if (conf.ordering == "random") {
			ordering = BCMOctree::RANDOM;
		} else if (conf.ordering == "PedigreeList") {
			ordering = BCMOctree::PEDIGREELIST;
		} else {
			exit(EX_READ_CONFIG);
		}
		PM_Stop(tm_Init_OrderBlock);

		tree = new BCMOctree(rootGrid, divider, ordering);
    tree->broadcast();
	} else {
    tree = BCMOctree::ReceiveFromMaster();
	}

  int numLeafNode = tree->getNumLeafNode();
	partition = new Partition(comm.Get_size(), numLeafNode);
	for(int n=0; n<comm.Get_size(); n++) {
		if( n==0 ) {
			PrintLog(2, "%-20s : #%05d [%04d:%04d] (%04d)",  "Partitions", n, partition->getStart(n), partition->getEnd(n)-1, partition->getEnd(n) - partition->getStart(n));
		} else {
			PrintLog(2, "%-20s : #%05d [%04d:%04d] (%04d)",  "", n, partition->getStart(n), partition->getEnd(n)-1, partition->getEnd(n) - partition->getStart(n));
		}
	}

  // �֥��å����Υ�����
	::Vec3i size(conf.size, conf.size, conf.size);

	Vec3r rootOrigin = conf.origin;
	REAL_TYPE rootLength = conf.rootLength;

  std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
  for (int id = partition->getStart(myrank); id < partition->getEnd(myrank); id++) {
    Node* node = leafNodeArray[id];
		int level = node->getLevel();
		Vec3r origin = tree->getOrigin(node) * rootLength + rootOrigin;
		Vec3r blockSize = node->getBlockSize() * rootLength;
		NeighborInfo* neighborInfo = tree->makeNeighborInfo(node, partition);

		for (int i = 0; i < NUM_FACE; ++i) {
			Face face = Face(i);
			switch(i) {
				case 0:
				case 1:
					if( conf.periodicX == true ) {
						continue;
					}
					break;
				case 2:
				case 3:
					if( conf.periodicY == true ) {
						continue;
					}
					break;
				case 4:
				case 5:
					if( conf.periodicZ == true ) {
						continue;
					}
					break;
			}
			bool bOuterBoundary = tree->checkOnOuterBoundary(node, face);
			neighborInfo[face].setOuterBoundary(bOuterBoundary);
		}

		BlockBase* block = new BlockBase(size, origin, blockSize, level, neighborInfo);
    blockManager.registerBlock(block);
  }
  blockManager.endRegisterBlock();
  blockManager.printBlockLayoutInfo();
  blockManager.printBlockLayoutInfo("data-block.txt");

	if( conf.GridGenerationMode ) {
		return EX_FAILURE;
		return EX_SUCCESS;
	}

	if( myrank == 0 ) {
		if (conf.polygonGroupList.size() > 0) {
			BlockBoundingBox bbb(tree, conf.origin, conf.rootLength, size, conf.vc);
			for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
				BoundingBox box;
				for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
					Node* node = leafNodeArray[id];
					box.addBox(bbb.getBoundingBox(node));
				}
				pl->set_bounding_box(iRank, box);
			}
			pl->send_to_all();
		}
	} else {
		if (conf.polygonGroupList.size() > 0) {
			pl->load_from_rank0();
		}
	}

	cutlib::RepairPolygonData(pl);

/*
	std::string file;
	pl->save_parallel(&file, "stl", "test");
*/

//  delete tree;
//  delete partition;
/////////////////////////////////////////////

/////////////////////////////////////////////
// Misc.
/////////////////////////////////////////////
	maxLevel = 0;
	minLevel = INT_MAX;
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		int level = block->getLevel();
		if( level >= maxLevel ) {
			maxLevel = level;
		}
		if( level <= minLevel ) {
			minLevel = level;
		}
	}
	diffLevel = maxLevel - conf.minLevel;
	vc = conf.vc;
/////////////////////////////////////////////



/////////////////////////////////////////////
// Init physical parameters
/////////////////////////////////////////////
	dt		= conf.dt;

	rhof	= conf.rhof;
	rhos	= conf.rhos;
	cpf		= conf.cpf;
	cps		= conf.cps;
	kf		= conf.kf;
	ks		= conf.ks;
	mu		= conf.mu;
/////////////////////////////////////////////


/////////////////////////////////////////////
// Init ILS
/////////////////////////////////////////////
	omegaU		= conf.omegaU;
	countMaxU	= conf.countMaxU;
	epsilonU	= conf.epsilonU;
	countPreConditionerU
						= conf.countPreConditionerU;
	countUX		= 0;
	residualUX= 0.0;
	countUY		= 0;
	residualUY= 0.0;
	countUZ		= 0;
	residualUZ= 0.0;

	omegaP		= conf.omegaP;
	countMaxP	= conf.countMaxP;
	epsilonP	= conf.epsilonP;
	countPreConditionerP
						= conf.countPreConditionerP;
	countP		= 0;
	residualP	= 0.0;

	omegaT		= conf.omegaT;
	countMaxT	= conf.countMaxT;
	epsilonT	= conf.epsilonT;
	countPreConditionerT
						= conf.countPreConditionerT;
	countT		= 0;
	residualT	= 0.0;

	pils = new ILS3D();

	int boundaryTypeNULL[NUM_FACE] = {
		1, 1, 1, 1, 1, 1,
	};

	real boundaryValueNULL[NUM_FACE] = {
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	};

	int boundaryValueNULLINT[NUM_FACE] = {
		0, 0, 0, 0, 0, 0,
	};

	plsr  = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsr0 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsp  = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsp_ = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsq_ = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plss  = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plss_ = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plst_ = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsr ->Fill(blockManager, 0.0);
	plsr0->Fill(blockManager, 0.0);
	plsp ->Fill(blockManager, 0.0);
	plsp_->Fill(blockManager, 0.0);
	plsq_->Fill(blockManager, 0.0);
	plss ->Fill(blockManager, 0.0);
	plss_->Fill(blockManager, 0.0);
	plst_->Fill(blockManager, 0.0);
/////////////////////////////////////////////


/////////////////////////////////////////////
// Init Force
/////////////////////////////////////////////
	plsFspx = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFspy = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFspz = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFsvx = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFsvy = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFsvz = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsFspx->Fill(blockManager, 0.0);
	plsFspy->Fill(blockManager, 0.0);
	plsFspz->Fill(blockManager, 0.0);
	plsFsvx->Fill(blockManager, 0.0);
	plsFsvy->Fill(blockManager, 0.0);
	plsFsvz->Fill(blockManager, 0.0);
/////////////////////////////////////////////


/////////////////////////////////////////////
// Init Cut
/////////////////////////////////////////////
	PrintLog(1, "Computing cuts");

	plsCut0 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut1 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut2 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut3 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut4 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCut5 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsCutId0 = new LocalScalar3D<int>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId1 = new LocalScalar3D<int>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId2 = new LocalScalar3D<int>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId3 = new LocalScalar3D<int>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId4 = new LocalScalar3D<int>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCutId5 = new LocalScalar3D<int>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsCut0->Fill(blockManager, 1.0);
	plsCut1->Fill(blockManager, 1.0);
	plsCut2->Fill(blockManager, 1.0);
	plsCut3->Fill(blockManager, 1.0);
	plsCut4->Fill(blockManager, 1.0);
	plsCut5->Fill(blockManager, 1.0);
	plsCutId0->Fill(blockManager, 0);
	plsCutId1->Fill(blockManager, 0);
	plsCutId2->Fill(blockManager, 0);
	plsCutId3->Fill(blockManager, 0);
	plsCutId4->Fill(blockManager, 0);
	plsCutId5->Fill(blockManager, 0);

	plsPhaseId = new LocalScalar3D<int>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULLINT, 2);
	plsPhaseId->Fill(blockManager, -1);


#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size      = block->getSize();
		::Vec3r origin    = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize  = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};

		REAL_TYPE bpos[3] = {origin.x, origin.y, origin.z};
		unsigned int bbsize[3] = {size.x, size.y, size.z};
		unsigned int gcsize[3] = {vc, vc, vc};
		REAL_TYPE dx[3] = {cellSize.x, cellSize.x, cellSize.x};
		size_t ncell[3];
		REAL_TYPE org[3];
		for(int i=0; i<3; i++) {
			ncell[i] = bbsize[i] + 2*gcsize[i];
			org[i] = bpos[i] - gcsize[i]*dx[i];
		}

		cutlib::GridAccessor* grid   = new cutlib::Cell(org, dx);
		cutlib::CutPosArray*  cutPos = new cutlib::CutPos32Array(ncell);
		cutlib::CutBidArray*  cutBid = new cutlib::CutBid5Array(ncell);

//		CutInfoCell(org, dx, pl, cutPos, cutBid);
		CalcCutInfo(grid, pl, cutPos, cutBid);

//		outputVtk("test_2_", n + blockManager.getNumBlock()*myrank, grid, cutPos, cutBid);

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);
#ifdef _BLOCK_IS_LARGE_
#pragma omp parallel for
#else
#endif
		for(int k=vc; k<vc+size.z; k++) {
			for(int j=vc; j<vc+size.y; j++) {
				for(int i=vc; i<vc+size.x; i++) {
					int m = i + (2*vc + size.x)*(j + (2*vc + size.y)*k);
					float cut0 = cutPos->getPos(i, j, k, 0);
					float cut1 = cutPos->getPos(i, j, k, 1);
					float cut2 = cutPos->getPos(i, j, k, 2);
					float cut3 = cutPos->getPos(i, j, k, 3);
					float cut4 = cutPos->getPos(i, j, k, 4);
					float cut5 = cutPos->getPos(i, j, k, 5);
					cutlib::BidType bid0 = cutBid->getBid(i, j, k, 0);
					cutlib::BidType bid1 = cutBid->getBid(i, j, k, 1);
					cutlib::BidType bid2 = cutBid->getBid(i, j, k, 2);
					cutlib::BidType bid3 = cutBid->getBid(i, j, k, 3);
					cutlib::BidType bid4 = cutBid->getBid(i, j, k, 4);
					cutlib::BidType bid5 = cutBid->getBid(i, j, k, 5);
					if( (int)bid0 != 0 ) {
						pCut0[m] = cut0;
					}
					if( (int)bid1 != 0 ) {
						pCut1[m] = cut1;
					}
					if( (int)bid2 != 0 ) {
						pCut2[m] = cut2;
					}
					if( (int)bid3 != 0 ) {
						pCut3[m] = cut3;
					}
					if( (int)bid4 != 0 ) {
						pCut4[m] = cut4;
					}
					if( (int)bid5 != 0 ) {
						pCut5[m] = cut5;
					}
					pCutId0[m] = bid0;
					pCutId1[m] = bid1;
					pCutId2[m] = bid2;
					pCutId3[m] = bid3;
					pCutId4[m] = bid4;
					pCutId5[m] = bid5;
				}
			}
		}
		delete grid;
		delete cutPos;
		delete cutBid;

		if( conf.cutoff ) {
			real eps[1] = {conf.cutoff_epsilon};
			bstl_cutoff_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							eps,
							sz, g);
		}

		if( conf.voxelization ) {
			bstl_voxelize_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							sz, g);
		}

		if( conf.symmetrization ) {
			bstl_symmetrize_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							sz, g);
		}
	}

	plsCut0->ImposeBoundaryCondition(blockManager);
	plsCut1->ImposeBoundaryCondition(blockManager);
	plsCut2->ImposeBoundaryCondition(blockManager);
	plsCut3->ImposeBoundaryCondition(blockManager);
	plsCut4->ImposeBoundaryCondition(blockManager);
	plsCut5->ImposeBoundaryCondition(blockManager);
	plsCutId0->ImposeBoundaryCondition(blockManager);
	plsCutId1->ImposeBoundaryCondition(blockManager);
	plsCutId2->ImposeBoundaryCondition(blockManager);
	plsCutId3->ImposeBoundaryCondition(blockManager);
	plsCutId4->ImposeBoundaryCondition(blockManager);
	plsCutId5->ImposeBoundaryCondition(blockManager);

	{
		int countLocal = 0;
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+: countLocal)
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i size = block->getSize();
			::Vec3r origin = block->getOrigin();
			::Vec3r blockSize = block->getBlockSize();
			::Vec3r cellSize = block->getCellSize();

			int sz[3] = {size.x, size.y, size.z};
			int g[1] = {vc};

			real* pCut0 = plsCut0->GetBlockData(block);
			real* pCut1 = plsCut1->GetBlockData(block);
			real* pCut2 = plsCut2->GetBlockData(block);
			real* pCut3 = plsCut3->GetBlockData(block);
			real* pCut4 = plsCut4->GetBlockData(block);
			real* pCut5 = plsCut5->GetBlockData(block);
			int* pCutId0 = plsCutId0->GetBlockData(block);
			int* pCutId1 = plsCutId1->GetBlockData(block);
			int* pCutId2 = plsCutId2->GetBlockData(block);
			int* pCutId3 = plsCutId3->GetBlockData(block);
			int* pCutId4 = plsCutId4->GetBlockData(block);
			int* pCutId5 = plsCutId5->GetBlockData(block);

			int count = 0;
			bstl_detect_zerocut_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							&count,
							sz, g);

			countLocal += count;
		}
		plsCut0->ImposeBoundaryCondition(blockManager);
		plsCut1->ImposeBoundaryCondition(blockManager);
		plsCut2->ImposeBoundaryCondition(blockManager);
		plsCut3->ImposeBoundaryCondition(blockManager);
		plsCut4->ImposeBoundaryCondition(blockManager);
		plsCut5->ImposeBoundaryCondition(blockManager);
		plsCutId0->ImposeBoundaryCondition(blockManager);
		plsCutId1->ImposeBoundaryCondition(blockManager);
		plsCutId2->ImposeBoundaryCondition(blockManager);
		plsCutId3->ImposeBoundaryCondition(blockManager);
		plsCutId4->ImposeBoundaryCondition(blockManager);
		plsCutId5->ImposeBoundaryCondition(blockManager);

		int countTmp = countLocal;
		MPI_Allreduce(&countTmp, &countLocal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		PrintLog(2, "%-20s : %d", "Zero distance", countLocal);
	}

	{
		int countLocal = 0;
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+: countLocal)
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i size = block->getSize();
			::Vec3r origin = block->getOrigin();
			::Vec3r blockSize = block->getBlockSize();
			::Vec3r cellSize = block->getCellSize();

			int sz[3] = {size.x, size.y, size.z};
			int g[1] = {vc};

			real* pCut0 = plsCut0->GetBlockData(block);
			real* pCut1 = plsCut1->GetBlockData(block);
			real* pCut2 = plsCut2->GetBlockData(block);
			real* pCut3 = plsCut3->GetBlockData(block);
			real* pCut4 = plsCut4->GetBlockData(block);
			real* pCut5 = plsCut5->GetBlockData(block);
			int* pCutId0 = plsCutId0->GetBlockData(block);
			int* pCutId1 = plsCutId1->GetBlockData(block);
			int* pCutId2 = plsCutId2->GetBlockData(block);
			int* pCutId3 = plsCutId3->GetBlockData(block);
			int* pCutId4 = plsCutId4->GetBlockData(block);
			int* pCutId5 = plsCutId5->GetBlockData(block);

			int count = 0;
			int bClose = 1;
			if( conf.holefilling == false ) {
				bClose = 0;
			}
			bstl_fill_holes_(
							pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
							pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
							&count,
							&bClose,
							sz, g);

			countLocal += count;
		}
		plsCut0->ImposeBoundaryCondition(blockManager);
		plsCut1->ImposeBoundaryCondition(blockManager);
		plsCut2->ImposeBoundaryCondition(blockManager);
		plsCut3->ImposeBoundaryCondition(blockManager);
		plsCut4->ImposeBoundaryCondition(blockManager);
		plsCut5->ImposeBoundaryCondition(blockManager);
		plsCutId0->ImposeBoundaryCondition(blockManager);
		plsCutId1->ImposeBoundaryCondition(blockManager);
		plsCutId2->ImposeBoundaryCondition(blockManager);
		plsCutId3->ImposeBoundaryCondition(blockManager);
		plsCutId4->ImposeBoundaryCondition(blockManager);
		plsCutId5->ImposeBoundaryCondition(blockManager);

		int countTmp = countLocal;
		MPI_Allreduce(&countTmp, &countLocal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		PrintLog(2, "%-20s : %d", "Hole faces(single)", countLocal);
	}

	{
		int nIterationCount = 0;
		int countLocal = 0;
		int countTotal = 0;
		do {
			countLocal = 0;

#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+: countLocal)
#endif
			for (int n=0; n<blockManager.getNumBlock(); ++n) {
				BlockBase* block = blockManager.getBlock(n);
				::Vec3i size = block->getSize();
				::Vec3r origin = block->getOrigin();
				::Vec3r blockSize = block->getBlockSize();
				::Vec3r cellSize = block->getCellSize();

				int sz[3] = {size.x, size.y, size.z};
				int g[1] = {vc};

				real* pCut0 = plsCut0->GetBlockData(block);
				real* pCut1 = plsCut1->GetBlockData(block);
				real* pCut2 = plsCut2->GetBlockData(block);
				real* pCut3 = plsCut3->GetBlockData(block);
				real* pCut4 = plsCut4->GetBlockData(block);
				real* pCut5 = plsCut5->GetBlockData(block);
				int* pCutId0 = plsCutId0->GetBlockData(block);
				int* pCutId1 = plsCutId1->GetBlockData(block);
				int* pCutId2 = plsCutId2->GetBlockData(block);
				int* pCutId3 = plsCutId3->GetBlockData(block);
				int* pCutId4 = plsCutId4->GetBlockData(block);
				int* pCutId5 = plsCutId5->GetBlockData(block);

				int count = 0;
				int bClose = 1;
				if( conf.holefilling2 == false ) {
					bClose = 0;
				}

				bstl_fill_holes_v2_(
								pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
								pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
								&count,
								&bClose,
								sz, g);

				countLocal += count;
			}
			plsCut0->ImposeBoundaryCondition(blockManager);
			plsCut1->ImposeBoundaryCondition(blockManager);
			plsCut2->ImposeBoundaryCondition(blockManager);
			plsCut3->ImposeBoundaryCondition(blockManager);
			plsCut4->ImposeBoundaryCondition(blockManager);
			plsCut5->ImposeBoundaryCondition(blockManager);
			plsCutId0->ImposeBoundaryCondition(blockManager);
			plsCutId1->ImposeBoundaryCondition(blockManager);
			plsCutId2->ImposeBoundaryCondition(blockManager);
			plsCutId3->ImposeBoundaryCondition(blockManager);
			plsCutId4->ImposeBoundaryCondition(blockManager);
			plsCutId5->ImposeBoundaryCondition(blockManager);

			int countTmp = countLocal;
			MPI_Allreduce(&countTmp, &countLocal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

			nIterationCount++;
			countTotal += countLocal;
		}while( countLocal > 0 && conf.holefilling2 == true);

		PrintLog(2, "%-20s : %d (#Ite.: %d)", "Hole faces(multi)", countTotal, nIterationCount);
	}

	PrintLog(1, "Filling");

	real xs = conf.seed.x;
	real ys = conf.seed.y;
	real zs = conf.seed.z;
	PrintLog(2, "%-20s : %f %f %f", "Seed for FLUID", xs, ys, zs);

#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;
		real org[3] = {origin.x, origin.y, origin.z};

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		bcut_set_fluidseed_(
			pPhaseId,
			&xs, &ys, &zs,
			&dx,
			org,
			sz, g);

	}
	plsPhaseId->ImposeBoundaryCondition(blockManager);

	{
		int nIterationCount = 0;
		long int nCellsChanged = 0;
		do {
			nCellsChanged = 0;
#ifdef _BLOCK_IS_LARGE_
#else
//#pragma omp parallel for reduction(+: nCellsChanged)
#endif
			for (int n=0; n<blockManager.getNumBlock(); ++n) {
				BlockBase* block = blockManager.getBlock(n);
				::Vec3i size = block->getSize();
				::Vec3r origin = block->getOrigin();
				::Vec3r blockSize = block->getBlockSize();
				::Vec3r cellSize = block->getCellSize();

				int sz[3] = {size.x, size.y, size.z};
				int g[1] = {vc};
				int nc[3] = {size.x + 2*vc, size.y + 2*vc, size.z + 2*vc};

				real* pCut0 = plsCut0->GetBlockData(block);
				real* pCut1 = plsCut1->GetBlockData(block);
				real* pCut2 = plsCut2->GetBlockData(block);
				real* pCut3 = plsCut3->GetBlockData(block);
				real* pCut4 = plsCut4->GetBlockData(block);
				real* pCut5 = plsCut5->GetBlockData(block);
				int* pCutId0 = plsCutId0->GetBlockData(block);
				int* pCutId1 = plsCutId1->GetBlockData(block);
				int* pCutId2 = plsCutId2->GetBlockData(block);
				int* pCutId3 = plsCutId3->GetBlockData(block);
				int* pCutId4 = plsCutId4->GetBlockData(block);
				int* pCutId5 = plsCutId5->GetBlockData(block);

				int* pPhaseId = plsPhaseId->GetBlockData(block);
#ifdef _BLOCK_IS_LARGE_
//#pragma omp parallel for reduction(+: nCellsChanged)
#else
#endif
				for(int k=vc; k<=size.z+vc-1; k++) {
					for(int j=vc; j<=size.y+vc-1; j++) {
						for(int i=vc; i<=size.x+vc-1; i++) {
							int mp = i + nc[0]*( j + nc[1]*k );
							int mw = i-1 + nc[0]*( j + nc[1]*k );
							int me = i+1 + nc[0]*( j + nc[1]*k );
							int ms = i + nc[0]*( j-1 + nc[1]*k );
							int mn = i + nc[0]*( j+1 + nc[1]*k );
							int mb = i + nc[0]*( j + nc[1]*(k-1) );
							int mt = i + nc[0]*( j + nc[1]*(k+1) );

							int cidp0 = pCutId0[mp];
							int cidp1 = pCutId1[mp];
							int cidp2 = pCutId2[mp];
							int cidp3 = pCutId3[mp];
							int cidp4 = pCutId4[mp];
							int cidp5 = pCutId5[mp];

							if( pPhaseId[mp] > 0 ) {
								continue;
							}

							if( (pPhaseId[mw] == 1 && cidp0 == 0) ||
									(pPhaseId[me] == 1 && cidp1 == 0) ||
									(pPhaseId[ms] == 1 && cidp2 == 0) ||
									(pPhaseId[mn] == 1 && cidp3 == 0) ||
									(pPhaseId[mb] == 1 && cidp4 == 0) ||
									(pPhaseId[mt] == 1 && cidp5 == 0) ) {
								pPhaseId[mp] = 1;
								nCellsChanged++;
							}
						}
					}
				}
			}
			plsPhaseId->ImposeBoundaryCondition(blockManager);

			long int nCellsChangedTmp = nCellsChanged;
			MPI_Allreduce(&nCellsChangedTmp, &nCellsChanged, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

			nIterationCount++;
		}while(nCellsChanged>0);

		long int count = 0;
		long int countS = 0;
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i size = block->getSize();
			::Vec3r origin = block->getOrigin();
			::Vec3r blockSize = block->getBlockSize();
			::Vec3r cellSize = block->getCellSize();

			int sz[3] = {size.x, size.y, size.z};
			int g[1] = {vc};
			int nc[3] = {size.x + 2*vc, size.y + 2*vc, size.z + 2*vc};

			int* pPhaseId = plsPhaseId->GetBlockData(block);
			for(int k=vc; k<=size.z+vc-1; k++) {
				for(int j=vc; j<=size.y+vc-1; j++) {
					for(int i=vc; i<=size.x+vc-1; i++) {
						int mp = i + nc[0]*( j + nc[1]*k );
						if( pPhaseId[mp] > 0 ) {
							count++;
						} else {
							countS++;
						}
					}
				}
			}
		}

		long int countTmp = count;
		MPI_Allreduce(&countTmp, &count, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

		countTmp = countS;
		MPI_Allreduce(&countTmp, &countS, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

		PrintLog(2, "%-20s : %d", "Iteration", nIterationCount);
		PrintLog(2, "%-20s : %d", "FLUID cells", count);
		PrintLog(2, "%-20s : %d", "SOLID cells", countS);
	}
/////////////////////////////////////////////


/////////////////////////////////////////////
// Init mask
/////////////////////////////////////////////
	plsMaskId = new LocalScalar3D<int>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULLINT);
	plsMaskId->Fill(blockManager, 0);

	plsM = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsM->Fill(blockManager, 0.0);

#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};

		real* pM = plsM->GetBlockData(block);
		int* pMaskId = plsMaskId->GetBlockData(block);

		setup_mask_(
			pM,
			pMaskId,
			sz, g);
	}
/////////////////////////////////////////////



/////////////////////////////////////////////
// Init Variables
/////////////////////////////////////////////
	int boundaryTypeUX[NUM_FACE] = {
		conf.boundaryTypeUX_X_M,
		conf.boundaryTypeUX_X_P,
		conf.boundaryTypeUX_Y_M,
		conf.boundaryTypeUX_Y_P,
		conf.boundaryTypeUX_Z_M,
		conf.boundaryTypeUX_Z_P,
	};
	real boundaryValueUX[NUM_FACE] = {
		conf.boundaryValueUX_X_M,
		conf.boundaryValueUX_X_P,
		conf.boundaryValueUX_Y_M,
		conf.boundaryValueUX_Y_P,
		conf.boundaryValueUX_Z_M,
		conf.boundaryValueUX_Z_P,
	};
	plsUX0 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeUX, boundaryValueUX, 1);
	plsUX1 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeUX, boundaryValueUX, 1);
	plsUXC = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUXCP= new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUXD = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUX0->Fill(blockManager, 0.0, conf.uxs0);
	plsUX1->Fill(blockManager, 0.0, conf.uxs0);
	plsUXC->Fill(blockManager, 0.0);
	plsUXCP->Fill(blockManager, 0.0);
	plsUXD->Fill(blockManager, 0.0);
	plsUX0->ImposeBoundaryCondition(blockManager);
	plsUX1->ImposeBoundaryCondition(blockManager);
	plsUXC->ImposeBoundaryCondition(blockManager);
	plsUXCP->ImposeBoundaryCondition(blockManager);
	plsUXD->ImposeBoundaryCondition(blockManager);

	int boundaryTypeUY[NUM_FACE] = {
		conf.boundaryTypeUY_X_M,
		conf.boundaryTypeUY_X_P,
		conf.boundaryTypeUY_Y_M,
		conf.boundaryTypeUY_Y_P,
		conf.boundaryTypeUY_Z_M,
		conf.boundaryTypeUY_Z_P,
	};
	real boundaryValueUY[NUM_FACE] = {
		conf.boundaryValueUY_X_M,
		conf.boundaryValueUY_X_P,
		conf.boundaryValueUY_Y_M,
		conf.boundaryValueUY_Y_P,
		conf.boundaryValueUY_Z_M,
		conf.boundaryValueUY_Z_P,
	};
	plsUY0 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeUY, boundaryValueUY, 1);
	plsUY1 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeUY, boundaryValueUY, 1);
	plsUYC = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUYCP= new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUYD = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUY0->Fill(blockManager, 0.0, conf.uys0);
	plsUY1->Fill(blockManager, 0.0, conf.uys0);
	plsUYC->Fill(blockManager, 0.0);
	plsUYCP->Fill(blockManager, 0.0);
	plsUYD->Fill(blockManager, 0.0);
	plsUY0->ImposeBoundaryCondition(blockManager);
	plsUY1->ImposeBoundaryCondition(blockManager);
	plsUYC->ImposeBoundaryCondition(blockManager);
	plsUYCP->ImposeBoundaryCondition(blockManager);
	plsUYD->ImposeBoundaryCondition(blockManager);

	int boundaryTypeUZ[NUM_FACE] = {
		conf.boundaryTypeUZ_X_M,
		conf.boundaryTypeUZ_X_P,
		conf.boundaryTypeUZ_Y_M,
		conf.boundaryTypeUZ_Y_P,
		conf.boundaryTypeUZ_Z_M,
		conf.boundaryTypeUZ_Z_P,
	};
	real boundaryValueUZ[NUM_FACE] = {
		conf.boundaryValueUZ_X_M,
		conf.boundaryValueUZ_X_P,
		conf.boundaryValueUZ_Y_M,
		conf.boundaryValueUZ_Y_P,
		conf.boundaryValueUZ_Z_M,
		conf.boundaryValueUZ_Z_P,
	};
	plsUZ0 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeUZ, boundaryValueUZ, 1);
	plsUZ1 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeUZ, boundaryValueUZ, 1);
	plsUZC = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUZCP= new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUZD = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsUZ0->Fill(blockManager, 0.0, conf.uzs0);
	plsUZ1->Fill(blockManager, 0.0, conf.uzs0);
	plsUZC->Fill(blockManager, 0.0);
	plsUZCP->Fill(blockManager, 0.0);
	plsUZD->Fill(blockManager, 0.0);
	plsUZ0->ImposeBoundaryCondition(blockManager);
	plsUZ1->ImposeBoundaryCondition(blockManager);
	plsUZC->ImposeBoundaryCondition(blockManager);
	plsUZCP->ImposeBoundaryCondition(blockManager);
	plsUZD->ImposeBoundaryCondition(blockManager);

	plsVw = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVe = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVs = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVn = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVb = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVt = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsVw->Fill(blockManager, 0.0);
	plsVe->Fill(blockManager, 0.0);
	plsVs->Fill(blockManager, 0.0);
	plsVn->Fill(blockManager, 0.0);
	plsVb->Fill(blockManager, 0.0);
	plsVt->Fill(blockManager, 0.0);
	plsVw->ImposeBoundaryCondition(blockManager);
	plsVe->ImposeBoundaryCondition(blockManager);
	plsVs->ImposeBoundaryCondition(blockManager);
	plsVn->ImposeBoundaryCondition(blockManager);
	plsVb->ImposeBoundaryCondition(blockManager);
	plsVt->ImposeBoundaryCondition(blockManager);

	int boundaryTypeT[NUM_FACE] = {
		conf.boundaryTypeT_X_M,
		conf.boundaryTypeT_X_P,
		conf.boundaryTypeT_Y_M,
		conf.boundaryTypeT_Y_P,
		conf.boundaryTypeT_Z_M,
		conf.boundaryTypeT_Z_P,
	};
	real boundaryValueT[NUM_FACE] = {
		conf.boundaryValueT_X_M,
		conf.boundaryValueT_X_P,
		conf.boundaryValueT_Y_M,
		conf.boundaryValueT_Y_P,
		conf.boundaryValueT_Z_M,
		conf.boundaryValueT_Z_P,
	};
	plsT0 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeT, boundaryValueT, 1);
	plsT1 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeT, boundaryValueT, 1);
	plsTC = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsTCP= new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsTD = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsT0->Fill(blockManager, 0.0);
	plsT1->Fill(blockManager, 0.0);
	plsTC->Fill(blockManager, 0.0);
	plsTCP->Fill(blockManager, 0.0);
	plsTD->Fill(blockManager, 0.0);
	plsT0->ImposeBoundaryCondition(blockManager);
	plsT1->ImposeBoundaryCondition(blockManager);
	plsTC->ImposeBoundaryCondition(blockManager);
	plsTCP->ImposeBoundaryCondition(blockManager);
	plsTD->ImposeBoundaryCondition(blockManager);

	int boundaryTypeP[NUM_FACE] = {
		conf.boundaryTypeP_X_M,
		conf.boundaryTypeP_X_P,
		conf.boundaryTypeP_Y_M,
		conf.boundaryTypeP_Y_P,
		conf.boundaryTypeP_Z_M,
		conf.boundaryTypeP_Z_P,
	};
	real boundaryValueP[NUM_FACE] = {
		conf.boundaryValueP_X_M,
		conf.boundaryValueP_X_P,
		conf.boundaryValueP_Y_M,
		conf.boundaryValueP_Y_P,
		conf.boundaryValueP_Z_M,
		conf.boundaryValueP_Z_P,
	};
	plsP0 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeP, boundaryValueP, 1);
	plsP1 = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeP, boundaryValueP, 1);
	plsLapP = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsP0->Fill(blockManager, 0.0);
	plsP1->Fill(blockManager, 0.0);
	plsLapP->Fill(blockManager, 0.0);
	plsP0->ImposeBoundaryCondition(blockManager);
	plsP1->ImposeBoundaryCondition(blockManager);
	plsLapP->ImposeBoundaryCondition(blockManager);
/////////////////////////////////////////////




/////////////////////////////////////////////
// Init A & b
/////////////////////////////////////////////
	plsAp = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAw = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAe = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAs = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAn = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAb = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAt = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsb  = new LocalScalar3D<real>(blockManager, vc, conf.updateMethod, boundaryTypeNULL, boundaryValueNULL);
	plsAp->Fill(blockManager, 0.0);
	plsAw->Fill(blockManager, 0.0);
	plsAe->Fill(blockManager, 0.0);
	plsAs->Fill(blockManager, 0.0);
	plsAn->Fill(blockManager, 0.0);
	plsAb->Fill(blockManager, 0.0);
	plsAt->Fill(blockManager, 0.0);
	plsb ->Fill(blockManager, 0.0);
	plsAp->ImposeBoundaryCondition(blockManager);
	plsAw->ImposeBoundaryCondition(blockManager);
	plsAe->ImposeBoundaryCondition(blockManager);
	plsAs->ImposeBoundaryCondition(blockManager);
	plsAn->ImposeBoundaryCondition(blockManager);
	plsAb->ImposeBoundaryCondition(blockManager);
	plsAt->ImposeBoundaryCondition(blockManager);
	plsb ->ImposeBoundaryCondition(blockManager);
/////////////////////////////////////////////


/////////////////////////////////////////////
	if( conf.BCMFileSave ) {
		BCMFileSaverInit(
					rootGrid,
					tree,
					partition,
					conf);
	}
	if( conf.PLOT3DSave ) {
		WriteXYZInPlot3DFormat(
					"xyz",
					diffLevel,
					rootGrid,
					tree,
					partition,
					conf);
	}
/////////////////////////////////////////////


	PrintLog(1, "Printing STL files for cut info");
	PrintCut();

	MPI_Barrier(MPI_COMM_WORLD);
	PrintLog(1, "Completed");

	return EX_SUCCESS;
}

int Solver::Loop() {
	int StepStart = conf.StepStart;
	int StepEnd   = conf.StepEnd;

	bRestart = false;
	if(StepStart == 0) {
		Print(0);
	} else {
		Load2(StepStart);
		bRestart = true;
	}

	for(int step=StepStart+1; step<=StepEnd; step++) {
		int nResult = Update(step);

		Print(step);

		if( step%conf.StepPrintBin == 0 ) {
			Dump2(step);
		}

		switch(nResult) {
			case EX_SUCCESS:
				break;
			default:
				break;
		}
	}

	return EX_SUCCESS;
}

int Solver::Print(int step) {
	if( step%conf.StepPrintTime == 0 ) {
		PrintTime(step);
	}

	if( step%conf.StepPrintStats == 0 ) {
		PrintILS(step);
		PrintStats(step);
		PrintForce(step);
	}

	if( step%conf.StepPrintData == 0 ) {
		PrintData(step);
		if( conf.BCMFileSave ) {
			BCMFileSaverPrint(step);
		}
	}

	return EX_SUCCESS;
}

void Solver::PrintTime(int step) {
	if( myrank != 0 ) {
		return;
	}

	ostringstream ossMessage;
	ossMessage.width(10);
	ossMessage.setf(std::ios::fixed);
	ossMessage.fill('0');
	ossMessage << step << " ";
	ossMessage.setf(std::ios::scientific, std::ios::floatfield);
	ossMessage.precision(5);
	ossMessage << times[0];
	PrintLog(0, ossMessage.str().c_str());

	std::string filename = "data-times.txt";

	std::ofstream ofs;
	if( step==0 ) {
		ofs.open(filename.c_str(), std::ios::out);
		ofs.close();
	}
	ofs.open(filename.c_str(), std::ios::out | std::ios::app);

	ofs.width(10);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << step << " ";

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(5);
	ofs << this->times[0] << " ";
	ofs << this->times[1] << " ";
	ofs << this->times[2] << " ";
	ofs << this->times[3] << " ";
	ofs << this->times[4] << " ";
	ofs << this->times[5] << " ";
	ofs << this->times[6] << " ";
	ofs << std::endl;
	ofs.close();
}

void Solver::PrintILS(int step) {
	if( myrank != 0 ) {
		return;
	}

	std::string filename = "data-ils.txt";

	std::ofstream ofs;
	if( step==0 ) {
		ofs.open(filename.c_str(), std::ios::out);
		ofs.close();
	}
	ofs.open(filename.c_str(), std::ios::out | std::ios::app);

	ofs.width(10);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << step << " ";

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(5);
	ofs << this->countUX << " ";
	ofs << this->residualUX << " ";
	ofs << this->countUY << " ";
	ofs << this->residualUY << " ";
	ofs << this->countUZ << " ";
	ofs << this->residualUZ << " ";
	ofs << this->countP << " ";
	ofs << this->residualP << " ";
	ofs << this->countT << " ";
	ofs << this->residualT << " ";
	ofs << std::endl;
	ofs.close();
}

void Solver::PrintStats(int step) {
	plsUX0->CalcStats(blockManager);
	plsUY0->CalcStats(blockManager);
	plsUZ0->CalcStats(blockManager);
	plsP0->CalcStats(blockManager);
	plsT0->CalcStats(blockManager);

	if( plsUX0->GetAbsMaxL() > 100.0 ) {
		std::cout << myrank << std::endl;
	}

	if( myrank != 0 ) {
		return;
	}

	std::string filename = "data-stats.txt";

	std::ofstream ofs;
	if( step==0 ) {
		ofs.open(filename.c_str(), std::ios::out);
		ofs.close();
	}
	ofs.open(filename.c_str(), std::ios::out | std::ios::app);

	ofs.width(10);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << step << " ";

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(5);
	ofs << plsUX0->GetSum() << " ";
	ofs << plsUY0->GetSum() << " ";
	ofs << plsUZ0->GetSum() << " ";
	ofs << plsP0->GetSum() << " ";
	ofs << plsT0->GetSum() << " ";
	ofs << plsUX0->GetMax() << " ";
	ofs << plsUY0->GetMax() << " ";
	ofs << plsUZ0->GetMax() << " ";
	ofs << plsP0->GetMax() << " ";
	ofs << plsT0->GetMax() << " ";
	ofs << plsUX0->GetMin() << " ";
	ofs << plsUY0->GetMin() << " ";
	ofs << plsUZ0->GetMin() << " ";
	ofs << plsP0->GetMin() << " ";
	ofs << plsT0->GetMin() << " ";
	ofs << plsUX0->GetAbsMax() << " ";
	ofs << plsUY0->GetAbsMax() << " ";
	ofs << plsUZ0->GetAbsMax() << " ";
	ofs << plsP0->GetAbsMax() << " ";
	ofs << plsT0->GetAbsMax() << " ";
	ofs << plsUX0->GetAbsMin() << " ";
	ofs << plsUY0->GetAbsMin() << " ";
	ofs << plsUZ0->GetAbsMin() << " ";
	ofs << plsP0->GetAbsMin() << " ";
	ofs << plsT0->GetAbsMin() << " ";
	ofs << std::endl;
	ofs.close();
}

void Solver::PrintForce(int step) {
	real fsp_local[3] = {0.0, 0.0, 0.0};
	real fsv_local[3] = {0.0, 0.0, 0.0};
#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;

		real* ux0  = plsUX0->GetBlockData(block);
		real* uy0  = plsUY0->GetBlockData(block);
		real* uz0  = plsUZ0->GetBlockData(block);
		real* p0   = plsP0 ->GetBlockData(block);

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real* fspx = plsFspx->GetBlockData(block);
		real* fspy = plsFspy->GetBlockData(block);
		real* fspz = plsFspz->GetBlockData(block);
		real* fsvx = plsFsvx->GetBlockData(block);
		real* fsvy = plsFsvy->GetBlockData(block);
		real* fsvz = plsFsvz->GetBlockData(block);

		real Us = 0.0;

		real fsp_block[3] = {0.0, 0.0, 0.0};
		real fsv_block[3] = {0.0, 0.0, 0.0};

		bcut_calc_f_p_(
				fspx,
				fspy,
				fspz,
				fsp_block,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&dx, &dt,
				sz, g);

		bcut_calc_f_v_(
				fsvx,
				fsvy,
				fsvz,
				fsv_block,
				ux0,
				uy0,
				uz0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&mu,
				&dx, &dt,
				&Us,
				sz, g);

		fsp_local[0] += fsp_block[0];
		fsp_local[1] += fsp_block[1];
		fsp_local[2] += fsp_block[2];
		fsv_local[0] += fsv_block[0];
		fsv_local[1] += fsv_block[1];
		fsv_local[2] += fsv_block[2];
	}

	real fsp_global[3] = {0.0, 0.0, 0.0};
	real fsv_global[3] = {0.0, 0.0, 0.0};

	real sum = fsp_local[0];
	real sum_tmp = sum;
	allreduce_(&sum, &sum_tmp);
//	MPI_Allreduce(&sum_tmp, &sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	fsp_global[0] = sum;

	sum = fsp_local[1];
	sum_tmp = sum;
	allreduce_(&sum, &sum_tmp);
//	MPI_Allreduce(&sum_tmp, &sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	fsp_global[1] = sum;

	sum = fsp_local[2];
	sum_tmp = sum;
	allreduce_(&sum, &sum_tmp);
//	MPI_Allreduce(&sum_tmp, &sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	fsp_global[2] = sum;

	sum = fsv_local[0];
	sum_tmp = sum;
	allreduce_(&sum, &sum_tmp);
//	MPI_Allreduce(&sum_tmp, &sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	fsv_global[0] = sum;

	sum = fsv_local[1];
	sum_tmp = sum;
	allreduce_(&sum, &sum_tmp);
//	MPI_Allreduce(&sum_tmp, &sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	fsv_global[1] = sum;

	sum = fsv_local[2];
	sum_tmp = sum;
	allreduce_(&sum, &sum_tmp);
//	MPI_Allreduce(&sum_tmp, &sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	fsv_global[2] = sum;

	if( myrank != 0 ) {
		return;
	}

	std::string filename = "data-force.txt";

	std::ofstream ofs;
	if( step==0 ) {
		ofs.open(filename.c_str(), std::ios::out);
		ofs.close();
	}
	ofs.open(filename.c_str(), std::ios::out | std::ios::app);

	ofs.width(10);
	ofs.setf(std::ios::fixed);
	ofs.fill('0');
	ofs << step << " ";

	ofs.setf(std::ios::scientific, std::ios::floatfield);
	ofs.precision(5);
	ofs << fsp_global[0] << " ";
	ofs << fsp_global[1] << " ";
	ofs << fsp_global[2] << " ";
	ofs << fsv_global[0] << " ";
	ofs << fsv_global[1] << " ";
	ofs << fsv_global[2] << " ";
	ofs << std::endl;
	ofs.close();
}

void Solver::PrintData(int step) {
	WriteDataInVTKFormat("flow", step, diffLevel, rootGrid, tree, partition, conf);
	if( conf.PLOT3DSave ) {
		WriteDataInPlot3DFormat("flow", step, diffLevel, rootGrid, tree, partition, conf);
	}
//	WriteDataInSiloFormat("flow", step, conf);
	plsT0->WriteDataInVTKFormat("t", step, diffLevel, rootGrid, tree, partition, conf);
	if( step == 0 ) {
		plsPhaseId->WriteDataInVTKFormat("solid", 0, diffLevel, rootGrid, tree, partition, conf);
	}
}

void Solver::PrintLog(int level, const char* format, ...) {
	if( myrank != 0 ) {
		return;
	}
	static int first = 1;

	va_list ap;
	va_start(ap, format);

	char* buffer;
	int size = vasprintf(&buffer, format, ap);
	va_end(ap);

	ostringstream ossMessage;
	switch(level) {
		case 1:
			ossMessage << "# ";
			break;
		case 2:
			ossMessage << "#   ";
			break;
		case 3:
			ossMessage << "#     ";
			break;
		case 0:
		default:
			break;
	}
	ossMessage << string(buffer);


	free(buffer);

	std::cout << ossMessage.str() << std::endl;

	std::string filename = "data-log.txt";
	std::ofstream ofs;
	if( first == 1 ) {
		ofs.open(filename.c_str(), std::ios::out);
		ofs.close();
	}
	ofs.open(filename.c_str(), std::ios::out | std::ios::app);
	ofs << ossMessage.str() << std::endl;
	ofs.close();

	first = 0;
}

int Solver::Post() {
	g_pPM->gather();

	if( this->myrank == 0 ) {
		char hostname[MPI_MAX_PROCESSOR_NAME] = {0};
		int nameLen;
		MPI_Get_processor_name(hostname, &nameLen);
		nameLen++;

		FILE* fp = fopen("data-pm.txt", "w");
		g_pPM->print(fp, hostname, conf.operatorname);
		g_pPM->printDetail(fp);
		fclose(fp);
	}

	return EX_SUCCESS;
}

int Solver::Update(int step) {
PM_Start(tm_Update, 0, 0, true);
double t0 = GetTime();

	if( conf.AccelDuration != 0 && step <= conf.STEP_ACCELDURATION ) {
		real vb = conf.boundaryValueUX_X_M*(real)step/(real)conf.STEP_ACCELDURATION;
		plsUX0->ResetBoundaryConditionValue(blockManager, 0, vb);
		plsUX1->ResetBoundaryConditionValue(blockManager, 0, vb);
	}

PM_Start(tm_UpdateT, 0, 0, true);
	if( conf.bHeat ) {
		UpdateT();
	}
double t1 = GetTime();
PM_Stop(tm_UpdateT);

PM_Start(tm_UpdateUX, 0, 0, true);
	UpdateUX();
PM_Stop(tm_UpdateUX);
double t2 = GetTime();

PM_Start(tm_UpdateUY, 0, 0, true);
	UpdateUY();
PM_Stop(tm_UpdateUY);
double t3 = GetTime();

PM_Start(tm_UpdateUZ, 0, 0, true);
	UpdateUZ();
PM_Stop(tm_UpdateUZ);
double t4 = GetTime();

PM_Start(tm_UpdateP, 0, 0, true);
	UpdateP();
PM_Stop(tm_UpdateP);
double t5 = GetTime();

PM_Start(tm_UpdateU, 0, 0, true);
	UpdateU();
PM_Stop(tm_UpdateU);
double t6 = GetTime();

PM_Stop(tm_Update);

	this->times[0] = t6 - t0;
	this->times[1] = t1 - t0;
	this->times[2] = t2 - t1;
	this->times[3] = t3 - t2;
	this->times[4] = t4 - t3;
	this->times[5] = t5 - t4;
	this->times[6] = t6 - t5;

	return EX_SUCCESS;
}

void Solver::UpdateUX() {
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;

		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* ux0  = plsUX0->GetBlockData(block);
		real* uxc0 = plsUXC->GetBlockData(block);
		real* uxcp = plsUXCP->GetBlockData(block);
		real* uxd0 = plsUXD->GetBlockData(block);
		real* p0   = plsP0->GetBlockData(block);

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real Us = 0.0;

		real* uy0  = plsUY0->GetBlockData(block);
		real* uz0  = plsUZ0->GetBlockData(block);
		real* uyd0 = plsUYD->GetBlockData(block);
		real* uzd0 = plsUZD->GetBlockData(block);

		bcut_calc_d_u_(
				uxd0, uyd0, uzd0,
				ux0, uy0, uz0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&mu,
				&dx, &dt,
				&Us,
				sz, g);
		if( conf.advection == "W3" ) {
			bcut_calc_c_f_w3_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					sz, g);
		} else if( conf.advection == "E3" ) {
			bcut_calc_c_f_e3_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					sz, g);
		} else if( conf.advection == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					uxc0,
					ux0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					sz, g);
		}
		int axis=0;
		real gx = conf.gx;
		real gy = conf.gy;
		real gz = conf.gz;
		bcut_calc_ab_u_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				ux0, uy0, uz0,
				uxc0, uxcp,
				uxd0,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&axis,
				&rhof,
				&mu,
				&dx, &dt,
				&Us,
				&gx, &gy, &gz,
				sz, g);
		copy_(
				uxcp,
				uxc0,
				sz, g);
	}
	plsUX0->ImposeBoundaryCondition(blockManager, plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt, plsb);

	if( g_pconf->masking == true ) {
		pils->BiCGSTAB_Mask(
							blockManager,
							plsUX0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUX1,
							plsMaskId,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUX,
							residualUX);
/*
		pils->Jacobi_Mask(
							blockManager,
							plsUX0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsUX1,
							plsMaskId,
							omegaU,
							countMaxU,
							epsilonU,
							countUX,
							residualUX);
*/
	} else {
		pils->BiCGSTAB(
							blockManager,
							plsUX0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUX1,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUX,
							residualUX);
	}
	plsUX0->ImposeBoundaryCondition(blockManager);
/////////////////////////////////////////////
}

void Solver::UpdateUY() {
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;

		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* uy0  = plsUY0->GetBlockData(block);
		real* uyc0 = plsUYC->GetBlockData(block);
		real* uycp = plsUYCP->GetBlockData(block);
		real* uyd0 = plsUYD->GetBlockData(block);
		real* p0   = plsP0->GetBlockData(block);

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real* ux0  = plsUX0->GetBlockData(block);
		real* uz0  = plsUZ0->GetBlockData(block);

		real Us = 0.0;

		if( conf.advection == "W3" ) {
			bcut_calc_c_f_w3_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					sz, g);
		} else if( conf.advection == "E3" ) {
			bcut_calc_c_f_e3_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					sz, g);
		} else if( conf.advection == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					uyc0,
					uy0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					sz, g);
		}
		int axis=1;
		real gx = conf.gx;
		real gy = conf.gy;
		real gz = conf.gz;
		bcut_calc_ab_u_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				ux0, uy0, uz0,
				uyc0, uycp,
				uyd0,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&axis,
				&rhof,
				&mu,
				&dx, &dt,
				&Us,
				&gx, &gy, &gz,
				sz, g);
		copy_(
				uycp,
				uyc0,
				sz, g);
	}
	plsUY0->ImposeBoundaryCondition(blockManager, plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt, plsb);

	if( g_pconf->masking == true ) {
		pils->BiCGSTAB_Mask(
							blockManager,
							plsUY0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUY1,
							plsMaskId,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUY,
							residualUY);
	} else {
		pils->BiCGSTAB(
							blockManager,
							plsUY0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUY1,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUY,
							residualUY);
	}
	plsUY0->ImposeBoundaryCondition(blockManager);
/////////////////////////////////////////////
}

void Solver::UpdateUZ() {
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;

		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* uz0  = plsUZ0->GetBlockData(block);
		real* uzc0 = plsUZC->GetBlockData(block);
		real* uzcp = plsUZCP->GetBlockData(block);
		real* uzd0 = plsUZD->GetBlockData(block);
		real* p0   = plsP0->GetBlockData(block);

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real* ux0  = plsUX0->GetBlockData(block);
		real* uy0  = plsUY0->GetBlockData(block);

		real Us = 0.0;

		if( conf.advection == "W3" ) {
			bcut_calc_c_f_w3_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					sz, g);
		} else if( conf.advection == "E3" ) {
			bcut_calc_c_f_e3_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					sz, g);
		} else if( conf.advection == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					uzc0,
					uz0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Us,
					sz, g);
		}
		int axis=2;
		real gx = conf.gx;
		real gy = conf.gy;
		real gz = conf.gz;
		bcut_calc_ab_u_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				ux0, uy0, uz0,
				uzc0, uzcp,
				uzd0,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&axis,
				&rhof,
				&mu,
				&dx, &dt,
				&Us,
				&gx, &gy, &gz,
				sz, g);
		copy_(
				uzcp,
				uzc0,
				sz, g);
	}
	plsUZ0->ImposeBoundaryCondition(blockManager, plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt, plsb);

	if( g_pconf->masking == true ) {
		pils->BiCGSTAB_Mask(
							blockManager,
							plsUZ0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUZ1,
							plsMaskId,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUZ,
							residualUZ);
	} else {
		pils->BiCGSTAB(
							blockManager,
							plsUZ0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsUZ1,
							omegaU,
							countPreConditionerU,
							countMaxU,
							epsilonU,
							countUZ,
							residualUZ);
	}
	plsUZ0->ImposeBoundaryCondition(blockManager);
/////////////////////////////////////////////
}

void Solver::UpdateP() {
/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;

		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* p0 = plsP0->GetBlockData(block);

		real* ux = plsUX0->GetBlockData(block);
		real* uy = plsUY0->GetBlockData(block);
		real* uz = plsUZ0->GetBlockData(block);

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		bcut_remove_p_(
				ux, uy, uz,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&dx, &dt,
				sz, g);
	}

	plsUX0->ImposeBoundaryCondition(blockManager);
	plsUY0->ImposeBoundaryCondition(blockManager);
	plsUZ0->ImposeBoundaryCondition(blockManager);

#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;

		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* p0 = plsP0->GetBlockData(block);

		real* ux = plsUX0->GetBlockData(block);
		real* uy = plsUY0->GetBlockData(block);
		real* uz = plsUZ0->GetBlockData(block);

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		bcut_calc_ab_p_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				vw, ve, vs, vn, vb, vt,
				p0,
				ux, uy, uz,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&dx, &dt,
				sz, g);
	}

	plsP0->ImposeBoundaryCondition(blockManager, plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt, plsb);

	if( g_pconf->masking == true ) {
		pils->BiCGSTAB_Mask(
							blockManager,
							plsP0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsP1,
							plsMaskId,
							omegaP,
							countPreConditionerP,
							countMaxP,
							epsilonP,
							countP,
							residualP);
	} else {
		pils->BiCGSTAB(
							blockManager,
							plsP0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsP1,
							omegaP,
							countPreConditionerP,
							countMaxP,
							epsilonP,
							countP,
							residualP);
	}
	plsP0->ImposeBoundaryCondition(blockManager);
/////////////////////////////////////////////
}

void Solver::UpdateU() {
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;

		real* ux0 = plsUX0->GetBlockData(block);
		real* uy0 = plsUY0->GetBlockData(block);
		real* uz0 = plsUZ0->GetBlockData(block);
		real* p0  = plsP0 ->GetBlockData(block);
		real* lapp= plsLapP->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		bcut_corr_u_(
				ux0, uy0, uz0,
				vw, ve, vs, vn, vb, vt,
				lapp,
				p0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof,
				&dx, &dt,
				sz, g);
	}
	plsUX0->ImposeBoundaryCondition(blockManager);
	plsUY0->ImposeBoundaryCondition(blockManager);
	plsUZ0->ImposeBoundaryCondition(blockManager);
	plsP0->ImposeBoundaryCondition(blockManager);
	plsLapP->ImposeBoundaryCondition(blockManager);
}

void Solver::UpdateT() {
/////////////////////////////////////////////
// Calc A & b
/////////////////////////////////////////////
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;

		real* Ap = plsAp->GetBlockData(block);
		real* Aw = plsAw->GetBlockData(block);
		real* Ae = plsAe->GetBlockData(block);
		real* As = plsAs->GetBlockData(block);
		real* An = plsAn->GetBlockData(block);
		real* Ab = plsAb->GetBlockData(block);
		real* At = plsAt->GetBlockData(block);
		real* b  = plsb ->GetBlockData(block);

		real* t0  = plsT0->GetBlockData(block);
		real* tc0 = plsTC->GetBlockData(block);
		real* tcp = plsTCP->GetBlockData(block);
		real* td0 = plsTD->GetBlockData(block);

		real* vw = plsVw->GetBlockData(block);
		real* ve = plsVe->GetBlockData(block);
		real* vs = plsVs->GetBlockData(block);
		real* vn = plsVn->GetBlockData(block);
		real* vb = plsVb->GetBlockData(block);
		real* vt = plsVt->GetBlockData(block);

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);

		int* pPhaseId = plsPhaseId->GetBlockData(block);

		real Tc = 1.0;

		bcut_calc_d_t_(
				td0,
				t0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof, &rhos,
				&cpf, &cps,
				&kf, &ks,
				&dx, &dt,
				&Tc,
				sz, g);
		if( conf.advection == "W3" ) {
			bcut_calc_c_f_w3_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		} else if( conf.advection == "E3" ) {
			bcut_calc_c_f_e3_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		} else if( conf.advection == "Blend" ) {
			real alpha = 0.95;
			bcut_calc_c_f_blend_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					&alpha,
					sz, g);
		} else {
			bcut_calc_c_f_c2_(
					tc0,
					t0,
					vw, ve, vs, vn, vb, vt,
					pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
					pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
					pPhaseId,
					&dx, &dt,
					&Tc,
					sz, g);
		}

/*
		bcut_calc_ab_t_1st_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				t0,
				tc0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof, &rhos,
				&cpf, &cps,
				&kf, &ks,
				&dx, &dt,
				&Tc,
				sz, g);
*/

		bcut_calc_ab_t_(
				Ap, Aw, Ae, As, An, Ab, At, b,
				t0,
				tc0, tcp,
				td0,
				pCut0, pCut1, pCut2, pCut3, pCut4, pCut5,
				pCutId0, pCutId1, pCutId2, pCutId3, pCutId4, pCutId5,
				pPhaseId,
				&rhof, &rhos,
				&cpf, &cps,
				&kf, &ks,
				&dx, &dt,
				&Tc,
				sz, g);
		copy_(
				tcp,
				tc0,
				sz, g);
	}
	plsT0->ImposeBoundaryCondition(blockManager, plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt, plsb);
/////////////////////////////////////////////

/////////////////////////////////////////////
// Implicit 1st-order Euler method
/////////////////////////////////////////////
	if( g_pconf->masking == true ) {
		pils->BiCGSTAB_Mask(
							blockManager,
							plsT0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsT1,
							plsMaskId,
							omegaT,
							countPreConditionerT,
							countMaxT,
							epsilonT,
							countT,
							residualT);
	} else {
		pils->BiCGSTAB(
							blockManager,
							plsT0,
							plsAp,
							plsAw,
							plsAe,
							plsAs,
							plsAn,
							plsAb,
							plsAt,
							plsb,
							plsr,
							plsr0,
							plsp,
							plsp_,
							plsq_,
							plss,
							plss_,
							plst_,
							plsT1,
							omegaT,
							countPreConditionerT,
							countMaxT,
							epsilonT,
							countT,
							residualT);
	}
	plsT0->ImposeBoundaryCondition(blockManager);
/////////////////////////////////////////////
}

#include <sys/time.h>
double Solver::GetTime() {
	struct timeval tp;
	int i = gettimeofday(&tp, 0);
	return ((double)(tp.tv_sec) + (double)(tp.tv_usec)*1.0e-6);
}

void Solver::WritePolygon(std::ofstream& ofs, float* pv) {

	ofs << "facet normal 0 0 0" << std::endl;
	ofs << "outer loop" << std::endl;
	ofs << "vertex " << pv[0] << " " << pv[1] << " " << pv[2] << std::endl;
	ofs << "vertex " << pv[3] << " " << pv[4] << " " << pv[5] << std::endl;
	ofs << "vertex " << pv[6] << " " << pv[7] << " " << pv[8] << std::endl;
	ofs << "endloop" << std::endl;
	ofs << "endfacet" << std::endl;

}

void Solver::PrintCut() {
	mkdir("STL", 0755);

	std::ostringstream ossFileName;
	ossFileName << "./STL/";
	ossFileName << "data-cut";
	ossFileName << "-";
	ossFileName.width(5);
	ossFileName.setf(std::ios::fixed);
	ossFileName.fill('0');
	ossFileName << myrank;
	ossFileName << ".stl";

	std::ofstream ofs;
	ofs.open(ossFileName.str().c_str(), std::ios::out);
	ofs.close();

#ifdef _BLOCK_IS_LARGE_
#else
#endif
	for (int n=0; n<blockManager.getNumBlock(); ++n) {
		BlockBase* block = blockManager.getBlock(n);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		int nc[3] = {size.x + 2*vc, size.y + 2*vc, size.z + 2*vc};
		real dx = cellSize.x;

		real* pCut0 = plsCut0->GetBlockData(block);
		real* pCut1 = plsCut1->GetBlockData(block);
		real* pCut2 = plsCut2->GetBlockData(block);
		real* pCut3 = plsCut3->GetBlockData(block);
		real* pCut4 = plsCut4->GetBlockData(block);
		real* pCut5 = plsCut5->GetBlockData(block);
		int* pCutId0 = plsCutId0->GetBlockData(block);
		int* pCutId1 = plsCutId1->GetBlockData(block);
		int* pCutId2 = plsCutId2->GetBlockData(block);
		int* pCutId3 = plsCutId3->GetBlockData(block);
		int* pCutId4 = plsCutId4->GetBlockData(block);
		int* pCutId5 = plsCutId5->GetBlockData(block);

		ofs.open(ossFileName.str().c_str(), std::ios::app);
		ofs << "solid" << std::endl;
		for(int k=vc; k<=size.z+vc-1; k++) {
			for(int j=vc; j<=size.y+vc-1; j++) {
				for(int i=vc; i<=size.x+vc-1; i++) {
					int mp = i + nc[0]*( j + nc[1]*k );

					float x0 = origin[0] + (i - vc)*dx;
					float y0 = origin[1] + (j - vc)*dx;
					float z0 = origin[2] + (k - vc)*dx;
					float x1 = origin[0] + (i + 1 - vc)*dx;
					float y1 = origin[1] + (j + 1 - vc)*dx;
					float z1 = origin[2] + (k + 1 - vc)*dx;

					float v[9];
					float x[2] = {x0, x1};
					float y[2] = {y0, y1};
					float z[2] = {z0, z1};

					int cidp0 = pCutId0[mp];
					int cidp1 = pCutId1[mp];
					int cidp2 = pCutId2[mp];
					int cidp3 = pCutId3[mp];
					int cidp4 = pCutId4[mp];
					int cidp5 = pCutId5[mp];

					if( cidp0 != 0 ) {
					}
					if( cidp1 != 0 ) {
						v[0] = x[1];
						v[1] = y[0];
						v[2] = z[0];
						v[3] = x[1];
						v[4] = y[0];
						v[5] = z[1];
						v[6] = x[1];
						v[7] = y[1];
						v[8] = z[0];
						WritePolygon(ofs, v);

						v[0] = x[1];
						v[1] = y[1];
						v[2] = z[1];
						v[6] = x[1];
						v[7] = y[0];
						v[8] = z[1];
						v[3] = x[1];
						v[4] = y[1];
						v[5] = z[0];
						WritePolygon(ofs, v);
					}
					if( cidp2 != 0 ) {
					}
					if( cidp3 != 0 ) {
						v[0] = x[0];
						v[1] = y[1];
						v[2] = z[0];
						v[6] = x[0];
						v[7] = y[1];
						v[8] = z[1];
						v[3] = x[1];
						v[4] = y[1];
						v[5] = z[0];
						WritePolygon(ofs, v);

						v[0] = x[1];
						v[1] = y[1];
						v[2] = z[1];
						v[6] = x[1];
						v[7] = y[1];
						v[8] = z[0];
						v[3] = x[0];
						v[4] = y[1];
						v[5] = z[1];
						WritePolygon(ofs, v);
					}
					if( cidp4 != 0 ) {
					}
					if( cidp5 != 0 ) {
						v[0] = x[0];
						v[1] = y[0];
						v[2] = z[1];
						v[6] = x[0];
						v[7] = y[1];
						v[8] = z[1];
						v[3] = x[1];
						v[4] = y[0];
						v[5] = z[1];
						WritePolygon(ofs, v);

						v[0] = x[1];
						v[1] = y[1];
						v[2] = z[1];
						v[3] = x[0];
						v[4] = y[1];
						v[5] = z[1];
						v[6] = x[1];
						v[7] = y[0];
						v[8] = z[1];
						WritePolygon(ofs, v);
					}

				}
			}
		}
		ofs << "endsolid" << std::endl;
		ofs.close();
	}
}

void Solver::Dump(const int step) {
	plsUX0->Dump(blockManager, step, "ux");
	plsUY0->Dump(blockManager, step, "uy");
	plsUZ0->Dump(blockManager, step, "uz");

	plsP0->Dump(blockManager, step, "p");

	plsVw->Dump(blockManager, step, "vw");
	plsVe->Dump(blockManager, step, "ve");
	plsVs->Dump(blockManager, step, "vs");
	plsVn->Dump(blockManager, step, "vn");
	plsVb->Dump(blockManager, step, "vb");
	plsVt->Dump(blockManager, step, "vt");

	plsT0->Dump(blockManager, step, "t");

	plsUXCP->Dump(blockManager, step, "uxcp");
	plsUYCP->Dump(blockManager, step, "uycp");
	plsUZCP->Dump(blockManager, step, "uzcp");
	plsTCP->Dump(blockManager, step, "tcp");
}

void Solver::Load(const int step) {
	plsUX0->Load(blockManager, step, "ux");
	plsUY0->Load(blockManager, step, "uy");
	plsUZ0->Load(blockManager, step, "uz");

	plsP0->Load(blockManager, step, "p");

	plsVw->Load(blockManager, step, "vw");
	plsVe->Load(blockManager, step, "ve");
	plsVs->Load(blockManager, step, "vs");
	plsVn->Load(blockManager, step, "vn");
	plsVb->Load(blockManager, step, "vb");
	plsVt->Load(blockManager, step, "vt");

	plsT0->Load(blockManager, step, "t");

	plsUXCP->Load(blockManager, step, "uxcp");
	plsUYCP->Load(blockManager, step, "uycp");
	plsUZCP->Load(blockManager, step, "uzcp");
	plsTCP->Load(blockManager, step, "tcp");
}

void Solver::Dump2(const int step) {
	plsUX0->Dump2(blockManager, step, "ux");
	plsUY0->Dump2(blockManager, step, "uy");
	plsUZ0->Dump2(blockManager, step, "uz");

	plsP0->Dump2(blockManager, step, "p");

	plsVw->Dump2(blockManager, step, "vw");
	plsVe->Dump2(blockManager, step, "ve");
	plsVs->Dump2(blockManager, step, "vs");
	plsVn->Dump2(blockManager, step, "vn");
	plsVb->Dump2(blockManager, step, "vb");
	plsVt->Dump2(blockManager, step, "vt");

	plsT0->Dump2(blockManager, step, "t");

	plsUXCP->Dump2(blockManager, step, "uxcp");
	plsUYCP->Dump2(blockManager, step, "uycp");
	plsUZCP->Dump2(blockManager, step, "uzcp");
	plsTCP->Dump2(blockManager, step, "tcp");
}

void Solver::Load2(const int step) {
	plsUX0->Load2(blockManager, step, "ux");
	plsUY0->Load2(blockManager, step, "uy");
	plsUZ0->Load2(blockManager, step, "uz");

	plsP0->Load2(blockManager, step, "p");

	plsVw->Load2(blockManager, step, "vw");
	plsVe->Load2(blockManager, step, "ve");
	plsVs->Load2(blockManager, step, "vs");
	plsVn->Load2(blockManager, step, "vn");
	plsVb->Load2(blockManager, step, "vb");
	plsVt->Load2(blockManager, step, "vt");

	plsT0->Load2(blockManager, step, "t");

	plsUXCP->Load2(blockManager, step, "uxcp");
	plsUYCP->Load2(blockManager, step, "uycp");
	plsUZCP->Load2(blockManager, step, "uzcp");
	plsTCP->Load2(blockManager, step, "tcp");
}

void Solver::Dump3(const int step) {
	plsUX0->Dump3(blockManager, step, "ux", partition, myrank);
	plsUY0->Dump3(blockManager, step, "uy", partition, myrank);
	plsUZ0->Dump3(blockManager, step, "uz", partition, myrank);

	plsP0->Dump3(blockManager, step, "p", partition, myrank);

	plsVw->Dump3(blockManager, step, "vw", partition, myrank);
	plsVe->Dump3(blockManager, step, "ve", partition, myrank);
	plsVs->Dump3(blockManager, step, "vs", partition, myrank);
	plsVn->Dump3(blockManager, step, "vn", partition, myrank);
	plsVb->Dump3(blockManager, step, "vb", partition, myrank);
	plsVt->Dump3(blockManager, step, "vt", partition, myrank);

	plsT0->Dump3(blockManager, step, "t", partition, myrank);

	plsUXCP->Dump3(blockManager, step, "uxcp", partition, myrank);
	plsUYCP->Dump3(blockManager, step, "uycp", partition, myrank);
	plsUZCP->Dump3(blockManager, step, "uzcp", partition, myrank);
	plsTCP->Dump3(blockManager, step, "tcp", partition, myrank);
}

void Solver::Load3(const int step) {
	plsUX0->Load3(blockManager, step, "ux", partition, myrank);
	plsUY0->Load3(blockManager, step, "uy", partition, myrank);
	plsUZ0->Load3(blockManager, step, "uz", partition, myrank);

	plsP0->Load3(blockManager, step, "p", partition, myrank);

	plsVw->Load3(blockManager, step, "vw", partition, myrank);
	plsVe->Load3(blockManager, step, "ve", partition, myrank);
	plsVs->Load3(blockManager, step, "vs", partition, myrank);
	plsVn->Load3(blockManager, step, "vn", partition, myrank);
	plsVb->Load3(blockManager, step, "vb", partition, myrank);
	plsVt->Load3(blockManager, step, "vt", partition, myrank);

	plsT0->Load3(blockManager, step, "t", partition, myrank);

	plsUXCP->Load3(blockManager, step, "uxcp", partition, myrank);
	plsUYCP->Load3(blockManager, step, "uycp", partition, myrank);
	plsUZCP->Load3(blockManager, step, "uzcp", partition, myrank);
	plsTCP->Load3(blockManager, step, "tcp", partition, myrank);
}
