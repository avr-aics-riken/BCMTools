/*
 * BCMTools
 *
 * Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <iostream>
#include "mpi.h"
#include "BCMTools.h"
#include "Config.h"
#include "RootGrid.h"
#include "BCMPolylib.h"
#include "BCMOctree.h"
#include "PolygonDivider.h"
#include "BoundaryConditionSetter.h"
#include "BlockFactory.h"
#include "Block.h"
#include "BlockManager.h"
#include "BlockBoundingBox.h"
#include "Scalar3D.h"
#include "Vector3D.h"

#include "BoundingBox.h"

#include "BCMFileSaver.h"
#include "BCMFileCommon.h"


template <typename T>
void print(const char* prefix, const int dcid, const int vc = 0){
	BlockManager& blockManager = BlockManager::getInstance();
	const MPI::Intracomm& comm = blockManager.getCommunicator();
	Vec3i sz = blockManager.getSize();

	char name[128]; sprintf(name, "%s_%03d.txt", prefix, comm.Get_rank());
	FILE *fp = fopen(name, "w");
	for(int id = 0; id < blockManager.getNumBlock(); ++id){
		BlockBase* block = blockManager.getBlock(id);

		Scalar3D<T> *mesh = dynamic_cast< Scalar3D<T>* >(block->getDataClass(dcid));
		T* data = mesh->getData();
		Index3DS idx = mesh->getIndex();

		fprintf(fp, "Data[%3d] (vc : %d) \n", id, vc);
		for(int z = -vc; z < sz.z + vc; z++){
			for(int y = -vc; y < sz.y + vc; y++){
				for(int x = -vc; x < sz.x + vc; x++){
					fprintf(fp, "%g ", data[idx(x, y, z)]);
				}
				fprintf(fp, "\n");
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

template <>
void print<unsigned char>(const char* prefix, const int dcid, const int vc){
	BlockManager& blockManager = BlockManager::getInstance();
	const MPI::Intracomm& comm = blockManager.getCommunicator();
	Vec3i sz = blockManager.getSize();

	char name[128]; sprintf(name, "%s_%03d.txt", prefix, comm.Get_rank());
	FILE *fp = fopen(name, "w");
	for(int id = 0; id < blockManager.getNumBlock(); ++id){
		BlockBase* block = blockManager.getBlock(id);

		Scalar3D<unsigned char> *mesh = dynamic_cast< Scalar3D<unsigned char>* >(block->getDataClass(dcid));
		unsigned char* data = mesh->getData();
		Index3DS idx = mesh->getIndex();

		fprintf(fp, "Data[%3d] (vc : %d) \n", id, vc);
		for(int z = -vc; z < sz.z + vc; z++){
			for(int y = -vc; y < sz.y + vc; y++){
				for(int x = -vc; x < sz.x + vc; x++){
					fprintf(fp, "%d ", data[idx(x, y, z)]);
				}
				fprintf(fp, "\n");
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

template <typename T>
bool SetScalarData(const int dataClassId, const int vc, const unsigned int step)
{
	BlockManager& blockManager = BlockManager::getInstance();
	Vec3i sz = blockManager.getSize();


	for(int id = 0; id < blockManager.getNumBlock(); ++id){
		BlockBase* block = blockManager.getBlock(id);
		Scalar3D<T>* mesh = dynamic_cast< Scalar3D<T>* >(block->getDataClass(dataClassId));
		T* data = mesh->getData();
		Index3DS idx = mesh->getIndex();

		const Vec3r& org = block->getOrigin();
		const Vec3r& pitch = block->getCellSize();

		for(int z = -vc; z < sz.z + vc; z++){
			for(int y = -vc; y < sz.y + vc; y++){
				for(int x = -vc; x < sz.x + vc; x++){
					T u = static_cast<T>( org.x + (x + 0.5) * pitch.x );
					T v = static_cast<T>( org.y + (y + 0.5) * pitch.y );
					T w = static_cast<T>( org.z + (z + 0.5) * pitch.z );
					data[idx(x, y, z)] = sqrt(u*u + v*v + w*w);
					//data[idx(x, y, z)] = id;
				}
			}
		}
	}
	return true;
}

template <typename T>
bool SetVectorData(const int dataClassId[3], const int vc, const unsigned int step)
{
	BlockManager& blockManager = BlockManager::getInstance();
	Vec3i sz = blockManager.getSize();

	for(int id = 0; id < blockManager.getNumBlock(); ++id)
	{
		BlockBase* block = blockManager.getBlock(id);
		Scalar3D<T>* meshU = dynamic_cast< Scalar3D<T>* >(block->getDataClass(dataClassId[0]));
		Scalar3D<T>* meshV = dynamic_cast< Scalar3D<T>* >(block->getDataClass(dataClassId[1]));
		Scalar3D<T>* meshW = dynamic_cast< Scalar3D<T>* >(block->getDataClass(dataClassId[2]));

		T* data[3]            = {  meshU->getData(),  meshV->getData(),  meshW->getData() };
		const Index3DS idx[3] = { meshU->getIndex(), meshV->getIndex(), meshW->getIndex() };

		const Vec3r& org = block->getOrigin();
		const Vec3r& pitch = block->getCellSize();

		for(int z = -vc; z < sz.z + vc; z++){
			for(int y = -vc; y < sz.y + vc; y++){
				for(int x = -vc; x < sz.x + vc; x++){
					T u = static_cast<T>( org.x + (x + 0.5) * pitch.x );
					T v = static_cast<T>( org.y + (y + 0.5) * pitch.y );
					T w = static_cast<T>( org.z + (z + 0.5) * pitch.z );
					data[0][idx[0](x, y, z)] = u;
					data[1][idx[1](x, y, z)] = v;
					data[2][idx[2](x, y, z)] = w;
					/*
					data[0][idx[0](x, y, z)] = id * 3 + 0;
					data[1][idx[1](x, y, z)] = id * 3 + 1;
					data[2][idx[2](x, y, z)] = id * 3 + 2;
					*/
				}
			}
		}
	}
	return true;
}

BoundingBox defineSearchRegion(const Pedigree& pedigree, int maxLevel, 
                               const Vec3r& origin, const double rootLength, const RootGrid* rootGrid)
{
	int level = pedigree.getLevel();
	int px = pedigree.getX();
	int py = pedigree.getY();
	int pz = pedigree.getZ();
	
	int rootID = pedigree.getRootID();
	int ix = rootGrid->rootID2indexX(rootID);
	int iy = rootGrid->rootID2indexY(rootID);
	int iz = rootGrid->rootID2indexZ(rootID);
	
	int max0 = pedigree.getUpperBound() - 1;  // 2^level - 1
	double d = 1.0 / (max0 + 1);  // ブロックサイズ
	
	double x0 = ix + px * d;
	double y0 = iy + py * d;
	double z0 = iz + pz * d;
	double x1 = ix + (px + 1) * d;
	double y1 = iy + (py + 1) * d;
	double z1 = iz + (pz + 1) * d;
	
	// リップル効果対策用のマージン設定
	if (level < maxLevel) {
	  int n = 1 << (maxLevel - level - 1);
	  double dd = d * (double)(n - 1) / n;
	  if (!(px == 0 && rootGrid->isOuterBoundary(rootID, X_M))) x0 -= dd;
	  if (!(py == 0 && rootGrid->isOuterBoundary(rootID, Y_M))) y0 -= dd;
	  if (!(pz == 0 && rootGrid->isOuterBoundary(rootID, Z_M))) z0 -= dd;
	  if (!(px == max0 && rootGrid->isOuterBoundary(rootID, X_P))) x1 += dd;
	  if (!(py == max0 && rootGrid->isOuterBoundary(rootID, Y_P))) y1 += dd;
	  if (!(pz == max0 && rootGrid->isOuterBoundary(rootID, Z_P))) z1 += dd;
	}
	
	// 原点移動，スケーリング
	x0 = origin.x + x0 * rootLength;
	y0 = origin.y + y0 * rootLength;
	z0 = origin.z + z0 * rootLength;
	x1 = origin.x + x1 * rootLength;
	y1 = origin.y + y1 * rootLength;
	z1 = origin.z + z1 * rootLength;
	
	return BoundingBox(x0, y0, z0, x1, y1, z1);
}


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

		std::cout <<	std::endl << "Configuration file: " << argv[1] << std::endl;

	}

	conf.load(argv[1]);

	RootGrid* rootGrid = new RootGrid(conf.rootN);

	if (myRank == 0) {
		conf.print();

		if (conf.polygonGroupList.size() > 0) {
			if (pl->load(conf.polylibConfig) != PLSTAT_OK) {
				comm.Abort(EX_OPEN_FILE);
			}
		}
		
		/*
		const std::string& polygonGroup = conf.polygonGroupList[0].polygonGroupName;
		PolylibNS::Vec3f min(-150.0, -13.4, -300.0);
		PolylibNS::Vec3f max( 300.0, 200.0,  600.0);
		std::vector<PolylibNS::Triangle*>* polygonList = pl->search_polygons(polygonGroup, min, max, false);
		*/

		// 仮想セル領域も考慮してポリゴン検索
		PolygonDivider* divider = new PolygonDivider(conf.origin, conf.rootLength,
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

	BlockFactory* blockFactory = new BlockFactory(tree, 
	                                              partition, 
	                                              boundaryConditionSetter,
												  conf.origin, 
												  conf.rootLength, 
												  size);


	BlockManager& blockManager = BlockManager::getInstance();

	for (int id = partition->getStart(myRank); id < partition->getEnd(myRank); id++) {
		Node* node = leafNodeArray[id];
		Block* block = blockFactory->makeBlock(node);
		blockManager.registerBlock(block);
	}

	blockManager.endRegisterBlock();

	blockManager.printBlockLayoutInfo();

	Vec3i sz = blockManager.getSize();
	int id_cellId = blockManager.setDataClass< Scalar3D<unsigned char> >(conf.vc);

	// Data Set for CellID
#if 0
	for (int id = 0; id < blockManager.getNumBlock(); ++id){
		BlockBase* block = blockManager.getBlock(id);

		Scalar3D<unsigned char>* mesh = dynamic_cast< Scalar3D<unsigned char>* >(block->getDataClass(id_cellId));
		unsigned char* data = mesh->getData();
		Index3DS idx = mesh->getIndex();

		for(int z = -conf.vc; z < sz.z + conf.vc; z++){
			for(int y = -conf.vc; y < sz.y + conf.vc; y++){
				for(int x = -conf.vc; x < sz.x + conf.vc; x++){
					//data[idx(x, y, z)] = (myRank % 32);
					data[idx(x, y, z)] = id % 32;
				}
			}
		}
	}
#else
	// リーフブロックの各セルとPolygonの交差判定 (Polylibを使って)
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

	for (int id = 0; id < blockManager.getNumBlock(); ++id){
		BlockBase* block = blockManager.getBlock(id);
		Scalar3D<unsigned char>* mesh = dynamic_cast< Scalar3D<unsigned char>* >(block->getDataClass(id_cellId));
		unsigned char* data = mesh->getData();
		Index3DS idx = mesh->getIndex();

		int lid = partition->getStart(myRank) + id;
		Pedigree p = leafNodeArray[lid]->getPedigree();
		int lv = p.getLevel();

		Vec3r rootRegion(conf.rootLength, conf.rootLength, conf.rootLength);
		Vec3r rootOrg(conf.origin.x + (rootRegion.x * rootGrid->rootID2indexX(p.getRootID())),
		              conf.origin.y + (rootRegion.y * rootGrid->rootID2indexY(p.getRootID())),
					  conf.origin.z + (rootRegion.z * rootGrid->rootID2indexZ(p.getRootID())));

		Vec3r rgn( rootRegion.x / (1 << lv),     rootRegion.y / (1 << lv),     rootRegion.z / (1 << lv));
		Vec3r org( rootOrg.x + p.getX() * rgn.x, rootOrg.y + p.getY() * rgn.y, rootOrg.z + p.getZ() * rgn.z);

		const std::string& polygonGroup = conf.polygonGroupList[0].polygonGroupName;
		int nPolygons = 0;
		std::vector<PolylibNS::Triangle*>* polygonList;

		PolylibNS::Vec3f bmin(org.x,         org.y,         org.z);
		PolylibNS::Vec3f bmax(org.x + rgn.x, org.y + rgn.y, org.z + rgn.z);
		polygonList = pl->search_polygons(polygonGroup, bmin, bmax, false);
		nPolygons = polygonList->size();
		delete polygonList;
		
		if( nPolygons > 0 ){
			double dx = rgn.x / (double)(sz.x);
			double dy = rgn.y / (double)(sz.y);
			double dz = rgn.z / (double)(sz.z);

			for(int z = 0; z < sz.z; z++){
				for(int y = 0; y < sz.y; y++){
					for(int x = 0; x < sz.x; x++){
						PolylibNS::Vec3f min(org.x + (dx*x), org.y + (dy*y), org.z + (dz*z));
						PolylibNS::Vec3f max(min[0] + dx, min[1] + dy, min[2] + dz);
						polygonList = pl->search_polygons(polygonGroup, min, max, false);
						if(polygonList->size() > 0){
							data[idx(x, y, z)] = 1;
						}else{
							data[idx(x, y, z)] = 0;
						}
						delete polygonList;
					}
				}
			}
		}else{
			for(int z = 0; z < sz.z; z++){
				for(int y = 0; y < sz.y; y++){
					for(int x = 0; x < sz.x; x++){
						data[idx(x, y, z)] = 0;
					}
				}
			}
		}
	}
#endif

	//print<unsigned char>("cid", id_cellId, conf.vc);

	// ファイル出力用にGlobalOrigin / GlobalRegionを設定
	Vec3r origin = conf.origin;
	Vec3r region( conf.rootLength * conf.rootN.x,
	              conf.rootLength * conf.rootN.y,
				  conf.rootLength * conf.rootN.z);
	
	BCMFileIO::IdxUnit unit;
	unit.length   = std::string("m");
	unit.L0_scale = 1.0;
	unit.velocity = std::string("Dimensional");
	unit.V0_scale = 0.1;

	//BCMFileIO::IdxStep tmpStep(0, 20, 4);
	//tmpStep.AddStep(2);
	//tmpStep.AddStep(10);
	//tmpStep.AddStep(15);
	//tmpStep.SubStep(16);
	BCMFileIO::IdxStep tmpStep(0, 1);

	// BCMFileSaver インスタンス生成
	BCMFileIO::BCMFileSaver saver(origin, region, tree, "out");
	// リーフブロック情報を登録 (CellID)
	saver.RegisterCellIDInformation(id_cellId, 5, conf.vc, "CellID", "cid", "lb", "cid");
	
	int id_s32 = blockManager.setDataClass< Scalar3D<float > >(conf.vc);
	int id_s64 = blockManager.setDataClass< Scalar3D<double> >(conf.vc);

	int id_v32[3] = {0}; // for Vector float
	id_v32[0] = blockManager.setDataClass< Scalar3D<float > >(conf.vc); // for u
	id_v32[1] = blockManager.setDataClass< Scalar3D<float > >(conf.vc); // for v
	id_v32[2] = blockManager.setDataClass< Scalar3D<float > >(conf.vc); // for w

	int id_v64[3] = {0}; // for Vector double
	id_v64[0] = blockManager.setDataClass< Scalar3D<double> >(conf.vc); // for u
	id_v64[1] = blockManager.setDataClass< Scalar3D<double> >(conf.vc); // for v
	id_v64[2] = blockManager.setDataClass< Scalar3D<double> >(conf.vc); // for w

	// リーフブロック情報を登録 (Scalar)
	saver.RegisterDataInformation(&id_s32, BCMFileIO::LB_SCALAR,  BCMFileIO::LB_FLOAT32, conf.vc, "Tmp32", "Tmp32", "lb", tmpStep, "TMP");
	saver.RegisterDataInformation(&id_s64, BCMFileIO::LB_SCALAR,  BCMFileIO::LB_FLOAT64, conf.vc, "Tmp64", "Tmp64", "lb", tmpStep, "TMP");
	// リーフブロック情報を登録 (Vector)
	saver.RegisterDataInformation( id_v32, BCMFileIO::LB_VECTOR3, BCMFileIO::LB_FLOAT32, conf.vc, "Vel32", "Vel32", "lb", tmpStep, "VEL");
	saver.RegisterDataInformation( id_v64, BCMFileIO::LB_VECTOR3, BCMFileIO::LB_FLOAT64, conf.vc, "Vel64", "Vel64", "lb", tmpStep, "VEL");

	// 単位系の設定
	saver.SetUnit(unit);
	// インデックスファイルとBCMOctreeをファイル出力 (cellid.bcm / tree.oct)
	saver.Save();
	// リーフブロックファイルを出力 (cid.lb)
	saver.SaveLeafBlock("CellID");
	
	const std::list<unsigned int>* tmpStepList = tmpStep.GetStepList();
	for(std::list<unsigned int>::const_iterator it = tmpStepList->begin(); it != tmpStepList->end(); ++it){
		// 各データに適当な値を書き込み
		SetScalarData<float >(id_s32, conf.vc, *it);
		SetScalarData<double>(id_s64, conf.vc, *it);
		SetVectorData<float >(id_v32, conf.vc, *it);
		SetVectorData<double>(id_v64, conf.vc, *it);
		
		/*
		// 各データの中身をテキストファイルに出力 (for debug)
		char prefix[128];
		sprintf(prefix, "tmp32_%03d", *it);   print<float >(prefix, id_s32, conf.vc);
		sprintf(prefix, "tmp64_%03d", *it);   print<double>(prefix, id_s64, conf.vc);

		sprintf(prefix, "vel32_u_%03d", *it); print<float >(prefix, id_v32[0], conf.vc);
		sprintf(prefix, "vel32_v_%03d", *it); print<float >(prefix, id_v32[1], conf.vc);
		sprintf(prefix, "vel32_w_%03d", *it); print<float >(prefix, id_v32[2], conf.vc);

		sprintf(prefix, "vel64_u_%03d", *it); print<double>(prefix, id_v64[0], conf.vc);
		sprintf(prefix, "vel64_v_%03d", *it); print<double>(prefix, id_v64[1], conf.vc);
		sprintf(prefix, "vel64_w_%03d", *it); print<double>(prefix, id_v64[2], conf.vc);
		*/

		// リーフブロックファイルを出力
		saver.SaveLeafBlock("Tmp32", *it);
		saver.SaveLeafBlock("Tmp64", *it);
		saver.SaveLeafBlock("Vel32", *it);
		saver.SaveLeafBlock("Vel64", *it);
	}

	delete tmpStepList;

	delete pl;
	delete tree;
	delete partition;
	delete boundaryConditionSetter;

	MPI::Finalize();

	return EX_SUCCESS;
}

