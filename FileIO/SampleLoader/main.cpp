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
#include "BoundaryConditionSetter.h"
#include "BlockFactory.h"
#include "Block.h"
#include "BlockManager.h"
#include "Scalar3D.h"

#include "BCMFileLoader.h"

#include "BCMFileSaver.h"

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


int main(int argc, char** argv)
{
	using namespace std;	

	if(argc != 3){
		printf("err : useage %s <input (cellid.bcm)> <input (data.bcm)>\n", argv[0]);
		return -1;
	}

	MPI::Init(argc, argv);

	BoundaryConditionSetter* bcsetter = new BoundaryConditionSetter;

	// argv[1] (cellid.bcm)の読み込み
	BCMFileIO::BCMFileLoader loader(argv[1], bcsetter);

	delete bcsetter;

	//	読み込んだOctreeのレイアウト表示
	BlockManager& blockManager = BlockManager::getInstance();
	blockManager.printBlockLayoutInfo();

	// argv[2] (data.bcm)の読み込み
	if( !loader.LoadAdditionalIndex(argv[2]) ){
		printf("err : Load File Error %s)\n", argv[2]);
		return -1;
	}
	
	// Temperatureのタイムステップを取得
	const BCMFileIO::IdxStep* tmpStep = loader.GetStep("Tmp32");
	
	// タイムステップのリストを取得
	const std::list<unsigned int>* tmpStepList = tmpStep->GetStepList();

	int vc = 2;
	// リーフブロック読込
	// CellIDは直接読み込む
	int id_cid = 0;
	int id_s32 = 0;
	int id_s64 = 0;
	int id_v32[3] = {0}; // 0:u, 1:v, 2:w
	int id_v64[3] = {0}; // 0:u, 1:v, 2:w

	loader.LoadLeafBlock(&id_cid, "CellID", vc);

	// Check
	//print<unsigned char>("cid", cellid_dcid, vc);

	// Temperatureは一度ブロックを生成する
	loader.CreateLeafBlock(&id_s32, "Tmp32", vc);
	loader.CreateLeafBlock(&id_s64, "Tmp64", vc);

	loader.CreateLeafBlock(id_v32, "Vel32", vc, true);
	loader.CreateLeafBlock(id_v64, "Vel64", vc, false);

	// 確認用ファイルセーバの準備
	BCMFileIO::BCMFileSaver saver(loader.GetGlobalOrigin(), loader.GetGlobalRegion(), loader.GetOctree(), "out");
	// 単位系をローダから取得し、セーバに渡す
	saver.SetUnit(loader.GetUnit());

	// CellIDの登録
	saver.RegisterCellIDInformation(id_cid, 5, vc, "CellID", "cid", "lb", "cid", true);
	// リーフブロック情報を登録 (Scalar)
	saver.RegisterDataInformation(&id_s32, BCMFileIO::LB_SCALAR,  BCMFileIO::LB_FLOAT32, vc, "Tmp32", "Tmp32", "lb", *tmpStep, "TMP");
	saver.RegisterDataInformation(&id_s64, BCMFileIO::LB_SCALAR,  BCMFileIO::LB_FLOAT64, vc, "Tmp64", "Tmp64", "lb", *tmpStep, "TMP");

	// リーフブロック情報を登録 (Vector)
	saver.RegisterDataInformation( id_v32, BCMFileIO::LB_VECTOR3, BCMFileIO::LB_FLOAT32, vc, "Vel32", "Vel32", "lb", *tmpStep, "VEL");
	saver.RegisterDataInformation( id_v64, BCMFileIO::LB_VECTOR3, BCMFileIO::LB_FLOAT64, vc, "Vel64", "Vel64", "lb", *tmpStep, "VEL");

	// インデックスファイルとOctreeを出力
	saver.Save();

	// CellIDのリーフブロックを出力
	saver.SaveLeafBlock("CellID");

	// Temperatureのステップでループ
	for(std::list<unsigned int>::const_iterator it = tmpStepList->begin(); it != tmpStepList->end(); ++it){
		// ステップに該当するリーフブロックを読み込み
		loader.LoadLeafBlock(&id_s32, "Tmp32", vc, *it);
		loader.LoadLeafBlock(&id_s64, "Tmp64", vc, *it);

		loader.LoadLeafBlock( id_v32, "Vel32", vc, *it);
		loader.LoadLeafBlock( id_v64, "Vel64", vc, *it);
		
		// 各データの中身をテキストファイルに出力 (for debug)
		char prefix[128];
		sprintf(prefix, "tmp32_%03d", *it);   print<float >(prefix, id_s32, vc);
		sprintf(prefix, "tmp64_%03d", *it);   print<double>(prefix, id_s64, vc);
                                             
		sprintf(prefix, "vel32_u_%03d", *it); print<float >(prefix, id_v32[0], vc);
		sprintf(prefix, "vel32_v_%03d", *it); print<float >(prefix, id_v32[1], vc);
		sprintf(prefix, "vel32_w_%03d", *it); print<float >(prefix, id_v32[2], vc);
                                             
		sprintf(prefix, "vel64_u_%03d", *it); print<double>(prefix, id_v64[0], vc);
		sprintf(prefix, "vel64_v_%03d", *it); print<double>(prefix, id_v64[1], vc);
		sprintf(prefix, "vel64_w_%03d", *it); print<double>(prefix, id_v64[2], vc);

		// リーフブロックファイルを出力
		//saver.SaveLeafBlock("Tmp32", *it);
		//saver.SaveLeafBlock("Tmp64", *it);
		//saver.SaveLeafBlock("Vel32", *it);
		//saver.SaveLeafBlock("Vel64", *it);
	}

	delete tmpStepList;

	MPI::Finalize();

	return EX_SUCCESS;
}

