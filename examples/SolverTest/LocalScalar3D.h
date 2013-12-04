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

#ifndef LOCALSCALAR3D_H
#define LOCALSCALAR3D_H

#include "BlockManager.h"

#include "BlockScalar3D.h"
#include "Scalar3DUpdater.h"
#include "Scalar3DUpdater1.h"
#include "Scalar3DUpdater2.h"
#include "Config.h"
#include "real.h"
#include "gv.h"

#include "VtkWriter.h"
//#include "SiloWriter.h"

template <typename T>
class LocalScalar3D {
public:
	LocalScalar3D(
					BlockManager& blockManager,
					int vc,
					std::string updateMethod,
					int* boundaryType,
					T* boundaryValue,
					int updaterType=1) {
		this->vc = vc;
		if( updaterType == 1 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater1<T> >(this->vc);
		} else if( updaterType == 10 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater1<T>, T >(this->vc);
		} else if( updaterType == 2 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater2<T> >(this->vc);
		} else {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater<T> >(this->vc);
		}

/*
		blockDataTable = new T* [blockManager.getNumBlock()];
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); n++) {
			BlockBase* block = blockManager.getBlock(n);
			T*    blockData = this->GetBlockData(block);
			this->blockDataTable[n] = blockData;
		}
*/

		this->updateMethod = VCUpdateMethod::AtOnce;
		if (updateMethod == "AtOnce") {
			this->updateMethod = VCUpdateMethod::AtOnce;
		} else if (updateMethod == "SeparateXYZ") {
			this->updateMethod = VCUpdateMethod::SeparateXYZ;
		} else if (updateMethod == "SeparateXYZ_SeparateLevelDiff") {
			this->updateMethod = VCUpdateMethod::SeparateXYZ_SeparateLevelDiff;
		}

		int tag = 100;
		blockManager.prepareForVCUpdate(this->id, tag, this->updateMethod);
		this->blockScalar3D = new BlockScalar3D<T> [blockManager.getNumBlock()];
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i size = block->getSize();
			::Vec3r origin = block->getOrigin();
			::Vec3r blockSize = block->getBlockSize();
			::Vec3r cellSize = block->getCellSize();
			T*    blockData = this->GetBlockData(block);
			const NeighborInfo* neighborInfo = block->getNeighborInfo();
			this->blockScalar3D[n].InitBoundaryCondition(
					size,
					origin,
					blockSize,
					cellSize,
					vc,
					blockData,
					boundaryType,
					boundaryValue,
					neighborInfo);
		}
	}
	~LocalScalar3D() {
		if( this->blockScalar3D ) {
			delete [] this->blockScalar3D;
			this->blockScalar3D = 0;
		}
	}

private:
	int vc;
	int id;
	VCUpdateMethod::Type updateMethod;
	BlockScalar3D<T>* blockScalar3D;
/*
	T** blockDataTable;
*/

private:
	T sum_l;
	T max_l;
	T min_l;
	T absmax_l;
	T absmin_l;

	T sum_g;
	T max_g;
	T min_g;
	T absmax_g;
	T absmin_g;

public:
	friend void LSSwap(LocalScalar3D<T>& x1, LocalScalar3D<T>& x0) {
		int id_tmp = x1.id;
		x1.id = x0.id;
		x0.id = id_tmp;
	}

public:
	int GetVC() {
		return vc;
	};

	int GetID() {
		return id;
	}

	T GetSum() {
		return sum_g;
	}

	T GetMax() {
		return max_g;
	}

	T GetMin() {
		return min_g;
	}

	T GetAbsMax() {
		return absmax_g;
	}

	T GetAbsMin() {
		return absmin_g;
	}

	T GetSumL() {
		return sum_l;
	}

	T GetMaxL() {
		return max_l;
	}

	T GetMinL() {
		return min_l;
	}

	T GetAbsMaxL() {
		return absmax_l;
	}

	T GetAbsMinL() {
		return absmin_l;
	}


	T* GetBlockData(BlockBase* block) {
		Scalar3D<T>* sdata = dynamic_cast<Scalar3D<T>*>(block->getDataClass(this->id));
		return sdata->getData();
	}

/*
	T* GetBlockData(int n) {
		return this->blockDataTable[n];
	}
*/

	void UpdateVirtualCells(
										BlockManager& blockManager) {
		switch(this->updateMethod) {
			case VCUpdateMethod::AtOnce:
				blockManager.beginUpdateVC(this->id);
				blockManager.endUpdateVC(this->id);
				break;
			case VCUpdateMethod::SeparateXYZ:
				blockManager.beginUpdateVC_X(this->id);
				blockManager.endUpdateVC_X(this->id);
				blockManager.beginUpdateVC_Y(this->id);
				blockManager.endUpdateVC_Y(this->id);
				blockManager.beginUpdateVC_Z(this->id);
				blockManager.endUpdateVC_Z(this->id);
				break;
			case VCUpdateMethod::SeparateXYZ_SeparateLevelDiff:
				blockManager.updateVC_X_F2C(this->id);
				blockManager.updateVC_Y_F2C(this->id);
				blockManager.updateVC_Z_F2C(this->id);
				blockManager.updateVC_X_Flat(this->id);
				blockManager.updateVC_Y_Flat(this->id);
				blockManager.updateVC_Z_Flat(this->id);
				blockManager.updateVC_X_C2F(this->id);
				blockManager.updateVC_Y_C2F(this->id);
				blockManager.updateVC_Z_C2F(this->id);
				break;
		}
	}

	void ImposeBoundaryCondition(
										BlockManager& blockManager) {
		switch(this->updateMethod) {
			case VCUpdateMethod::AtOnce:
				blockManager.beginUpdateVC(this->id);
				blockManager.endUpdateVC(this->id);
				ImposeBlockBoundaryCondition(blockManager);
				break;
			case VCUpdateMethod::SeparateXYZ:
				blockManager.beginUpdateVC_X(this->id);
				blockManager.endUpdateVC_X(this->id);
				blockManager.beginUpdateVC_Y(this->id);
				blockManager.endUpdateVC_Y(this->id);
				blockManager.beginUpdateVC_Z(this->id);
				blockManager.endUpdateVC_Z(this->id);
				ImposeBlockBoundaryCondition(blockManager);
				break;
			case VCUpdateMethod::SeparateXYZ_SeparateLevelDiff:
				blockManager.updateVC_X_F2C(this->id);
				blockManager.updateVC_Y_F2C(this->id);
				blockManager.updateVC_Z_F2C(this->id);
				blockManager.updateVC_X_Flat(this->id);
				blockManager.updateVC_Y_Flat(this->id);
				blockManager.updateVC_Z_Flat(this->id);
				blockManager.updateVC_X_C2F(this->id);
				blockManager.updateVC_Y_C2F(this->id);
				blockManager.updateVC_Z_C2F(this->id);
				ImposeBlockBoundaryCondition(blockManager);
				break;
		}
	}

	void ImposeBoundaryCondition(
										BlockManager& blockManager,
										LocalScalar3D<T>* plsAp,
										LocalScalar3D<T>* plsAw,
										LocalScalar3D<T>* plsAe,
										LocalScalar3D<T>* plsAs,
										LocalScalar3D<T>* plsAn,
										LocalScalar3D<T>* plsAb,
										LocalScalar3D<T>* plsAt,
										LocalScalar3D<T>* plsb) {
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* b  = plsb ->GetBlockData(block);
			this->blockScalar3D[n].ImposeBoundaryCondition(Ap, Aw, Ae, As, An, Ab, At, b);
		}
	}

	void ResetBoundaryConditionValue(
									BlockManager& blockManager,
									int i_face,
									T boundaryValue) {
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			this->blockScalar3D[n].ResetBoundaryConditionValue(i_face, boundaryValue);
		}
	}

private:
	void ImposeBlockBoundaryCondition(BlockManager& blockManager) {
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			this->blockScalar3D[n].ImposeBoundaryCondition();
		}
	}

public:
	void WriteDataInVTKFormat(
						const char* dataname,
						int step,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Config& conf) {
		VtkWriter writer;
		writer.writeScalar<T>(
						this->id,
						this->vc,
						string(dataname),
						step,
						difflevel,
						rootGrid,
						tree,
						partition,
						conf.origin,
						conf.rootLength);
	}

/*
	void WriteDataInSILOFormat(const char* dataname, int step) {
		std::ostringstream ossFileName;
		ossFileName << "data-";
		ossFileName << dataname;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".silo";

		std::string outputFile(ossFileName.str().c_str());
		SiloWriter writer(outputFile, "mesh", false);
		writer.writeDomain("block_mesh", "domain");
		writer.writeScalar<T>(this->id, "result");
	}
*/

	void Fill(BlockManager& blockManager, T value, T deviation=0.0) {
#ifdef _BLOCK_IS_LARGE_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i size = block->getSize();
			T*      blockData = this->GetBlockData(block);
			int ix = size.x + 2*vc;
			int jx = size.y + 2*vc;
			int kx = size.z + 2*vc;
			for(int k=0; k<kx; k++) {
				for(int j=0; j<jx; j++) {
					for(int i=0; i<ix; i++) {
						int mp = i + ix*(j + jx*k);
//						blockData[mp] = RandomNormal(value, deviation);
						blockData[mp] = value;
					}
				}
			}
		}
	}

	void CalcStats(BlockManager& BlockManager) {
	}

	void Dump(BlockManager& blockManager, const int step, const char* label) {
	}

	void Load(BlockManager& blockManager, const int step, const char* label) {
	}

	void Dump2(BlockManager& blockManager, const int step, const char* label) {
	}

	void Load2(BlockManager& blockManager, const int step, const char* label) {
	}

	void Dump3(BlockManager& blockManager, const int step, const char* label, Partition* partition, int myrank) {
	}

	void Load3(BlockManager& blockManager, const int step, const char* label, Partition* partition, int myrank) {
	}

};

template <>
void LocalScalar3D<real>::CalcStats(BlockManager& blockManager);

template <>
void LocalScalar3D<real>::Dump(BlockManager& blockManager, const int step, const char* label);

template <>
void LocalScalar3D<real>::Load(BlockManager& blockManager, const int step, const char* label);

template <>
void LocalScalar3D<real>::Dump2(BlockManager& blockManager, const int step, const char* label);

template <>
void LocalScalar3D<real>::Load2(BlockManager& blockManager, const int step, const char* label);

template <>
void LocalScalar3D<real>::Dump3(BlockManager& blockManager, const int step, const char* label, Partition* partition, int myrank);

template <>
void LocalScalar3D<real>::Load3(BlockManager& blockManager, const int step, const char* label, Partition* partition, int myrank);


#endif

