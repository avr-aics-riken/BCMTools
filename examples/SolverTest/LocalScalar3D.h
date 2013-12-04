#ifndef LOCALSCALAR3D_H
#define LOCALSCALAR3D_H

#include "BlockManager.h"

#include "BlockScalar3D.h"
#include "Scalar3DUpdater.h"
#include "Scalar3DUpdater2.h"
#include "Config.h"
#include "real.h"

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
		if( updaterType == 2 ) {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater2<T> >(this->vc);
		} else {
			this->id = blockManager.setDataClass<Scalar3D<T>, Scalar3DUpdater<T> >(this->vc);
		}

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
#ifdef _LARGE_BLOCK_
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

private:
	T sum;
	T max;
	T min;
	T absmax;
	T absmin;

public:
/*
	friend void LSSwap(LocalScalar3D<T>& x1, LocalScalar3D<T>& x0) {
		int vc_tmp = x1.vc;
		x1.vc = x0.vc;
		x0.vc = vc_tmp;

		int id_tmp = x1.id;
		x1.id = x0.id;
		x0.id = id_tmp;

		bool separateVCUpdate_tmp = x1.separateVCUpdate;
		x1.separateVCUpdate = x0.separateVCUpdate;
		x0.separateVCUpdate = separateVCUpdate_tmp;

		BlockScalar3D<T>* blockScalar3D_tmp = x1.blockScalar3D;
		x1.blockScalar3D = x0.blockScalar3D;
		x0.blockScalar3D = blockScalar3D_tmp;
	}
*/

public:
	int GetVC() {
		return vc;
	};

	int GetID() {
		return id;
	}

	T GetSum() {
		return sum;
	}

	T GetMax() {
		return max;
	}

	T GetMin() {
		return min;
	}

	T GetAbsMax() {
		return absmax;
	}

	T GetAbsMin() {
		return absmin;
	}

	T* GetBlockData(BlockBase* block) {
		Scalar3D<T>* sdata = dynamic_cast<Scalar3D<T>*>(block->getDataClass(this->id));
		return sdata->getData();
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
				blockManager.updateVC_X_Flat(this->id);
				blockManager.updateVC_Y_Flat(this->id);
				blockManager.updateVC_Z_Flat(this->id);
				blockManager.updateVC_X_F2C(this->id);
				blockManager.updateVC_Y_F2C(this->id);
				blockManager.updateVC_Z_F2C(this->id);
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
#ifdef _LARGE_BLOCK_
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
#ifdef _LARGE_BLOCK_
#else
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			this->blockScalar3D[n].ResetBoundaryConditionValue(i_face, boundaryValue);
		}
	}

private:
	void ImposeBlockBoundaryCondition(BlockManager& blockManager) {
#ifdef _LARGE_BLOCK_
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

	void Fill(BlockManager& blockManager, T value) {
#ifdef _LARGE_BLOCK_
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

