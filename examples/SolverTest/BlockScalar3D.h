/*
 * BCMTools
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef BLOCKSCALAR3D_H
#define BLOCKSCALAR3D_H

#include "real.h"
#include "BlockManager.h"
#include "Scalar3D.h"

#include "bcx.h"
#include "bca.h"
#include "bcax.h"

template <typename T>
class BlockScalar3D {
public:
  enum BoundaryType {
    DIRICHLET,  ///< Dirichlet
    NEUMANN,    ///< Neumann
    PERIODIC,   ///< ¼þ´ü¶­³¦
    INNER,      ///< ÆâÉô¶­³¦
		POISEUILLE_U = 100,
		POISEUILLE_P = 101,
  };

public:
	BlockScalar3D() {
	}
	~BlockScalar3D() {
	}

private:
	int  size[3];
	double origin[3];
	double blockSize[3];
	double cellSize[3];
	int  vc;
	T*   blockData;
	int* blockBoundaryType;
	T*   blockBoundaryValue;

public:
	void ResetBoundaryConditionValue(
					int i_face,
					T boundaryValue) {
		Face face = Face(i_face);
		if( this->blockBoundaryType[face] != INNER ) {
			this->blockBoundaryValue[face] = boundaryValue;
		}
	}

	void InitBoundaryCondition(
					::Vec3i size,
					::Vec3d origin,
					::Vec3d blockSize,
					::Vec3d cellSize,
					int vc,
					T* blockData,
					int* boundaryType,
					T* boundaryValue,
					const NeighborInfo* neighborInfo) {
		this->size[0]       = size.x;
		this->size[1]       = size.y;
		this->size[2]       = size.z;
		this->origin[0]     = origin.x;
		this->origin[1]     = origin.y;
		this->origin[2]     = origin.z;
		this->blockSize[0]  = blockSize.x;
		this->blockSize[1]  = blockSize.y;
		this->blockSize[2]  = blockSize.z;
		this->cellSize[0]   = cellSize.x;
		this->cellSize[1]   = cellSize.y;
		this->cellSize[2]   = cellSize.z;
		this->vc                 = vc;
		this->blockData          = blockData;
		this->blockBoundaryType  = new int [NUM_FACE];
		this->blockBoundaryValue = new T   [NUM_FACE];

		for (int i=0; i<NUM_FACE; i++) {
			Face face = Face(i);
			if( neighborInfo[face].isOuterBoundary() ) {
				this->blockBoundaryType[face]  = boundaryType[face];
				this->blockBoundaryValue[face] = boundaryValue[face];
				if( this->blockBoundaryType[face] == NEUMANN ) {
					this->blockBoundaryValue[face] = boundaryValue[face]*this->cellSize[0];
				}
			} else {
				this->blockBoundaryType[face]  = INNER;
				this->blockBoundaryValue[face] = 0.0;
			}
		}

		if( this->blockBoundaryType[X_M] == DIRICHLET ) {
			this->pfImposeBlockBoundaryCondition_X_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_X_M_D;
			this->pfImposeBlockBoundaryCondition_Aw  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Aw_D;
		} else if( this->blockBoundaryType[X_M] == NEUMANN ) {
			this->pfImposeBlockBoundaryCondition_X_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_X_M_N;
			this->pfImposeBlockBoundaryCondition_Aw  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Aw_N;
		} else if( this->blockBoundaryType[X_M] == PERIODIC ) {
			this->pfImposeBlockBoundaryCondition_X_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_X_M_P;
			this->pfImposeBlockBoundaryCondition_Aw  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Aw_P;
		} else if( this->blockBoundaryType[X_M] == POISEUILLE_U ) {
			this->pfImposeBlockBoundaryCondition_X_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_X_M_POISEUILLE_U;
			this->pfImposeBlockBoundaryCondition_Aw  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Aw_POISEUILLE_U;
		} else if( this->blockBoundaryType[X_M] == POISEUILLE_P ) {
			this->pfImposeBlockBoundaryCondition_X_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_X_M_POISEUILLE_P;
			this->pfImposeBlockBoundaryCondition_Aw  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Aw_POISEUILLE_P;
		} else {
			this->pfImposeBlockBoundaryCondition_X_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
			this->pfImposeBlockBoundaryCondition_Aw  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
		}

		if( this->blockBoundaryType[X_P] == DIRICHLET ) {
			this->pfImposeBlockBoundaryCondition_X_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_X_P_D;
			this->pfImposeBlockBoundaryCondition_Ae  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Ae_D;
		} else if( this->blockBoundaryType[X_P] == NEUMANN ) {
			this->pfImposeBlockBoundaryCondition_X_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_X_P_N;
			this->pfImposeBlockBoundaryCondition_Ae  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Ae_N;
		} else if( this->blockBoundaryType[X_P] == PERIODIC ) {
			this->pfImposeBlockBoundaryCondition_X_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_X_P_P;
			this->pfImposeBlockBoundaryCondition_Ae  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Ae_P;
		} else {
			this->pfImposeBlockBoundaryCondition_X_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
			this->pfImposeBlockBoundaryCondition_Ae  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
		}

		if( this->blockBoundaryType[Y_M] == DIRICHLET ) {
			this->pfImposeBlockBoundaryCondition_Y_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Y_M_D;
			this->pfImposeBlockBoundaryCondition_As  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_As_D;
		} else if( this->blockBoundaryType[Y_M] == NEUMANN ) {
			this->pfImposeBlockBoundaryCondition_Y_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Y_M_N;
			this->pfImposeBlockBoundaryCondition_As  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_As_N;
		} else if( this->blockBoundaryType[Y_M] == PERIODIC ) {
			this->pfImposeBlockBoundaryCondition_Y_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Y_M_P;
			this->pfImposeBlockBoundaryCondition_As  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_As_P;
		} else {
			this->pfImposeBlockBoundaryCondition_Y_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
			this->pfImposeBlockBoundaryCondition_As  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
		}

		if( this->blockBoundaryType[Y_P] == DIRICHLET ) {
			this->pfImposeBlockBoundaryCondition_Y_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Y_P_D;
			this->pfImposeBlockBoundaryCondition_An  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_An_D;
		} else if( this->blockBoundaryType[Y_P] == NEUMANN ) {
			this->pfImposeBlockBoundaryCondition_Y_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Y_P_N;
			this->pfImposeBlockBoundaryCondition_An  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_An_N;
		} else if( this->blockBoundaryType[Y_P] == PERIODIC ) {
			this->pfImposeBlockBoundaryCondition_Y_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Y_P_P;
			this->pfImposeBlockBoundaryCondition_An  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_An_P;
		} else {
			this->pfImposeBlockBoundaryCondition_Y_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
			this->pfImposeBlockBoundaryCondition_An  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
		}

		if( this->blockBoundaryType[Z_M] == DIRICHLET ) {
			this->pfImposeBlockBoundaryCondition_Z_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Z_M_D;
			this->pfImposeBlockBoundaryCondition_Ab  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Ab_D;
		} else if( this->blockBoundaryType[Z_M] == NEUMANN ) {
			this->pfImposeBlockBoundaryCondition_Z_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Z_M_N;
			this->pfImposeBlockBoundaryCondition_Ab  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Ab_N;
		} else if( this->blockBoundaryType[Z_M] == PERIODIC ) {
			this->pfImposeBlockBoundaryCondition_Z_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Z_M_P;
			this->pfImposeBlockBoundaryCondition_Ab  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Ab_P;
		} else {
			this->pfImposeBlockBoundaryCondition_Z_M = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
			this->pfImposeBlockBoundaryCondition_Ab  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
		}

		if( this->blockBoundaryType[Z_P] == DIRICHLET ) {
			this->pfImposeBlockBoundaryCondition_Z_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Z_P_D;
			this->pfImposeBlockBoundaryCondition_At  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_At_D;
		} else if( this->blockBoundaryType[Z_P] == NEUMANN ) {
			this->pfImposeBlockBoundaryCondition_Z_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Z_P_N;
			this->pfImposeBlockBoundaryCondition_At  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_At_N;
		} else if( this->blockBoundaryType[Z_P] == PERIODIC ) {
			this->pfImposeBlockBoundaryCondition_Z_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Z_P_P;
			this->pfImposeBlockBoundaryCondition_At  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_At_P;
		} else {
			this->pfImposeBlockBoundaryCondition_Z_P = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
			this->pfImposeBlockBoundaryCondition_At  = &BlockScalar3D<T>::ImposeBlockBoundaryCondition_Dummy;
		}
	}

	void ImposeBoundaryCondition() {
		(this->*pfImposeBlockBoundaryCondition_X_M)();
		(this->*pfImposeBlockBoundaryCondition_X_P)();
		(this->*pfImposeBlockBoundaryCondition_Y_M)();
		(this->*pfImposeBlockBoundaryCondition_Y_P)();
		(this->*pfImposeBlockBoundaryCondition_Z_M)();
		(this->*pfImposeBlockBoundaryCondition_Z_P)();
	}

	void ImposeBoundaryCondition(real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At, real* b) {
		(this->*pfImposeBlockBoundaryCondition_Aw)(Ap, Aw, Ae, b);
		(this->*pfImposeBlockBoundaryCondition_Ae)(Ap, Aw, Ae, b);
		(this->*pfImposeBlockBoundaryCondition_As)(Ap, As, An, b);
		(this->*pfImposeBlockBoundaryCondition_An)(Ap, As, An, b);
		(this->*pfImposeBlockBoundaryCondition_Ab)(Ap, Ab, At, b);
		(this->*pfImposeBlockBoundaryCondition_At)(Ap, Ab, At, b);
	}

private:
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_X_M)();
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_X_P)();
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_Y_M)();
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_Y_P)();
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_Z_M)();
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_Z_P)();
	void ImposeBlockBoundaryCondition_Dummy() {
	}
	void ImposeBlockBoundaryCondition_X_M_D() {
		bc_x3_d_(this->blockData, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_X_M_N() {
		bc_x3_n_(this->blockData, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_X_M_P() {
//		bc_x3_p_(this->blockData, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_X_P_D() {
		bc_x1_d_(this->blockData, &(this->blockBoundaryValue[X_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_X_P_N() {
		bc_x1_n_(this->blockData, &(this->blockBoundaryValue[X_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_X_P_P() {
		bc_x1_p_(this->blockData, &(this->blockBoundaryValue[X_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Y_M_D() {
		bc_x4_d_(this->blockData, &(this->blockBoundaryValue[Y_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Y_M_N() {
		bc_x4_n_(this->blockData, &(this->blockBoundaryValue[Y_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Y_M_P() {
//		bc_x4_p_(this->blockData, &(this->blockBoundaryValue[Y_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Y_P_D() {
		bc_x2_d_(this->blockData, &(this->blockBoundaryValue[Y_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Y_P_N() {
		bc_x2_n_(this->blockData, &(this->blockBoundaryValue[Y_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Y_P_P() {
		bc_x2_p_(this->blockData, &(this->blockBoundaryValue[Y_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Z_M_D() {
		bc_x6_d_(this->blockData, &(this->blockBoundaryValue[Z_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Z_M_N() {
		bc_x6_n_(this->blockData, &(this->blockBoundaryValue[Z_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Z_M_P() {
//		bc_x6_p_(this->blockData, &(this->blockBoundaryValue[Z_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Z_P_D() {
		bc_x5_d_(this->blockData, &(this->blockBoundaryValue[Z_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Z_P_N() {
		bc_x5_n_(this->blockData, &(this->blockBoundaryValue[Z_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Z_P_P() {
		bc_x5_p_(this->blockData, &(this->blockBoundaryValue[Z_P]), this->size, (int*)&(this->vc));
	}

	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_Aw)(real* Ap, real* Aw, real* Ae, real* b);
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_Ae)(real* Ap, real* Aw, real* Ae, real* b);
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_As)(real* Ap, real* As, real* An, real* b);
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_An)(real* Ap, real* As, real* An, real* b);
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_Ab)(real* Ap, real* Ab, real* At, real* b);
	void (BlockScalar3D<T>::*pfImposeBlockBoundaryCondition_At)(real* Ap, real* Ab, real* At, real* b);
	void ImposeBlockBoundaryCondition_Dummy(real* Ap, real* Aw, real* Ae, real* b) {
	}
	void ImposeBlockBoundaryCondition_Aw_D(real* Ap, real* Aw, real* Ae, real* b) {
		bc_aw_d_(Ap, Aw, b, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Aw_N(real* Ap, real* Aw, real* Ae, real* b) {
		bc_aw_n_(Ap, Aw, b, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Aw_P(real* Ap, real* Aw, real* Ae, real* b) {
	}

	void ImposeBlockBoundaryCondition_Ae_D(real* Ap, real* Aw, real* Ae, real* b) {
		bc_ae_d_(Ap, Ae, b, &(this->blockBoundaryValue[X_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Ae_N(real* Ap, real* Aw, real* Ae, real* b) {
		bc_ae_n_(Ap, Ae, b, &(this->blockBoundaryValue[X_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Ae_P(real* Ap, real* Aw, real* Ae, real* b) {
	}
	void ImposeBlockBoundaryCondition_As_D(real* Ap, real* As, real* An, real* b) {
		bc_as_d_(Ap, As, b, &(this->blockBoundaryValue[Y_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_As_N(real* Ap, real* As, real* An, real* b) {
		bc_as_n_(Ap, As, b, &(this->blockBoundaryValue[Y_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_As_P(real* Ap, real* As, real* An, real* b) {
	}
	void ImposeBlockBoundaryCondition_An_D(real* Ap, real* As, real* An, real* b) {
		bc_an_d_(Ap, An, b, &(this->blockBoundaryValue[Y_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_An_N(real* Ap, real* As, real* An, real* b) {
		bc_an_n_(Ap, An, b, &(this->blockBoundaryValue[Y_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_An_P(real* Ap, real* As, real* An, real* b) {
	}
	void ImposeBlockBoundaryCondition_Ab_D(real* Ap, real* Ab, real* At, real* b) {
		bc_ab_d_(Ap, Ab, b, &(this->blockBoundaryValue[Z_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Ab_N(real* Ap, real* Ab, real* At, real* b) {
		bc_ab_n_(Ap, Ab, b, &(this->blockBoundaryValue[Z_M]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_Ab_P(real* Ap, real* Ab, real* At, real* b) {
	}
	void ImposeBlockBoundaryCondition_At_D(real* Ap, real* Ab, real* At, real* b) {
		bc_at_d_(Ap, At, b, &(this->blockBoundaryValue[Z_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_At_N(real* Ap, real* Ab, real* At, real* b) {
		bc_at_n_(Ap, At, b, &(this->blockBoundaryValue[Z_P]), this->size, (int*)&(this->vc));
	}
	void ImposeBlockBoundaryCondition_At_P(real* Ap, real* Ab, real* At, real* b) {
	}


	void ImposeBlockBoundaryCondition_X_M_POISEUILLE_U() {
//		bc_x3_poiseuille_u_(this->blockData, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc), this->origin, this->blockSize, this->cellSize);
	}
	void ImposeBlockBoundaryCondition_Aw_POISEUILLE_U(real* Ap, real* Aw, real* Ae, real* b) {
//		bc_aw_poiseuille_u_(Ap, Aw, b, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc), this->origin, this->blockSize, this->cellSize);
	}
	void ImposeBlockBoundaryCondition_X_M_POISEUILLE_P() {
//		bc_x3_poiseuille_p_(this->blockData, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc), this->origin, this->blockSize, this->cellSize);
	}
	void ImposeBlockBoundaryCondition_Aw_POISEUILLE_P(real* Ap, real* Aw, real* Ae, real* b) {
//		bc_aw_poiseuille_p_(Ap, Aw, b, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc), this->origin, this->blockSize, this->cellSize);
	}
};

template <>
void BlockScalar3D<real>::ImposeBlockBoundaryCondition_X_M_POISEUILLE_U();
template <>
void BlockScalar3D<real>::ImposeBlockBoundaryCondition_Aw_POISEUILLE_U(real* Ap, real* Aw, real* Ae, real* b);
template <>
void BlockScalar3D<real>::ImposeBlockBoundaryCondition_X_M_POISEUILLE_P();
template <>
void BlockScalar3D<real>::ImposeBlockBoundaryCondition_Aw_POISEUILLE_P(real* Ap, real* Aw, real* Ae, real* b);


template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Dummy();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_M_D();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_P_D();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_M_D();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_P_D();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_M_D();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_P_D();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_M_N();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_P_N();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_M_N();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_P_N();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_M_N();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_P_N();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_M_P();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_P_P();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_M_P();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_P_P();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_M_P();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_P_P();

template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Dummy(real* Ap, real* Aw, real* Ae, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Aw_D(real* Ap, real* Aw, real* Ae, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ae_D(real* Ap, real* Aw, real* Ae, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_As_D(real* Ap, real* As, real* An, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_An_D(real* Ap, real* As, real* An, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ab_D(real* Ap, real* Ab, real* At, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_At_D(real* Ap, real* Ab, real* At, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Aw_N(real* Ap, real* Aw, real* Ae, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ae_N(real* Ap, real* Aw, real* Ae, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_As_N(real* Ap, real* As, real* An, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_An_N(real* Ap, real* As, real* An, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ab_N(real* Ap, real* Ab, real* At, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_At_N(real* Ap, real* Ab, real* At, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Aw_P(real* Ap, real* Aw, real* Ae, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ae_P(real* Ap, real* Aw, real* Ae, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_As_P(real* Ap, real* As, real* An, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_An_P(real* Ap, real* As, real* An, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ab_P(real* Ap, real* Ab, real* At, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_At_P(real* Ap, real* Ab, real* At, real* b);

template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_M_POISEUILLE_U();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Aw_POISEUILLE_U(real* Ap, real* Aw, real* Ae, real* b);
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_M_POISEUILLE_P();
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Aw_POISEUILLE_P(real* Ap, real* Aw, real* Ae, real* b);

#endif

