/*
 * BCMTools
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef ILS3D_H
#define ILS3D_H

#include <cmath>
#include <cfloat>

#include "BlockManager.h"
#include "LocalScalar3D.h"

#include "real.h"
#include "blas.h"
#include "comm.h"
#include "PM.h"

class ILS3D {
public:
	ILS3D() {
	}

	~ILS3D() {
	}

private:
	bool IsConverged(
						BlockManager& blockManager,
						real& residual,
						real& rr,
						real& bb,
						real& epsilon,
						int count,
						int countMax) {
		residual = (fabs(bb)>FLT_MIN) ? fabs(rr)/fabs(bb) : fabs(rr);
		if( residual <= epsilon*epsilon ) {
			return true;
		}
		if( isnan(residual) ) {
			std::cout << bb << " " << rr << std::endl;
			Exit(0);
		}
		if( isinf(residual) ) {
			Exit(0);
		}
		if( count == countMax ) {
		}
		return false;
	}

	void Jacobi_Smoother(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsb,
						LocalScalar3D<real>* plsx0,
						real omega) {
		int vc = plsx->GetVC();
PM_Start(tm_JacobiSmoother, 0, 0, true);

		int NB = blockManager.getNumBlock();
		BlockBase* block0 = blockManager.getBlock(0);
		::Vec3i size = block0->getSize();
		int NX = size.x;
		int NY = size.y;
		int NZ = size.z;
		int sz[3] = {NX, NY, NZ};

PM_Start(tm_JacobiSmoother_Calc, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<NB; ++n) {
			BlockBase* block = blockManager.getBlock(n);
			real* x  = plsx ->GetBlockData(block);
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* b  = plsb ->GetBlockData(block);
			real* x0 = plsx0->GetBlockData(block);
			real pomega = omega;

			jacobi_smoother_(
							x0, x,
							Ap, Aw, Ae, As, An, Ab, At,
							b,
							&pomega,
							sz, &vc);
		}
PM_Stop(tm_JacobiSmoother_Calc, 0, 0, 16.0*NX*NY*NZ, NB);

PM_Start(tm_JacobiSmoother_Swap, 0, 0, true);
		LSSwap(*plsx, *plsx0);
PM_Stop(tm_JacobiSmoother_Swap);

PM_Start(tm_JacobiSmoother_Comm, 0, 0, true);
		plsx->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_JacobiSmoother_Comm);

PM_Stop(tm_JacobiSmoother);
	}

	void CalcAx(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsAx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsx) {
		int vc = plsx->GetVC();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* Ax = plsAx->GetBlockData(block);
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* x  = plsx ->GetBlockData(block);

			calc_ax_(
							Ax,
							Ap, Aw, Ae, As, An, Ab, At,
							x,
							sz, &vc);
		}
	}

	void CalcR(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsr,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsb) {
		int vc = plsx->GetVC();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* r  = plsr ->GetBlockData(block);
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* x  = plsx ->GetBlockData(block);
			real* b  = plsb ->GetBlockData(block);

			calc_r_(
							r,
							Ap, Aw, Ae, As, An, Ab, At,
							x,
							b,
							sz, &vc);
		}
	}

	void CalcR2(
						BlockManager& blockManager,
						real& rr,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsb) {
		int vc = plsx->GetVC();
		double rr_local = 0.0;
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+: rr_local)
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* x  = plsx ->GetBlockData(block);
			real* b  = plsb ->GetBlockData(block);

			real rr_block = 0.0;
			calc_r2_(
							&rr_block,
							Ap, Aw, Ae, As, An, Ab, At,
							x,
							b,
							sz, &vc);

			rr_local += rr_block;
//			std::cout << rr_block << std::endl;
		}

		const MPI::Intracomm& comm = blockManager.getCommunicator();

		double rr_global = 0.0;
		allreduce_(&rr_global, &rr_local);
/*
#ifdef _REAL_IS_DOUBLE_		
		comm.Allreduce(&rr_local, &rr_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM);
#else
		comm.Allreduce(&rr_local, &rr_global, 1, MPI_FLOAT, MPI_SUM);
#endif
*/

		rr = rr_global;
	}

	void DOT(
						BlockManager& blockManager,
						real& xy,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsy) {
		int vc = plsx->GetVC();
		double xy_local = 0.0;
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+: xy_local)
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* x = plsx->GetBlockData(block);
			real* y = plsy->GetBlockData(block);

			real xy_block = 0.0;
			dot_(&xy_block, x, y, sz, &vc);

			xy_local += xy_block;
		}

		const MPI::Intracomm& comm = blockManager.getCommunicator();

		double xy_global = 0.0;
		allreduce_(&xy_global, &xy_local);
/*
#ifdef _REAL_IS_DOUBLE_		
		comm.Allreduce(&xy_local, &xy_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM);
#else
		comm.Allreduce(&xy_local, &xy_global, 1, MPI_FLOAT, MPI_SUM);
#endif
*/

		xy = xy_global;
	}

	void AXPY(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsy,
						LocalScalar3D<real>* plsx,
						real a) {
		int vc = plsx->GetVC();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* y = plsy->GetBlockData(block);
			real* x = plsx->GetBlockData(block);

			axpy_(y, x, &a, sz, &vc);
		}
	}
	void XPAY(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsy,
						LocalScalar3D<real>* plsx,
						real a) {
		int vc = plsx->GetVC();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* y = plsy->GetBlockData(block);
			real* x = plsx->GetBlockData(block);

			xpay_(y, x, &a, sz, &vc);
		}
	}
	void AXPYZ(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsz,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsy,
						real a) {
		int vc = plsx->GetVC();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* z = plsz->GetBlockData(block);
			real* x = plsx->GetBlockData(block);
			real* y = plsy->GetBlockData(block);

			axpyz_(z, x, y, &a, sz, &vc);
		}
	}
	void AXPBYPZ(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsz,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsy,
						real a,
						real b) {
		int vc = plsx->GetVC();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* z = plsz->GetBlockData(block);
			real* x = plsx->GetBlockData(block);
			real* y = plsy->GetBlockData(block);

			axpbypz_(z, x, y, &a, &b, sz, &vc);
		}
	}

public:
	void Fill(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						real value) {
		int vc = plsx->GetVC();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* x = plsx->GetBlockData(block);

			real v = value;
			fill_(x, &v, sz, &vc);
		}
	}
	
	void Copy(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsy,
						LocalScalar3D<real>* plsx) {
		int vc = plsx->GetVC();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* x = plsx->GetBlockData(block);
			real* y = plsy->GetBlockData(block);

			copy_(y, x, sz, &vc);
		}
	}

	void Jacobi_PreConditioner(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsb,
						LocalScalar3D<real>* plsx0,
						real omega,
						int countPreConditioner ) {
		for(int count=0; count<countPreConditioner; count++) {
			Jacobi_Smoother(
							blockManager,
							plsx,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsb,
							plsx0,
							omega);
		}
	}


	void Jacobi(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsb,
						LocalScalar3D<real>* plsx0,
						real omega,
						int countMax,
						real epsilon,
						int& count,
						real& residual) {
		real bb = 0.0;
		DOT(blockManager, bb, plsb, plsb);

		for(count=1; count<=countMax; ++count) {
			Jacobi_Smoother(
							blockManager,
							plsx, 
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsb,
							plsx0,
							omega);

			real rr = 0.0;
			CalcR2(
							blockManager,
							rr,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsx,
							plsb);

			bool bResult = IsConverged(
							blockManager,
							residual,
							rr,
							bb,
							epsilon,
							count,
							countMax);
			if( bResult == true ) {
				break;
			}

		}
		plsx->ImposeBoundaryCondition(blockManager);
	}

	void CG(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsb,
						LocalScalar3D<real>* plsr,
						LocalScalar3D<real>* plsp,
						LocalScalar3D<real>* plsq,
						LocalScalar3D<real>* plsz,
						LocalScalar3D<real>* plsx0,
						real omega,
						int countPreConditioner,
						int countMax,
						real epsilon,
						int& count,
						real& residual) {
		real bb = 0.0;
		DOT(blockManager, bb, plsb, plsb);

		if( fabs(bb) < FLT_MIN ) {
			residual = 0.0;
			count = 0;
			return;
		}

		CalcR(
						blockManager,
						plsr,
						plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
						plsx,
						plsb);
		plsr->ImposeBoundaryCondition(blockManager);

		real rr0 = 1.0;
		real rr1 = 0.0;
		for(count=1; count<=countMax; ++count) {
			Fill(blockManager, plsz, 0.0);
			Jacobi_PreConditioner(
							blockManager,
							plsz,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsr,
							plsx0,
							omega,
							countPreConditioner);

			rr1 = 0.0;
			DOT(
							blockManager,
							rr1,
							plsr,
							plsz);

			if( fabs(rr1) < FLT_MIN ) {
				residual = rr1;
				count = 0;
				break;
			}

			real beta = rr1/rr0;

			if( count==1 ) {
				Copy(
							blockManager,
							plsp,
							plsz);
			} else {
				XPAY(
							blockManager,
							plsp,
							plsz,
							beta);
			}
			plsp->ImposeBoundaryCondition(blockManager);

			CalcAx(
							blockManager,
							plsq,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsp);

			real qp = 0.0;
			DOT(
							blockManager,
							qp,
							plsq,
							plsp);

			real alpha = rr1/qp;

			AXPY(
							blockManager,
							plsx,
							plsp,
							alpha);
			AXPY(
							blockManager,
							plsr,
							plsq,
							-alpha);
			plsr->ImposeBoundaryCondition(blockManager);

			real rr = 0.0;
			DOT(
							blockManager,
							rr,
							plsr,
							plsr);

			bool bResult = IsConverged(
							blockManager,
							residual,
							rr,
							bb,
							epsilon,
							count,
							countMax);
			if( bResult == true ) {
				break;
			}

			rr0 = rr1;
		}
		plsx->ImposeBoundaryCondition(blockManager);
	}

	void BiCGSTAB(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsb,
						LocalScalar3D<real>* plsr,
						LocalScalar3D<real>* plsr0,
						LocalScalar3D<real>* plsp,
						LocalScalar3D<real>* plsp_,
						LocalScalar3D<real>* plsq_,
						LocalScalar3D<real>* plss,
						LocalScalar3D<real>* plss_,
						LocalScalar3D<real>* plst_,
						LocalScalar3D<real>* plsx0,
						real omega,
						int countPreConditioner,
						int countMax,
						real epsilon,
						int& count,
						real& residual) {
		real bb = 0.0;
		DOT(blockManager, bb, plsb, plsb);

		if( fabs(bb) < FLT_MIN ) {
			residual = 0.0;
			count = 0;
			return;
		}

		CalcR(
						blockManager,
						plsr,
						plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
						plsx,
						plsb);

		Copy(
						blockManager,
						plsr0,
						plsr);

		real rr0 = 1.0;
		real alpha = 0.0;
		real gamma = 1.0;
		for(count=1; count<=countMax; ++count) {
			real rr1 = 0.0;
			DOT(blockManager, rr1, plsr, plsr0);

			if( fabs(rr1) < FLT_MIN ) {
				residual = rr1;
				count = 0;
				break;
			}

			if( count == 1 ) {
				Copy(
							blockManager,
							plsp,
							plsr);
			} else {
				real beta = rr1/rr0*alpha/gamma;
				AXPY(
							blockManager,
							plsp,
							plsq_,
							-gamma);
				XPAY(
							blockManager,
							plsp,
							plsr,
							beta);
			}
			plsp->ImposeBoundaryCondition(blockManager);

			Fill(blockManager, plsp_, 0.0);
			Jacobi_PreConditioner(
							blockManager,
							plsp_,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsp,
							plsx0,
							omega,
							countPreConditioner);

			CalcAx(
							blockManager,
							plsq_,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsp_);

			real q_r0 = 0.0;
			DOT(
							blockManager,
							q_r0,
							plsq_,
							plsr0);

			alpha = rr1/q_r0;

			AXPYZ(
							blockManager,
							plss,
							plsq_,
							plsr,
							-alpha);
			plss->ImposeBoundaryCondition(blockManager);

			Fill(blockManager, plss_, 0.0);
			Jacobi_PreConditioner(
							blockManager,
							plss_,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plss,
							plsx0,
							omega,
							countPreConditioner);

			CalcAx(
							blockManager,
							plst_,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plss_);

			real t_s = 0.0;
			DOT(
							blockManager,
							t_s,
							plst_,
							plss);

			real t_t_ = 0.0;
			DOT(
							blockManager,
							t_t_,
							plst_,
							plst_);

			gamma = t_s/t_t_;

			AXPBYPZ(
							blockManager,
							plsx,
							plsp_,
							plss_,
							alpha,
							gamma);
			AXPYZ(
							blockManager,
							plsr,
							plst_,
							plss,
							-gamma);

			real rr = 0.0;
			DOT(
							blockManager,
							rr,
							plsr,
							plsr);

			bool bResult = IsConverged(
							blockManager,
							residual,
							rr,
							bb,
							epsilon,
							count,
							countMax);
			if( bResult == true ) {
				break;
			}

			rr0 = rr1;
		}
		plsx->ImposeBoundaryCondition(blockManager);
	}

	void Jacobi_Mask(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsb,
						LocalScalar3D<real>* plsx0,
						LocalScalar3D<int>* plsMaskId,
						real omega,
						int countMax,
						real epsilon,
						int& count,
						real& residual) {
		real bb = 0.0;
		DOT(blockManager, bb, plsb, plsb);

		for(count=1; count<=countMax; ++count) {
			Jacobi_Smoother_Mask(
							blockManager,
							plsx, 
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsb,
							plsx0,
							plsMaskId,
							omega);

			real rr = 0.0;
			CalcR2(
							blockManager,
							rr,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsx,
							plsb);

			bool bResult = IsConverged(
							blockManager,
							residual,
							rr,
							bb,
							epsilon,
							count,
							countMax);
			if( bResult == true ) {
				break;
			}

		}
		plsx->ImposeBoundaryCondition(blockManager);
	}

	void BiCGSTAB_Mask(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsb,
						LocalScalar3D<real>* plsr,
						LocalScalar3D<real>* plsr0,
						LocalScalar3D<real>* plsp,
						LocalScalar3D<real>* plsp_,
						LocalScalar3D<real>* plsq_,
						LocalScalar3D<real>* plss,
						LocalScalar3D<real>* plss_,
						LocalScalar3D<real>* plst_,
						LocalScalar3D<real>* plsx0,
						LocalScalar3D<int>* plsMaskId,
						real omega,
						int countPreConditioner,
						int countMax,
						real epsilon,
						int& count,
						real& residual) {

		real bb = 0.0;
		DOT_Mask(blockManager, bb, plsb, plsb, plsMaskId);

		if( fabs(bb) < FLT_MIN ) {
			residual = 0.0;
			count = 0;
			return;
		}

		CalcR(
						blockManager,
						plsr,
						plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
						plsx,
						plsb);

		Copy(
						blockManager,
						plsr0,
						plsr);

		real rr0 = 1.0;
		real alpha = 0.0;
		real gamma = 1.0;
		for(count=1; count<=countMax; ++count) {
			real rr1 = 0.0;
			DOT_Mask(blockManager, rr1, plsr, plsr0, plsMaskId);

			if( fabs(rr1) < FLT_MIN ) {
				residual = rr1;
				count = 0;
				break;
			}

			if( count == 1 ) {
				Copy(
							blockManager,
							plsp,
							plsr);
			} else {
				real beta = rr1/rr0*alpha/gamma;
				AXPY(
							blockManager,
							plsp,
							plsq_,
							-gamma);
				XPAY(
							blockManager,
							plsp,
							plsr,
							beta);
			}
			plsp->ImposeBoundaryCondition(blockManager);

			Fill(blockManager, plsp_, 0.0);

			Jacobi_PreConditioner_Mask(
							blockManager,
							plsp_,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsp,
							plsx0,
							plsMaskId,
							omega,
							countPreConditioner);

			CalcAx_Mask(
							blockManager,
							plsq_,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsp_,
							plsMaskId);

			real q_r0 = 0.0;
			DOT_Mask(
							blockManager,
							q_r0,
							plsq_,
							plsr0,
							plsMaskId);

			alpha = rr1/q_r0;

			AXPYZ(
							blockManager,
							plss,
							plsq_,
							plsr,
							-alpha);

			plss->ImposeBoundaryCondition(blockManager);

			Fill(blockManager, plss_, 0.0);

			Jacobi_PreConditioner_Mask(
							blockManager,
							plss_,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plss,
							plsx0,
							plsMaskId,
							omega,
							countPreConditioner);

			CalcAx_Mask(
							blockManager,
							plst_,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plss_,
							plsMaskId);

			real t_s = 0.0;
			DOT_Mask(
							blockManager,
							t_s,
							plst_,
							plss,
							plsMaskId);

			real t_t_ = 0.0;
			DOT_Mask(
							blockManager,
							t_t_,
							plst_,
							plst_,
							plsMaskId);

			gamma = t_s/t_t_;

			AXPBYPZ(
							blockManager,
							plsx,
							plsp_,
							plss_,
							alpha,
							gamma);
			AXPYZ(
							blockManager,
							plsr,
							plst_,
							plss,
							-gamma);

			real rr = 0.0;
			DOT_Mask(
							blockManager,
							rr,
							plsr,
							plsr,
							plsMaskId);

			bool bResult = IsConverged(
							blockManager,
							residual,
							rr,
							bb,
							epsilon,
							count,
							countMax);
			if( bResult == true ) {
				break;
			}

			rr0 = rr1;
		}

		plsx->ImposeBoundaryCondition(blockManager);

	}

	void DOT_Mask(
						BlockManager& blockManager,
						real& xy,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsy,
						LocalScalar3D<int>* plsMaskId) {
PM_Start(tm_DOT, 0, 0, true);

		int vc = plsx->GetVC();
		double xy_local = 0.0;

		int NB = blockManager.getNumBlock();
		BlockBase* block0 = blockManager.getBlock(0);
		::Vec3i size = block0->getSize();
		int NX = size.x;
		int NY = size.y;
		int NZ = size.z;
		int sz[3] = {size.x, size.y, size.z};

PM_Start(tm_DOT_Calc, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for reduction(+: xy_local)
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			real* x = plsx->GetBlockData(block);
			real* y = plsy->GetBlockData(block);
			int* mask = plsMaskId->GetBlockData(block);

			real xy_block = 0.0;
			dot_mask_(&xy_block, x, y, mask, sz, &vc);

			xy_local += xy_block;
		}
PM_Stop(tm_DOT_Calc, 0, 0, 2.0*NX*NY*NZ, NB);

		const MPI::Intracomm& comm = blockManager.getCommunicator();
		double xy_global = 0.0;

PM_Start(tm_DOT_Comm, 0, 0, true);
		allreduce_(&xy_global, &xy_local);
PM_Stop(tm_DOT_Comm);

		xy = xy_global;

PM_Stop(tm_DOT);
	}

	void Jacobi_PreConditioner_Mask(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsb,
						LocalScalar3D<real>* plsx0,
						LocalScalar3D<int>* plsMaskId,
						real omega,
						int countPreConditioner ) {
		for(int count=0; count<countPreConditioner; count++) {
			Jacobi_Smoother_Mask(
							blockManager,
							plsx,
							plsAp, plsAw, plsAe, plsAs, plsAn, plsAb, plsAt,
							plsb,
							plsx0,
							plsMaskId,
							omega);
		}
	}

	void Jacobi_Smoother_Mask(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsb,
						LocalScalar3D<real>* plsx0,
						LocalScalar3D<int>* plsMaskId,
						real omega) {
		int vc = plsx->GetVC();
PM_Start(tm_JacobiSmoother, 0, 0, true);

		int NB = blockManager.getNumBlock();
		BlockBase* block0 = blockManager.getBlock(0);
		::Vec3i size = block0->getSize();
		int NX = size.x;
		int NY = size.y;
		int NZ = size.z;
		int sz[3] = {size.x, size.y, size.z};

PM_Start(tm_JacobiSmoother_Calc, 0, 0, true);
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			real* x  = plsx ->GetBlockData(block);
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* b  = plsb ->GetBlockData(block);
			real* x0 = plsx0->GetBlockData(block);
			int* mask = plsMaskId->GetBlockData(block);
			real pomega = omega;

			jacobi_smoother_mask_(
							x0, x,
							Ap, Aw, Ae, As, An, Ab, At,
							b,
							mask,
							&pomega,
							sz, &vc);
		}
PM_Stop(tm_JacobiSmoother_Calc, 0, 0, 16.0*NX*NY*NZ, NB);

PM_Start(tm_JacobiSmoother_Swap, 0, 0, true);
		LSSwap(*plsx, *plsx0);
PM_Stop(tm_JacobiSmoother_Swap);

PM_Start(tm_JacobiSmoother_Comm, 0, 0, true);
		plsx->ImposeBoundaryCondition(blockManager);
PM_Stop(tm_JacobiSmoother_Comm);

PM_Stop(tm_JacobiSmoother);
	}

	void Jacobi_Smoother_Mask_2(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsb,
						LocalScalar3D<real>* plsx0,
						LocalScalar3D<int>* plsMaskId,
						real omega) {
		int vc = plsx->GetVC();

#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* x  = plsx ->GetBlockData(block);
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* b  = plsb ->GetBlockData(block);
			real* x0 = plsx0->GetBlockData(block);
			int* mask = plsMaskId->GetBlockData(block);
			real pomega = omega;

			jacobi_smoother_mask_(
							x0, x,
							Ap, Aw, Ae, As, An, Ab, At,
							b,
							mask,
							&pomega,
							sz, &vc);
		}

		plsx0->UpdateVirtualCells(blockManager);

#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* x  = plsx ->GetBlockData(block);
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* b  = plsb ->GetBlockData(block);
			real* x0 = plsx0->GetBlockData(block);
			int* mask = plsMaskId->GetBlockData(block);
			real pomega = omega;

			jacobi_smoother_mask_(
							x, x0,
							Ap, Aw, Ae, As, An, Ab, At,
							b,
							mask,
							&pomega,
							sz, &vc);
		}

		plsx->ImposeBoundaryCondition(blockManager);

	}

	void CalcAx_Mask(
						BlockManager& blockManager,
						LocalScalar3D<real>* plsAx,
						LocalScalar3D<real>* plsAp,
						LocalScalar3D<real>* plsAw,
						LocalScalar3D<real>* plsAe,
						LocalScalar3D<real>* plsAs,
						LocalScalar3D<real>* plsAn,
						LocalScalar3D<real>* plsAb,
						LocalScalar3D<real>* plsAt,
						LocalScalar3D<real>* plsx,
						LocalScalar3D<int>* plsMaskId) {
		int vc = plsx->GetVC();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for (int n=0; n<blockManager.getNumBlock(); ++n) {
			BlockBase* block = blockManager.getBlock(n);
			::Vec3i blockSize = block->getSize();
			::Vec3d cellSize  = block->getCellSize();
			int sz[3] = {blockSize.x, blockSize.y, blockSize.z};

			real* Ax = plsAx->GetBlockData(block);
			real* Ap = plsAp->GetBlockData(block);
			real* Aw = plsAw->GetBlockData(block);
			real* Ae = plsAe->GetBlockData(block);
			real* As = plsAs->GetBlockData(block);
			real* An = plsAn->GetBlockData(block);
			real* Ab = plsAb->GetBlockData(block);
			real* At = plsAt->GetBlockData(block);
			real* x  = plsx ->GetBlockData(block);
			int* mask = plsMaskId->GetBlockData(block);

			calc_ax_mask_(
							Ax,
							Ap, Aw, Ae, As, An, Ab, At,
							x,
							mask,
							sz, &vc);
		}
	}

};

#endif

