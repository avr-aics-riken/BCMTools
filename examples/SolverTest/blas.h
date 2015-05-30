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

#ifndef __BLAS_H__
#define __BLAS_H__

#include "real.h"

extern "C" {
	void fill_(
					real *x, real *a
				, int *sz, int *g);
	void fill_2_(
					real *x, real *a
				, int *sz, int *g);

	void fill_vf3d_(
					real *x, real *a
				, int *sz, int *g, int *ne);

	void copy_(
					real *y, real *x
				, int *sz, int *g);
	void copy_integer_(
					int *y, int *x
				, int *sz, int *g);
	void add_(
					real *C, real *A, real *B
				, int *sz, int *g);
	void triad_(
					real *C, real *A, real *B, real *d
				, int *sz, int *g);
	void scal_(
					real *y, real *a
				, int *sz, int *g);
	void axpy_(
					real *y, real *x, real *a
				, int *sz, int *g);
	void xpay_(
					real *y, real *x, real *a
				, int *sz, int *g);
	void axpyz_(
					real *z, real *x, real *y, real *a
				, int *sz, int *g);
	void axpbypz_(
					real *z, real *x, real *y, real *a, real *b
				, int *sz, int *g);
	void dot_(
					real *xy, real *y, real *x
				, int *sz, int *g);
	void dotx2_(
					real *xy, real *xz, real *x, real *y, real *z
				, int *sz, int *g);


	void jacobi_smoother_(
					real* x1, real* x0,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* b,
					real* omega,
					int* sz, int* g);
	void jacobi_smoother2_(
					real* x1, real* x0,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* b,
					real* omega,
					int* sz, int* g);
	void calc_ax_( 
					real* Ax,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* x,
					int* sz, int* g);
	void calc_r_( 
					real* r,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* x,
					real* b,
					int* sz, int* g);
	void calc_r2_( 
					real* r,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* x,
					real* b,
					int* sz, int* g);


	void setup_mask_(
					real* m,
					int* mask,
					int* sz, int* g);
	void copy_mask_(
					real *y, real *x,
					int* mask,
					int *sz, int *g);
	void dot_mask_(
					real *xy, real *y, real *x,
					int* mask,
					int *sz, int *g);
	void calc_ax_mask_( 
					real* Ax,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* x,
					int* mask,
					int* sz, int* g);
	void jacobi_smoother_mask_(
					real* x1, real* x0,
					real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At,
					real* b,
					int* mask,
					real* omega,
					int* sz, int* g);
}

#endif

