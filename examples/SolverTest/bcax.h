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

#ifndef __BCAX_H__
#define __BCAX_H__

#include "real.h"

extern "C" {
	void bc_x3_poiseuille_u_(
					real *x, real *xc
				, int *sz, int *g
				, double *center
				, double *org, double *blockSize, double *cellSize);

	void bc_aw_poiseuille_u_(
					real* Ap,
					real* Aw,
					real* b,
					real* xc,
					int* sz, int* g,
					double *center,
					double *org, double *blockSize, double *cellSize);

	void bc_x3_poiseuille_p_(
					real *x, real *xc
				, int *sz, int *g
				, double *center
				, double *org, double *blockSize, double *cellSize);

	void bc_aw_poiseuille_p_(
					real* Ap,
					real* Aw,
					real* b,
					real* xc,
					int* sz, int* g,
					double *center,
					double *org, double *blockSize, double *cellSize);
}

#endif

