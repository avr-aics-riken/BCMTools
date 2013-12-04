#ifndef __BCAX_H__
#define __BCAX_H__

#include "real.h"

extern "C" {
	void bc_x3_poiseuille_u_(
					real *x, real *xc
				, int *sz, int *g
				, double *org, double *blockSize, double *cellSize);

	void bc_aw_poiseuille_u_(
					real* Ap,
					real* Aw,
					real* b,
					real* xc,
					int* sz, int* g,
					double *org, double *blockSize, double *cellSize);

	void bc_x3_poiseuille_p_(
					real *x, real *xc
				, int *sz, int *g
				, double *org, double *blockSize, double *cellSize);

	void bc_aw_poiseuille_p_(
					real* Ap,
					real* Aw,
					real* b,
					real* xc,
					int* sz, int* g,
					double *org, double *blockSize, double *cellSize);
}

#endif

