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

#ifndef __BCA7_H__
#define __BCA7_H__

#include "real.h"

extern "C" {
	void bc_aw_d_(
					real* Ap,
					real* Aw,
					real* b,
					real* xc,
					int* sz, int* g);
	void bc_ae_d_(
					real* Ap,
					real* Ae,
					real* b,
					real* xc,
					int* sz, int* g);
	void bc_as_d_(
					real* Ap,
					real* As,
					real* b,
					real* xc,
					int* sz, int* g);
	void bc_an_d_(
					real* Ap,
					real* An,
					real* b,
					real* xc,
					int* sz, int* g);
	void bc_ab_d_(
					real* Ap,
					real* Ab,
					real* b,
					real* xc,
					int* sz, int* g);
	void bc_at_d_(
					real* Ap,
					real* At,
					real* b,
					real* xc,
					int* sz, int* g);

	void bc_aw_n_(
					real* Ap,
					real* Aw,
					real* b,
					real* xc,
					int* sz, int* g);
	void bc_ae_n_(
					real* Ap,
					real* Ae,
					real* b,
					real* xc,
					int* sz, int* g);
	void bc_as_n_(
					real* Ap,
					real* As,
					real* b,
					real* xc,
					int* sz, int* g);
	void bc_an_n_(
					real* Ap,
					real* An,
					real* b,
					real* xc,
					int* sz, int* g);
	void bc_ab_n_(
					real* Ap,
					real* Ab,
					real* b,
					real* xc,
					int* sz, int* g);
	void bc_at_n_(
					real* Ap,
					real* At,
					real* b,
					real* xc,
					int* sz, int* g);
}

#endif

