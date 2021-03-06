/*
###################################################################################
#
# BCMTools
#
# Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2016 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
*/

#ifndef __BCX_H__
#define __BCX_H__

#include "real.h"

extern "C" {
	void bc_x1_d_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x3_d_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x2_d_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x4_d_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x5_d_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x6_d_(
					real *x, real *xc
				, int *sz, int *g);

	void bc_x1_d_f_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x3_d_f_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x2_d_f_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x4_d_f_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x5_d_f_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x6_d_f_(
					real *x, real *xc
				, int *sz, int *g);

	void bc_x1_n_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x3_n_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x2_n_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x4_n_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x5_n_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x6_n_(
					real *x, real *xc
				, int *sz, int *g);

	void bc_x1_p_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x2_p_(
					real *x, real *xc
				, int *sz, int *g);
	void bc_x5_p_(
					real *x, real *xc
				, int *sz, int *g);
}

#endif
