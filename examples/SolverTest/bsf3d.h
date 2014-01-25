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

#ifndef __BSF3D_H__
#define __BSF3D_H__

#include "real.h"

extern "C" {
	void sf3d_copy_x2_(
					real *x
				, real *xc
				, int *sz, int *g);

	void sf3d_calc_stats_(
					real *sum
				, real *max
				, real *min
				, real *absmax
				, real *absmin
				, real *data
				, int *sz, int *g);
}

#endif


