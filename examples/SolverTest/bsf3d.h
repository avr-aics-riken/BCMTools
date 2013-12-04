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


