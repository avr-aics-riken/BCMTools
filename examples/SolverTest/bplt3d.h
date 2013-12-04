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

#ifndef __BPLT3D_H__
#define __BPLT3D_H__

extern "C" {
	void bplt3d_open_file_(char* filename, int* filenamelength, int* unit);
	void bplt3d_close_file_(int* unit);
	void bplt3d_write_xyz_header_(int* ix, int* jx, int* kx, int* ngrid, int* unit);
	void bplt3d_write_xyz_block_(float* x, float* y, float* z, int* ix, int* jx, int* kx, int* unit);
	void bplt3d_write_func_header_(int* ix, int* jx, int* kx, int* nvar, int* ngrid, int* unit);
	void bplt3d_write_func_block_(float* p, float* ux, float* uy, float* uz, int* ix, int* jx, int* kx, int* unit);
}

#endif

