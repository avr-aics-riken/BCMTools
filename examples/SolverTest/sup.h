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

#ifndef __SUP_H__
#define __SUP_H__

#include "real.h"

extern "C" {
	void sup_copy_from_neighbor_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup_copy_from_neighbor_c2f_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup_copy_from_neighbor_f2c_(
				real* data_dst,
				int* i1_dst,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup_copy_to_buffer_(
				real* buffer,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup_copy_to_buffer_c2f_(
				real* buffer,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup_copy_to_buffer_f2c_(
				real* buffer,
				real* data_src,
				int* i1_src,
				int* sz_c,
				int* sz, int* g);

	void sup_copy_from_buffer_(
				real* data_dst,
				int* i1_dst,
				real* buffer,
				int* sz_c,
				int* sz, int* g);

	void sup_copy_from_buffer_c2f_(
				real* data_dst,
				int* i1_dst,
				real* buffer,
				int* sz_c,
				int* sz, int* g);

	void sup_copy_from_buffer_f2c_(
				real* data_dst,
				int* i1_dst,
				real* buffer,
				int* sz_c,
				int* sz, int* g);
}

#endif

