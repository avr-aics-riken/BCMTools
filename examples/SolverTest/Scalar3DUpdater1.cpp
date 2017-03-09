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

#include "Scalar3DUpdater1.h"

#include "sup.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

/// 隣接データクラスから仮想セルデータをコピー(同レベル間).
template <>
void Scalar3DUpdater1<real>::copyFromNeighbor(Face face)
{
  Scalar3D<real>* dc_src = neighborDataClass[face][0];
  real* data_src = dc_src->getData();

	if( !data_src ) {
		return;
	}

	real* data_dst = dataClass->getData();
	int i1_dst[3] = {1, 1, 1};
	int i1_src[3] = {1, 1, 1};
	int sz_c[3] = {nx, ny, nz};
	int sz[3] = {nx, ny, nz};
	int g[1] = {vc};

	switch( face ) {
		case X_M:
			i1_dst[0] = 1 - vc;
			i1_dst[1] = 1;
			i1_dst[2] = 1;
			i1_src[0] = 1 + nx - vc;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = vc;
			sz_c[1] = ny;
			sz_c[2] = nz;
			break;
		case X_P:
			i1_dst[0] = 1 + nx;
			i1_dst[1] = 1;
			i1_dst[2] = 1;
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = vc;
			sz_c[1] = ny;
			sz_c[2] = nz;
			break;
		case Y_M:
			i1_dst[0] = 1;
			i1_dst[1] = 1 - vc;
			i1_dst[2] = 1;
			i1_src[0] = 1;
			i1_src[1] = 1 + ny - vc;
			i1_src[2] = 1;
			sz_c[0] = nx;
			sz_c[1] = vc;
			sz_c[2] = nz;
			break;
		case Y_P:
			i1_dst[0] = 1;
			i1_dst[1] = 1 + ny;
			i1_dst[2] = 1;
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = nx;
			sz_c[1] = vc;
			sz_c[2] = nz;
			break;
		case Z_M:
			i1_dst[0] = 1;
			i1_dst[1] = 1;
			i1_dst[2] = 1 - vc;
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1 + nz - vc;
			sz_c[0] = nx;
			sz_c[1] = ny;
			sz_c[2] = vc;
			break;
		case Z_P:
			i1_dst[0] = 1;
			i1_dst[1] = 1;
			i1_dst[2] = 1 + nz;
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = nx;
			sz_c[1] = ny;
			sz_c[2] = vc;
			break;
		default:
			break;
	}

	sup_copy_from_neighbor_(data_dst, i1_dst, data_src, i1_src, sz_c, sz, g);
}

/// 送信バッファに仮想セルデータをコピー(同レベル間).
template <>
void Scalar3DUpdater1<real>::copyToCommBuffer(Face face)
{
  real* buffer = sendBuffer[face][0];
  if ( !buffer ) {
		return;
	}

	real* data_src = dataClass->getData();
	int i1_src[3] = {1, 1, 1};
	int sz_c[3] = {nx, ny, nz};
	int sz[3] = {nx, ny, nz};
	int g[1] = {vc};

  switch (face) {
    case X_M:
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = vc;
			sz_c[1] = ny;
			sz_c[2] = nz;
      break;
    case X_P:
			i1_src[0] = 1 + nx - vc;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = vc;
			sz_c[1] = ny;
			sz_c[2] = nz;
      break;
    case Y_M:
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = nx;
			sz_c[1] = vc;
			sz_c[2] = nz;
      break;
    case Y_P:
			i1_src[0] = 1;
			i1_src[1] = 1 + ny - vc;
			i1_src[2] = 1;
			sz_c[0] = nx;
			sz_c[1] = vc;
			sz_c[2] = nz;
      break;
    case Z_M:
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = nx;
			sz_c[1] = ny;
			sz_c[2] = vc;
      break;
    case Z_P:
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1 + nz - vc;
			sz_c[0] = nx;
			sz_c[1] = ny;
			sz_c[2] = vc;
      break;
    default:
      break;
  }

	sup_copy_to_buffer_(buffer, data_src, i1_src, sz_c, sz, g);
}


/// 受信バッファから仮想セルデータをコピー(同レベル間).
template <>
void Scalar3DUpdater1<real>::copyFromCommBuffer(Face face)
{
  real* buffer = recvBuffer[face][0];
  if ( !buffer ) {
		return;
	}

	real* data_dst = dataClass->getData();
	int i1_dst[3] = {1, 1, 1};
	int sz_c[3] = {nx, ny, nz};
	int sz[3] = {nx, ny, nz};
	int g[1] = {vc};

  switch (face) {
    case X_M:
			i1_dst[0] = 1 - vc;
			i1_dst[1] = 1;
			i1_dst[2] = 1;
			sz_c[0] = vc;
			sz_c[1] = ny;
			sz_c[2] = nz;
      break;
    case X_P:
			i1_dst[0] = 1 + nx;
			i1_dst[1] = 1;
			i1_dst[2] = 1;
			sz_c[0] = vc;
			sz_c[1] = ny;
			sz_c[2] = nz;
      break;
    case Y_M:
			i1_dst[0] = 1;
			i1_dst[1] = 1 - vc;
			i1_dst[2] = 1;
			sz_c[0] = nx;
			sz_c[1] = vc;
			sz_c[2] = nz;
      break;
    case Y_P:
			i1_dst[0] = 1;
			i1_dst[1] = 1 + ny;
			i1_dst[2] = 1;
			sz_c[0] = nx;
			sz_c[1] = vc;
			sz_c[2] = nz;
      break;
    case Z_M:
			i1_dst[0] = 1;
			i1_dst[1] = 1;
			i1_dst[2] = 1 - vc;
			sz_c[0] = nx;
			sz_c[1] = ny;
			sz_c[2] = vc;
      break;
    case Z_P:
			i1_dst[0] = 1;
			i1_dst[1] = 1;
			i1_dst[2] = 1 + nz;
			sz_c[0] = nx;
			sz_c[1] = ny;
			sz_c[2] = vc;
      break;
    default:
      break;
  }

	sup_copy_from_buffer_(data_dst, i1_dst, buffer, sz_c, sz, g);
}


/// 隣接データクラスから仮想セルデータをコピー(レベルL→L+1).
template <>
void Scalar3DUpdater1<real>::copyFromNeighborC2F(Face face, Subface subface)
{
  Scalar3D<real>* dc_src = neighborDataClass[face][0];
  real* data_src = dc_src->getData();

	if( !data_src ) {
		return;
	}

	real* data_dst = dataClass->getData();
	int i1_dst[3] = {1, 1, 1};
	int i1_src[3] = {1, 1, 1};
	int sz_c[3] = {nx, ny, nz};
	int sz[3] = {nx, ny, nz};
	int g[1] = {vc};

  switch (face) {
    case X_M:
    {
			i1_dst[0] = 1 - vc;
			i1_dst[1] = 1;
			i1_dst[2] = 1;
			i1_src[0] = 1 + 2*nx - vc;
			i1_src[1] = 1 + ny * subfaceOrigin0(subface);
			i1_src[2] = 1 + nz * subfaceOrigin1(subface);
			sz_c[0] = vc;
			sz_c[1] = ny;
			sz_c[2] = nz;
      break;
    }
    case X_P:
    {
			i1_dst[0] = 1 + nx;
			i1_dst[1] = 1;
			i1_dst[2] = 1;
			i1_src[0] = 1;
			i1_src[1] = 1 + ny * subfaceOrigin0(subface);
			i1_src[2] = 1 + nz * subfaceOrigin1(subface);
			sz_c[0] = vc;
			sz_c[1] = ny;
			sz_c[2] = nz;
      break;
    }
    case Y_M:
    {
			i1_dst[0] = 1;
			i1_dst[1] = 1 - vc;
			i1_dst[2] = 1;
			i1_src[0] = 1 + nx * subfaceOrigin1(subface);
			i1_src[1] = 1 + 2*ny - vc;
			i1_src[2] = 1 + nz * subfaceOrigin0(subface);
			sz_c[0] = nx;
			sz_c[1] = vc;
			sz_c[2] = nz;
      break;
    }
    case Y_P:
    {
			i1_dst[0] = 1;
			i1_dst[1] = 1 + ny;
			i1_dst[2] = 1;
			i1_src[0] = 1 + nx * subfaceOrigin1(subface);
			i1_src[1] = 1;
			i1_src[2] = 1 + nz * subfaceOrigin0(subface);
			sz_c[0] = nx;
			sz_c[1] = vc;
			sz_c[2] = nz;
      break;
    }
    case Z_M:
    {
			i1_dst[0] = 1;
			i1_dst[1] = 1;
			i1_dst[2] = 1 - vc;
			i1_src[0] = 1 + nx * subfaceOrigin0(subface);
			i1_src[1] = 1 + ny * subfaceOrigin1(subface);
			i1_src[2] = 1 + 2*nz - vc;
			sz_c[0] = nx;
			sz_c[1] = ny;
			sz_c[2] = vc;
      break;
    }
    case Z_P:
    {
			i1_dst[0] = 1;
			i1_dst[1] = 1;
			i1_dst[2] = 1 + nz;
			i1_src[0] = 1 + nx * subfaceOrigin0(subface);
			i1_src[1] = 1 + ny * subfaceOrigin1(subface);
			i1_src[2] = 1;
			sz_c[0] = nx;
			sz_c[1] = ny;
			sz_c[2] = vc;
      break;
    }
    default:
      break;
  }

	sup_copy_from_neighbor_c2f_(data_dst, i1_dst, data_src, i1_src, sz_c, sz, g);
}


/// 送信バッファに仮想セルデータをコピー(レベルL→L+1).
template <>
void Scalar3DUpdater1<real>::copyToCommBufferC2F(Face face, Subface subface)
{
  real* buffer = sendBuffer[face][subface];
  if ( !buffer ) {
		return;
	}

	real* data_src = dataClass->getData();
	int i1_src[3] = {1, 1, 1};
	int sz_c[3] = {nx, ny, nz};
	int sz[3] = {nx, ny, nz};
	int g[1] = {vc};

  switch (face) {
    case X_M:
			i1_src[0] = 1;
			i1_src[1] = 1 + ny * subfaceOrigin0(subface);
			i1_src[2] = 1 + nz * subfaceOrigin1(subface);
			sz_c[0] = vc;
			sz_c[1] = ny;
			sz_c[2] = nz;
      break;
    case X_P:
			i1_src[0] = 1 + 2*nx - vc;
			i1_src[1] = 1 + ny * subfaceOrigin0(subface);
			i1_src[2] = 1 + nz * subfaceOrigin1(subface);
			sz_c[0] = vc;
			sz_c[1] = ny;
			sz_c[2] = nz;
      break;
    case Y_M:
			i1_src[0] = 1 + nx * subfaceOrigin1(subface);
			i1_src[1] = 1;
			i1_src[2] = 1 + nz * subfaceOrigin0(subface);
			sz_c[0] = nx;
			sz_c[1] = vc;
			sz_c[2] = nz;
      break;
    case Y_P:
			i1_src[0] = 1 + nx * subfaceOrigin1(subface);
			i1_src[1] = 1 + 2*ny - vc;
			i1_src[2] = 1 + nz * subfaceOrigin0(subface);
			sz_c[0] = nx;
			sz_c[1] = vc;
			sz_c[2] = nz;
      break;
    case Z_M:
			i1_src[0] = 1 + nx * subfaceOrigin0(subface);
			i1_src[1] = 1 + ny * subfaceOrigin1(subface);
			i1_src[2] = 1;
			sz_c[0] = nx;
			sz_c[1] = ny;
			sz_c[2] = vc;
      break;
    case Z_P:
			i1_src[0] = 1 + nx * subfaceOrigin0(subface);
			i1_src[1] = 1 + ny * subfaceOrigin1(subface);
			i1_src[2] = 1 + 2*nz - vc;
			sz_c[0] = nx;
			sz_c[1] = ny;
			sz_c[2] = vc;
      break;
    default:
      break;
  }

	sup_copy_to_buffer_c2f_(buffer, data_src, i1_src, sz_c, sz, g);
}


/// 受信バッファから仮想セルデータをコピー(レベルL→L+1).
template <>
void Scalar3DUpdater1<real>::copyFromCommBufferC2F(Face face, Subface subface)
{
	Scalar3DUpdater1<real>::copyFromCommBuffer(face);
}


/// 隣接データクラスから仮想セルデータをコピー(レベルL+1→L).
template <>
void Scalar3DUpdater1<real>::copyFromNeighborF2C(Face face, Subface subface)
{
  Scalar3D<real>* dc_src = neighborDataClass[face][subface];
  real* data_src = dc_src->getData();

	if( !data_src ) {
		return;
	}

	real* data_dst = dataClass->getData();
	int i1_dst[3] = {1, 1, 1};
	int i1_src[3] = {1, 1, 1};
	int sz_c[3] = {nx, ny, nz};
	int sz[3] = {nx, ny, nz};
	int g[1] = {vc};

	switch( face ) {
		case X_M:
			i1_dst[0] = 1 - vc;
			i1_dst[1] = 1 + (ny/2) * subfaceOrigin0(subface);
			i1_dst[2] = 1 + (nz/2) * subfaceOrigin1(subface);
			i1_src[0] = 1 + nx/2 - vc;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = vc;
			sz_c[1] = ny/2;
			sz_c[2] = nz/2;
			break;
		case X_P:
			i1_dst[0] = 1 + nx;
			i1_dst[1] = 1 + (ny/2) * subfaceOrigin0(subface);
			i1_dst[2] = 1 + (nz/2) * subfaceOrigin1(subface);
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = vc;
			sz_c[1] = ny/2;
			sz_c[2] = nz/2;
			break;
		case Y_M:
			i1_dst[0] = 1 + (nx/2) * subfaceOrigin1(subface);
			i1_dst[1] = 1 - vc;
			i1_dst[2] = 1 + (nz/2) * subfaceOrigin0(subface);
			i1_src[0] = 1;
			i1_src[1] = 1 + ny/2 - vc;
			i1_src[2] = 1;
			sz_c[0] = nx/2;
			sz_c[1] = vc;
			sz_c[2] = nz/2;
			break;
		case Y_P:
			i1_dst[0] = 1 + (nx/2) * subfaceOrigin1(subface);
			i1_dst[1] = 1 + ny;
			i1_dst[2] = 1 + (nz/2) * subfaceOrigin0(subface);
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = nx/2;
			sz_c[1] = vc;
			sz_c[2] = nz/2;
			break;
		case Z_M:
			i1_dst[0] = 1 + (nx/2) * subfaceOrigin0(subface);
			i1_dst[1] = 1 + (ny/2) * subfaceOrigin1(subface);
			i1_dst[2] = 1 - vc;
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1 + nz/2 - vc;
			sz_c[0] = nx/2;
			sz_c[1] = ny/2;
			sz_c[2] = vc;
			break;
		case Z_P:
			i1_dst[0] = 1 + (nx/2) * subfaceOrigin0(subface);
			i1_dst[1] = 1 + (ny/2) * subfaceOrigin1(subface);
			i1_dst[2] = 1 + nz;
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = nx/2;
			sz_c[1] = ny/2;
			sz_c[2] = vc;
			break;
		default:
			break;
	}

	sup_copy_from_neighbor_f2c_(data_dst, i1_dst, data_src, i1_src, sz_c, sz, g);
}


/// 送信バッファに仮想セルデータをコピー(レベルL+1→L).
template <>
void Scalar3DUpdater1<real>::copyToCommBufferF2C(Face face, Subface subface)
{
  real* buffer = sendBuffer[face][0];

	real* data_src = dataClass->getData();
	int i1_src[3] = {1, 1, 1};
	int sz_c[3] = {nx, ny, nz};
	int sz[3] = {nx, ny, nz};
	int g[1] = {vc};

  switch (face) {
    case X_M:
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = vc;
			sz_c[1] = ny/2;
			sz_c[2] = nz/2;
      break;
    case X_P:
			i1_src[0] = 1 + nx/2 - vc;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = vc;
			sz_c[1] = ny/2;
			sz_c[2] = nz/2;
      break;
    case Y_M:
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = nx/2;
			sz_c[1] = vc;
			sz_c[2] = nz/2;
      break;
    case Y_P:
			i1_src[0] = 1;
			i1_src[1] = 1 + ny/2 - vc;
			i1_src[2] = 1;
			sz_c[0] = nx/2;
			sz_c[1] = vc;
			sz_c[2] = nz/2;
      break;
    case Z_M:
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1;
			sz_c[0] = nx/2;
			sz_c[1] = ny/2;
			sz_c[2] = vc;
      break;
    case Z_P:
			i1_src[0] = 1;
			i1_src[1] = 1;
			i1_src[2] = 1 + nz/2 - vc;
			sz_c[0] = nx/2;
			sz_c[1] = ny/2;
			sz_c[2] = vc;
      break;
    default:
      break;
  }

	sup_copy_to_buffer_f2c_(buffer, data_src, i1_src, sz_c, sz, g);
}


/// 受信バッファから仮想セルデータをコピー(レベルL+1→L).
template <>
void Scalar3DUpdater1<real>::copyFromCommBufferF2C(Face face, Subface subface)
{
  real* buffer = recvBuffer[face][subface];
  if ( !buffer ) {
		return;
	}

	real* data_dst = dataClass->getData();
	int i1_dst[3] = {1, 1, 1};
	int sz_c[3] = {nx, ny, nz};
	int sz[3] = {nx, ny, nz};
	int g[1] = {vc};

  switch (face) {
    case X_M:
			i1_dst[0] = 1 - vc;
			i1_dst[1] = 1 + (ny/2) * subfaceOrigin0(subface);
			i1_dst[2] = 1 + (nz/2) * subfaceOrigin1(subface);
			sz_c[0] = vc;
			sz_c[1] = ny/2;
			sz_c[2] = nz/2;
      break;
    case X_P:
			i1_dst[0] = 1 + nx;
			i1_dst[1] = 1 + (ny/2) * subfaceOrigin0(subface);
			i1_dst[2] = 1 + (nz/2) * subfaceOrigin1(subface);
			sz_c[0] = vc;
			sz_c[1] = ny/2;
			sz_c[2] = nz/2;
      break;
    case Y_M:
			i1_dst[0] = 1 + (nx/2) * subfaceOrigin1(subface);
			i1_dst[1] = 1 - vc;
			i1_dst[2] = 1 + (nz/2) * subfaceOrigin0(subface);
			sz_c[0] = nx/2;
			sz_c[1] = vc;
			sz_c[2] = nz/2;
      break;
    case Y_P:
			i1_dst[0] = 1 + (nx/2) * subfaceOrigin1(subface);
			i1_dst[1] = 1 + ny;
			i1_dst[2] = 1 + (nz/2) * subfaceOrigin0(subface);
			sz_c[0] = nx/2;
			sz_c[1] = vc;
			sz_c[2] = nz/2;
      break;
    case Z_M:
			i1_dst[0] = 1 + (nx/2) * subfaceOrigin0(subface);
			i1_dst[1] = 1 + (ny/2) * subfaceOrigin1(subface);
			i1_dst[2] = 1 - vc;
			sz_c[0] = nx/2;
			sz_c[1] = ny/2;
			sz_c[2] = vc;
      break;
    case Z_P:
			i1_dst[0] = 1 + (nx/2) * subfaceOrigin0(subface);
			i1_dst[1] = 1 + (ny/2) * subfaceOrigin1(subface);
			i1_dst[2] = 1 + nz;
			sz_c[0] = nx/2;
			sz_c[1] = ny/2;
			sz_c[2] = vc;
      break;
    default:
      break;
  }

	sup_copy_from_buffer_(data_dst, i1_dst, buffer, sz_c, sz, g);
}


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif
