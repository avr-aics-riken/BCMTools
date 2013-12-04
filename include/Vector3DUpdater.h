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

///
/// @file Vector3DUpdater.h
/// @brief ベクトルデータクラス仮想セルアップデータ(未実装)
///

#ifndef VECTOR_3D_UPDATER_H
#define VECTOR_3D_UPDATER_H

#include "BCMTools.h"
#include "VCUpdater.h"
#include "Vector3D.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif


/// ベクトルデータクラス仮想セルアップデータ(未実装).
template <typename T>
class Vector3DUpdater : public VCUpdater {

private:

  Vector3DUpdater(const NeighborInfo* neighborInfo,
                  const MPI::Comm& comm = MPI::COMM_WORLD);

  ~Vector3DUpdater();
};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // VECTOR_3D_UPDATER_H
