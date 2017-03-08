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
