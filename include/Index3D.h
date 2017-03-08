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
/// @file Index3D.h
/// @brief データクラス用インデックスファンクタクラス
///

#ifndef INDEX_3D_H
#define INDEX_3D_H

#include "Vec3.h"

using namespace Vec3class;

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

/// スカラデータクラス用インデックスファンクタ.
class Index3DS {

  int c_j, c_k, c_0;

public:

  /// コンストラクタ.
  ///
  ///  @param[in] size セル分割数
  ///  @param[in] vc 仮想セル幅
  ///
  Index3DS(const Vec3i& size, int vc) {
    int nx0 = size[0] + 2*vc;
    int ny0 = size[1] + 2*vc;
  //int nz0 = size[2] + 2*vc;
    c_j = nx0;
    c_k = nx0 * ny0;
    c_0 = (1 + nx0 + nx0 * ny0) * vc;
  }

  /// 3次元→1次元 インデックス変換.
  int operator() (int i, int j, int k) const {
    return i + c_j * j + c_k * k + c_0;
  }
};


/// ベクトルデータクラス(3成分連続格納)用インデックスファンクタ.
class Index3DV {

  int c_i, c_j, c_k, c_0;

public:

  /// コンストラクタ.
  ///
  ///  @param[in] size セル分割数
  ///  @param[in] vc 仮想セル幅
  ///
  Index3DV(const Vec3i& size, int vc) {
    int nx0 = size[0] + 2*vc;
    int ny0 = size[1] + 2*vc;
  //int nz0 = size[2] + 2*vc;

    c_i = 3;
    c_j = 3 * nx0;
    c_k = 3 * nx0 * ny0;
    c_0 = 3 * (1 + nx0 + nx0 * ny0) * vc;
  }

  /// 3次元→1次元 インデックス変換.
  int operator() (int i, int j, int k) const {
    return c_i * i + c_j * j + c_k * k + c_0;
  }

  /// 3次元＋ベクトル成分→1次元 インデックス変換.
  int operator() (int i, int j, int k, int l) const {
    return l + c_i * i + c_j * j + c_k * k + c_0;
  }

};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // INDEX_3D_H
