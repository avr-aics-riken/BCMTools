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
/// @file MultiRootDivider.h
/// @brief ブロック分割判定クラス(マルチルート用派生クラス)
///

#ifndef MULTI_ROOT_DIVIDER_H
#define MULTI_ROOT_DIVIDER_H

#include "Divider.h"
#include "BoundingBox.h"

/// ブロック分割判定クラス(マルチルート用派生クラス).
class MultiRootDivider : public Divider {

  const Vec3r& origin;      ///< 最初のルートブロックの原点位置
  const REAL_TYPE rootLength;  ///< ルートブロックの辺長

  const RootGrid* rootGrid; ///< ルートブロック配置情報

public:

  /// コンストラクタ.
  ///
  ///  @param[in] origin 最初のルートブロックの原点位置
  ///  @param[in] rootLength ルートブロックの辺長
  ///  @param[in] rootGrid ルートブロック配置情報
  ///
  MultiRootDivider(const Vec3r& origin, REAL_TYPE rootLength,
                   const RootGrid* rootGrid)
    : origin(origin), rootLength(rootLength), rootGrid(rootGrid) {
  }

protected:

  /// 境界探査領域の決定.
  ///
  ///  @param[in] pedigree 現ブロックのPedigree
  ///  @param[in] maxLevel 最大分割レベル
  ///  @return 探査領域
  ///
  BoundingBox defineSearchRegion(const Pedigree& pedigree, int maxLevel) {
    int level = pedigree.getLevel();
    int px = pedigree.getX();
    int py = pedigree.getY();
    int pz = pedigree.getZ();

    int rootID = pedigree.getRootID();
    int ix = rootGrid->rootID2indexX(rootID);
    int iy = rootGrid->rootID2indexY(rootID);
    int iz = rootGrid->rootID2indexZ(rootID);

    int max0 = pedigree.getUpperBound() - 1;  // 2^level - 1
    REAL_TYPE d = 1.0 / (max0 + 1);  // ブロックサイズ

    REAL_TYPE x0 = ix + px * d;
    REAL_TYPE y0 = iy + py * d;
    REAL_TYPE z0 = iz + pz * d;
    REAL_TYPE x1 = ix + (px + 1) * d;
    REAL_TYPE y1 = iy + (py + 1) * d;
    REAL_TYPE z1 = iz + (pz + 1) * d;

    // リップル効果対策用のマージン設定
    if (level < maxLevel) {
      int n = 1 << (maxLevel - level - 1);
      REAL_TYPE dd = d * (REAL_TYPE)(n - 1) / n;
      if (!(px == 0 && rootGrid->isOuterBoundary(rootID, X_M))) x0 -= dd;
      if (!(py == 0 && rootGrid->isOuterBoundary(rootID, Y_M))) y0 -= dd;
      if (!(pz == 0 && rootGrid->isOuterBoundary(rootID, Z_M))) z0 -= dd;
      if (!(px == max0 && rootGrid->isOuterBoundary(rootID, X_P))) x1 += dd;
      if (!(py == max0 && rootGrid->isOuterBoundary(rootID, Y_P))) y1 += dd;
      if (!(pz == max0 && rootGrid->isOuterBoundary(rootID, Z_P))) z1 += dd;
    }

    // 原点移動，スケーリング
    x0 = origin.x + x0 * rootLength;
    y0 = origin.y + y0 * rootLength;
    z0 = origin.z + z0 * rootLength;
    x1 = origin.x + x1 * rootLength;
    y1 = origin.y + y1 * rootLength;
    z1 = origin.z + z1 * rootLength;

    return BoundingBox(x0, y0, z0, x1, y1, z1);
  }



};

#endif // MULTI_ROOT_DIVIDER_H
