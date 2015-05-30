/*
 * BCMTools
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file CylinderDivider.h
/// @brief ブロック分割判定クラス(CylinderDivider)
/// 

#ifndef CYLINDER_DIVIDER_H
#define CYLINDER_DIVIDER_H

#include <cmath>
#include <algorithm> // for max, min
#include "MultiRootDivider.h"
#include "BoundingBox.h"

/// Z方向に延びる円筒内部をブロック分割するDivider.
///
/// ルートブロックの大きさと位置は固定
///   - 大きさ: 単位立方体
///   - 位置: 最初のルートブロックの原点が(0,0,0)
///
class CylinderDivider : public MultiRootDivider {

  const int minLevel;  ///< 最小分割レベル
  const int maxLevel;  ///< 最大分割レベル

  const double ox;  ///< 円柱中心座標
  const double oy;  ///< 円柱中心座標
  const double r;   ///< 円柱半径

  const double r2;  ///< r^2

  Vec3d margin;  ///< 追加マージン

public:

  /// コンストラクタ.
  ///
  ///  @param[in] rootGrid ルートブロック配置情報
  ///  @param[in] minLevel 最小分割レベル
  ///  @param[in] maxLevel 最大分割レベル
  ///  @param[in] ox,oy    円柱中心座標
  ///  @param[in] r        円柱半径
  ///  @param[in] extraMarginRatio 追加マージンの幅
  ///
  ///  @note 境界探査領域にマージンを追加する場合には，
  ///  追加マージン幅の最大分割レベルブロック辺長に対する比率を
  ///  extraMarginRatioに指定する．
  ///  例えば，仮想セル領域を境界面探査領域に追加するには
  ///  「(double)仮想セル数/ブロック内分割数」を指定する．
  ///
  CylinderDivider(const RootGrid* rootGrid, int minLevel, int maxLevel,
                  double ox, double oy, double r,
                  double extraMarginRatio = 0.0)
    : MultiRootDivider(Vec3d(0.0,0.0,0.0), 1.0, rootGrid),
      minLevel(minLevel), maxLevel(maxLevel),
      ox(ox), oy(oy), r(r), r2(r*r) {
    assert(minLevel >= 0);
    assert(maxLevel >= minLevel);

    margin = Vec3d(extraMarginRatio / (1 << maxLevel));
  }

  /// デストラクタ.
  ~CylinderDivider() {}

  /// ブロックを分割するかどうかを判定.
  ///
  ///   @param[in] pedigree ブロックのPedigree
  ///   @return ブロックタイプ
  ///
  NodeType operator() (const Pedigree& pedigree) {
    int level = pedigree.getLevel();

    if  (level < minLevel) return BRANCH;

    BoundingBox region = defineSearchRegion(pedigree, maxLevel);

    Vec3d min = region.getMin() - margin;
    Vec3d max = region.getMax() + margin;

  //std::cout << "(x0,y0,z0) = " << min << std::endl;
  //std::cout << "(x1,y1,z1) = " << max << std::endl;

    if (checkInner(min.x, min.y, min.z, max.x, max.y, max.z)) return LEAF_ACTIVE;

    if (checkOuter(min.x, min.y, min.z, max.x, max.y, max.z)) return LEAF_NO_ACTIVE;

    if  (level == maxLevel) return LEAF_ACTIVE;

    return BRANCH;
  }


private:

  /// 矩形領域が円筒内に完全に含まれるかどうかチェック.
  ///
  /// 矩形領域内で球中心から最も遠い点を求め，その距離と球半径を比較する.
  ///
  bool checkInner(double x0, double y0, double z0,
                  double x1, double y1, double z1) {

    double x = std::max(fabs(x0-ox), fabs(x1-ox));
    double r2max = x * x;
    if (r2max >= r2) return false;

    double y = std::max(fabs(y0-oy), fabs(y1-oy));
    r2max += y * y;
    if (r2max >= r2) return false;

    return true;
  }


  /// 矩形領域が球面の外部にあるかどうかチェック.
  ///
  /// 矩形領域内で球中心に最も近い点を求め，その距離と球半径を比較する.
  ///
  bool checkOuter(double x0, double y0, double z0,
                  double x1, double y1, double z1) {

    double x = 0.0;
    if      (ox < x0) x = x0 - ox;
    else if (ox > x1) x = ox - x1;
    double r2min = x * x;
    if (r2min >= r2) return true;

    double y = 0.0;
    if      (oy < y0) y = y0 - oy;
    else if (oy > y1) y = oy - y1;
    r2min += y * y;
    if (r2min >= r2) return true;

    return false;
  }

};



#endif // CYLINDER_DIVIDER_H
