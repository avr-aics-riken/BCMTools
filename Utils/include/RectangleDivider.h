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

///
/// @file RectangleDivider.h
/// @brief ブロック分割判定クラス(RectangleDivider)
/// 

#ifndef RECTANGLE_DIVIDER_H
#define RECTANGLE_DIVIDER_H

#include <cmath>
#include <algorithm> // for max, min
#include "BCMTools.h"
#include "MultiRootDivider.h"
#include "BoundingBox.h"

/// 球面外部をブロック分割するDivider.
///
/// ルートブロックの大きさと位置は固定
///   - 大きさ: 単位立方体
///   - 位置: 最初のルートブロックの原点が(0,0,0)
///
class RectangleDivider : public MultiRootDivider {

  const int minLevel;  ///< 最小分割レベル
  const int maxLevel;  ///< 最大分割レベル

	const double rx0;
	const double ry0;
	const double rz0;
	const double rx1;
	const double ry1;
	const double rz1;

  Vec3r margin;  ///< 追加マージン

public:

  /// コンストラクタ.
  ///
  ///  @param[in] rootGrid ルートブロック配置情報
  ///  @param[in] minLevel 最小分割レベル
  ///  @param[in] maxLevel 最大分割レベル
  ///  @param[in] x0,y0,z0,x1,y1,z1 BoundingBox
  ///  @param[in] extraMarginRatio 追加マージンの幅
  ///
  ///  @note 境界探査領域にマージンを追加する場合には，
  ///  追加マージン幅の最大分割レベルブロック辺長に対する比率を
  ///  extraMarginRatioに指定する．
  ///  例えば，仮想セル領域を境界面探査領域に追加するには
  ///  「(double)仮想セル数/ブロック内分割数」を指定する．
  ///
  RectangleDivider(const RootGrid* rootGrid, int minLevel, int maxLevel,
								double rx0, double ry0, double rz0, double rx1, double ry1, double rz1,
                double extraMarginRatio = 0.0)
    : MultiRootDivider(Vec3r(0.0,0.0,0.0), 1.0, rootGrid),
      minLevel(minLevel), maxLevel(maxLevel),
      rx0(rx0), ry0(ry0), rz0(rz0), rx1(rx1), ry1(ry1), rz1(rz1) {
    assert(minLevel >= 0);
    assert(maxLevel >= minLevel);

    margin = Vec3r(extraMarginRatio / (1 << maxLevel));
  }

  /// デストラクタ.
  ~RectangleDivider() {}

  /// ブロックを分割するかどうかを判定.
  ///
  ///   @param[in] pedigree ブロックのPedigree
  ///   @return ブロックタイプ
  ///
  NodeType operator() (const Pedigree& pedigree) {
    int level = pedigree.getLevel();

    if  (level < minLevel) return BRANCH;

    BoundingBox region = defineSearchRegion(pedigree, maxLevel);

    Vec3r min = region.getMin() - margin;
    Vec3r max = region.getMax() + margin;

  //std::cout << "(x0,y0,z0) = " << min << std::endl;
  //std::cout << "(x1,y1,z1) = " << max << std::endl;

//    if (checkInner(min.x, min.y, min.z, max.x, max.y, max.z)) return LEAF_NO_ACTIVE;
//    if (checkInner(min.x, min.y, min.z, max.x, max.y, max.z)) return LEAF_ACTIVE;

    if (checkOuter(min.x, min.y, min.z, max.x, max.y, max.z)) return LEAF_ACTIVE;

    if  (level == maxLevel) return LEAF_ACTIVE;

    return BRANCH;
  }


private:

  /// 矩形領域が球面の外部にあるかどうかチェック.
  ///
  /// 矩形領域内で球中心に最も近い点を求め，その距離と球半径を比較する.
  ///
  bool checkOuter(double x0, double y0, double z0,
                  double x1, double y1, double z1) {

		if( x0 > rx1 ) {
			return true;
		}

		if( y0 > ry1 ) {
			return true;
		}

		if( z0 > rz1 ) {
			return true;
		}

		if( x1 < rx0 ) {
			return true;
		}

		if( y1 < ry0 ) {
			return true;
		}

		if( z1 < rz0 ) {
			return true;
		}

    return false;
  }

};



#endif // RECTANGLE_DIVIDER_H

