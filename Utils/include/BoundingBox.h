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
/// @file BoundingBox.h
/// @brief バウンディングボックスクラス
/// 

#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include <algorithm>
#include <iostream>
#include "Vec3.h"

using namespace Vec3class;

/// バウンディングボックスクラス.
class BoundingBox {

  Vec3d min;  ///< 最小位置
  Vec3d max;  ///< 最大位置

public:

  /// コンストラクタ(値設定なし).
  BoundingBox()
    : min(HUGE_VAL), max(-HUGE_VAL) {}

  /// コンストラクタ(コピー).
  BoundingBox(const BoundingBox& bb)
    : min(bb.min), max(bb.max) {}

  /// コンストラクタ(値設定).
  BoundingBox(const Vec3d& min, const Vec3d& max)
    : min(min), max(max) {}

  /// コンストラクタ(値設定).
  BoundingBox(double minX, double minY, double minZ,
              double maxX, double maxY, double maxZ)
    : min(minX, minY, minZ), max(maxX, maxY, maxZ) {}

  /// デストラクタ.
  ~BoundingBox() {}

  /// 最小位置設定.
  void setMin(const Vec3d& min) { this->min = min; }

  /// 最大位置設定.
  void setMax(const Vec3d& max) { this->max = max; }

  /// 点を追加してバウンディングボックスを更新.
  void addPoint(const Vec3d& p) {
    min.x = std::min(min.x, p.x);
    min.y = std::min(min.y, p.y);
    min.z = std::min(min.z, p.z);
    max.x = std::max(max.x, p.x);
    max.y = std::max(max.y, p.y);
    max.z = std::max(max.z, p.z);
  }

  /// 直方体領域を指定してバウンディングボックスを更新.
  void addBox(const BoundingBox& box) {
    const Vec3d& boxMin = box.getMin();
    min.x = std::min(min.x, boxMin.x);
    min.y = std::min(min.y, boxMin.y);
    min.z = std::min(min.z, boxMin.z);
    const Vec3d& boxMax = box.getMax();
    max.x = std::max(max.x, boxMax.x);
    max.y = std::max(max.y, boxMax.y);
    max.z = std::max(max.z, boxMax.z);
  }

  /// 指定長だけバウンディングボックスを拡張.
  void setMargin(double margin) {
    min -= Vec3d(margin);
    max += Vec3d(margin);
  }

  /// 最小位置を取得.
  const Vec3d& getMin() const { return min; }

  /// 最大位置を取得.
  const Vec3d& getMax() const { return max; }

  /// 他のバウンディングボックスとの交わりがあるか確認.
  bool intersects(const BoundingBox& bbox) const {
    if (max.x < bbox.min.x || bbox.max.x < min.x ||
        max.y < bbox.min.y || bbox.max.y < min.y ||
        max.z < bbox.min.z || bbox.max.z < min.z) return false;
    return true;
  }

  bool inside(const BoundingBox& bbox) const {
    if (min.x < bbox.min.x && bbox.max.x < max.x &&
        min.y < bbox.min.y && bbox.max.y < max.y &&
        min.z < bbox.min.z && bbox.max.z < max.z) return true;
    return false;
	}

  bool outside(const BoundingBox& bbox) const {
		if( inside(bbox) ){
			return false;
		}
		if( intersects(bbox) ) {
			return false;
		}
		return true;
	}

  /// バウンディングボックス情報をストリームに出力.
  friend std::ostream& operator<<(std::ostream& os, const BoundingBox& bb) {
    return os << bb.min << "-" << bb.max;
  }

  /// バウンディングボックス情報をストリームから入力.
  friend std::istream& operator>>(std::istream& is, BoundingBox& bb) {
    return is >> bb.min >> bb.max;
  }

};



#endif // BOUNDING_BOX_H
