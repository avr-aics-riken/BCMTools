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
/// @file SphereDivider3.h
/// @brief ブロック分割判定クラス(PolygonBBoxDivider)
/// 

#ifndef SPHERE_DIVIDER3_H
#define SPHERE_DIVIDER3_H

#include <string>
#include <vector>
#include "BCMPolylib.h"
#include "BCMTools.h"
#include "MultiRootDivider.h"
#include "BoundingBox.h"
#include "PolygonBBoxDivider.h"

class SphereDivider3 : public MultiRootDivider {

  const int minLevel;  ///< 最小分割レベル

  const PolylibNS::BCMPolylib* pl;  ///< BCMPolylibオブジェクト

  /// <ポリゴングループ名, 分割レベル>ペアのリスト
  const std::vector<PolygonGroupSpec>& polygonGroupList;

  /// <バウンディングボックス, 分割レベル>ペアのリスト
  const std::vector<BoundingBoxSpec>& boundingBoxList; 

  const std::vector<BoundingBoxSpec>& sphericalBoxList; 

  double extraMarginRatio;  ///< 追加マージン幅の最大分割ブロック辺長に対する比

public:

  /// コンストラクタ.
  ///
  ///  @param[in] origin 最初のルートブロックの原点位置
  ///  @param[in] rootLength ルートブロックの辺長
  ///  @param[in] rootGrid ルートブロック配置情報
  ///  @param[in] minLevel 最小分割レベル
  ///  @param[in] pl       Polylibオブジェクト
  ///  @param[in] polygonGroupList <ポリゴングループ名, 分割レベル>ペアのリスト
  ///  @param[in] boundingBoxList <バウンディングボックス, 分割レベル>ペアのリスト
  ///  @param[in] extraMarginRatio 追加マージンの幅
  ///
  ///  @note 境界探査領域にマージンを追加する場合には，
  ///  追加マージン幅の最大分割レベルブロック辺長に対する比率を
  ///  extraMarginRatioに指定する．
  ///  例えば，仮想セル領域を境界面探査領域に追加するには
  ///  「(double)仮想セル数/ブロック内分割数」を指定する．
  ///
  SphereDivider3(const Vec3r& origin, double rootLength, const RootGrid* rootGrid, 
                 int minLevel, const PolylibNS::BCMPolylib* pl,
                 const std::vector<PolygonGroupSpec>& polygonGroupList,
                 const std::vector<BoundingBoxSpec>& boundingBoxList, 
                 const std::vector<BoundingBoxSpec>& sphericalBoxList, 
                 double extraMarginRatio = 0.0)
    : MultiRootDivider(origin, rootLength, rootGrid),
      minLevel(minLevel), pl(pl),
      polygonGroupList(polygonGroupList), boundingBoxList(boundingBoxList), sphericalBoxList(sphericalBoxList),
      extraMarginRatio(extraMarginRatio) {
  }

  /// デストラクタ.
  ~SphereDivider3() {}

  /// ブロックを分割するかどうかを判定.
  ///
  ///   @param[in] pedigree ブロックのPedigree
  ///   @return ブロックタイプ
  ///
  NodeType operator() (const Pedigree& pedigree) {
    int level = pedigree.getLevel();

    if  (level < minLevel) return BRANCH;

    NodeType ret = LEAF_ACTIVE;

		BoundingBox sphere = sphericalBoxList.begin()[0].boundingBox;
		int maxLevel       = sphericalBoxList.begin()[0].level;

	  BoundingBox region = defineSearchRegion(pedigree, maxLevel);
    region.setMargin(extraMarginRatio / (1 << maxLevel));

		if( checkInner(sphere, region) ) {
			return LEAF_NO_ACTIVE;
		}

		if( checkOuter(sphere, region) ) {
			ret = LEAF_ACTIVE;
			for (std::vector<BoundingBoxSpec>::const_iterator it = sphericalBoxList.begin();
					 it != sphericalBoxList.end(); ++it) {
				int maxLevel = it->level;
				if (level < maxLevel) {
					BoundingBox box = it->boundingBox;
					if (!checkOuter(box, region)) ret = BRANCH;
				}
			}
			for (std::vector<BoundingBoxSpec>::const_iterator it = boundingBoxList.begin();
					 it != boundingBoxList.end(); ++it) {
				int maxLevel = it->level;
				if (level < maxLevel) {
					BoundingBox box = it->boundingBox;
					if (box.intersects(region)) ret = BRANCH;
				}
			}
			return ret;
		}

		if( level == maxLevel ) {
			return LEAF_ACTIVE;
		}

		return BRANCH;
  }

	bool checkIntersect(BoundingBox &box, BoundingBox &rgn) {
		if( checkInner(box, rgn) ) return false;
		if( checkOuter(box, rgn) ) return false;
		return true;
	}

	bool checkInner(BoundingBox &box, BoundingBox &rgn) {
		double x0 = rgn.getMin().x;
		double y0 = rgn.getMin().y;
		double z0 = rgn.getMin().z;
		double x1 = rgn.getMax().x;
		double y1 = rgn.getMax().y;
		double z1 = rgn.getMax().z;
		double ox = 0.5*(box.getMin().x + box.getMax().x);
		double oy = 0.5*(box.getMin().y + box.getMax().y);
		double oz = 0.5*(box.getMin().z + box.getMax().z);
		double r = 0.5*(box.getMax().x - box.getMin().x);
    double x = std::max(fabs(x0-ox), fabs(x1-ox));
    double r2max = x * x;
		double r2 = r*r;
    if (r2max >= r2) return false;

    double y = std::max(fabs(y0-oy), fabs(y1-oy));
    r2max += y * y;
    if (r2max >= r2) return false;

    double z = std::max(fabs(z0-oz), fabs(z1-oz));
    r2max += z * z;
    if (r2max >= r2) return false;

    return true;
	}

	bool checkOuter(BoundingBox &box, BoundingBox &rgn) {
		double x0 = rgn.getMin().x;
		double y0 = rgn.getMin().y;
		double z0 = rgn.getMin().z;
		double x1 = rgn.getMax().x;
		double y1 = rgn.getMax().y;
		double z1 = rgn.getMax().z;
		double ox = 0.5*(box.getMin().x + box.getMax().x);
		double oy = 0.5*(box.getMin().y + box.getMax().y);
		double oz = 0.5*(box.getMin().z + box.getMax().z);
		double r = 0.5*(box.getMax().x - box.getMin().x);
		double r2 = r*r;
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

    double z = 0.0;
    if      (oz < z0) z = z0 - oz;
    else if (oz > z1) z = oz - z1;
    r2min += z * z;
    if (r2min >= r2) return true;

    return false;
  }
};

#endif // SPHERE_DIVIDER3_H
