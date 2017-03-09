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
/// @file PolygonBBoxDivider.h
/// @brief ブロック分割判定クラス(PolygonBBoxDivider)
///

#ifndef BB_DIVIDER_H
#define BB_DIVIDER_H

#include <string>
#include <vector>
#include "BCMPolylib.h"
#include "BCMTools.h"
#include "MultiRootDivider.h"
#include "BoundingBox.h"
#include "PolygonBBoxDivider.h"

class BBDivider : public MultiRootDivider {

  const int minLevel;  ///< 最小分割レベル

  const PolylibNS::BCMPolylib* pl;  ///< BCMPolylibオブジェクト

  /// <ポリゴングループ名, 分割レベル>ペアのリスト
  const std::vector<PolygonGroupSpec>& polygonGroupList;

  /// <バウンディングボックス, 分割レベル>ペアのリスト
  const std::vector<BoundingBoxSpec>& boundingBoxList;

  const std::vector<BoundingBoxSpec>& sphericalBoxList;

  REAL_TYPE extraMarginRatio;  ///< 追加マージン幅の最大分割ブロック辺長に対する比

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
  ///  「(REAL_TYPE)仮想セル数/ブロック内分割数」を指定する．
  ///
  BBDivider(const Vec3r& origin, REAL_TYPE rootLength, const RootGrid* rootGrid,
                 int minLevel, const PolylibNS::BCMPolylib* pl,
                 const std::vector<PolygonGroupSpec>& polygonGroupList,
                 const std::vector<BoundingBoxSpec>& boundingBoxList,
                 const std::vector<BoundingBoxSpec>& sphericalBoxList,
                 REAL_TYPE extraMarginRatio = 0.0)
    : MultiRootDivider(origin, rootLength, rootGrid),
      minLevel(minLevel), pl(pl),
      polygonGroupList(polygonGroupList), boundingBoxList(boundingBoxList), sphericalBoxList(sphericalBoxList),
      extraMarginRatio(extraMarginRatio) {
  }

  /// デストラクタ.
  ~BBDivider() {}

  /// ブロックを分割するかどうかを判定.
  ///
  ///   @param[in] pedigree ブロックのPedigree
  ///   @return ブロックタイプ
  ///
  NodeType operator() (const Pedigree& pedigree) {
    int level = pedigree.getLevel();

    if  (level < minLevel) return BRANCH;

    NodeType ret = LEAF_ACTIVE;

    for (std::vector<BoundingBoxSpec>::const_iterator it = sphericalBoxList.begin();
         it != sphericalBoxList.end(); ++it) {
			if( it->level > 0 ) {
				continue;
			}
      int maxLevel = abs(it->level);
      BoundingBox box = it->boundingBox;
      BoundingBox region = defineSearchRegion(pedigree, maxLevel);
      region.setMargin(extraMarginRatio / (1 << maxLevel));
      if (level < maxLevel) {
				if (checkInner(box, region)) ret = LEAF_NO_ACTIVE;
      }
    }

    for (std::vector<BoundingBoxSpec>::const_iterator it = boundingBoxList.begin();
         it != boundingBoxList.end(); ++it) {
			if( it->level < 0 ) {
				continue;
			}
      int maxLevel = abs(it->level);
      BoundingBox box = it->boundingBox;
      BoundingBox region = defineSearchRegion(pedigree, maxLevel);
      region.setMargin(extraMarginRatio / (1 << maxLevel));
      if (level < maxLevel) {
        if (box.intersects(region)) ret = BRANCH;
      }
    }

    for (std::vector<BoundingBoxSpec>::const_iterator it = sphericalBoxList.begin();
         it != sphericalBoxList.end(); ++it) {
			if( it->level < 0 ) {
				continue;
			}
      int maxLevel = abs(it->level);
      BoundingBox box = it->boundingBox;
      BoundingBox region = defineSearchRegion(pedigree, maxLevel);
      region.setMargin(extraMarginRatio / (1 << maxLevel));
      if (level < maxLevel) {
        if (checkIntersect(box, region)) ret = BRANCH;
      }
    }

    for (std::vector<BoundingBoxSpec>::const_iterator it = sphericalBoxList.begin();
         it != sphericalBoxList.end(); ++it) {
			if( it->level > 0 ) {
				continue;
			}
      int maxLevel = abs(it->level);
      BoundingBox box = it->boundingBox;
      BoundingBox region = defineSearchRegion(pedigree, maxLevel);
      region.setMargin(extraMarginRatio / (1 << maxLevel));
      if (level < maxLevel) {
				if (checkInner(box, region)) ret = LEAF_NO_ACTIVE;
      }
    }


    for (std::vector<PolygonGroupSpec>::const_iterator it = polygonGroupList.begin();
         it != polygonGroupList.end(); ++it) {
      int maxLevel = it->level;
      if (level < maxLevel) {
        const std::string& polygonGroup = it->polygonGroupName;
        BoundingBox region = defineSearchRegion(pedigree, maxLevel);
        region.setMargin(extraMarginRatio / (1 << maxLevel));
        Vec3f min(region.getMin().x, region.getMin().y, region.getMin().z);
        Vec3f max(region.getMax().x, region.getMax().y, region.getMax().z);
        std::vector<PolylibNS::Triangle*>* polygonList
                        = pl->search_polygons(polygonGroup, min, max, false);
        int nPolygon = polygonList->size();
        delete polygonList;
        if (nPolygon > 0) ret = BRANCH;
      }
    }

    return ret;
  }

	bool checkIntersect(BoundingBox &box, BoundingBox &rgn) {
		if( checkInner(box, rgn) ) return false;
		if( checkOuter(box, rgn) ) return false;
		return true;
	}

	bool checkInner(BoundingBox &box, BoundingBox &rgn) {
		REAL_TYPE x0 = rgn.getMin().x;
		REAL_TYPE y0 = rgn.getMin().y;
		REAL_TYPE z0 = rgn.getMin().z;
		REAL_TYPE x1 = rgn.getMax().x;
		REAL_TYPE y1 = rgn.getMax().y;
		REAL_TYPE z1 = rgn.getMax().z;
		REAL_TYPE ox = 0.5*(box.getMin().x + box.getMax().x);
		REAL_TYPE oy = 0.5*(box.getMin().y + box.getMax().y);
		REAL_TYPE oz = 0.5*(box.getMin().z + box.getMax().z);
		REAL_TYPE r = 0.5*(box.getMax().x - box.getMin().x);
    REAL_TYPE x = std::max(fabs(x0-ox), fabs(x1-ox));
    REAL_TYPE r2max = x * x;
		REAL_TYPE r2 = r*r;
    if (r2max >= r2) return false;

    REAL_TYPE y = std::max(fabs(y0-oy), fabs(y1-oy));
    r2max += y * y;
    if (r2max >= r2) return false;

    REAL_TYPE z = std::max(fabs(z0-oz), fabs(z1-oz));
    r2max += z * z;
    if (r2max >= r2) return false;

    return true;
	}

	bool checkOuter(BoundingBox &box, BoundingBox &rgn) {
		REAL_TYPE x0 = rgn.getMin().x;
		REAL_TYPE y0 = rgn.getMin().y;
		REAL_TYPE z0 = rgn.getMin().z;
		REAL_TYPE x1 = rgn.getMax().x;
		REAL_TYPE y1 = rgn.getMax().y;
		REAL_TYPE z1 = rgn.getMax().z;
		REAL_TYPE ox = 0.5*(box.getMin().x + box.getMax().x);
		REAL_TYPE oy = 0.5*(box.getMin().y + box.getMax().y);
		REAL_TYPE oz = 0.5*(box.getMin().z + box.getMax().z);
		REAL_TYPE r = 0.5*(box.getMax().x - box.getMin().x);
		REAL_TYPE r2 = r*r;
    REAL_TYPE x = 0.0;
    if      (ox < x0) x = x0 - ox;
    else if (ox > x1) x = ox - x1;
    REAL_TYPE r2min = x * x;
    if (r2min >= r2) return true;

    REAL_TYPE y = 0.0;
    if      (oy < y0) y = y0 - oy;
    else if (oy > y1) y = oy - y1;
    r2min += y * y;
    if (r2min >= r2) return true;

    REAL_TYPE z = 0.0;
    if      (oz < z0) z = z0 - oz;
    else if (oz > z1) z = oz - z1;
    r2min += z * z;
    if (r2min >= r2) return true;

    return false;
  }
};

#endif // BB_DIVIDER_H
