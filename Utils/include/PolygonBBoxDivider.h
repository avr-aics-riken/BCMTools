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
/// @file PolygonBBoxDivider.h
/// @brief ブロック分割判定クラス(PolygonBBoxDivider)
/// 

#ifndef POLYGON_DIVIDER_H
#define POLYGON_DIVIDER_H

#include <string>
#include <vector>
#include "BCMPolylib.h"
#include "BCMTools.h"
#include "MultiRootDivider.h"
#include "BoundingBox.h"

/// <ポリゴングループ名, 分割レベル>ペア.
struct PolygonGroupSpec {
  std::string polygonGroupName;  ///< ポリゴングループ名
  int level;  ///< 最大分割レベル

  /// ストリーム出力.
  friend std::ostream& operator<<(std::ostream& os, const PolygonGroupSpec& pgs) {
    return os << pgs.polygonGroupName << " -> " << pgs.level;
  }

  /// ストリーム入力.
  friend std::istream& operator>>(std::istream& is, PolygonGroupSpec& pgs) {
    return is >> pgs.polygonGroupName >> pgs.level;
  }
};


/// <バウンディングボックス, 分割レベル>ペアのリスト.
struct BoundingBoxSpec {
  BoundingBox boundingBox;  ///< バウンディングボックス
  int level;  ///< 最大分割レベル

  /// ストリーム出力.
  friend std::ostream& operator<<(std::ostream& os, const BoundingBoxSpec& bbs) {
    return os << bbs.boundingBox << " -> " << bbs.level;
  }

  /// ストリーム入力.
  friend std::istream& operator>>(std::istream& is, BoundingBoxSpec& bbs) {
    return is >> bbs.boundingBox >> bbs.level;
  }
};


/// 複数のポリゴンデータ境界指定によりブロック分割するDivider.
///
/// 複数のバウンディングボックス指定により，ブロック分割の調整が可能.
///
class PolygonBBoxDivider : public MultiRootDivider {

  const int minLevel;  ///< 最小分割レベル

  const PolylibNS::BCMPolylib* pl;  ///< BCMPolylibオブジェクト

  /// <ポリゴングループ名, 分割レベル>ペアのリスト
  const std::vector<PolygonGroupSpec>& polygonGroupList;

  /// <バウンディングボックス, 分割レベル>ペアのリスト
  const std::vector<BoundingBoxSpec>& boundingBoxList; 

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
  PolygonBBoxDivider(const Vec3d& origin, double rootLength, const RootGrid* rootGrid, 
                 int minLevel, const PolylibNS::BCMPolylib* pl,
                 const std::vector<PolygonGroupSpec>& polygonGroupList,
                 const std::vector<BoundingBoxSpec>& boundingBoxList, 
                 double extraMarginRatio = 0.0)
    : MultiRootDivider(origin, rootLength, rootGrid),
      minLevel(minLevel), pl(pl),
      polygonGroupList(polygonGroupList), boundingBoxList(boundingBoxList),
      extraMarginRatio(extraMarginRatio) {
  }

  /// デストラクタ.
  ~PolygonBBoxDivider() {}

  /// ブロックを分割するかどうかを判定.
  ///
  ///   @param[in] pedigree ブロックのPedigree
  ///   @return ブロックタイプ
  ///
  NodeType operator() (const Pedigree& pedigree) {

    int level = pedigree.getLevel();

    if  (level < minLevel) return BRANCH;

    NodeType ret = LEAF_ACTIVE;

    for (std::vector<BoundingBoxSpec>::const_iterator it = boundingBoxList.begin();
         it != boundingBoxList.end(); ++it) {
      int maxLevel = it->level;
      if (level < maxLevel) {
        BoundingBox box = it->boundingBox;
        BoundingBox region = defineSearchRegion(pedigree, maxLevel);
        region.setMargin(extraMarginRatio / (1 << maxLevel));
        if (box.intersects(region)) ret = BRANCH;
      }
    }


    for (std::vector<PolygonGroupSpec>::const_iterator it = polygonGroupList.begin();
         it != polygonGroupList.end(); ++it) {
      int maxLevel = it->level;
      if (level < maxLevel) {
        const std::string& polygonGroup = it->polygonGroupName;
        BoundingBox region = defineSearchRegion(pedigree, maxLevel);
        region.setMargin(extraMarginRatio / (1 << maxLevel));
        Vec3r min(region.getMin().x, region.getMin().y, region.getMin().z);
        Vec3r max(region.getMax().x, region.getMax().y, region.getMax().z);
        std::vector<PolylibNS::Triangle*>* polygonList 
                        = pl->search_polygons(polygonGroup, min, max, false);
        int nPolygon = polygonList->size();
        delete polygonList;
        if (nPolygon > 0) ret = BRANCH;
      }
    }

    return ret;
  }

};

#endif // POLYGON_DIVIDER_H
