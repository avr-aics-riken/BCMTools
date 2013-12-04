///
/// @file PolygonDivider.h
/// @brief ブロック分割判定クラス(PolygonDivider)
/// 

#ifndef POLYGON_DIVIDER_H
#define POLYGON_DIVIDER_H

#include <string>
#include <vector>
#include "Polylib.h"

#include "BCMTools.h"
#include "MultiRootDivider.h"
#include "BoundingBox.h"

/// ポリゴンデータ境界指定によりブロック分割するDivider.
///
/// @note 使用するPolylibにはsearch nearestパッチを当てておくことが必要．
///
/// @note 使用するSTLファイル内のポリゴンデータは，
/// 法線ベクトルの向きで定義される境界面の表裏に整合性があることが必要．
/// 本Dividerは，境界の裏面領域にあるブロックは非アクティブと判定する．
///
class PolygonDivider : public MultiRootDivider {

  const int minLevel;  ///< 最小分割レベル
  const int maxLevel;  ///< 最大分割レベル

  const PolylibNS::Polylib* polylib;  ///< Polylibオブジェクト

  const std::string& polygonGroup;  ///< 境界のポリゴングループ名

  bool polygonInsideOut;  ///< 境界面の表裏反転フラグ

  Vec3r margin;  ///< 追加マージン

public:

  /// コンストラクタ.
  ///
  ///  @param[in] origin 最初のルートブロックの原点位置
  ///  @param[in] rootLength ルートブロックの辺長
  ///  @param[in] rootGrid ルートブロック配置情報
  ///  @param[in] minLevel 最小分割レベル
  ///  @param[in] maxLevel 最大分割レベル
  ///  @param[in] polylib  Polylibオブジェクト
  ///  @param[in] polygonGroup ポリゴングループ名
  ///  @param[in] polygonInsideOut 境界面の表裏反転フラグ
  ///  @param[in] extraMarginRatio 追加マージンの幅
  ///
  ///  @note 境界探査領域にマージンを追加する場合には，
  ///  追加マージン幅の最大分割レベルブロック辺長に対する比率を
  ///  extraMarginRatioに指定する．
  ///  例えば，仮想セル領域を境界面探査領域に追加するには
  ///  「(double)仮想セル数/ブロック内分割数」を指定する．
  ///
  PolygonDivider(const Vec3r& origin, double rootLength, const RootGrid* rootGrid, 
                 int minLevel, int maxLevel,
                 const PolylibNS::Polylib* polylib, std::string& polygonGroup,
                 bool polygonInsideOut,
                 double extraMarginRatio = 0.0)
    : MultiRootDivider(origin, rootLength, rootGrid),
      minLevel(minLevel), maxLevel(maxLevel), polylib(polylib),
      polygonGroup(polygonGroup), polygonInsideOut(polygonInsideOut) {
    assert(minLevel >= 0);
    assert(maxLevel >= minLevel);

    margin = Vec3r(extraMarginRatio / (1 << maxLevel));
  }

  /// デストラクタ.
  ~PolygonDivider() {}

  /// ブロックを分割するかどうかを判定.
  ///
  ///   @param[in] pedigree ブロックのPedigree
  ///   @return ブロックタイプ
  ///
  NodeType operator() (const Pedigree& pedigree) {

    int level = pedigree.getLevel();

    if  (level < minLevel) return BRANCH;

    BoundingBox region = defineSearchRegion(pedigree, maxLevel);

    Vec3r rmin = region.getMin() - margin;
    Vec3r rmax = region.getMax() + margin;

  //std::cout << "(x0,y0,z0) = " << rmin << std::endl;
  //std::cout << "(x1,y1,z1) = " << rmax << std::endl;

    PolylibNS::Vec3f min(rmin.x, rmin.y, rmin.z);
    PolylibNS::Vec3f max(rmax.x, rmax.y, rmax.z);

    std::vector<PolylibNS::Triangle*>* polygonList 
                                       = polylib->search_polygons(polygonGroup, min, max, false);
    int nPolygon = polygonList->size();
    delete polygonList;

    if (nPolygon > 0) {
      if (level == maxLevel) return LEAF_ACTIVE;
      else                   return BRANCH;
    } else {
      PolylibNS::Vec3f pos(0.5*(rmin.x+rmax.x),
                           0.5*(rmin.y+rmax.y),
                           0.5*(rmin.z+rmax.z));
      const PolylibNS::Triangle* t = polylib->search_nearest_polygon(polygonGroup, pos);
      assert(t);
      PolylibNS::Vec3f n = t->get_normal();
      PolylibNS::Vec3f* v = t->get_vertex();
      PolylibNS::Vec3f c((v[0][0]+v[1][0]+v[2][0])/3.0,
                         (v[0][1]+v[1][1]+v[2][1])/3.0,
                         (v[0][2]+v[1][2]+v[2][2])/3.0);
      if (polygonInsideOut) {
        return PolylibNS::dot(pos-c, n) < 0 ? LEAF_ACTIVE : LEAF_NO_ACTIVE;
      } else {
        return PolylibNS::dot(pos-c, n) > 0 ? LEAF_ACTIVE : LEAF_NO_ACTIVE;
      }
    }

  }



};



#endif // POLYGON_DIVIDER_H
