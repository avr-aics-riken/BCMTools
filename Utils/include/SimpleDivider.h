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
/// @file SimpleDivider.h
/// @brief ブロック分割判定クラス(FlatDivider, SimpleDivider)
///

#ifndef SIMPLE_DIVIDIER_H
#define SIMPLE_DIVIDIER_H

#include "Divider.h"
#include "BCMTools.h"

/// 均一にブロック分割するDivider.
///
/// ルートブロックの大きさと位置は固定
///   - 大きさ: 単位立方体
///   - 位置: 最初のルートブロックの原点が(0,0,0)
///
class FlatDivider : public Divider {

  const RootGrid* rootGrid;  ///< ルートブロック配置
  int maxLevel;  ///< 分割レベル

public:

  /// コンストラクタ.
  ///
  ///  @param[in] rootGrid ルートブロック配置情報
  ///  @param[in] maxLevel 分割レベル
  ///
  FlatDivider(const RootGrid* rootGrid, int maxLevel)
    : rootGrid(rootGrid), maxLevel(maxLevel) {}

  /// デストラクタ.
  ~FlatDivider() {}

  /// ブロックを分割するかどうかを判定.
  ///
  ///   @param[in] pedigree ブロックのPedigree
  ///   @return ブロックタイプ
  ///
  NodeType operator() (const Pedigree& pedigree) {
    return  (pedigree.getLevel() >= maxLevel) ? LEAF_ACTIVE : BRANCH;
  }
};


/// 指定された外部境界面に接するブロックを分割するDivider.
///
/// ルートブロックの大きさと位置は固定
///   @li 大きさ: 単位立方体
///   @li 位置: 最初のルートブロックの原点が(0,0,0)
///
class SimpleDivider : public Divider {

public:

  /// 外部境界面を指定する型(ビット和による合成が可能)
  typedef unsigned char Direction;

  static const Direction NONE = 0;  ///< 指定面なし
  static const Direction XM = 1;   ///< -X面
  static const Direction XP = 2;   ///< +X面
  static const Direction YM = 4;   ///< -Y面
  static const Direction YP = 8;   ///< +Y面
  static const Direction ZM = 16;  ///< -Z面
  static const Direction ZP = 32;  ///< +Z面
  static const Direction ALL = XM | XP | YM | YP | ZM | ZP; ///< 全外部境界面

private:

  const RootGrid* rootGrid;  ///< ルートブロック配置
  const int minLevel;  ///< 最小分割レベル
  const int maxLevel;  ///< 最小分割レベル

  bool divide_xm;  ///< -X方向分割フラグ
  bool divide_xp;  ///< +X方向分割フラグ
  bool divide_ym;  ///< -Y方向分割フラグ
  bool divide_yp;  ///< +Y方向分割フラグ
  bool divide_zm;  ///< -Z方向分割フラグ
  bool divide_zp;  ///< +Z方向分割フラグ

public:

  /// コンストラクタ.
  ///
  ///  @param[in] rootGrid ルートブロック配置情報
  ///  @param[in] minLevel 最小分割レベル
  ///  @param[in] maxLevel 最大分割レベル
  ///  @param[in] direction 分割を進めていく向き(外部境界面)
  ///
  SimpleDivider(const RootGrid* rootGrid, int minLevel, int maxLevel,
                Direction direction = ALL)
    : rootGrid(rootGrid), minLevel(minLevel), maxLevel(maxLevel) {
    assert(minLevel >= 0);
    assert(maxLevel >= minLevel);
    divide_xm = direction & XM;
    divide_xp = direction & XP;
    divide_ym = direction & YM;
    divide_yp = direction & YP;
    divide_zm = direction & ZM;
    divide_zp = direction & ZP;
  }

  /// デストラクタ.
  ~SimpleDivider() {}

  /// ブロックを分割するかどうかを判定.
  ///
  ///   @param[in] pedigree ブロックのPedigree
  ///   @return ブロックタイプ
  ///
  NodeType operator() (const Pedigree& pedigree) {
    int level = pedigree.getLevel();

    if  (level < minLevel) return BRANCH;
    if  (level >= maxLevel) return LEAF_ACTIVE;

    int x = pedigree.getX();
    int y = pedigree.getY();
    int z = pedigree.getZ();
    int rootID = pedigree.getRootID();

    int max0 = pedigree.getUpperBound() - 1;  // 2^level - 1

    if (divide_xm && x == 0 && rootGrid->isOuterBoundary(rootID, X_M)) return BRANCH;
    if (divide_ym && y == 0 && rootGrid->isOuterBoundary(rootID, Y_M)) return BRANCH;
    if (divide_zm && z == 0 && rootGrid->isOuterBoundary(rootID, Z_M)) return BRANCH;
    if (divide_xp && x == max0 && rootGrid->isOuterBoundary(rootID, X_P)) return BRANCH;
    if (divide_yp && y == max0 && rootGrid->isOuterBoundary(rootID, Y_P)) return BRANCH;
    if (divide_zp && z == max0 && rootGrid->isOuterBoundary(rootID, Z_P)) return BRANCH;

    return LEAF_ACTIVE;
  }
};

#endif // SIMPLE_DIVIDER_H
