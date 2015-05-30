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
/// @file Block.h
/// @brief ブロック(サンプル)サンプル
/// 

#ifndef BLOCK_H
#define BLOCK_H

#include "BlockBase.h"
#include "BoundaryInfo.h"


/// ブロック(サンプル)クラス.
class Block : public BlockBase {

  BoundaryInfo* boundaryInfo;  ///< 境界条件情報配列

public:

  /// コンストラクタ.
  ///
  ///  @param[in] size セル分割数
  ///  @param[in] origin 原点座標
  ///  @param[in] blockSize ブロックサイズ
  ///  @param[in] level ツリーレベル
  ///  @param[in] neighborInfo 隣接情報配列
  ///  @param[in] boundaryInfo 境界条件情報配列
  ///
  ///  @note セル分割数は偶数であること
  ///
  Block(const Vec3i& size, const Vec3d& origin, const Vec3d& blockSize,
        int level, NeighborInfo* neighborInfo, BoundaryInfo* boundaryInfo) 
   : BlockBase(size, origin, blockSize, level, neighborInfo),
     boundaryInfo(boundaryInfo) {}

  /// デストラクタ.
  virtual ~Block() {
    delete[] boundaryInfo;
  }

  /// 境界条件情報配列を取得.
  const BoundaryInfo* getBoundaryInfo() const { return boundaryInfo; }

};

#endif // BLOCK_H
