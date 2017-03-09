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
/// @file BlockBase.h
/// @brief ブロック基底クラス
///

#ifndef BLOCK_BASE_H
#define BLOCK_BASE_H

#include <vector>

#include "Vec3.h"
#include "BCMTools.h"
#include "NeighborInfo.h"
#include "DataClass.h"

using namespace Vec3class;

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

/// ブロック基底クラス.
class BlockBase{

  Vec3i size;       ///< セル分割数
  Vec3r origin;     ///< 原点座標
  Vec3r blockSize;  ///< ブロックサイズ
  Vec3r cellSize;   ///< セルサイズ

  int level;  ///< ツリーレベル

  NeighborInfo* neighborInfo;  ///< 隣接情報配列

  /// データクラスへのポインターテーブル型
  typedef std::vector<DataClass*> DataClassTable;

  DataClassTable dataClassTable;  ///< データクラスへのポインタテーブル

public:

  /// コンストラクタ.
  ///
  ///  @param[in] size セル分割数
  ///  @param[in] origin 原点座標
  ///  @param[in] blockSize ブロックサイズ
  ///  @param[in] level ツリーレベル
  ///  @param[in] neighborInfo 隣接情報配列
  ///
  ///  @note セル分割数は偶数であること
  ///
  BlockBase(const Vec3i& size, const Vec3r& origin, const Vec3r& blockSize,
            int level, NeighborInfo* neighborInfo)
   : size(size), origin(origin), blockSize(blockSize), level(level),
     neighborInfo(neighborInfo) {
    cellSize.x = blockSize.x / size.x;
    cellSize.y = blockSize.y / size.y;
    cellSize.z = blockSize.z / size.z;
    if (size.x % 2 != 0 || size.y % 2 != 0 || size.z % 2 != 0) {
      std::cout << "*** error: cell size must be even number" << std::endl;
      Exit(EX_FAILURE);
    }
  }

  /// デストラクタ.
  ///
  /// @note 隣接情報配列，登録されているデータクラスのメモリも解放する
  ///
  virtual ~BlockBase() {
    delete[] neighborInfo;
    DataClassTable::iterator it = dataClassTable.begin();
    for (; it != dataClassTable.end(); ++it) delete *it;
  }

  /// データクラスを登録.
  ///
  ///  @return データクラスID(登録済完了データクラス数-1)
  ///
  int setDataClass(DataClass* dataClass) {
    dataClassTable.push_back(dataClass);
    return dataClassTable.size() - 1;
  }

  /// セル分割を取得.
  const Vec3i& getSize() const { return size; }

  /// 座標原点を取得.
  const Vec3r& getOrigin() const { return origin; }

  /// ブロックサイズを取得.
  const Vec3r& getBlockSize() const { return blockSize; }

  /// セルサイズを取得.
  const Vec3r& getCellSize() const { return cellSize; }

  /// 登録されているデータクラスを取得.
  ///
  ///  @param[in] id データクラスID
  ///  @return データクラスへのポインタ
  ///
  ///  @note ID値が不正な場合，0が返る
  ///
  DataClass* getDataClass(int id) {
    if (id < 0 || id >= dataClassTable.size()) return 0;
    return dataClassTable[id];
  }

  /// ツリーレベルを取得.
  int getLevel() const { return level; }

  /// 隣接情報配列を取得.
  const NeighborInfo* getNeighborInfo() const { return neighborInfo; }

private:

  /// コピーコンストラクタ(コピー禁止).
  BlockBase(const BlockBase& rhs);

  /// 代入演算子(コピー禁止).
  BlockBase& operator=(const BlockBase& rhs);

};

#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // BLOCK_BASE_H
