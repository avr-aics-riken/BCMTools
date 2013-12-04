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
/// @file Pedigree.h
/// @brief Octree用Pedigreeクラス
/// 

#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <stdint.h>  /* for uint64_t */
#include <iostream>
#include "BCMTools.h"


/// Pedigreeクラス.
///
///  @code
///   64ビットによる実装.
///   x:      [63:48] (16bit)
///   y:      [47:32] (16bit)
///   z:      [31:16] (16bit)
///   rootID: [15:4]  (12bit)
///   level:  [3:0]   ( 4bit)
///  @endcode
///
class Pedigree {

public:

  const static unsigned MaxLevel = 0xf; ///< 最大レベル(4ビット)
  const static unsigned MaxRootID = 0xfff; ///< 最大ルート数(12ビット)
  const static unsigned MaxCoord = 0xffff; ///< 最大座標値(16ビット)

private:

  uint64_t p;  ///< Pedigre格納用内部64ビットデータ

  /// 内部データにPedigree値を設定.
  ///
  ///  @param[in] level ツリーレベル
  ///  @param[in] x x位置
  ///  @param[in] y y位置
  ///  @param[in] z z位置
  ///  @param[in] rootID ルートID
  ///
  void setPedigree(unsigned level, unsigned x, unsigned y, unsigned z,
                   unsigned rootID) {
    assert(level <= MaxLevel);
    assert(rootID <= MaxRootID);
    assert(x <= MaxCoord);
    assert(y <= MaxCoord);
    assert(z <= MaxCoord);
    uint64_t xx = x;
    uint64_t yy = y;
    uint64_t zz = z;
    p = (xx << 48) + (yy << 32) + (zz << 16) + (rootID << 4) + level;
  }


public:

  /// コンストラクタ(ルートノード).
  Pedigree(unsigned rootID = 0) {
    setPedigree(0, 0, 0, 0, rootID);
  }

  /// コンストラクタ.
  ///
  ///  @param[in] level ツリーレベル
  ///  @param[in] x x位置
  ///  @param[in] y y位置
  ///  @param[in] z z位置
  ///  @param[in] rootID ルートID
  ///
  Pedigree(unsigned level, unsigned x, unsigned y, unsigned z,
           unsigned rootID = 0) {
    if (level >= 0) {
      assert(x < (1 << level));
      assert(y < (1 << level));
      assert(z < (1 << level));
    }
    setPedigree(level, x, y, z, rootID);
  }

  /// コンストラクタ(子ノード).
  ///
  ///  @param[in] parent 親ノードのPedigree
  ///  @param[in] ijk 子ノード番号(0〜7)
  ///
  Pedigree(const Pedigree& parent, unsigned ijk) {
    assert(ijk < 8);
    unsigned i = ijk & 0x01;
    unsigned j = (ijk >> 1) & 0x01;
    unsigned k = (ijk >> 2) & 0x01;
    unsigned x = parent.getX() * 2 + i;
    unsigned y = parent.getY() * 2 + j;
    unsigned z = parent.getZ() * 2 + k;
    unsigned level = parent.getLevel() + 1;
    unsigned rootID = parent.getRootID();
    setPedigree(level, x, y, z, rootID);
  }


  /// デストラクタ.
  ~Pedigree() {}

  /// ツリーレベルを取得.
  ///
  ///  @return ツリーレベル
  ///
  unsigned getLevel() const { return p & 0xf; }

  /// X方向位置を取得.
  ///
  ///   @return X方向位置
  ///
  unsigned getX() const { return (p >> 48) & 0xffff; }

  /// Y方向位置を取得.
  ///
  ///   @return Y方向位置
  ///
  unsigned getY() const { return (p >> 32) & 0xffff; }

  /// Z方向位置を取得.
  ///
  ///   @return Z方向位置
  ///
  unsigned getZ() const { return (p >> 16) & 0xffff; }

  /// ルートIDを取得.
  ///
  ///   @return ルートID
  ///
  unsigned getRootID() const { return (p >> 4) & 0xfff; }

  /// そのレベルでの最大座標値を取得
  ///
  ///   @return 最大座標値(= 2のレベル値乗)
  ///
  unsigned getUpperBound() const { return 1 << getLevel(); }  // return 2^level

  /// 指定されたレベルのX座標値を取得.
  ///
  ///   @param[in] level レベル
  ///   @return X座標値(0 or 1)
  ///
  unsigned getX(unsigned level) const {
    assert(level <= MaxLevel);
    return (p >> (getLevel() - level + 48)) & 0x01;
  }

  /// 指定されたレベルのY座標値を取得.
  ///
  ///   @param[in] level レベル
  ///   @return Y座標値(0 or 1)
  ///
  unsigned getY(unsigned level) const {
    assert(level <= MaxLevel);
    return (p >> (getLevel() - level + 32)) & 0x01;
  }

  /// 指定されたレベルのZ座標値を取得.
  ///
  ///   @param[in] level レベル
  ///   @return Z座標値(0 or 1)
  ///
  unsigned getZ(unsigned level) const {
    assert(level <= MaxLevel);
    return (p >> (getLevel() - level + 16)) & 0x01;
  }

  /// 指定されたレベルでの子ノード番号を取得.
  ///
  ///   @param[in] level レベル
  ///   @return 子ノード番号(0〜7)
  ///
  unsigned getChildId(unsigned level) const {
    return getX(level) + getY(level) * 2 + getZ(level) * 4;
  }

  /// シリアライズに必要なバイト数を取得.
  ///
  ///   @return バイト数
  ///
//  static size_t GetSerializeSize() { return sizeof(p); }
  static size_t GetSerializeSize() { return sizeof(uint64_t); }

  /// シリアライズ.
  ///
  ///   @param[out] buf シリアライズデータの出力先領域
  ///
  void serialize(void *buf) const {
    *static_cast<uint64_t*>(buf) = p;
  }

  /// デシリアライズ.
  ///
  ///   @param[in] buf シリアライズデータ
  ///
  void deserialize(const void *buf) {
    p = *static_cast<const uint64_t*>(buf);
  }

};


/// Pedigree情報のストリームへの出力.
inline std::ostream& operator<<(std::ostream& os, const Pedigree& p) {
  return os << "(" << p.getRootID() << ":" << p.getLevel() << "("
            << p.getX() << "," << p.getY() << "," << p.getZ() << "))";
}

#endif // PEDIGREE_H
