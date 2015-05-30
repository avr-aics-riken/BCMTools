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
/// @file DataClass.h
/// @brief データクラス基底クラス
///

#ifndef DATA_CLASS_H
#define DATA_CLASS_H

#include "Vec3.h"
#include "VCUpdater.h"
#include "PointerSetter.h"
#include "BCMTools.h"

using namespace Vec3class;

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

/// データクラス基底クラス.
class DataClass {

protected:

  Vec3i size;  ///< セル分割数
  int vc;      ///< 仮想セル幅

public:

  /// コンストラクタ.
  ///
  ///  @param[in] size 分割数
  ///  @param[in] vc 仮想セル幅
  ///
  DataClass(const Vec3i& size, int vc) : size(size), vc(vc) {
    if (size.x / 2 < vc || size.y / 2 < vc || size.z / 2 < vc) {
      std::cout << "*** error: cell size must be > 2*vc" << std::endl;
      Exit(EX_FAILURE);
    }
  }

  /// デストラクタ.
  virtual ~DataClass() {}

  /// セル分割数を取得.
  const Vec3i& getSize() const { return size; }

  /// x方向セル分割数を取得.
  int getSizeX() const { return size[0]; }

  /// y方向セル分割数を取得.
  int getSizeY() const { return size[1]; }

  /// z方向セル分割数を取得.
  int getSizeZ() const { return size[2]; }

  /// 仮想セル幅を取得.
  int getVCSize() const { return vc; }

};


/// 仮想セル同期可能なデータクラス(仮想クラス).
///
/// @note 仮想セル同期計算は仮想セルアップデータに丸投げ
///
class UpdatableDataClass : public DataClass {

  VCUpdater* updater;   ///< 仮想セルアップデータ

public:

  /// コンストラクタ.
  ///
  ///  @param[in] size 分割数
  ///  @param[in] vc 仮想セル幅
  ///
  ///  @note upaterは0に初期化する(仮想セル同期をしない使い方もOK)
  ///
  UpdatableDataClass(const Vec3i& size, int vc)
    : DataClass(size, vc), updater(0) {}

  /// デストラクタ.
  ///
  ///   @note 仮想セルアップデータも(登録されていたなら)解放する
  ///
  virtual ~UpdatableDataClass() { delete updater; }

  /// 仮想セルアップデータの登録.
  ///
  ///  @param[in] updater 仮想セルアップデータ
  ///
  void setUpdater(VCUpdater* updater) {
    this->updater = updater;
    this->updater->setDataClass(this);
  }

  /// 同並列計算ノード内の隣接データクラスを登録.
  void setNeighbor(Face face, Subface subface, DataClass* dataClass) {
    checkVCUpdater();
    updater->setNeighbor(face, subface, dataClass);
  }

  /// 同並列計算ノード内の隣接データクラスを登録.
  void setNeighbor(Face face, DataClass* dataClass) {
    setNeighbor(face, Subface(0), dataClass);
  }

  /// 仮想セル同期データ送信に必要なバッファサイズを取得.
  size_t getSendBufferByteSize(Face face, Subface subface) const {
    checkVCUpdater();
    return updater->getSendBufferByteSize(face, subface);
  }

  /// 仮想セル同期データ受信に必要なバッファサイズを取得.
  size_t getRecvBufferByteSize(Face face, Subface subface) const {
    checkVCUpdater();
    return updater->getRecvBufferByteSize(face, subface);
  }

  /// 仮想セル同期データ送信バッファ用PointerSetterオブジェクトを取得.
  PointerSetterBase* getSendBufferPointerSetter(Face face, Subface subface) {
    checkVCUpdater();
    return updater->getSendBufferPointerSetter(face, subface);
  }

  /// 仮想セル同期データ受信バッファ用PointerSetterオブジェクトを取得.
  PointerSetterBase* getRecvBufferPointerSetter(Face face, Subface subface) {
    checkVCUpdater();
    return updater->getRecvBufferPointerSetter(face, subface);
  }

  /// 隣接データクラスから仮想セルデータをコピー.
  void copyFromNeighbor(Face face, Subface subface) {
    checkVCUpdater();
    updater->copyFromNeighbor(face, subface);
  }

  /// 送信バッファに仮想セルデータをコピー.
  void copyToCommBuffer(Face face, Subface subface) {
    checkVCUpdater();
    updater->copyToCommBuffer(face, subface);
  }

  /// 受信バッファから仮想セルデータをコピー.
  void copyFromCommBuffer(Face face, Subface subface) {
    checkVCUpdater();
    updater->copyFromCommBuffer(face, subface);
  }

private:

  /// 仮想セルアップデータが登録されているか確認.
  void checkVCUpdater() const {
    if (updater == 0) {
      std::cout << "*** error: UpdatableDataClass has no VCUpdater" << std::endl;
      Exit(EX_FAILURE);
    }
  }

};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // DATA_CLASS_H
