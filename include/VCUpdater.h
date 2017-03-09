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
/// @file VCUpdater.h
/// @brief 仮想セルアップデータ基底クラス.
///

#ifndef VC_UPDATER_H
#define VC_UPDATER_H

#include "mpi.h"
#include "BCMTools.h"
#include "PointerSetter.h"
#include "NeighborInfo.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

class DataClass;

/// 仮想セルアップデータ基底クラス.
class VCUpdater {

  const NeighborInfo* neighborInfo;  ///< 隣接情報配列

  const MPI::Comm& comm;  ///< MPIコミュニケータ

  int myrank;  ///< 所属ランク番号

public:

  /// コンストラクタ.
  ///
  ///  @param[in] neighborInfo 隣接情報配列
  ///  @param[in] comm MPIコミュニケータ(ディフォルトMPI::COMM_WORLD)
  ///
  VCUpdater(const NeighborInfo* neighborInfo, const MPI::Comm& comm = MPI::COMM_WORLD)
    : neighborInfo(neighborInfo), comm(comm) {
    myrank = comm.Get_rank();
  }

  /// デストラクタ.
  virtual ~VCUpdater() {}

  /// 仮想セル同期を行うデータクラスの登録.
  virtual void setDataClass(DataClass* dataClass) = 0;

  /// 隣接データクラスから仮想セルデータをコピー.
  void copyFromNeighbor(Face face, Subface subface) {
    assert(neighborInfo[face].getRank(subface) == myrank);
    if (neighborInfo[face].getLevelDifference() == 0) {
      assert(subface == Subface(0));
      copyFromNeighbor(face);
    }
    else if (neighborInfo[face].getLevelDifference() == 1) {
      copyFromNeighborF2C(face, subface);
    }
    else if (neighborInfo[face].getLevelDifference() == -1) {
      assert(subface == Subface(0));
      Subface subface = neighborInfo[face].getNeighborSubface();
      copyFromNeighborC2F(face, subface);
    }
    else {
      Exit(EX_FAILURE);
    }
  }

  /// 送信バッファに仮想セルデータをコピー.
  void copyToCommBuffer(Face face, Subface subface) {
    assert(neighborInfo[face].getRank(subface) != myrank);
    assert(neighborInfo[face].getRank(subface) != MPI::PROC_NULL);
    if (neighborInfo[face].getLevelDifference() == 0) {
      assert(subface == Subface(0));
      copyToCommBuffer(face);
    }
    else if (neighborInfo[face].getLevelDifference() == 1) {
      copyToCommBufferC2F(face, subface);
    }
    else if (neighborInfo[face].getLevelDifference() == -1) {
      assert(subface == Subface(0));
      Subface subface = neighborInfo[face].getNeighborSubface();
      copyToCommBufferF2C(face, subface);
    }
    else {
      Exit(EX_FAILURE);
    }
  }

  /// 受信バッファから仮想セルデータをコピー.
  void copyFromCommBuffer(Face face, Subface subface) {
    assert(neighborInfo[face].getRank(subface) != myrank);
    assert(neighborInfo[face].getRank(subface) != MPI::PROC_NULL);
    if (neighborInfo[face].getLevelDifference() == 0) {
      assert(subface == Subface(0));
      copyFromCommBuffer(face);
    }
    else if (neighborInfo[face].getLevelDifference() == 1) {
      copyFromCommBufferF2C(face, subface);
    }
    else if (neighborInfo[face].getLevelDifference() == -1) {
      assert(subface == Subface(0));
      Subface subface = neighborInfo[face].getNeighborSubface();
      copyFromCommBufferC2F(face, subface);
    }
    else {
      Exit(EX_FAILURE);
    }
  }

  /// 隣接データクラスの登録.
  void setNeighbor(Face face, DataClass* dataClass) {
    setNeighbor(face, Subface(0), dataClass);
  }

  /// 隣接データクラスの登録.
  virtual void setNeighbor(Face face, Subface subface, DataClass* dataClass) = 0;

  /// 隣接データクラスの登録解除.
  void clearNeighbor() {
    for (int i = 0; i < NUM_FACE; i++) {
      Face face = Face(i);
      for (int j = 0; j < NUM_SUBFACE; j++) {
        Subface subface = Subface(j);
        clearNeighbor(face, subface);
      }
    }
  }

  /// 隣接データクラスの登録解除.
  virtual void clearNeighbor(Face face, Subface subface) = 0;

  /// 仮想セル同期データ送信に必要なバッファサイズを取得.
  size_t getSendBufferByteSize(Face face, Subface subface) const {
    assert(neighborInfo[face].getRank(subface) != myrank);
    assert(neighborInfo[face].getRank(subface) != MPI::PROC_NULL);
    if (neighborInfo[face].getLevelDifference() == 0) {
      assert(subface == Subface(0));
      return getSendBufferByteSize(face);
    }
    else if (neighborInfo[face].getLevelDifference() == 1) {
      return getSendBufferByteSizeC2F(face, subface);
    }
    else if (neighborInfo[face].getLevelDifference() == -1) {
      assert(subface == Subface(0));
      Subface subface = neighborInfo[face].getNeighborSubface();
      return getSendBufferByteSizeF2C(face, subface);
    }
    else {
      Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }

  /// 仮想セル同期データ受信に必要なバッファサイズを取得.
  size_t getRecvBufferByteSize(Face face, Subface subface) const {
    assert(neighborInfo[face].getRank(subface) != myrank);
    assert(neighborInfo[face].getRank(subface) != MPI::PROC_NULL);
    if (neighborInfo[face].getLevelDifference() == 0) {
      assert(subface == Subface(0));
      return getRecvBufferByteSize(face);
    }
    else if (neighborInfo[face].getLevelDifference() == 1) {
      return getRecvBufferByteSizeF2C(face, subface);
    }
    else if (neighborInfo[face].getLevelDifference() == -1) {
      assert(subface == Subface(0));
      Subface subface = neighborInfo[face].getNeighborSubface();
      return getRecvBufferByteSizeC2F(face, subface);
    }
    else {
      Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }

  /// 仮想セル同期データ送信バッファ用PointerSetterオブジェクトを取得.
  virtual PointerSetterBase* getSendBufferPointerSetter(Face face, Subface subface) = 0;

  /// 仮想セル同期データ受信バッファ用PointerSetterオブジェクトを取得.
  virtual PointerSetterBase* getRecvBufferPointerSetter(Face face, Subface subface) = 0;

private:

  /// 隣接データクラスから仮想セルデータをコピー(同レベル間).
  virtual void copyFromNeighbor(Face face) = 0;

  /// 隣接データクラスから仮想セルデータをコピー(レベルL+1→L).
  virtual void copyFromNeighborF2C(Face face, Subface subface) = 0;

  /// 隣接データクラスから仮想セルデータをコピー(レベルL→L+1).
  virtual void copyFromNeighborC2F(Face face, Subface subface) = 0;


  /// 送信バッファに仮想セルデータをコピー(同レベル間).
  virtual void copyToCommBuffer(Face face) = 0;

  /// 送信バッファに仮想セルデータをコピー(レベルL+1→L).
  virtual void copyToCommBufferF2C(Face face, Subface subface) = 0;

  /// 送信バッファに仮想セルデータをコピー(レベルL→L+1).
  virtual void copyToCommBufferC2F(Face face, Subface subface) = 0;


  /// 受信バッファから仮想セルデータをコピー(同レベル間).
  virtual void copyFromCommBuffer(Face face) = 0;

  /// 受信バッファから仮想セルデータをコピー(レベルL+1→L).
  virtual void copyFromCommBufferF2C(Face face, Subface subface) = 0;

  /// 受信バッファから仮想セルデータをコピー(レベルL→L+1).
  virtual void copyFromCommBufferC2F(Face face, Subface subface) = 0;


  /// 仮想セル同期データ送信に必要なバッファサイズを取得(同レベル間).
  virtual size_t getSendBufferByteSize(Face face) const = 0;

  /// 仮想セル同期データ送信に必要なバッファサイズを取得(レベルL+1→L).
  virtual size_t getSendBufferByteSizeF2C(Face face, Subface subface) const = 0;

  /// 仮想セル同期データ送信に必要なバッファサイズを取得(レベルL→L+1).
  virtual size_t getSendBufferByteSizeC2F(Face face, Subface subface) const = 0;


  /// 仮想セル同期データ受信に必要なバッファサイズを取得(同レベル間).
  virtual size_t getRecvBufferByteSize(Face face) const = 0;

  /// 仮想セル同期データ受信に必要なバッファサイズを取得(レベルL+1→L).
  virtual size_t getRecvBufferByteSizeF2C(Face face, Subface subface) const = 0;

  /// 仮想セル同期データ受信に必要なバッファサイズを取得(レベルL→L+1).
  virtual size_t getRecvBufferByteSizeC2F(Face face, Subface subface) const = 0;

protected:

  /// サブフェイスの第1座標の起点(返り値が1の場合は中点から).
  static int subfaceOrigin0(Subface subface) {
    switch (subface) {
      case SF_00: return 0;
      case SF_01: return 1;
      case SF_10: return 0;
      case SF_11: return 1;
      default: Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }

  /// サブフェイスの第2座標の起点(返り値が1の場合は中点から).
  static int subfaceOrigin1(Subface subface) {
    switch (subface) {
      case SF_00: return 0;
      case SF_01: return 0;
      case SF_10: return 1;
      case SF_11: return 1;
      default: Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }

};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // VC_UPDATER_H
