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
/// @file NeighborInfo.h
/// @brief 隣接情報クラス
///

#ifndef NEIGHBOR_INFO_H
#define NEIGHBOR_INFO_H

#include "BCMTools.h"
#include "mpi.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

/// 隣接情報クラス.
class NeighborInfo {
  
  /// 隣接ブロックID.
  /// (隣接ブロックが存在しない場合は-1を入れる)
  int neighborID[NUM_SUBFACE];

  /// 隣接ブロック所属ランク.
  /// (周期境界以外の外部境界面の場合はMPI::PROC_NULLを入れる)
  int neighborRank[NUM_SUBFACE];

  /// 外部境界フラグ.
  /// (周期境界の場合もtrue)
  bool outerBoundary;

  /// 隣接ブロックとのレベル差(-1, 0, +1).
  /// (隣のレベル - 自分のレベル)
  int levelDiffarence;

  /// 隣接する相手のサブブロック番号.
  /// (levelDiffarence = -1 以外の場合は0を入れる)
  int neighborSubface;

public:

  /// コンストラクタ.
  ///
  ///  @note 初期値は, neighborID[i]=-1(隣接ブロックなし),
  ///        neighborRank[i]=MPI::PROC_NULL(隣接ブロックなし),
  ///        outBoundary=false(内部境界),
  ///        levelDiffarence=0(レベル差なし),
  ///        neighborSubface=0(サブブロック0と隣接)
  ///
  NeighborInfo() {
    for (int i = 0; i < NUM_SUBFACE; i++) {
      neighborID[i] = -1;
      neighborRank[i] = MPI::PROC_NULL;
    }
    outerBoundary = false;
    levelDiffarence = 0;
    neighborSubface = 0;
  }

  /// デストラクタ.
  ~NeighborInfo() {}

  /// レベル差を設定.
  void setLevelDifference(int dLevel) {
    assert(-1 <= dLevel && dLevel <= 1);
    levelDiffarence = dLevel;
  }

  /// レベル差を取得.
  int getLevelDifference() const { return levelDiffarence; }

  /// 隣接ブロックIDを設定.
  void setID(int id) { 
    assert(levelDiffarence == 0 || levelDiffarence == -1);
    neighborID[0] = id;
  }

  /// 隣接ブロックIDを設定.
  void setID(Subface subface, int id) { 
    assert(levelDiffarence == 1);
    neighborID[subface] = id;
  }

  /// 隣接ブロックIDを取得.
  int getID() const {
    assert(levelDiffarence == 0 || levelDiffarence == -1);
    return neighborID[0];
  }

  /// 隣接ブロックIDを取得.
  int getID(Subface subface) const {
    return neighborID[subface];
  }

  /// 隣接ブロックランクを設定.
  void setRank(int rank) { 
    assert(levelDiffarence == 0 || levelDiffarence == -1);
    neighborRank[0] = rank;
  }

  /// 隣接ブロックランクを設定.
  void setRank(Subface subface, int rank) { 
    assert(levelDiffarence == 1);
    neighborRank[subface] = rank;
  }

  /// 隣接ブロックランクを取得.
  int getRank() const {
    assert(levelDiffarence == 0 || levelDiffarence == -1);
    return neighborRank[0];
  }

  /// 隣接ブロックランクを取得.
  int getRank(Subface subface) const {
    return neighborRank[subface];
  }

  /// 隣接ブロックのサブフェイス番号を設定.
  void setNeighborSubface(Subface subface) {
    assert(levelDiffarence == -1);
    neighborSubface = subface;
  }

  /// 隣接ブロックのサブフェイス番号を取得.
  Subface getNeighborSubface() const {
  //assert(levelDiffarence == -1);
    return Subface(neighborSubface);
  }

  /// 外部境界フラグを(オンに)設定.
  void setOuterBoundary(bool flag = true) {
    outerBoundary = flag;
  }

  /// 外部境界であるか確認.
  bool isOuterBoundary() const {
    return outerBoundary; 
  }

  /// 隣接ブロックが存在するか(内部境界or周期境界)確認.
  bool exists() const {
    return neighborID[0] >= 0;
  }

  /// デバッグ情報出力.
  void print() const {
    std::cout << "levelDIff=" << levelDiffarence 
              << ", neighborSubface=" << neighborSubface << std::endl;
    std::cout << "  (ID,rank)=";
    for (int i = 0; i < NUM_SUBFACE; i++) {
      std::cout << "(" << neighborID[i] << "," << neighborRank[i] << ")";
    }
    std::cout << std:: endl;
  }

  /// 指定した隣接面を含む隣接ブロックにおける子ブロック番号を取得.
//static int getNeighborChild(Face face, Subface subface) {
  static int getNeighborChildId(Face face, Subface subface) {
    const int neighborChild_XM[] = { 1, 3, 5, 7 };
    const int neighborChild_XP[] = { 0, 2, 4, 6 };
    const int neighborChild_YM[] = { 2, 6, 3, 7 };
    const int neighborChild_YP[] = { 0, 4, 1, 5 };
    const int neighborChild_ZM[] = { 4, 5, 6, 7 };
    const int neighborChild_ZP[] = { 0, 1, 2, 3 };
    switch (face) {
      case X_M: return neighborChild_XM[subface];
      case X_P: return neighborChild_XP[subface];
      case Y_M: return neighborChild_YM[subface];
      case Y_P: return neighborChild_YP[subface];
      case Z_M: return neighborChild_ZM[subface];
      case Z_P: return neighborChild_ZP[subface];
      default: Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }

  /// 指定した接触面における子ブロックのSubface番号を取得.
  static Subface childIdToSubface(Face face, int childId) {
    const int subface_X[] = { 0, 0, 1, 1, 2, 2, 3, 3 };
    const int subface_Y[] = { 0, 2, 0, 2, 1, 3, 1, 3 };
    const int subface_Z[] = { 0, 1, 2, 3, 0, 1, 2, 3 };
    assert(0 <= childId && childId < 8);
    switch (face) {
      case X_M: return Subface(subface_X[childId]);
      case X_P: return Subface(subface_X[childId]);
      case Y_M: return Subface(subface_Y[childId]);
      case Y_P: return Subface(subface_Y[childId]);
      case Z_M: return Subface(subface_Z[childId]);
      case Z_P: return Subface(subface_Z[childId]);
      default: Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }

  /// Subface番号を対面ブロックのものに変換.
  static Face reverseFace(Face face) {
    switch (face) {
      case X_M: return X_P;
      case X_P: return X_M;
      case Y_M: return Y_P;
      case Y_P: return Y_M;
      case Z_M: return Z_P;
      case Z_P: return Z_M;
      default: Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }


};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // NEIGHBOR_INFO_H
