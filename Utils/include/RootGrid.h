/*
 * BCMTools
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file RootGrid.h
/// @brief マルチルートOctree用のルートブロック配置管理クラス
/// 

#ifndef ROOT_GRID_H
#define ROOT_GRID_H

#include "BCMTools.h"
#include "Vec3.h"
#include "mpi.h"

using namespace Vec3class;

/// マルチルートOctree用のルートブロック配置管理クラス.
///
///  @note 位置(i,j,k)のルートブロックのIDは，i + nx*j + nx*ny*k
///
class RootGrid {

  int nx;  ///< X方向ルート数
  int ny;  ///< Y方向ルート数
  int nz;  ///< Z方向ルート数

  bool periodicX;  ///< X方向周期境界条件フラグ
  bool periodicY;  ///< Y方向周期境界条件フラグ
  bool periodicZ;  ///< Z方向周期境界条件フラグ

public:

  /// コンストラクタ.
  ///
  ///  @param[in] nx X方向ルート数
  ///  @param[in] ny Y方向ルート数
  ///  @param[in] nz Z方向ルート数
  ///
  RootGrid(int nx, int ny, int nz) : nx(nx), ny(ny), nz(nz),
    periodicX(false), periodicY(false), periodicZ(false) {}

  /// コンストラクタ.
  ///
  ///  @param[in] n ルート数ベクトル
  ///
  RootGrid(const ::Vec3i& n) : nx(n.x), ny(n.y), nz(n.z),
    periodicX(false), periodicY(false), periodicZ(false) {}

  /// デストラクタ.
  ~RootGrid() {}

  /// ルート総数を取得.
  ///
  ///   @return ルート総数
  ///
  int getSize() const { return nx*ny*nz; }

  /// X方向ルート数を取得.
  ///
  ///   @return X方向ルート数
  ///
  int getSizeX() const { return nx; }

  /// Y方向ルート数を取得.
  ///
  ///   @return Y方向ルート数
  ///
  int getSizeY() const { return ny; }

  /// Z方向ルート数を取得.
  ///
  ///   @return Z方向ルート数
  ///
  int getSizeZ() const { return nz; }

  /// X方向に周期境界条件を設定.
  void setPeriodicX() { periodicX = true; }

  /// Y方向に周期境界条件を設定.
  void setPeriodicY() { periodicY = true; }

  /// Z方向に周期境界条件を設定.
  void setPeriodicZ() { periodicZ = true; }

  /// X方向の周期境界条件を解除.
  void clearPeriodicX() { periodicX = false; }

  /// Y方向の周期境界条件を解除.
  void clearPeriodicY() { periodicY = false; }

  /// Z方向の周期境界条件を解除.
  void clearPeriodicZ() { periodicZ = false; }

  /// X方向インデクスを取得.
  ///
  ///  @param[in] rootID ルートID
  ///  @return X方向インデクス
  ///
  int rootID2indexX(int rootID) const { return rootID % nx; }

  /// Y方向インデクスを取得.
  ///
  ///  @param[in] rootID ルートID
  ///  @return Y方向インデクス
  ///
  int rootID2indexY(int rootID) const { return (rootID / nx) % ny; }

  /// Z方向インデクスを取得.
  ///
  ///  @param[in] rootID ルートID
  ///  @return Z方向インデクス
  ///
  int rootID2indexZ(int rootID) const { return (rootID / nx) / ny; }

  /// 位置インデクスをルートIDに変換.
  ///
  ///  @param[in] ix X方向インデクス
  ///  @param[in] iy Y方向インデクス
  ///  @param[in] iz Z方向インデクス
  ///  @return ルートID
  ///
  int index2rootID(int ix, int iy, int iz) const { return ix + nx*iy + nx*ny*iz; }

  /// 隣接するルートのIDを返す.
  ///
  ///  @param[in] rootID ルートID
  ///  @param[in] face 隣接面
  ///  @return 隣接するルートのルートID
  ///
  ///  @note 隣接ルートが存在しない場合は-1を返す.
  ///
  int getNeighborRoot(int rootID, Face face) const {
    assert(0 <= rootID && rootID < nx*ny*nz);
    int ix = rootID2indexX(rootID);
    int iy = rootID2indexY(rootID);
    int iz = rootID2indexZ(rootID);

    switch (face) {
      case X_M:
        ix--;
        if (ix < 0) {
          if (periodicX) ix += nx;
          else return -1;
        }
        return index2rootID(ix, iy, iz);
      case X_P:
        ix++;
        if (ix >= nx) {
          if (periodicX) ix -= nx;
          else return -1;
        }
        return index2rootID(ix, iy, iz);
      case Y_M:
        iy--;
        if (iy < 0) {
          if (periodicY) iy += ny;
          else return -1;
        }
        return index2rootID(ix, iy, iz);
      case Y_P:
        iy++;
        if (iy >= ny) {
          if (periodicY) iy -= ny;
          else return -1;
        }
        return index2rootID(ix, iy, iz);
        break;
      case Z_M:
        iz--;
        if (iz < 0) {
          if (periodicZ) iz += nz;
          else return -1;
        }
        return index2rootID(ix, iy, iz);
      case Z_P:
        iz++;
        if (iz >= nz) {
          if (periodicZ) iz -= nz;
          else return -1;
        }
        return index2rootID(ix, iy, iz);
      default:
        Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }

  /// 指定した面が外部境界かどうかチェック.
  ///
  ///  @param[in] rootID ルートID
  ///  @param[in] face 隣接面
  ///  @return 指定した面が外部境界ならtrue
  ///
  bool isOuterBoundary(int rootID, Face face) const {
    assert(0 <= rootID && rootID < nx*ny*nz);
    switch (face) {
      case X_M:
        return rootID2indexX(rootID) == 0 ? true : false;
      case X_P:
        return rootID2indexX(rootID) == nx - 1 ? true : false;
      case Y_M:
        return rootID2indexY(rootID) == 0 ? true : false;
      case Y_P:
        return rootID2indexY(rootID) == ny - 1 ? true : false;
      case Z_M:
        return rootID2indexZ(rootID) == 0 ? true : false;
      case Z_P:
        return rootID2indexZ(rootID) == nz - 1 ? true : false;
      default:
        Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }

  /// ルート配置情報を他プロセスにブロードキャスト.
  ///
  ///  @param[in] comm MPIコミュニケータ
  ///
  void broadcast(MPI::Intracomm& comm = MPI::COMM_WORLD) {
    int buf[4];
    buf[0] = nx;
    buf[1] = ny;
    buf[2] = nz;
    buf[3] = 0;
    if (periodicX) buf[3] += 1;
    if (periodicY) buf[3] += 2;
    if (periodicZ) buf[3] += 4;
    comm.Bcast(buf, 4, MPI::INT, 0);
  }

  /// ランク0からルート配置情報を受信.
  ///
  ///  @param[in] comm MPIコミュニケータ
  ///
  static RootGrid* ReceiveFromMaster(MPI::Intracomm& comm = MPI::COMM_WORLD) {
    int buf[4];
    comm.Bcast(buf, 4, MPI::INT, 0);
    RootGrid* rootGrid = new RootGrid(buf[0], buf[1], buf[2]);
    if (buf[3] & 1) rootGrid->setPeriodicX();
    if (buf[3] & 2) rootGrid->setPeriodicY();
    if (buf[3] & 4) rootGrid->setPeriodicZ();
    return rootGrid;
  }

};

#endif // ROOT_GRID_H
