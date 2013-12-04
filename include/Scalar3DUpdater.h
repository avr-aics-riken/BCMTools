///
/// @file Scalar3DUpdater.h
/// @brief スカラデータクラス仮想セルアップデータ
///

#ifndef SCALAR_3D_UPDATER_H
#define SCALAR_3D_UPDATER_H

#include "BCMTools.h"
#include "VCUpdater.h"
#include "Scalar3D.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif


/// スカラデータクラス仮想セルアップデータ.
///
///  @note 通信と補間の順序は，簡単のためL→L+1もL+1→Lも，
///        送信元で補間を行なってから通信．
///
///  @todo 補間計算部分をFortranで実装
///
///
template <typename T>
class Scalar3DUpdater : public VCUpdater {

private:

  Scalar3D<T>* dataClass;   ///< 仮想セル同期対象データクラス

  T* sendBuffer[NUM_FACE][NUM_SUBFACE];  ///< 送信データバッファテーブル
  T* recvBuffer[NUM_FACE][NUM_SUBFACE];  ///< 受信データバッファテーブル

  Scalar3D<T>* neighborDataClass[NUM_FACE][NUM_SUBFACE];  ///< 隣接データクラステーブル

  int nx, ny, nz, vc;

public:

  /// コンストラクタ.
  ///
  ///  @param[in] neighborInfo 隣接情報配列
  ///  @param[in] comm MPIコミュニケータ(ディフォルトMPI::COMM_WORLD)
  ///
  Scalar3DUpdater(const NeighborInfo* neighborInfo,
                  const MPI::Comm& comm = MPI::COMM_WORLD)
    : VCUpdater(neighborInfo, comm) {
    clearCommBufferPointer();
    clearNeighbor();
  }

  /// デストラクタ.
  ~Scalar3DUpdater() {}

  /// 仮想セル同期対象データクラスを登録.
  void setDataClass(DataClass* dc) {
    dataClass = dynamic_cast<Scalar3D<T>*>(dc);
    nx = dataClass->getSizeX();
    ny = dataClass->getSizeY();
    nz = dataClass->getSizeZ();
    vc = dataClass->getVCSize();
  }

  /// 仮想セル同期データ送信に必要なバッファサイズを取得(同レベル間).
  size_t getSendBufferByteSize(Face face) const {
    return sizeof(T) * getCommBufferSize(face);
  }

  /// 仮想セル同期データ送信に必要なバッファサイズを取得(レベルL+1→L).
  size_t getSendBufferByteSizeF2C(Face face, Subface subface) const {
    return sizeof(T) * getCommBufferSize(face) / 4;
  }

  /// 仮想セル同期データ送信に必要なバッファサイズを取得(レベルL→L+1).
  size_t getSendBufferByteSizeC2F(Face face, Subface subface) const {
    return sizeof(T) * getCommBufferSize(face);
  }

  /// 仮想セル同期データ受信に必要なバッファサイズを取得(同レベル間).
  size_t getRecvBufferByteSize(Face face) const {
    return sizeof(T) * getCommBufferSize(face);
  }

  /// 仮想セル同期データ受信に必要なバッファサイズを取得(レベルL+1→L).
  size_t getRecvBufferByteSizeF2C(Face face, Subface subface) const {
    return sizeof(T) * getCommBufferSize(face) / 4;
  }

  /// 仮想セル同期データ受信に必要なバッファサイズを取得(レベルL→L+1).
  size_t getRecvBufferByteSizeC2F(Face face, Subface subface) const {
    return sizeof(T) * getCommBufferSize(face);
  }

  /// 仮想セル同期データ送信バッファ用PointerSetterオブジェクトを取得.
  PointerSetterBase* getSendBufferPointerSetter(Face face, Subface subface) {
    return new PointerSetter<T>(&sendBuffer[face][subface]);
  }

  /// 仮想セル同期データ受信バッファ用PointerSetterオブジェクトを取得.
  PointerSetterBase* getRecvBufferPointerSetter(Face face, Subface subface) {
    return new PointerSetter<T>(&recvBuffer[face][subface]);
  }


public:

  /// 同並列計算ノード内の隣接データクラスを登録.
  void setNeighbor(Face face, Subface subface, DataClass* dataClass) {
    neighborDataClass[face][subface] = dynamic_cast<Scalar3D<T>*>(dataClass);
  }

  /// 隣接データクラスの登録解除.
  void clearNeighbor(Face face, Subface subface) {
    neighborDataClass[face][subface] = 0;
  }

  /// 隣接データクラスの登録解除.
  void clearNeighbor() {
    for (int i = 0; i < NUM_FACE; ++i) {
      for (int j = 0; j < NUM_SUBFACE; ++j) {
        clearNeighbor(Face(i), Subface(j));
      }
    }
  }

  /// 通信バッファテーブルのエントリをクリア.
  void clearCommBufferPointer(Face face, Subface subface) {
    sendBuffer[face][subface] = recvBuffer[face][subface] = 0;
  }

  /// 通信バッファテーブルをクリア.
  void clearCommBufferPointer() {
    for (int i = 0; i < NUM_FACE; ++i) {
      for (int j = 0; j < NUM_SUBFACE; ++j) {
        clearCommBufferPointer(Face(i), Subface(j));
      }
    }
  }

private:

  /// 通信バッファサイズを計算.
  size_t getCommBufferSize(Face face) const {
    switch (face) {
      case X_M:
      case X_P:
        return ny * nz * vc;
      case Y_M:
      case Y_P:
        return nz * nx * vc;
      case Z_M:
      case Z_P:
        return nx * ny * vc;
      default:
        Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }


  /// レベルL+1→Lの線形補間 (細f(i,j,k) → 粗c(I,J,K)).
  T interpolateF2C(const Scalar3D<T>& f, int I, int J, int K) {
    int i = 2 * I;
    int j = 2 * J;
    int k = 2 * K;
    return 0.125 * (f(i,j,k)   + f(i+1,j,k)   + f(i,j+1,k)   + f(i+1,j+1,k)
                  + f(i,j,k+1) + f(i+1,j,k+1) + f(i,j+1,k+1) + f(i+1,j+1,k+1));
  }

  /// レベルL+1→Lの線形補間 (細f(i,j,k) → 粗c(I,J,K)).
  T interpolateF2C(const T* fData, const Index3DS& fIndex, int I, int J, int K) {
    int i = 2 * I;
    int j = 2 * J;
    int k = 2 * K;
    return 0.125 * (fData[fIndex(i  ,j  ,k  )] + fData[fIndex(i+1,j  ,k  )]
                  + fData[fIndex(i  ,j+1,k  )] + fData[fIndex(i+1,j+1,k  )]
                  + fData[fIndex(i  ,j  ,k+1)] + fData[fIndex(i+1,j  ,k+1)]
                  + fData[fIndex(i  ,j+1,k+1)] + fData[fIndex(i+1,j+1,k+1)]);
  }



  /// レベルL→L+1の線形補間 (粗c(I,J,K) → 細f(i,j,k)).
  T interpolateC2F(const Scalar3D<T>& c, int i, int j, int k) {
    int I, J, K;
    double r, s, t;
    linearInterpolate(i, nx, I, r);
    linearInterpolate(j, ny, J, s);
    linearInterpolate(k, nz, K, t);

    return (1.0-t)*( 
                     (1.0-s)*( (1.0-r)*c(I  ,J  ,K  ) + r*c(I+1,J  ,K  ) )
                         + s*( (1.0-r)*c(I  ,J+1,K  ) + r*c(I+1,J+1,K  ) )
                    )
                +t*(
                     (1.0-s)*( (1.0-r)*c(I  ,J  ,K+1) + r*c(I+1,J  ,K+1) )
                         + s*( (1.0-r)*c(I  ,J+1,K+1) + r*c(I+1,J+1,K+1) )
                   );
  }

  /// レベルL→L+1の線形補間 (粗c(I,J,K) → 細f(i,j,k)).
  T interpolateC2F(const T* cData, const Index3DS& cIndex, int i, int j, int k) {
    int I, J, K;
    double r, s, t;
    linearInterpolate(i, nx, I, r);
    linearInterpolate(j, ny, J, s);
    linearInterpolate(k, nz, K, t);

    return (1.0-t)*( 
         (1.0-s)*( (1.0-r)*cData[cIndex(I  ,J  ,K  )] + r*cData[cIndex(I+1,J  ,K  )] )
             + s*( (1.0-r)*cData[cIndex(I  ,J+1,K  )] + r*cData[cIndex(I+1,J+1,K  )] )
        )
       +t*(
         (1.0-s)*( (1.0-r)*cData[cIndex(I  ,J  ,K+1)] + r*cData[cIndex(I+1,J  ,K+1)] )
             + s*( (1.0-r)*cData[cIndex(I  ,J+1,K+1)] + r*cData[cIndex(I+1,J+1,K+1)] )
        );
  }

  /// C2F補間における補間パラメータの計算.
  ///
  ///  @note 端点では，内挿ではなく外挿
  ///
  void linearInterpolate(int i, int n, int& I, double& r) {
#if 1
    I = std::min(std::max(i/2 - 1 + i%2, 0), n - 2);
    r = -0.25 + 0.5 * i - double(I);
#else
    if (i == 0) {
      // 外挿
      I = 0;
      r = -0.25;
    }
    else if (i == 2*n-1) {
      // 外挿
      I = n - 2;
      r = 1.25;
    }
    else if (i%2 == 0) {
      I = i/2 - 1;
      r = 0.75;
    }
    else {
      I = i/2;
      r = 0.25;
    }
#endif
  }


  /// 隣接データクラスから仮想セルデータをコピー(同レベル間).
  void copyFromNeighbor(Face face);

  /// 隣接データクラスから仮想セルデータをコピー(レベルL+1→L).
  void copyFromNeighborF2C(Face face, Subface subface);

  /// 隣接データクラスから仮想セルデータをコピー(レベルL→L+1).
  void copyFromNeighborC2F(Face face, Subface subface);

  /// 送信バッファに仮想セルデータをコピー(同レベル間).
  void copyToCommBuffer(Face face);

  /// 送信バッファに仮想セルデータをコピー(レベルL+1→L).
  void copyToCommBufferF2C(Face face, Subface subface);

  /// 送信バッファに仮想セルデータをコピー(レベルL→L+1).
  void copyToCommBufferC2F(Face face, Subface subface);

  /// 受信バッファから仮想セルデータをコピー(同レベル間).
  void copyFromCommBuffer(Face face);

  /// 受信バッファから仮想セルデータをコピー(レベルL+1→L).
  void copyFromCommBufferF2C(Face face, Subface subface);

  /// 受信バッファから仮想セルデータをコピー(レベルL→L+1).
  void copyFromCommBufferC2F(Face face, Subface subface);



  void copyFromNeighborF2C_0(int nx, int ny, int nz, int vc,
                             Face face, Subface subface,
                             const T* fData, Index3DS fIndex,
                             T* cData, Index3DS cIndex);

  void copyFromNeighborC2F_0(int nx, int ny, int nz, int vc,
                             Face face, Subface subface,
                             const T* cData, Index3DS cIndex,
                             T* fData, Index3DS fIndex);

  void copyToCommBufferC2F_0(int nx, int ny, int nz, int vc,
                             Face face, Subface subface,
                             const T* cData, Index3DS cIndex,
                             T* buffer);

  void copyToCommBufferF2C_0(int nx, int ny, int nz, int vc,
                             Face face, Subface subface,
                             const T* fData, Index3DS fIndex,
                             T* buffer);

};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // SCALAR_3D_UPDATER_H
