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
/// @file Scalar3D.h
/// @brief スカラーデータクラス
///

#ifndef SCALAR_3D_H
#define SCALAR_3D_H

#include "VCUpdater.h"
#include "DataClass.h"
#include "Index3D.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif


/// スカラーデータクラス.
template <typename T>
class Scalar3D : public UpdatableDataClass {

private:

  int nx0, ny0, nz0;  ///< 内部配列サイズ(仮想セルを含めたセル分割数)

  int c_0, c_j, c_k;

  T* data;   ///< 1次元データ配列

  Index3DS index;  ///< インデックスファンクタ

	int vc0;

public:

  /// コンストラクタ.
  ///
  ///  @param[in] size 分割数
  ///  @param[in] vc 仮想セル幅
  ///
  Scalar3D(const ::Vec3i& size, int vc) : UpdatableDataClass(size, vc), index(size, vc) {
    nx0 = size[0] + 2*vc;
    ny0 = size[1] + 2*vc;
    nz0 = size[2] + 2*vc;
    data = new T[nx0*ny0*nz0];
		vc0 = vc;

    c_j = nx0;
    c_k = nx0 * ny0;
    c_0 = (1 + nx0 + nx0 * ny0) * vc;
  }

  /// コンストラクタ(for contiguous memory access).
  ///
  ///  @param[in] size 分割数
  ///  @param[in] vc 仮想セル幅
  ///
  Scalar3D(const ::Vec3i& size, int vc, T* data0) : UpdatableDataClass(size, vc), index(size, vc) {
    nx0 = size[0] + 2*vc;
    ny0 = size[1] + 2*vc;
    nz0 = size[2] + 2*vc;
    data = data0;
		vc0 = vc;

    c_j = nx0;
    c_k = nx0 * ny0;
    c_0 = (1 + nx0 + nx0 * ny0) * vc;
  }

  /// デストラクタ.
  ~Scalar3D() { 
    delete[] data;
  }

  /// データ領域の取得.
  T* getData() const { return data; }

	int getVCsize() { return vc0; }

  /// インデックスファンクタの取得.
  Index3DS getIndex() const { return index; }

  /// 3次元添字によるデータアクセス.
  T& operator() (int i, int j, int k) {
    return data[i + c_j * j + c_k * k + c_0];
  }

  /// 3次元添字によるデータアクセス.
  const T& operator() (int i, int j, int k) const {
    return data[i + c_j * j + c_k * k + c_0];
  }


/*
  /// 直方体領域からバッファへのデータコピー(シリアライズ).
  ///
  ///  @param[in] i0,j0,k0 コピー元の直方体領域の起点
  ///  @param[in] nx,ny,nz コピー元の直方体領域のサイズ
  ///  @param[out] buffer  コピー先バッファのアドレス
  ///
  void copyToBuffer(int i0, int j0, int k0, int nx, int ny, int nz, T* buffer) const;

  /// バッファから直方体領域へのデータコピー(デシリアライズ).
  ///
  ///  @param[in] i0,j0,k0 コピー先の直方体領域の起点
  ///  @param[in] nx,ny,nz コピー先の直方体領域のサイズ
  ///  @param[in] buffer  コピー元バッファのアドレス
  ///
  void copyFromBuffer(int i0, int j0, int k0, int nx, int ny, int nz, const T* buffer);

  /// 他データクラスの直方体領域から直方体領域へのデータコピー.
  ///
  ///  @param[in] i0,j0,k0 コピー元の直方体領域の起点
  ///  @param[in] i1,j1,k1 コピー先の直方体領域の起点
  ///  @param[in] nx,ny,nz 直方体領域のサイズ(コピー元/コピー先で共通)
  ///  @param[in] dataClass コピー元データクラス
  ///
  ///  @todo (i0,j0,k0)と(i1,j1,k1)を逆した方が分かりやすい？
  ///
  void copyFromDataClass(int i0, int j0, int k0, int i1, int j1, int k1,
                         int nx, int ny, int nz, const DataClass* dataClass);
*/

private:

  /// コピーコンストラクタ(コピー禁止).
  Scalar3D(const Scalar3D<T>& rhs);

  /// 代入演算子(コピー禁止).
  Scalar3D& operator=(const Scalar3D<T>& rhs);

/*
  void copyToBuffer_0(int i0, int j0, int k0, int nx, int ny, int nz,
                      const T* data, Index3DS index, T* buffer) const;

  void copyFromBuffer_0(int i0, int j0, int k0, int nx, int ny, int nz,
                        const T* buffer, T* data, Index3DS index);

  void copyFromDataClass_0(int i0, int j0, int k0, int i1, int j1, int k1,
                           int nx, int ny, int nz,
                           const T* sData, Index3DS sIndex,
                           T* dData, Index3DS dIndex);
*/


#define USE_PRIVATE_METHODS

public:
/// 直方体領域からバッファへのデータコピー(シリアライズ).
void copyToBuffer(int i0, int j0, int k0, int nx, int ny, int nz, T* buffer) const
{
#ifdef USE_PRIVATE_METHODS
  copyToBuffer_0(i0, j0, k0, nx, ny, nz, data, index, buffer);
#else
  for (int k = k0; k < k0 + nz; ++k) {
    for (int j = j0; j < j0 + ny; ++j) {
      for (int i = i0; i < i0 + nx; ++i) {
        buffer[i-i0 + nx*(j-j0) + (nx*ny)*(k-k0)] = data[index(i,j,k)];
      }
    }
  }
#endif
}


/// バッファから直方体領域へのデータコピー(デシリアライズ).
void copyFromBuffer(int i0, int j0, int k0, int nx, int ny, int nz, const T* buffer)
{
#ifdef USE_PRIVATE_METHODS
  copyFromBuffer_0(i0, j0, k0, nx, ny, nz, buffer, data, index);
#else
  for (int k = k0; k < k0 + nz; ++k) {
    for (int j = j0; j < j0 + ny; ++j) {
      for (int i = i0; i < i0 + nx; ++i) {
        data[index(i,j,k)] = buffer[i-i0 + nx*(j-j0) + (nx*ny)*(k-k0)];
      }
    }
  }
#endif
}


/// 他データクラスの直方体領域から直方体領域へのデータコピー.
void copyFromDataClass(int i0, int j0, int k0, int i1, int j1, int k1,
                                    int nx, int ny, int nz, const DataClass* dataClass)
{
  const Scalar3D<T>* s = dynamic_cast<const Scalar3D<T>*>(dataClass);
  T* sData = s->getData();
  Index3DS sIndex = s->getIndex();
#ifdef USE_PRIVATE_METHODS
  copyFromDataClass_0(i0, j0, k0, i1, j1, k1, nx, ny, nz, sData, sIndex, data, index);
#else
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        data[index(i0+i,j0+j,k0+k)] = sData[sIndex(i1+i,j1+j,k1+k)];
      }
    }
  }
#endif
}


private:
void copyToBuffer_0(int i0, int j0, int k0, int nx, int ny, int nz,
                                 const T* data, Index3DS index, T* buffer) const
{
#ifdef _BLOCK_IS_LARGE_
#pragma omp parallel for
#else
#endif
  for (int k = k0; k < k0 + nz; ++k) {
    for (int j = j0; j < j0 + ny; ++j) {
      for (int i = i0; i < i0 + nx; ++i) {
        buffer[i-i0 + nx*(j-j0) + (nx*ny)*(k-k0)] = data[index(i,j,k)];
      }
    }
  }
}


void copyFromBuffer_0(int i0, int j0, int k0, int nx, int ny, int nz,
                                   const T* buffer, T* data, Index3DS index)
{
#ifdef _BLOCK_IS_LARGE_
#pragma omp parallel for
#else
#endif
  for (int k = k0; k < k0 + nz; ++k) {
    for (int j = j0; j < j0 + ny; ++j) {
      for (int i = i0; i < i0 + nx; ++i) {
        data[index(i,j,k)] = buffer[i-i0 + nx*(j-j0) + (nx*ny)*(k-k0)];
      }
    }
  }
} 


void copyFromDataClass_0(int i0, int j0, int k0, int i1, int j1, int k1,
                                       int nx, int ny, int nz,
                                       const T* sData, Index3DS sIndex,
                                       T* dData, Index3DS dIndex)
{
#ifdef _BLOCK_IS_LARGE_
#pragma omp parallel for
#else
#endif
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
//      dData[dIndex(i0+i,j0+j,k0+k)] = sData[sIndex(i1+i,j1+j,k1+k)];
        int ii = dIndex(i0+i,j0+j,k0+k);
        int jj = sIndex(i1+i,j1+j,k1+k);
        dData[ii] = sData[jj];
      }
    }
  }
}

};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // SCALAR_3D_H
