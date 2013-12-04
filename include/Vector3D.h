///
/// @file Vector3D.h
/// @brief ベクトルデータクラス
///

#ifndef VECTOR_3D_H
#define VECTOR_3D_H

#include "DataClass.h"
#include "VCUpdater.h"
#include "Index3D.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

/// ベクトルデータクラス(3成分連続格納).
template <typename T>
class Vector3D : public UpdatableDataClass {

private:

  int nx0, ny0, nz0;  ///< 内部配列サイズ(仮想セルを含めたセル分割数)

  int c_0, c_i, c_j, c_k;

  T* data;   ///< 1次元データ配列

  Index3DV index;  ///< インデックスファンクタ

public:

  /// コンストラクタ.
  ///
  ///  @param[in] size 分割数
  ///  @param[in] vc 仮想セル幅
  ///
  Vector3D(const Vec3i& size, int vc) : UpdatableDataClass(size, vc), index(size, vc) {
    nx0 = size[0] + 2*vc;
    ny0 = size[1] + 2*vc;
    nz0 = size[2] + 2*vc;
    data = new T[3*nx0*ny0*nz0];

    c_i = 3;
    c_j = 3 * nx0;
    c_k = 3 * nx0 * ny0;
    c_0 = 3 * (1 + nx0 + nx0 * ny0) * vc;
  }

  /// デストラクタ.
  ~Vector3D() { 
    delete[] data;
  }

  /// データ領域の取得.
  T* getData() const { return data; }

  /// インデックスファンクタの取得.
  Index3DV getIndex() const { return index; }

  /// 3次元添字によるデータアクセス.
  T& operator() (int i, int j, int k) {
    return data[c_i * i + c_j * j + c_k * k + c_0];
  }

  /// 3次元添字によるデータアクセス.
  const T& operator() (int i, int j, int k) const {
    return data[c_i * i + c_j * j + c_k * k + c_0];
  }

  /// 3次元添字によるデータアクセス.
  T& operator() (int i, int j, int k, int l) {
    return data[l + c_i * i + c_j * j + c_k * k + c_0];
  }

  /// 3次元添字によるデータアクセス.
  const T& operator() (int i, int j, int k, int l) const {
    return data[l + c_i * i + c_j * j + c_k * k + c_0];
  }

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


private:

  /// コピーコンストラクタ(コピー禁止).
  Vector3D(const Vector3D<T>& rhs);

  /// 代入演算子(コピー禁止).
  Vector3D& operator=(const Vector3D<T>& rhs);


  void copyToBuffer_0(int i0, int j0, int k0, int nx, int ny, int nz,
                      const T* data, const Index3DV& index, T* buffer) const;

  void copyFromBuffer_0(int i0, int j0, int k0, int nx, int ny, int nz,
                        const T* buffer, T* data, const Index3DV& index);

  void copyFromDataClass_0(int i0, int j0, int k0, int i1, int j1, int k1,
                           int nx, int ny, int nz,
                           const T* sData, const Index3DV& sIndex,
                           T* dData, const Index3DV& dIndex);


};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // VECTOR_3D_H
