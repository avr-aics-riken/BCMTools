#include "Vector3D.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

#define USE_PRIVATE_METHODS


/// 直方体領域からバッファへのデータコピー(シリアライズ).
template <typename T>
void Vector3D<T>::copyToBuffer(int i0, int j0, int k0, int nx, int ny, int nz, T* buffer) const
{
#ifdef USE_PRIVATE_METHODS
  copyToBuffer_0(i0, j0, k0, nx, ny, nz, data, index, buffer);
#else
  for (int k = k0; k < k0 + nz; ++k) {
    for (int j = j0; j < j0 + ny; ++j) {
      for (int i = i0; i < i0 + nx; ++i) {
        int ii = 3*(i-i0) + (3*nx)*(j-j0) + (3*nx*ny)*(k-k0);
        int jj = index(i,j,k);
        buffer[ii  ] = data[jj  ];
        buffer[ii+1] = data[jj+1];
        buffer[ii+2] = data[jj+2];
      }
    }
  }
#endif
}


/// バッファから直方体領域へのデータコピー(デシリアライズ).
template <typename T>
void Vector3D<T>::copyFromBuffer(int i0, int j0, int k0, int nx, int ny, int nz, const T* buffer)
{
#ifdef USE_PRIVATE_METHODS
  copyFromBuffer_0(i0, j0, k0, nx, ny, nz, buffer, data, index);
#else
  for (int k = k0; k < k0 + nz; ++k) {
    for (int j = j0; j < j0 + ny; ++j) {
      for (int i = i0; i < i0 + nx; ++i) {
        int ii = 3*(i-i0) + (3*nx)*(j-j0) + (3*nx*ny)*(k-k0);
        int jj = index(i,j,k);
        data[jj  ] = buffer[ii  ];
        data[jj+1] = buffer[ii+1];
        data[jj+2] = buffer[ii+2];
      }
    }
  }
#endif
}


/// 他データクラスの直方体領域から直方体領域へのデータコピー.
template <typename T>
void Vector3D<T>::copyFromDataClass(int i0, int j0, int k0, int i1, int j1, int k1,
                       int nx, int ny, int nz, const DataClass* dataClass)
{
  const Vector3D<T>* s = dynamic_cast<const Vector3D<T>*>(dataClass);
  T* sData = s->getData();
  Index3DV sIndex = s->getIndex();
#ifdef USE_PRIVATE_METHODS
  copyFromDataClass_0(i0, j0, k0, i1, j1, k1, nx, ny, nz, sData, sIndex, data, index);
#else
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int ii = index(i0+i,j0+j,k0+k);
        int jj = sIndex(i1+i,j1+j,k1+k);
        data[ii  ] = sData[jj  ];
        data[ii+1] = sData[jj+1];
        data[ii+2] = sData[jj+2];
      }
    }
  }
#endif
}



template <typename T>
void Vector3D<T>::copyToBuffer_0(int i0, int j0, int k0, int nx, int ny, int nz,
                                 const T* data, const Index3DV& index, T* buffer) const
{
#pragma loop noalias
  for (int k = k0; k < k0 + nz; ++k) {
    for (int j = j0; j < j0 + ny; ++j) {
      for (int i = i0; i < i0 + nx; ++i) {
        int ii = 3*(i-i0) + (3*nx)*(j-j0) + (3*nx*ny)*(k-k0);
        int jj = index(i,j,k);
        buffer[ii  ] = data[jj  ];
        buffer[ii+1] = data[jj+1];
        buffer[ii+2] = data[jj+2];
      }
    }
  }
}


template <typename T>
void Vector3D<T>::copyFromBuffer_0(int i0, int j0, int k0, int nx, int ny, int nz,
                                   const T* buffer, T* data, const Index3DV& index)
{
#pragma loop noalias
  for (int k = k0; k < k0 + nz; ++k) {
    for (int j = j0; j < j0 + ny; ++j) {
      for (int i = i0; i < i0 + nx; ++i) {
        int ii = 3*(i-i0) + (3*nx)*(j-j0) + (3*nx*ny)*(k-k0);
        int jj = index(i,j,k);
        data[jj  ] = buffer[ii  ];
        data[jj+1] = buffer[ii+1];
        data[jj+2] = buffer[ii+2];
      }
    }
  }
} 


template <typename T>
void Vector3D<T>:: copyFromDataClass_0(int i0, int j0, int k0, int i1, int j1, int k1,
                                       int nx, int ny, int nz,
                                       const T* sData, const Index3DV& sIndex,
                                       T* dData, const Index3DV& dIndex)
{
//#pragma omp parallel for
#pragma loop noalias
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int ii = dIndex(i0+i,j0+j,k0+k);
        int jj = sIndex(i1+i,j1+j,k1+k);
        dData[ii  ] = sData[jj  ];
        dData[ii+1] = sData[jj+1];
        dData[ii+2] = sData[jj+2];
      }
    }
  }
}



#ifdef USE_DATACLASS_VECTOR_DOUBLE
template class Vector3D<double>;
#endif

#ifdef USE_DATACLASS_VECTOR_FLOAT
template class Vector3D<float>;
#endif

#ifdef USE_DATACLASS_VECTOR_INT
template class Vector3D<int>;
#endif

#if (defined USE_DATACLASS_VECTOR_UNSIGNED_INT || defined USE_DATACLASS_VECTOR_UNSIGNED)
template class Vector3D<unsigned int>;
#endif

#ifdef USE_DATACLASS_VECTOR_CHAR
template class Vector3D<char>;
#endif

#ifdef USE_DATACLASS_VECTOR_UNSIGNED_CHAR
template class Vector3D<unsigned char>;
#endif

#ifdef USE_DATACLASS_VECTOR_LONG
template class Vector3D<long>;
#endif

#ifdef USE_DATACLASS_VECTOR_UNSIGNED_LONG
template class Vector3D<unsigned long>;
#endif


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif
