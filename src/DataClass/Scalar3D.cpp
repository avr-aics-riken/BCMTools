#include "Scalar3D.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

#define USE_PRIVATE_METHODS


/// 直方体領域からバッファへのデータコピー(シリアライズ).
template <typename T>
void Scalar3D<T>::copyToBuffer(int i0, int j0, int k0, int nx, int ny, int nz, T* buffer) const
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
template <typename T>
void Scalar3D<T>::copyFromBuffer(int i0, int j0, int k0, int nx, int ny, int nz, const T* buffer)
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
template <typename T>
void Scalar3D<T>::copyFromDataClass(int i0, int j0, int k0, int i1, int j1, int k1,
                                    int nx, int ny, int nz, const DataClass* dataClass)
{
  const Scalar3D<T>* s = dynamic_cast<const Scalar3D<T>*>(dataClass);
  T* sData = s->getData();
  Index3DS sIndex = s->getIndex();
#if USE_PRIVATE_METHODs
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



template <typename T>
void Scalar3D<T>::copyToBuffer_0(int i0, int j0, int k0, int nx, int ny, int nz,
                                 const T* data, Index3DS index, T* buffer) const
{
#pragma loop noalias
#pragma omp parallel for if(nz >= 16)
  for (int k = k0; k < k0 + nz; ++k) {
    for (int j = j0; j < j0 + ny; ++j) {
      for (int i = i0; i < i0 + nx; ++i) {
        buffer[i-i0 + nx*(j-j0) + (nx*ny)*(k-k0)] = data[index(i,j,k)];
      }
    }
  }
}


template <typename T>
void Scalar3D<T>::copyFromBuffer_0(int i0, int j0, int k0, int nx, int ny, int nz,
                                   const T* buffer, T* data, Index3DS index)
{
#pragma loop noalias
#pragma omp parallel for if(nz >= 16)
  for (int k = k0; k < k0 + nz; ++k) {
    for (int j = j0; j < j0 + ny; ++j) {
      for (int i = i0; i < i0 + nx; ++i) {
        data[index(i,j,k)] = buffer[i-i0 + nx*(j-j0) + (nx*ny)*(k-k0)];
      }
    }
  }
} 


template <typename T>
void Scalar3D<T>:: copyFromDataClass_0(int i0, int j0, int k0, int i1, int j1, int k1,
                                       int nx, int ny, int nz,
                                       const T* sData, Index3DS sIndex,
                                       T* dData, Index3DS dIndex)
{
#pragma loop noalias
#pragma omp parallel for schedule(static,1) if(nz >= 16)
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



#ifdef USE_DATACLASS_SCALAR_DOUBLE
template class Scalar3D<double>;
#endif

#ifdef USE_DATACLASS_SCALAR_FLOAT
template class Scalar3D<float>;
#endif

#ifdef USE_DATACLASS_SCALAR_INT
template class Scalar3D<int>;
#endif

#if (defined USE_DATACLASS_SCALAR_UNSIGNED_INT || defined USE_DATACLASS_SCALAR_UNSIGNED)
template class Scalar3D<unsigned int>;
#endif

#ifdef USE_DATACLASS_SCALAR_CHAR
template class Scalar3D<char>;
#endif

#ifdef USE_DATACLASS_SCALAR_UNSIGNED_CHAR
template class Scalar3D<unsigned char>;
#endif

#ifdef USE_DATACLASS_SCALAR_LONG
template class Scalar3D<long>;
#endif

#ifdef USE_DATACLASS_SCALAR_UNSIGNED_LONG
template class Scalar3D<unsigned long>;
#endif


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif
