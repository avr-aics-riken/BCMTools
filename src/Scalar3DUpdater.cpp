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

#include "Scalar3DUpdater.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

/*
/// 隣接データクラスから仮想セルデータをコピー(同レベル間).
template <typename T>
void Scalar3DUpdater<T>::copyFromNeighbor(Face face)
{
  Scalar3D<T>* dc = neighborDataClass[face][0];
  if (!dc) return;
  switch (face) {
    case X_M:
      dataClass->copyFromDataClass(-vc, 0, 0,  dc->getSizeX()-vc, 0, 0, vc, ny, nz,  dc);
      break;
    case X_P:
      dataClass->copyFromDataClass(nx, 0, 0,  0, 0, 0,  vc, ny, nz,  dc);
      break;
    case Y_M:
      dataClass->copyFromDataClass(0, -vc, 0,  0, dc->getSizeY()-vc, 0, nx, vc, nz,  dc);
      break;
    case Y_P:
      dataClass->copyFromDataClass(0, ny, 0,  0, 0, 0,  nx, vc, nz,  dc);
      break;
    case Z_M:
      dataClass->copyFromDataClass(0, 0, -vc,  0, 0, dc->getSizeZ()-vc, nx, ny, vc,  dc);
      break;
    case Z_P:
      dataClass->copyFromDataClass(0, 0, nz,  0, 0, 0,  nx, ny, vc,  dc);
      break;
    default:
      break;
  }
}


/// 隣接データクラスから仮想セルデータをコピー(レベルL+1→L).
template <typename T>
void Scalar3DUpdater<T>::copyFromNeighborF2C(Face face, Subface subface)
{
  T* cData = dataClass->getData();
  Index3DS cIndex = dataClass->getIndex();
  Scalar3D<T>* f = neighborDataClass[face][subface];
  T* fData = f->getData();
  Index3DS fIndex = f->getIndex();

  copyFromNeighborF2C_0(nx, ny, nz, vc, face, subface, fData, fIndex, cData, cIndex);
}


/// 隣接データクラスから仮想セルデータをコピー(レベルL→L+1).
template <typename T>
void Scalar3DUpdater<T>::copyFromNeighborC2F(Face face, Subface subface)
{
  T* fData = dataClass->getData();
  Index3DS fIndex = dataClass->getIndex();
  Scalar3D<T>* c = neighborDataClass[face][0];
  T* cData = c->getData();
  Index3DS cIndex = c->getIndex();

  copyFromNeighborC2F_0(nx, ny, nz, vc, face, subface, cData, cIndex, fData, fIndex);
}


/// 送信バッファに仮想セルデータをコピー(同レベル間).
template <typename T>
void Scalar3DUpdater<T>::copyToCommBuffer(Face face)
{
  T* buffer = sendBuffer[face][0];
  if (!buffer) return;
  switch (face) {
    case X_M:
      dataClass->copyToBuffer(0, 0, 0,  vc, ny, nz,  buffer);
      break;
    case X_P:
      dataClass->copyToBuffer(nx-vc, 0, 0,  vc, ny, nz,  buffer);
      break;
    case Y_M:
      dataClass->copyToBuffer(0, 0, 0,  nx, vc, nz,  buffer);
      break;
    case Y_P:
      dataClass->copyToBuffer(0, ny-vc, 0,  nx, vc, nz,  buffer);
      break;
    case Z_M:
      dataClass->copyToBuffer(0, 0, 0,  nx, ny, vc,  buffer);
      break;
    case Z_P:
      dataClass->copyToBuffer(0, 0, nz-vc,  nx, ny, vc,  buffer);
      break;
    default:
      break;
  }
}


/// 送信バッファに仮想セルデータをコピー(レベルL+1→L).
template <typename T>
void Scalar3DUpdater<T>::copyToCommBufferF2C(Face face, Subface subface)
{
  T* buffer = sendBuffer[face][0];
  T* fData = dataClass->getData();
  Index3DS fIndex = dataClass->getIndex();

  copyToCommBufferF2C_0(nx, ny, nz, vc, face, subface, fData, fIndex, buffer);
}


/// 送信バッファに仮想セルデータをコピー(レベルL→L+1).
template <typename T>
void Scalar3DUpdater<T>::copyToCommBufferC2F(Face face, Subface subface)
{
  T* cData = dataClass->getData();
  Index3DS cIndex = dataClass->getIndex();
  T* buffer = sendBuffer[face][subface];

  copyToCommBufferC2F_0(nx, ny, nz, vc, face, subface, cData, cIndex, buffer);
}


/// 受信バッファから仮想セルデータをコピー(同レベル間).
template <typename T>
void Scalar3DUpdater<T>::copyFromCommBuffer(Face face)
{
  T* buffer = recvBuffer[face][0];
  if (!buffer) return;
  switch (face) {
    case X_M:
      dataClass->copyFromBuffer(-vc, 0, 0,  vc, ny, nz,  buffer);
      break;
    case X_P:
      dataClass->copyFromBuffer(nx, 0, 0,  vc, ny, nz,  buffer);
      break;
    case Y_M:
      dataClass->copyFromBuffer(0, -vc, 0,  nx, vc, nz,  buffer);
      break;
    case Y_P:
      dataClass->copyFromBuffer(0, ny, 0,  nx, vc, nz,  buffer);
      break;
    case Z_M:
      dataClass->copyFromBuffer(0, 0, -vc,  nx, ny, vc,  buffer);
      break;
    case Z_P:
      dataClass->copyFromBuffer(0, 0, nz,  nx, ny, vc,  buffer);
      break;
    default:
      break;
  }
}


/// 受信バッファから仮想セルデータをコピー(レベルL+1→L).
template <typename T>
void Scalar3DUpdater<T>::copyFromCommBufferF2C(Face face, Subface subface)
{
  T* buffer = recvBuffer[face][subface];
  switch (face) {
    case X_M:
    {
      int j0 = (ny/2) * subfaceOrigin0(subface);
      int k0 = (nz/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(-vc, j0, k0, vc, ny/2, nz/2, buffer);
      break;
    }
    case X_P:
    {
      int j0 = (ny/2) * subfaceOrigin0(subface);
      int k0 = (nz/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(nx, j0, k0, vc, ny/2, nz/2, buffer);
      break;
    }
    case Y_M:
    {
      int k0 = (nz/2) * subfaceOrigin0(subface);
      int i0 = (nx/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(i0, -vc, k0,  nx/2, vc, nz/2,  buffer);
      break;
    }
    case Y_P:
    {
      int k0 = (nz/2) * subfaceOrigin0(subface);
      int i0 = (nx/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(i0, ny, k0,  nx/2, vc, nz/2,  buffer);
      break;
    }
    case Z_M:
    {
      int i0 = (nx/2) * subfaceOrigin0(subface);
      int j0 = (ny/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(i0, j0, -vc,  nx/2, ny/2, vc,  buffer);
      break;
    }
    case Z_P:
    {
      int i0 = (nx/2) * subfaceOrigin0(subface);
      int j0 = (ny/2) * subfaceOrigin1(subface);
      dataClass->copyFromBuffer(i0, j0, nz,  nx/2, ny/2, vc,  buffer);
      break;
    }
    default:
      break;
  }
}


/// 受信バッファから仮想セルデータをコピー(レベルL→L+1).
template <typename T>
void Scalar3DUpdater<T>::copyFromCommBufferC2F(Face face, Subface subface)
{
  copyFromCommBuffer(face);
}



template <typename T>
void Scalar3DUpdater<T>::copyFromNeighborF2C_0(int nx, int ny, int nz, int vc,
                                               Face face, Subface subface,
                                               const T* fData, Index3DS fIndex,
                                               T* cData, Index3DS cIndex)
{
  switch (face) {
    case X_M:
    {
      int j0 = (ny/2) * subfaceOrigin0(subface);
      int k0 = (nz/2) * subfaceOrigin1(subface);
#pragma loop noalias
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < ny/2; j++) {
          for (int i = 0; i < vc; i++) {
            cData[cIndex(i-vc, j+j0, k+k0)] = interpolateF2C(fData, fIndex, i+nx/2-vc, j, k);
          }
        }
      }
      break;
    }
    case X_P:
    {
      int j0 = (ny/2) * subfaceOrigin0(subface);
      int k0 = (nz/2) * subfaceOrigin1(subface);
#pragma loop noalias
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < ny/2; j++) {
          for (int i = 0; i < vc; i++) {
            cData[cIndex(i+nx, j+j0, k+k0)] = interpolateF2C(fData, fIndex, i, j, k);
          }
        }
      }
      break;
    }
    case Y_M:
    {
      int k0 = (nz/2) * subfaceOrigin0(subface);
      int i0 = (nx/2) * subfaceOrigin1(subface);
#pragma loop noalias
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < vc; j++) {
          for (int i = 0; i < nx/2; i++) {
            cData[cIndex(i+i0, j-vc, k+k0)] = interpolateF2C(fData, fIndex, i, j+ny/2-vc, k);
          }
        }
      }
      break;
    }
    case Y_P:
    {
      int k0 = (nz/2) * subfaceOrigin0(subface);
      int i0 = (nx/2) * subfaceOrigin1(subface);
#pragma loop noalias
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < vc; j++) {
          for (int i = 0; i < nx/2; i++) {
            cData[cIndex(i+i0, j+ny, k+k0)] = interpolateF2C(fData, fIndex, i, j, k);
          }
        }
      }
      break;
    }
    case Z_M:
    {
      int i0 = (nx/2) * subfaceOrigin0(subface);
      int j0 = (ny/2) * subfaceOrigin1(subface);
#pragma loop noalias
      for (int k = 0; k < vc; k++) {
        for (int j = 0; j < ny/2; j++) {
          for (int i = 0; i < nx/2; i++) {
            cData[cIndex(i+i0, j+j0, k-vc)] = interpolateF2C(fData, fIndex, i, j, k+nz/2-vc);
          }
        }
      }
      break;
    }
    case Z_P:
    {
      int i0 = (nx/2) * subfaceOrigin0(subface);
      int j0 = (ny/2) * subfaceOrigin1(subface);
#pragma loop noalias
      for (int k = 0; k < vc; k++) {
        for (int j = 0; j < ny/2; j++) {
          for (int i = 0; i < nx/2; i++) {
            cData[cIndex(i+i0, j+j0, k+nz)] = interpolateF2C(fData, fIndex, i, j, k);
          }
        }
      }
      break;
    }
    default:
      break;
  }
}


template <typename T>
void Scalar3DUpdater<T>::copyFromNeighborC2F_0(int nx, int ny, int nz, int vc,
                                               Face face, Subface subface,
                                               const T* cData, Index3DS cIndex,
                                               T* fData, Index3DS fIndex)
{
  switch (face) {
    case X_M:
    {
      int J0 = ny * subfaceOrigin0(subface);
      int K0 = nz * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < nz; K++) {
        for (int J = 0; J < ny; J++) {
          for (int I = 0; I < vc; I++) {
            fData[fIndex(I-vc, J, K)] = interpolateC2F(cData, cIndex, I+2*nx-vc, J+J0, K+K0);
          }
        }
      }
      break;
    }
    case X_P:
    {
      int J0 = ny * subfaceOrigin0(subface);
      int K0 = nz * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < nz; K++) {
        for (int J = 0; J < ny; J++) {
          for (int I = 0; I < vc; I++) {
            fData[fIndex(I+nx, J, K)] = interpolateC2F(cData, cIndex, I, J+J0, K+K0);
          }
        }
      }
      break;
    }
    case Y_M:
    {
      int K0 = nz * subfaceOrigin0(subface);
      int I0 = nx * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < nz; K++) {
        for (int J = 0; J < vc; J++) {
          for (int I = 0; I < nx; I++) {
            fData[fIndex(I, J-vc, K)] = interpolateC2F(cData, cIndex, I+I0, J+2*ny-vc, K+K0);
          }
        }
      }
      break;
    }
    case Y_P:
    {
      int K0 = nz * subfaceOrigin0(subface);
      int I0 = nx * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < nz; K++) {
        for (int J = 0; J < vc; J++) {
          for (int I = 0; I < nx; I++) {
            fData[fIndex(I, J+ny, K)] = interpolateC2F(cData, cIndex, I+I0, J, K+K0);
          }
        }
      }
      break;
    }
    case Z_M:
    {
      int I0 = nx * subfaceOrigin0(subface);
      int J0 = ny * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < vc; K++) {
        for (int J = 0; J < ny; J++) {
          for (int I = 0; I < nx; I++) {
            fData[fIndex(I, J, K-vc)] = interpolateC2F(cData, cIndex, I+I0, J+J0, K+2*nz-vc);
          }
        }
      }
      break;
    }
    case Z_P:
    {
      int I0 = nx * subfaceOrigin0(subface);
      int J0 = ny * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < vc; K++) {
        for (int J = 0; J < ny; J++) {
          for (int I = 0; I < nx; I++) {
            fData[fIndex(I, J, K+nz)] = interpolateC2F(cData, cIndex, I+I0, J+J0, K);
          }
        }
      }
      break;
    }
    default:
      break;
  }
}


template <typename T>
void Scalar3DUpdater<T>::copyToCommBufferC2F_0(int nx, int ny, int nz, int vc,
                                               Face face, Subface subface,
                                               const T* cData, Index3DS cIndex,
                                               T* buffer)
{
  int ii = 0;
  switch (face) {
    case X_M:
    {
      int J0 = ny * subfaceOrigin0(subface);
      int K0 = nz * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < nz; K++) {
        for (int J = 0; J < ny; J++) {
          for (int I = 0; I < vc; I++) {
            buffer[ii++] = interpolateC2F(cData, cIndex, I, J+J0, K+K0);
          }
        }
      }
      break;
    }
    case X_P:
    {
      int J0 = ny * subfaceOrigin0(subface);
      int K0 = nz * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < nz; K++) {
        for (int J = 0; J < ny; J++) {
          for (int I = 0; I < vc; I++) {
            buffer[ii++] = interpolateC2F(cData, cIndex, I+2*nx-vc, J+J0, K+K0);
          }
        }
      }
      break;
    }
    case Y_M:
    {
      int K0 = nz * subfaceOrigin0(subface);
      int I0 = nx * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < nz; K++) {
        for (int J = 0; J < vc; J++) {
          for (int I = 0; I < nx; I++) {
            buffer[ii++] = interpolateC2F(cData, cIndex, I+I0, J, K+K0);
          }
        }
      }
      break;
    }
    case Y_P:
    {
      int K0 = nz * subfaceOrigin0(subface);
      int I0 = nx * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < nz; K++) {
        for (int J = 0; J < vc; J++) {
          for (int I = 0; I < nx; I++) {
            buffer[ii++] = interpolateC2F(cData, cIndex, I+I0, J+2*ny-vc, K+K0);
          }
        }
      }
      break;
    }
    case Z_M:
    {
      int I0 = nx * subfaceOrigin0(subface);
      int J0 = ny * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < vc; K++) {
        for (int J = 0; J < ny; J++) {
          for (int I = 0; I < nx; I++) {
            buffer[ii++] = interpolateC2F(cData, cIndex, I+I0, J+J0, K);
          }
        }
      }
      break;
    }
    case Z_P:
    {
      int I0 = nx * subfaceOrigin0(subface);
      int J0 = ny * subfaceOrigin1(subface);
#pragma loop noalias
      for (int K = 0; K < vc; K++) {
        for (int J = 0; J < ny; J++) {
          for (int I = 0; I < nx; I++) {
            buffer[ii++] = interpolateC2F(cData, cIndex, I+I0, J+J0, K+2*nz-vc);
          }
        }
      }
      break;
    }
    default:
      break;
  }
}

template <typename T>
void Scalar3DUpdater<T>::copyToCommBufferF2C_0(int nx, int ny, int nz, int vc,
                                               Face face, Subface subface,
                                               const T* fData, Index3DS fIndex,
                                               T* buffer)
{
  int ii = 0;
  switch (face) {
    case X_M:
    {
#pragma loop noalias
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < ny/2; j++) {
          for (int i = 0; i < vc; i++) {
            buffer[ii++] = interpolateF2C(fData, fIndex, i, j, k);
          }
        }
      }
      break;
    }
    case X_P:
    {
#pragma loop noalias
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < ny/2; j++) {
          for (int i = 0; i < vc; i++) {
            buffer[ii++] = interpolateF2C(fData, fIndex, i+nx/2-vc, j, k);
          }
        }
      }
      break;
    }
    case Y_M:
    {
#pragma loop noalias
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < vc; j++) {
          for (int i = 0; i < nx/2; i++) {
            buffer[ii++] = interpolateF2C(fData, fIndex, i, j, k);
          }
        }
      }
      break;
    }
    case Y_P:
    {
#pragma loop noalias
      for (int k = 0; k < nz/2; k++) {
        for (int j = 0; j < vc; j++) {
          for (int i = 0; i < nx/2; i++) {
            buffer[ii++] =  interpolateF2C(fData, fIndex, i, j+ny/2-vc, k);
          }
        }
      }
      break;
    }
    case Z_M:
    {
#pragma loop noalias
      for (int k = 0; k < vc; k++) {
        for (int j = 0; j < ny/2; j++) {
          for (int i = 0; i < nx/2; i++) {
            buffer[ii++] = interpolateF2C(fData, fIndex, i, j, k);
          }
        }
      }
      break;
    }
    case Z_P:
    {
#pragma loop noalias
      for (int k = 0; k < vc; k++) {
        for (int j = 0; j < ny/2; j++) {
          for (int i = 0; i < nx/2; i++) {
            buffer[ii++] = interpolateF2C(fData, fIndex, i, j, k+nz/2-vc);
          }
        }
      }
      break;
    }
    default:
      break;
  }
}


#ifdef USE_UPDATER_SCALAR_DOUBLE
template class Scalar3DUpdater<double>;
#endif

#ifdef USE_UPDATER_SCALAR_FLOAT
template class Scalar3DUpdater<float>;
#endif

#ifdef USE_UPDATER_SCALAR_INT
//#error Scalar3DUpdater not implemented for integer type yet.
template class Scalar3DUpdater<int>;
#endif

#if (defined USE_UPDATER_SCALAR_UNSIGNED_INT || defined USE_UPDATER_SCALAR_UNSIGNED)
#error Scalar3DUpdater not implemented for integer type yet.
//template class Scalar3DUpdater<unsigned int>;
#endif
*/

#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif
