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
/// @file VtkWriter.h
/// @brief VTKフォーマットによるBCMデータ出力クラス
/// 
/// @note 
///  
///

#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>

#include "mpi.h"
#include "BCMTools.h"
#include "BlockManager.h"
#include "Scalar3D.h"
#include "Vector3D.h"
#include "BCMOctree.h"
#include "Partition.h"

using namespace std;

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif


class VtkWriter {
  BlockManager& blockManager;
  const MPI::Intracomm& comm;

public:
	VtkWriter()
			: blockManager(BlockManager::getInstance()),
				comm(blockManager.getCommunicator()) {
    int myrank = comm.Get_rank();

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
    }

    if (myrank == 0) {
    } else {
    }
  
  }

  /// デストラクタ.
  ~VtkWriter() {
  }

  template <typename T>
	void printVTIC(T* sData, 
								const char* label, int step, int rank, int block,
								int NX, int NY, int NZ,
								int vc,
								double ox, double oy, double oz,
								double dx) {
		ostringstream ossFileName;
		ossFileName << "./VTK/";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << "/";
		ossFileName << "data-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << rank;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << block;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".vti";

		int iNX1 = 0;
		int iNY1 = 0;
		int iNZ1 = 0;
		int iNXN = NX;
		int iNYN = NY;
		int iNZN = NZ;

		unsigned int nSizeX = NX;
		unsigned int nSizeY = NY;
		unsigned int nSizeZ = NZ;
		unsigned int nSize = nSizeX*nSizeY*nSizeZ;

		ofstream ofs;
		ofs.open(ossFileName.str().c_str(), ios::out);
		ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"";
#ifdef __FUJITSU
		ofs << "BigEndian";
#else
		ofs << "LittleEndian";
#endif
		ofs << "\">" << endl;
		ofs << "<ImageData WholeExtent=\"";
		ofs << iNX1 << " ";
		ofs << iNXN << " ";
		ofs << iNY1 << " ";
		ofs << iNYN << " ";
		ofs << iNZ1 << " ";
		ofs << iNZN << "\" ";
		ofs << "Origin=\"";
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(16);
		ofs << ox << " ";
		ofs << oy << " ";
		ofs << oz << "\" ";
		ofs << "Spacing=\"";
		ofs << dx << " ";
		ofs << dx << " ";
		ofs << dx << "\">" << std::endl;
		ofs << "<Piece Extent=\"";
		ofs << iNX1 << " ";
		ofs << iNXN << " ";
		ofs << iNY1 << " ";
		ofs << iNYN << " ";
		ofs << iNZ1 << " ";
		ofs << iNZN << "\">" << std::endl;
		ofs << "<PointData>" << endl;
		ofs << "</PointData>" << endl;
		ofs << "<CellData>" << endl;

		ofs << "<DataArray type=\"Float32\" Name=\"";
		ofs << label;
		ofs << "\" format=\"appended\" offset=\"";
		ofs << (sizeof(int) + sizeof(float)*nSize)*0;
		ofs << "\"/>" << endl;

		ofs << "</CellData>" << endl;
		ofs << "<Coordinates>" << endl;
		ofs << "</Coordinates>" << endl;
		ofs << "</Piece>" << endl;
		ofs << "</ImageData>" << endl;
		ofs << "<AppendedData encoding=\"raw\">" << endl;
		ofs << "_";
		ofs.close();

		ofs.open(ossFileName.str().c_str(), ios::out | ios::app | ios::binary);

		float* pScalar = new float[nSize];
		int nBytes = sizeof(float)*nSize;

		for(int k=0; k<NZ; k++) {
			for(int j=0; j<NY; j++) {
				for(int i=0; i<NX; i++) {
					int i0 = i + vc;
					int j0 = j + vc;
					int k0 = k + vc;
					int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
					float data = sData[m0];
					int m = i + NX*( j + NY*k );
					pScalar[m] = data;
				}
			}
		}

		ofs.write((char*)&nBytes, sizeof(int));
		ofs.write((char*)pScalar, sizeof(float)*(nSize));

		ofs.close();

		ofs.open(ossFileName.str().c_str(), ios::out | ios::app);
		ofs << endl;
		ofs << "</AppendedData>" << endl;
		ofs << "</VTKFile>" << endl;
		ofs.close();

		delete [] pScalar;
	}

  template <typename T>
	void printVTIP(T* sData, 
								const char* label, int step, int rank, int block,
								int NX, int NY, int NZ,
								int vc,
								double ox, double oy, double oz,
								double dx) {
		ostringstream ossFileName;
		ossFileName << "./VTK/";
		ossFileName << "data-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << rank;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << block;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".vti";

		int iNX1 = 0;
		int iNY1 = 0;
		int iNZ1 = 0;
		int iNXN = NX;
		int iNYN = NY;
		int iNZN = NZ;

		unsigned int nSizeX = NX + 1;
		unsigned int nSizeY = NY + 1;
		unsigned int nSizeZ = NZ + 1;
		unsigned int nSize = nSizeX*nSizeY*nSizeZ;

		ofstream ofs;
		ofs.open(ossFileName.str().c_str(), ios::out);
		ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"";
#ifdef __FUJITSU
		ofs << "BigEndian";
#else
		ofs << "LittleEndian";
#endif
		ofs << "\">" << endl;
		ofs << "<ImageData WholeExtent=\"";
		ofs << iNX1 << " ";
		ofs << iNXN << " ";
		ofs << iNY1 << " ";
		ofs << iNYN << " ";
		ofs << iNZ1 << " ";
		ofs << iNZN << "\" ";
		ofs << "Origin=\"";
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(16);
		ofs << ox << " ";
		ofs << oy << " ";
		ofs << oz << "\" ";
		ofs << "Spacing=\"";
		ofs << dx << " ";
		ofs << dx << " ";
		ofs << dx << "\">" << std::endl;
		ofs << "<Piece Extent=\"";
		ofs << iNX1 << " ";
		ofs << iNXN << " ";
		ofs << iNY1 << " ";
		ofs << iNYN << " ";
		ofs << iNZ1 << " ";
		ofs << iNZN << "\">" << std::endl;
		ofs << "<PointData>" << endl;

		ofs << "<DataArray type=\"Float32\" Name=\"";
		ofs << label;
		ofs << "\" format=\"appended\" offset=\"";
		ofs << (sizeof(int) + sizeof(float)*nSize)*0;
		ofs << "\"/>" << endl;

		ofs << "</PointData>" << endl;
		ofs << "<CellData>" << endl;
		ofs << "</CellData>" << endl;
		ofs << "<Coordinates>" << endl;
		ofs << "</Coordinates>" << endl;
		ofs << "</Piece>" << endl;
		ofs << "</ImageData>" << endl;
		ofs << "<AppendedData encoding=\"raw\">" << endl;
		ofs << "_";
		ofs.close();

		ofs.open(ossFileName.str().c_str(), ios::out | ios::app | ios::binary);

		float* pScalar = new float[nSize];
		int nBytes = sizeof(float)*nSize;

		for(int k=0; k<=NZ; k++) {
			for(int j=0; j<=NY; j++) {
				for(int i=0; i<=NX; i++) {
/*
					int il = i + vc - 1>=0 ? i + vc - 1 : 0;
					int jl = j + vc - 1>=0 ? j + vc - 1 : 0;
					int kl = k + vc - 1>=0 ? k + vc - 1 : 0;
					int iu = i + vc<NX ? i + vc : NX-1;
					int ju = j + vc<NY ? j + vc : NY-1;
					int ku = k + vc<NZ ? k + vc : NZ-1;
*/
					int il = i + vc - 1;
					int jl = j + vc - 1;
					int kl = k + vc - 1;
					int iu = i + vc;
					int ju = j + vc;
					int ku = k + vc;

					int i0 = il;
					int j0 = jl;
					int k0 = kl;
					int i1 = iu;
					int j1 = jl;
					int k1 = kl;
					int i2 = il;
					int j2 = ju;
					int k2 = kl;
					int i3 = il;
					int j3 = jl;
					int k3 = ku;
					int i4 = iu;
					int j4 = ju;
					int k4 = kl;
					int i5 = il;
					int j5 = ju;
					int k5 = ku;
					int i6 = iu;
					int j6 = jl;
					int k6 = ku;
					int i7 = iu;
					int j7 = ju;
					int k7 = ku;

					int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
					int m1 = i1 + (NX+2*vc)*( j1 + (NY+2*vc)*k1 );
					int m2 = i2 + (NX+2*vc)*( j2 + (NY+2*vc)*k2 );
					int m3 = i3 + (NX+2*vc)*( j3 + (NY+2*vc)*k3 );
					int m4 = i4 + (NX+2*vc)*( j4 + (NY+2*vc)*k4 );
					int m5 = i5 + (NX+2*vc)*( j5 + (NY+2*vc)*k5 );
					int m6 = i6 + (NX+2*vc)*( j6 + (NY+2*vc)*k6 );
					int m7 = i7 + (NX+2*vc)*( j7 + (NY+2*vc)*k7 );

					float data0 = sData[m0];
					float data1 = sData[m1];
					float data2 = sData[m2];
					float data3 = sData[m3];
					float data4 = sData[m4];
					float data5 = sData[m5];
					float data6 = sData[m6];
					float data7 = sData[m7];

					float data = 0.125*(data0 + data1 + data2 + data3 + data4 + data5 + data6 + data7);

					int m = i + (NX+1)*( j + (NY+1)*k );
					pScalar[m] = data;
				}
			}
		}

		ofs.write((char*)&nBytes, sizeof(int));
		ofs.write((char*)pScalar, sizeof(float)*(nSize));

		ofs.close();

		ofs.open(ossFileName.str().c_str(), ios::out | ios::app);
		ofs << endl;
		ofs << "</AppendedData>" << endl;
		ofs << "</VTKFile>" << endl;
		ofs.close();

		delete [] pScalar;
	}

  template <typename T>
	void printVTIC(
								T* sDataP, 
								T* sDataUX, 
								T* sDataUY, 
								T* sDataUZ, 
								const char* label, int step, int rank, int block,
								int NX, int NY, int NZ,
								int vc,
								double ox, double oy, double oz,
								double dx) {
		ostringstream ossFileName;
		ossFileName << "./VTK/";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << "/";
		ossFileName << "data-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << rank;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << block;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".vti";

		int iNX1 = 0;
		int iNY1 = 0;
		int iNZ1 = 0;
		int iNXN = NX;
		int iNYN = NY;
		int iNZN = NZ;

		unsigned int nSizeX = NX;
		unsigned int nSizeY = NY;
		unsigned int nSizeZ = NZ;
		unsigned int nSize = nSizeX*nSizeY*nSizeZ;

		ofstream ofs;
		ofs.open(ossFileName.str().c_str(), ios::out);
		ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"";
#ifdef __FUJITSU
		ofs << "BigEndian";
#else
		ofs << "LittleEndian";
#endif
		ofs << "\">" << endl;
		ofs << "<ImageData WholeExtent=\"";
		ofs << iNX1 << " ";
		ofs << iNXN << " ";
		ofs << iNY1 << " ";
		ofs << iNYN << " ";
		ofs << iNZ1 << " ";
		ofs << iNZN << "\" ";
		ofs << "Origin=\"";
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(16);
		ofs << ox << " ";
		ofs << oy << " ";
		ofs << oz << "\" ";
		ofs << "Spacing=\"";
		ofs << dx << " ";
		ofs << dx << " ";
		ofs << dx << "\">" << std::endl;
		ofs << "<Piece Extent=\"";
		ofs << iNX1 << " ";
		ofs << iNXN << " ";
		ofs << iNY1 << " ";
		ofs << iNYN << " ";
		ofs << iNZ1 << " ";
		ofs << iNZN << "\">" << std::endl;
		ofs << "<PointData>" << endl;
		ofs << "</PointData>" << endl;
		ofs << "<CellData>" << endl;

		ofs << "<DataArray type=\"Float32\" Name=\"";
		ofs << "p";
		ofs << "\" format=\"appended\" offset=\"";
		ofs << (sizeof(int) + sizeof(float)*nSize)*0;
		ofs << "\"/>" << endl;

		ofs << "<DataArray type=\"Float32\" Name=\"";
		ofs << "u";
		ofs << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"";
		ofs << (sizeof(int) + sizeof(float)*nSize)*1;
		ofs << "\"/>" << endl;

		ofs << "</CellData>" << endl;
		ofs << "<Coordinates>" << endl;
		ofs << "</Coordinates>" << endl;
		ofs << "</Piece>" << endl;
		ofs << "</ImageData>" << endl;
		ofs << "<AppendedData encoding=\"raw\">" << endl;
		ofs << "_";
		ofs.close();

		ofs.open(ossFileName.str().c_str(), ios::out | ios::app | ios::binary);

		float* pScalar = new float[nSize];
		float* pVector = new float[nSize*3];
		int nBytesS = sizeof(float)*nSize;
		int nBytesV = sizeof(float)*nSize*3;

		for(int k=0; k<NZ; k++) {
			for(int j=0; j<NY; j++) {
				for(int i=0; i<NX; i++) {
					int i0 = i + vc;
					int j0 = j + vc;
					int k0 = k + vc;
					int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
					float data = sDataP[m0];
					int m = i + NX*( j + NY*k );
					pScalar[m] = data;
				}
			}
		}
		ofs.write((char*)&nBytesS, sizeof(int));
		ofs.write((char*)pScalar, sizeof(float)*(nSize));

		for(int k=0; k<NZ; k++) {
			for(int j=0; j<NY; j++) {
				for(int i=0; i<NX; i++) {
					int i0 = i + vc;
					int j0 = j + vc;
					int k0 = k + vc;
					int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
					float data0 = sDataUX[m0];
					float data1 = sDataUY[m0];
					float data2 = sDataUZ[m0];
					int m = i + NX*( j + NY*k );
					pVector[3*m + 0] = data0;
					pVector[3*m + 1] = data1;
					pVector[3*m + 2] = data2;
				}
			}
		}
		ofs.write((char*)&nBytesV, sizeof(int));
		ofs.write((char*)pVector, sizeof(float)*(nSize*3));

		ofs.close();

		ofs.open(ossFileName.str().c_str(), ios::out | ios::app);
		ofs << endl;
		ofs << "</AppendedData>" << endl;
		ofs << "</VTKFile>" << endl;
		ofs.close();

		delete [] pScalar;
		delete [] pVector;
	}

  template <typename T>
	void printVTIC(
								T* sDataP, 
								T* sDataUX, 
								T* sDataUY, 
								T* sDataUZ, 
								T* sDataT, 
								const char* label, int step, int rank, int block,
								int NX, int NY, int NZ,
								int vc,
								double ox, double oy, double oz,
								double dx) {
		ostringstream ossFileName;
		ossFileName << "./VTK/";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << "/";
		ossFileName << "data-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << rank;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << block;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".vti";

		int iNX1 = 0;
		int iNY1 = 0;
		int iNZ1 = 0;
		int iNXN = NX;
		int iNYN = NY;
		int iNZN = NZ;

		unsigned int nSizeX = NX;
		unsigned int nSizeY = NY;
		unsigned int nSizeZ = NZ;
		unsigned int nSize = nSizeX*nSizeY*nSizeZ;

		ofstream ofs;
		ofs.open(ossFileName.str().c_str(), ios::out);
		ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"";
#ifdef __FUJITSU
		ofs << "BigEndian";
#else
		ofs << "LittleEndian";
#endif
		ofs << "\">" << endl;
		ofs << "<ImageData WholeExtent=\"";
		ofs << iNX1 << " ";
		ofs << iNXN << " ";
		ofs << iNY1 << " ";
		ofs << iNYN << " ";
		ofs << iNZ1 << " ";
		ofs << iNZN << "\" ";
		ofs << "Origin=\"";
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(16);
		ofs << ox << " ";
		ofs << oy << " ";
		ofs << oz << "\" ";
		ofs << "Spacing=\"";
		ofs << dx << " ";
		ofs << dx << " ";
		ofs << dx << "\">" << std::endl;
		ofs << "<Piece Extent=\"";
		ofs << iNX1 << " ";
		ofs << iNXN << " ";
		ofs << iNY1 << " ";
		ofs << iNYN << " ";
		ofs << iNZ1 << " ";
		ofs << iNZN << "\">" << std::endl;
		ofs << "<PointData>" << endl;
		ofs << "</PointData>" << endl;
		ofs << "<CellData>" << endl;

		ofs << "<DataArray type=\"Float32\" Name=\"";
		ofs << "p";
		ofs << "\" format=\"appended\" offset=\"";
		ofs << sizeof(int)*0 + sizeof(float)*nSize*0;
		ofs << "\"/>" << endl;

		ofs << "<DataArray type=\"Float32\" Name=\"";
		ofs << "u";
		ofs << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"";
		ofs << sizeof(int)*1 + sizeof(float)*nSize*1;
		ofs << "\"/>" << endl;

		ofs << "<DataArray type=\"Float32\" Name=\"";
		ofs << "t";
		ofs << "\" format=\"appended\" offset=\"";
		ofs << sizeof(int)*2 + sizeof(float)*nSize*4;
		ofs << "\"/>" << endl;

		ofs << "</CellData>" << endl;
		ofs << "<Coordinates>" << endl;
		ofs << "</Coordinates>" << endl;
		ofs << "</Piece>" << endl;
		ofs << "</ImageData>" << endl;
		ofs << "<AppendedData encoding=\"raw\">" << endl;
		ofs << "_";
		ofs.close();

		ofs.open(ossFileName.str().c_str(), ios::out | ios::app | ios::binary);

		float* pScalar = new float[nSize];
		float* pVector = new float[nSize*3];
		int nBytesS = sizeof(float)*nSize;
		int nBytesV = sizeof(float)*nSize*3;

		for(int k=0; k<NZ; k++) {
			for(int j=0; j<NY; j++) {
				for(int i=0; i<NX; i++) {
					int i0 = i + vc;
					int j0 = j + vc;
					int k0 = k + vc;
					int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
					float data = sDataP[m0];
					int m = i + NX*( j + NY*k );
					pScalar[m] = data;
				}
			}
		}
		ofs.write((char*)&nBytesS, sizeof(int));
		ofs.write((char*)pScalar, sizeof(float)*(nSize));

		for(int k=0; k<NZ; k++) {
			for(int j=0; j<NY; j++) {
				for(int i=0; i<NX; i++) {
					int i0 = i + vc;
					int j0 = j + vc;
					int k0 = k + vc;
					int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
					float data0 = sDataUX[m0];
					float data1 = sDataUY[m0];
					float data2 = sDataUZ[m0];
					int m = i + NX*( j + NY*k );
					pVector[3*m + 0] = data0;
					pVector[3*m + 1] = data1;
					pVector[3*m + 2] = data2;
				}
			}
		}
		ofs.write((char*)&nBytesV, sizeof(int));
		ofs.write((char*)pVector, sizeof(float)*(nSize*3));

		for(int k=0; k<NZ; k++) {
			for(int j=0; j<NY; j++) {
				for(int i=0; i<NX; i++) {
					int i0 = i + vc;
					int j0 = j + vc;
					int k0 = k + vc;
					int m0 = i0 + (NX+2*vc)*( j0 + (NY+2*vc)*k0 );
					float data = sDataT[m0];
					int m = i + NX*( j + NY*k );
					pScalar[m] = data;
				}
			}
		}
		ofs.write((char*)&nBytesS, sizeof(int));
		ofs.write((char*)pScalar, sizeof(float)*(nSize));

		ofs.close();

		ofs.open(ossFileName.str().c_str(), ios::out | ios::app);
		ofs << endl;
		ofs << "</AppendedData>" << endl;
		ofs << "</VTKFile>" << endl;
		ofs.close();

		delete [] pScalar;
		delete [] pVector;
	}


	template <typename T>
  void writeScalar(
						int dataClassID,
						int vc,
						const std::string& name,
						int step,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Vec3d rootOrigin,
						double rootLength) {
		ostringstream ossFileNameTime;
		ossFileNameTime << "./VTK/";
		mkdir(ossFileNameTime.str().c_str(), 0755);
		ossFileNameTime.width(10);
		ossFileNameTime.setf(ios::fixed);
		ossFileNameTime.fill('0');
		ossFileNameTime << step;
		mkdir(ossFileNameTime.str().c_str(), 0755);

    const ::Vec3i& size = blockManager.getSize();
    int myrank = comm.Get_rank();

    float* data = new float[(size.x) * (size.y) * (size.z)];

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
			::Vec3i size = block->getSize();
			Vec3d origin = block->getOrigin();
			Vec3d blockSize = block->getBlockSize();
			Vec3d cellSize = block->getCellSize();
			int level = block->getLevel();

      Scalar3D<T>* s = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID));
      T* sData = s->getData();
			int vc = s->getVCsize();

			printVTIC(sData, name.c_str(), step, myrank, id, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
//			printVTIP(sData, name.c_str(), step, myrank, id, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);

    }

    delete[] data;

		ostringstream ossFileName;
		ossFileName << "./VTK/";
		ossFileName << "data-";
		ossFileName << name.c_str();
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".vthb";

		if( myrank == 0 ) {
			ofstream ofs;
			ofs.open(ossFileName.str().c_str(), ios::out);
			ofs << "<VTKFile type=\"vtkHierarchicalBoxDataSet\" version=\"0.1\">" << endl;
			ofs << "<vtkHierarchicalBoxDataSet>" << endl;

			for(int n=0; n<difflevel+1; n++) {
				ofs << "<RefinementRatio level=\"";
				ofs << n;
				ofs << "\" refinement=\"2\"/>";
				ofs << endl;
			}
			ofs << endl;

			std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
			for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
				for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
					Node* node = leafNodeArray[id];
//					Vec3d origin = tree->getOrigin(node) * rootLength + rootOrigin;
					Vec3d origin = tree->getOrigin(node) * rootLength;
					Vec3d blockSize = node->getBlockSize() * rootLength;
					Vec3d cellSize;
					cellSize.x = blockSize.x / size.x;
					cellSize.y = blockSize.y / size.y;
					cellSize.z = blockSize.z / size.z;
					int level = node->getLevel();

					ostringstream ossFileName2;
					ossFileName2 << "./";
					ossFileName2.width(10);
					ossFileName2.setf(ios::fixed);
					ossFileName2.fill('0');
					ossFileName2 << step;
					ossFileName2 << "/";
					ossFileName2 << "data-";
					ossFileName2 << name.c_str();
					ossFileName2 << "-";
					ossFileName2.width(5);
					ossFileName2.setf(ios::fixed);
					ossFileName2.fill('0');
					ossFileName2 << iRank;
					ossFileName2 << "-";
					ossFileName2.width(5);
					ossFileName2.setf(ios::fixed);
					ossFileName2.fill('0');
					ossFileName2 << id - partition->getStart(iRank);
					ossFileName2 << "-";
					ossFileName2.width(10);
					ossFileName2.setf(ios::fixed);
					ossFileName2.fill('0');
					ossFileName2 << step;
					ossFileName2 << ".vti";

/*
					int nx0 = origin.x/cellSize.x;
					int ny0 = origin.y/cellSize.y;
					int nz0 = origin.z/cellSize.z;
*/

					int lx = size.x*(1 << level);
					int ly = size.y*(1 << level);
					int lz = size.z*(1 << level);

					Vec3d origin2 = tree->getOrigin(node);

					double nx0 = (origin2.x)*lx;
					double ny0 = (origin2.y)*ly;
					double nz0 = (origin2.z)*lz;

					double nx1 = nx0 + size.x - 1;
					double ny1 = ny0 + size.y - 1;
					double nz1 = nz0 + size.z - 1;

					ofs << "<DataSet group=\"";
					ofs << level;
					ofs << "\" dataset=\"";
					ofs << id + blockManager.getStartID();
					ofs << "\" amr_box=\"";
					ofs << nx0;
					ofs << " ";
					ofs << nx1;
					ofs << " ";
					ofs << ny0;
					ofs << " ";
					ofs << ny1;
					ofs << " ";
					ofs << nz0;
					ofs << " ";
					ofs << nz1;
					ofs << "\" file=\"";
					ofs << ossFileName2.str().c_str();
					ofs << "\"/>";
					ofs << endl;
				}
			}

			ofs << "</vtkHierarchicalBoxDataSet>" << endl;
			ofs << "</VTKFile>" << endl;
			ofs.close();
		}

  }

  template <typename T>
  void writePU(
						int dataClassID_P,
						int dataClassID_UX, int dataClassID_UY, int dataClassID_UZ,
						int vc,
						const std::string& name,
						int step,
						int difflevel) {
		ostringstream ossFileNameTime;
		ossFileNameTime << "./VTK/";
		mkdir(ossFileNameTime.str().c_str(), 0755);
		ossFileNameTime.width(10);
		ossFileNameTime.setf(ios::fixed);
		ossFileNameTime.fill('0');
		ossFileNameTime << step;
		mkdir(ossFileNameTime.str().c_str(), 0755);

    const ::Vec3i& size = blockManager.getSize();
    int myrank = comm.Get_rank();

    float* dataP  = new float[(size.x) * (size.y) * (size.z)];
    float* dataUX = new float[(size.x) * (size.y) * (size.z)];
    float* dataUY = new float[(size.x) * (size.y) * (size.z)];
    float* dataUZ = new float[(size.x) * (size.y) * (size.z)];

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
			::Vec3i size = block->getSize();
			Vec3d origin = block->getOrigin();
			Vec3d blockSize = block->getBlockSize();
			Vec3d cellSize = block->getCellSize();
			int level = block->getLevel();

      Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
      T* sDataP = sp->getData();
      Scalar3D<T>* sux = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UX));
      T* sDataUX = sux->getData();
      Scalar3D<T>* suy = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UY));
      T* sDataUY = suy->getData();
      Scalar3D<T>* suz = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UZ));
      T* sDataUZ = suz->getData();

			printVTIC(sDataP, sDataUX, sDataUY, sDataUZ, name.c_str(), step, myrank, id, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
//			printVTIP(sData, name.c_str(), step, myrank, id, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
    }

    delete[] dataP;
    delete[] dataUX;
    delete[] dataUY;
    delete[] dataUZ;

		ostringstream ossFileName;
		ossFileName << "./VTK/";
		ossFileName << "data-";
		ossFileName << name.c_str();
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".vthb";

		ofstream ofs;
		if( myrank == 0 ) {
			ofs.open(ossFileName.str().c_str(), ios::out);
			ofs << "<VTKFile type=\"vtkHierarchicalBoxDataSet\" version=\"0.1\">" << endl;
			ofs << "<vtkHierarchicalBoxDataSet>" << endl;

			for(int n=0; n<difflevel+1; n++) {
				ofs << "<RefinementRatio level=\"";
				ofs << n;
				ofs << "\" refinement=\"2\"/>";
				ofs << endl;
			}
			ofs << endl;
			ofs.close();
		}

    int Nprocs = comm.Get_size();
		for(int n=0; n<Nprocs; n++) {
			MPI_Barrier(MPI_COMM_WORLD);
			if( n != myrank ) {
				continue;
			}

			ofs.open(ossFileName.str().c_str(), ios::out|ios::app);
			for (int id = 0; id < blockManager.getNumBlock(); ++id) {
				BlockBase* block = blockManager.getBlock(id);
				::Vec3i size = block->getSize();
				Vec3d origin = block->getOrigin();
				Vec3d blockSize = block->getBlockSize();
				Vec3d cellSize = block->getCellSize();
				int level = block->getLevel();

				ostringstream ossFileName2;
				ossFileName2 << "./";
				ossFileName2.width(10);
				ossFileName2.setf(ios::fixed);
				ossFileName2.fill('0');
				ossFileName2 << step;
				ossFileName2 << "/";
				ossFileName2 << "data-";
				ossFileName2 << name.c_str();
				ossFileName2 << "-";
				ossFileName2.width(5);
				ossFileName2.setf(ios::fixed);
				ossFileName2.fill('0');
				ossFileName2 << myrank;
				ossFileName2 << "-";
				ossFileName2.width(5);
				ossFileName2.setf(ios::fixed);
				ossFileName2.fill('0');
				ossFileName2 << id;
				ossFileName2 << "-";
				ossFileName2.width(10);
				ossFileName2.setf(ios::fixed);
				ossFileName2.fill('0');
				ossFileName2 << step;
				ossFileName2 << ".vti";

				int nx0 = origin.x/cellSize.x;
				int ny0 = origin.y/cellSize.y;
				int nz0 = origin.z/cellSize.z;
				ofs << "<DataSet group=\"";
				ofs << level;
				ofs << "\" dataset=\"";
				ofs << id + blockManager.getStartID();
				ofs << "\" amr_box=\"";
				ofs << nx0;
				ofs << " ";
				ofs << nx0 + size[0] - 1;
				ofs << " ";
				ofs << ny0;
				ofs << " ";
				ofs << ny0 + size[1] - 1;
				ofs << " ";
				ofs << nz0;
				ofs << " ";
				ofs << nz0 + size[2] - 1;
				ofs << "\" file=\"";
				ofs << ossFileName2.str().c_str();
				ofs << "\"/>";
				ofs << endl;
			}
			ofs.close();
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if( myrank == 0 ) {
			ofs.open(ossFileName.str().c_str(), ios::out|ios::app);
			ofs << "</vtkHierarchicalBoxDataSet>" << endl;
			ofs << "</VTKFile>" << endl;
			ofs.close();
		}
  }

  template <typename T>
  void writePUT(
						int dataClassID_P,
						int dataClassID_UX, int dataClassID_UY, int dataClassID_UZ,
						int dataClassID_T,
						int vc,
						const std::string& name,
						int step,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Vec3d rootOrigin,
						double rootLength) {
		ostringstream ossFileNameTime;
		ossFileNameTime << "./VTK/";
		mkdir(ossFileNameTime.str().c_str(), 0755);
		ossFileNameTime.width(10);
		ossFileNameTime.setf(ios::fixed);
		ossFileNameTime.fill('0');
		ossFileNameTime << step;
		mkdir(ossFileNameTime.str().c_str(), 0755);

    const ::Vec3i& size = blockManager.getSize();
    int myrank = comm.Get_rank();

    float* dataP  = new float[(size.x) * (size.y) * (size.z)];
    float* dataUX = new float[(size.x) * (size.y) * (size.z)];
    float* dataUY = new float[(size.x) * (size.y) * (size.z)];
    float* dataUZ = new float[(size.x) * (size.y) * (size.z)];
    float* dataT  = new float[(size.x) * (size.y) * (size.z)];

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
			::Vec3i size = block->getSize();
			Vec3d origin = block->getOrigin();
			Vec3d blockSize = block->getBlockSize();
			Vec3d cellSize = block->getCellSize();
			int level = block->getLevel();

      Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
      T* sDataP = sp->getData();
      Scalar3D<T>* sux = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UX));
      T* sDataUX = sux->getData();
      Scalar3D<T>* suy = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UY));
      T* sDataUY = suy->getData();
      Scalar3D<T>* suz = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UZ));
      T* sDataUZ = suz->getData();
      Scalar3D<T>* st = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_T));
      T* sDataT = st->getData();

			printVTIC(sDataP, sDataUX, sDataUY, sDataUZ, sDataT, name.c_str(), step, myrank, id, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
//			printVTIP(sData, name.c_str(), step, myrank, id, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
    }

    delete[] dataP;
    delete[] dataUX;
    delete[] dataUY;
    delete[] dataUZ;
    delete[] dataT;

		ostringstream ossFileName;
		ossFileName << "./VTK/";
		ossFileName << "data-";
		ossFileName << name.c_str();
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".vthb";

		if( myrank == 0 ) {
			ofstream ofs;
			ofs.open(ossFileName.str().c_str(), ios::out);
			ofs << "<VTKFile type=\"vtkHierarchicalBoxDataSet\" version=\"0.1\">" << endl;
			ofs << "<vtkHierarchicalBoxDataSet>" << endl;

			for(int n=0; n<difflevel+1; n++) {
				ofs << "<RefinementRatio level=\"";
				ofs << n;
				ofs << "\" refinement=\"2\"/>";
				ofs << endl;
			}
			ofs << endl;

			std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
			for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
				for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
					Node* node = leafNodeArray[id];
//					Vec3d origin = tree->getOrigin(node) * rootLength + rootOrigin;
					Vec3d origin = tree->getOrigin(node) * rootLength;
					Vec3d blockSize = node->getBlockSize() * rootLength;
					Vec3d cellSize;
					cellSize.x = blockSize.x / size.x;
					cellSize.y = blockSize.y / size.y;
					cellSize.z = blockSize.z / size.z;
					int level = node->getLevel();

					ostringstream ossFileName2;
					ossFileName2 << "./";
					ossFileName2.width(10);
					ossFileName2.setf(ios::fixed);
					ossFileName2.fill('0');
					ossFileName2 << step;
					ossFileName2 << "/";
					ossFileName2 << "data-";
					ossFileName2 << name.c_str();
					ossFileName2 << "-";
					ossFileName2.width(5);
					ossFileName2.setf(ios::fixed);
					ossFileName2.fill('0');
					ossFileName2 << iRank;
					ossFileName2 << "-";
					ossFileName2.width(5);
					ossFileName2.setf(ios::fixed);
					ossFileName2.fill('0');
					ossFileName2 << id - partition->getStart(iRank);
					ossFileName2 << "-";
					ossFileName2.width(10);
					ossFileName2.setf(ios::fixed);
					ossFileName2.fill('0');
					ossFileName2 << step;
					ossFileName2 << ".vti";

/*
					int nx0 = origin.x/cellSize.x;
					int ny0 = origin.y/cellSize.y;
					int nz0 = origin.z/cellSize.z;
*/

					int lx = size.x*(1 << level);
					int ly = size.y*(1 << level);
					int lz = size.z*(1 << level);

					Vec3d origin2 = tree->getOrigin(node);

					double nx0 = (origin2.x)*lx;
					double ny0 = (origin2.y)*ly;
					double nz0 = (origin2.z)*lz;

					double nx1 = nx0 + size.x - 1;
					double ny1 = ny0 + size.y - 1;
					double nz1 = nz0 + size.z - 1;

					ofs << "<DataSet group=\"";
					ofs << level;
					ofs << "\" dataset=\"";
					ofs << id + blockManager.getStartID();
					ofs << "\" amr_box=\"";
					ofs << nx0;
					ofs << " ";
					ofs << nx1;
					ofs << " ";
					ofs << ny0;
					ofs << " ";
					ofs << ny1;
					ofs << " ";
					ofs << nz0;
					ofs << " ";
					ofs << nz1;
//					ofs << " " << lx << " " << ly << " " << lz << " ";
					ofs << "\" file=\"";
					ofs << ossFileName2.str().c_str();
					ofs << "\"/>";
					ofs << endl;
				}
			}

			ofs << "</vtkHierarchicalBoxDataSet>" << endl;
			ofs << "</VTKFile>" << endl;
			ofs.close();
		}
  }

  template <typename T>
  void writeVtkOverlappingAMR_PUT(
						int dataClassID_P,
						int dataClassID_UX, int dataClassID_UY, int dataClassID_UZ,
						int dataClassID_T,
						int vc,
						const std::string& name,
						int step,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Vec3d rootOrigin,
						double rootLength) {
		ostringstream ossFileNameTime;
		ossFileNameTime << "./VTK/";
		mkdir(ossFileNameTime.str().c_str(), 0755);
		ossFileNameTime.width(10);
		ossFileNameTime.setf(ios::fixed);
		ossFileNameTime.fill('0');
		ossFileNameTime << step;
		mkdir(ossFileNameTime.str().c_str(), 0755);

    const ::Vec3i& size = blockManager.getSize();
    int myrank = comm.Get_rank();

    float* dataP  = new float[(size.x) * (size.y) * (size.z)];
    float* dataUX = new float[(size.x) * (size.y) * (size.z)];
    float* dataUY = new float[(size.x) * (size.y) * (size.z)];
    float* dataUZ = new float[(size.x) * (size.y) * (size.z)];
    float* dataT  = new float[(size.x) * (size.y) * (size.z)];

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
			::Vec3i size = block->getSize();
			Vec3d origin = block->getOrigin();
			Vec3d blockSize = block->getBlockSize();
			Vec3d cellSize = block->getCellSize();
			int level = block->getLevel();

      Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
      T* sDataP = sp->getData();
      Scalar3D<T>* sux = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UX));
      T* sDataUX = sux->getData();
      Scalar3D<T>* suy = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UY));
      T* sDataUY = suy->getData();
      Scalar3D<T>* suz = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UZ));
      T* sDataUZ = suz->getData();
      Scalar3D<T>* st = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_T));
      T* sDataT = st->getData();

//			printVTIC(sDataP, sDataUX, sDataUY, sDataUZ, sDataT, name.c_str(), step, myrank, id, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
//			printVTIP(sData, name.c_str(), step, myrank, id, size[0], size[1], size[2], vc, origin[0], origin[1], origin[2], cellSize[0]);
    }

    delete[] dataP;
    delete[] dataUX;
    delete[] dataUY;
    delete[] dataUZ;
    delete[] dataT;

		ostringstream ossFileName;
		ossFileName << "./VTK/";
		ossFileName << "data-";
		ossFileName << name.c_str();
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".vth";

		if( myrank == 0 ) {
			ofstream ofs;
			ofs.open(ossFileName.str().c_str(), ios::out);
			ofs << "<VTKFile type=\"vtkOverlappingAMR\" version=\"1.1\">" << endl;
			ofs << "<vtkOverlappingAMR origin=\"0 0 0\" grid_description=\"XYZ\">" << endl;

			std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
			int lmax = difflevel+1;
			for(int n=0; n<lmax; n++) {
				double dx = rootLength/(1 << n);
				ofs << "\t<Block level=\"";
				ofs << n;
				ofs << "\" spacing=\"";
				ofs << dx;
				ofs << " ";
				ofs << dx;
				ofs << " ";
				ofs << dx;
				ofs << "\">";
				ofs << endl;
				for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
					for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
						Node* node = leafNodeArray[id];
						Vec3d origin = tree->getOrigin(node) * rootLength;
						Vec3d blockSize = node->getBlockSize() * rootLength;
						Vec3d cellSize;
						cellSize.x = blockSize.x / size.x;
						cellSize.y = blockSize.y / size.y;
						cellSize.z = blockSize.z / size.z;
						int level = node->getLevel();

						ostringstream ossFileName2;
						ossFileName2 << "./";
						ossFileName2.width(10);
						ossFileName2.setf(ios::fixed);
						ossFileName2.fill('0');
						ossFileName2 << step;
						ossFileName2 << "/";
						ossFileName2 << "data-";
						ossFileName2 << name.c_str();
						ossFileName2 << "-";
						ossFileName2.width(5);
						ossFileName2.setf(ios::fixed);
						ossFileName2.fill('0');
						ossFileName2 << iRank;
						ossFileName2 << "-";
						ossFileName2.width(5);
						ossFileName2.setf(ios::fixed);
						ossFileName2.fill('0');
						ossFileName2 << id - partition->getStart(iRank);
						ossFileName2 << "-";
						ossFileName2.width(10);
						ossFileName2.setf(ios::fixed);
						ossFileName2.fill('0');
						ossFileName2 << step;
						ossFileName2 << ".vti";

						int lx = size.x*(1 << level);
						int ly = size.y*(1 << level);
						int lz = size.z*(1 << level);

						Vec3d origin2 = tree->getOrigin(node);

						double nx0 = (origin2.x)*lx;
						double ny0 = (origin2.y)*ly;
						double nz0 = (origin2.z)*lz;

						double nx1 = nx0 + size.x - 1;
						double ny1 = ny0 + size.y - 1;
						double nz1 = nz0 + size.z - 1;

						if( level == n ) {
							ofs << "\t\t<DataSet index=\"";
							ofs << id + blockManager.getStartID();
							ofs << "\" amr_box=\"";
							ofs << nx0;
							ofs << " ";
							ofs << nx1;
							ofs << " ";
							ofs << ny0;
							ofs << " ";
							ofs << ny1;
							ofs << " ";
							ofs << nz0;
							ofs << " ";
							ofs << nz1;
							ofs << "\" file=\"";
							ofs << ossFileName2.str().c_str();
							ofs << "\">";
							ofs << endl;
							ofs << "\t\t</DataSet>" << endl;	
						}
					}
				}
				ofs << "\t</Block>";
				ofs << endl;
				ofs << endl;
			}

			ofs << "</vtkOverlappingAMR>" << endl;
			ofs << "</VTKFile>" << endl;
			ofs.close();
		}
  }

  template <typename T>
  void writeVtkOverlappingAMR_PUT_0(
						const std::string& name,
						int step,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Vec3d rootOrigin,
						double rootLength) {
    const ::Vec3i& size = blockManager.getSize();
    int myrank = comm.Get_rank();

		ostringstream ossFileName;
		ossFileName << "data-";
		ossFileName << name.c_str();
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".vth";

		if( myrank == 0 ) {
			ofstream ofs;
			ofs.open(ossFileName.str().c_str(), ios::out);
			ofs << "<VTKFile type=\"vtkOverlappingAMR\" version=\"1.1\">" << endl;
			ofs << "<vtkOverlappingAMR origin=\"0 0 0\" grid_description=\"XYZ\">" << endl;

			std::vector<Node*>& leafNodeArray = tree->getLeafNodeArray();
			int lmax = difflevel+1;
			for(int n=0; n<lmax; n++) {
				double dx = rootLength/(1 << n);
				ofs << "\t<Block level=\"";
				ofs << n;
				ofs << "\" spacing=\"";
				ofs << dx;
				ofs << " ";
				ofs << dx;
				ofs << " ";
				ofs << dx;
				ofs << "\">";
				ofs << endl;
				for (int iRank = 0; iRank < comm.Get_size(); iRank++) {
					for (int id = partition->getStart(iRank); id < partition->getEnd(iRank); id++) {
						Node* node = leafNodeArray[id];
						Vec3d origin = tree->getOrigin(node) * rootLength;
						Vec3d blockSize = node->getBlockSize() * rootLength;
						Vec3d cellSize;
						cellSize.x = blockSize.x / size.x;
						cellSize.y = blockSize.y / size.y;
						cellSize.z = blockSize.z / size.z;
						int level = node->getLevel();

						ostringstream ossFileName2;
						ossFileName2 << "./VTK/";
						ossFileName2.width(10);
						ossFileName2.setf(ios::fixed);
						ossFileName2.fill('0');
						ossFileName2 << step;
						ossFileName2 << "/";
						ossFileName2 << "data-";
						ossFileName2 << name.c_str();
						ossFileName2 << "-";
						ossFileName2.width(5);
						ossFileName2.setf(ios::fixed);
						ossFileName2.fill('0');
						ossFileName2 << iRank;
						ossFileName2 << "-";
						ossFileName2.width(5);
						ossFileName2.setf(ios::fixed);
						ossFileName2.fill('0');
						ossFileName2 << id - partition->getStart(iRank);
						ossFileName2 << "-";
						ossFileName2.width(10);
						ossFileName2.setf(ios::fixed);
						ossFileName2.fill('0');
						ossFileName2 << step;
						ossFileName2 << ".vti";

						int lx = size.x*(1 << level);
						int ly = size.y*(1 << level);
						int lz = size.z*(1 << level);

						Vec3d origin2 = tree->getOrigin(node);

						double nx0 = (origin2.x)*lx;
						double ny0 = (origin2.y)*ly;
						double nz0 = (origin2.z)*lz;

						double nx1 = nx0 + size.x - 1;
						double ny1 = ny0 + size.y - 1;
						double nz1 = nz0 + size.z - 1;

						if( level == n ) {
							ofs << "\t\t<DataSet index=\"";
							ofs << id + blockManager.getStartID();
							ofs << "\" amr_box=\"";
							ofs << nx0;
							ofs << " ";
							ofs << nx1;
							ofs << " ";
							ofs << ny0;
							ofs << " ";
							ofs << ny1;
							ofs << " ";
							ofs << nz0;
							ofs << " ";
							ofs << nz1;
							ofs << "\" file=\"";
							ofs << ossFileName2.str().c_str();
							ofs << "\">";
							ofs << endl;
							ofs << "\t\t</DataSet>" << endl;	
						}
					}
				}
				ofs << "\t</Block>";
				ofs << endl;
				ofs << endl;
			}

			ofs << "</vtkOverlappingAMR>" << endl;
			ofs << "</VTKFile>" << endl;
			ofs.close();
		}
  }


private:
};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // VTK_WRITER_H

