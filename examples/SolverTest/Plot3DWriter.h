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
/// @file Plot3DWriter.h
/// @brief Plot3DフォーマットによるBCMデータ出力クラス
///
/// @note
///
///

#ifndef PLOT3D_WRITER_H
#define PLOT3D_WRITER_H

#include "mpi.h"
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>

#include "BCMTools.h"
#include "BlockManager.h"
#include "Scalar3D.h"
#include "Vector3D.h"
#include "BCMOctree.h"
#include "Partition.h"

#include "bplt3d.h"

using namespace std;

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif


class Plot3DWriter {
  BlockManager& blockManager;
  const MPI::Intracomm& comm;
//  const std::string& fileName;

public:
	Plot3DWriter()
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
  ~Plot3DWriter() {
  }

  void writeXYZ(
						int dataClassID_Blank,
						int vc,
						const std::string& name,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Vec3r rootOrigin,
						REAL_TYPE rootLength) {
		ostringstream ossFileNameTime;
		ossFileNameTime << "./PLOT3D/";
		mkdir(ossFileNameTime.str().c_str(), 0755);
		ossFileNameTime << "./XYZ/";
		mkdir(ossFileNameTime.str().c_str(), 0755);

    const ::Vec3i& size = blockManager.getSize();
    int myrank = comm.Get_rank();

		ostringstream ossFileName;
		ossFileName << "./PLOT3D/";
		ossFileName << "./XYZ/";
		ossFileName << "data-";
		ossFileName << name.c_str();
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << myrank;
		ossFileName << ".xyz";

/*
		std::ofstream ofs;
		ofs.open(ossFileName.str().c_str(), std::ios::out);
*/

		int unit = 55;
		int filenamelength = strlen(ossFileName.str().c_str());
		bplt3d_open_file_((char*)ossFileName.str().c_str(), &filenamelength, &unit);

		BlockBase* block0 = blockManager.getBlock(0);
		::Vec3i size0 = block0->getSize();
		int ix = size0.x;
		int jx = size0.y;
		int kx = size0.z;
		int ngrid = blockManager.getNumBlock();

/*
		ofs.write((char*)&ngrid, sizeof(int));
    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
			ofs.write((char*)&ix, sizeof(int));
			ofs.write((char*)&jx, sizeof(int));
			ofs.write((char*)&kx, sizeof(int));
		}
*/

		bplt3d_write_xyz_header_(&ix, &jx, &kx, &ngrid, &unit);

		float* px = new float [ix*jx*kx];
		float* py = new float [ix*jx*kx];
		float* pz = new float [ix*jx*kx];
		int* blank = new int [ix*jx*kx];
    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
			::Vec3i size = block->getSize();
			Vec3r origin = block->getOrigin();
			Vec3r blockSize = block->getBlockSize();
			Vec3r cellSize = block->getCellSize();
			int level = block->getLevel();
			REAL_TYPE dx = cellSize.x;

      Scalar3D<int>* sblank = dynamic_cast<Scalar3D<int>*>(block->getDataClass(dataClassID_Blank));
      int* sDataBlank = sblank->getData();

			int ix0 = ix + 2*vc;
			int jx0 = jx + 2*vc;
			int kx0 = kx + 2*vc;
			for(int k=0; k<kx; k++) {
				for(int j=0; j<jx; j++) {
					for(int i=0; i<ix; i++) {
						int m = i + ix*(j + jx*k);
						int m0 = i + vc + ix0*(j + vc + jx0*(k + vc));
						px[m] = origin.x + (i + 0.5)*dx;
						py[m] = origin.y + (j + 0.5)*dx;
						pz[m] = origin.z + (k + 0.5)*dx;
						blank[m] = (sDataBlank[m0] == 1 ? 1 : 0);
					}
				}
			}
/*
			ofs.write((char*)px, sizeof(float)*ix*jx*kx);
			ofs.write((char*)py, sizeof(float)*ix*jx*kx);
			ofs.write((char*)pz, sizeof(float)*ix*jx*kx);
*/
			bplt3d_write_xyz_block_(px, py, pz, &ix, &jx, &kx, &unit);
    }

/*
		ofs.close();
*/
		bplt3d_close_file_(&unit);

		delete [] px;
		delete [] py;
		delete [] pz;
		delete [] blank;
  }

  void writeXYZ_binary(
						int dataClassID_Blank,
						int vc,
						const std::string& name,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Vec3r rootOrigin,
						REAL_TYPE rootLength) {
		ostringstream ossFileNameTime;
		ossFileNameTime << "./PLOT3D/";
		mkdir(ossFileNameTime.str().c_str(), 0755);
		ossFileNameTime << "./XYZ/";
		mkdir(ossFileNameTime.str().c_str(), 0755);

    const ::Vec3i& size = blockManager.getSize();
    int myrank = comm.Get_rank();

		ostringstream ossFileName;
		ossFileName << "./PLOT3D/";
		ossFileName << "./XYZ/";
		ossFileName << "data-";
		ossFileName << name.c_str();
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << myrank;
		ossFileName << ".xyz";

		std::ofstream ofs;
		ofs.open(ossFileName.str().c_str(), std::ios::out);

/*
		int unit = 55;
		int filenamelength = strlen(ossFileName.str().c_str());
		bplt3d_open_file_((char*)ossFileName.str().c_str(), &filenamelength, &unit);
*/

		BlockBase* block0 = blockManager.getBlock(0);
		::Vec3i size0 = block0->getSize();
		int ix = size0.x;
		int jx = size0.y;
		int kx = size0.z;
		int ngrid = blockManager.getNumBlock();

		ofs.write((char*)&ngrid, sizeof(int));
    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
			ofs.write((char*)&ix, sizeof(int));
			ofs.write((char*)&jx, sizeof(int));
			ofs.write((char*)&kx, sizeof(int));
		}

/*
		bplt3d_write_xyz_header_(&ix, &jx, &kx, &ngrid, &unit);
*/

		float* px = new float [ix*jx*kx];
		float* py = new float [ix*jx*kx];
		float* pz = new float [ix*jx*kx];
		int* blank = new int [ix*jx*kx];
    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
			::Vec3i size = block->getSize();
			Vec3r origin = block->getOrigin();
			Vec3r blockSize = block->getBlockSize();
			Vec3r cellSize = block->getCellSize();
			int level = block->getLevel();
			REAL_TYPE dx = cellSize.x;

      Scalar3D<int>* sblank = dynamic_cast<Scalar3D<int>*>(block->getDataClass(dataClassID_Blank));
      int* sDataBlank = sblank->getData();

			int ix0 = ix + 2*vc;
			int jx0 = jx + 2*vc;
			int kx0 = kx + 2*vc;
			for(int k=0; k<kx; k++) {
				for(int j=0; j<jx; j++) {
					for(int i=0; i<ix; i++) {
						int m = i + ix*(j + jx*k);
						int m0 = i + vc + ix0*(j + vc + jx0*(k + vc));
						px[m] = origin.x + (i + 0.5)*dx;
						py[m] = origin.y + (j + 0.5)*dx;
						pz[m] = origin.z + (k + 0.5)*dx;
						blank[m] = (sDataBlank[m0] == 1 ? 1 : 0);
					}
				}
			}
			ofs.write((char*)px, sizeof(float)*ix*jx*kx);
			ofs.write((char*)py, sizeof(float)*ix*jx*kx);
			ofs.write((char*)pz, sizeof(float)*ix*jx*kx);
//			ofs.write((char*)blank, sizeof(int)*ix*jx*kx);
/*
			bplt3d_write_xyz_block_(px, py, pz, &ix, &jx, &kx, &unit);
*/
    }

		ofs.close();
/*
		bplt3d_close_file_(&unit);
*/

		delete [] px;
		delete [] py;
		delete [] pz;
		delete [] blank;
  }

  template <typename T>
  void writeData(
						int dataClassID_P,
						int dataClassID_UX, int dataClassID_UY, int dataClassID_UZ,
						int vc,
						const std::string& name,
						int step,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Vec3r rootOrigin,
						REAL_TYPE rootLength) {
		ostringstream ossFileNameTime;
		ossFileNameTime << "./PLOT3D/";
		mkdir(ossFileNameTime.str().c_str(), 0755);
		ossFileNameTime.width(10);
		ossFileNameTime.setf(ios::fixed);
		ossFileNameTime.fill('0');
		ossFileNameTime << step;
		mkdir(ossFileNameTime.str().c_str(), 0755);

    const ::Vec3i& size = blockManager.getSize();
    int myrank = comm.Get_rank();

		ostringstream ossFileName;
		ossFileName << "./PLOT3D/";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << "/";
		ossFileName << "data-";
		ossFileName << name.c_str();
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << myrank;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".func";

/*
		std::ofstream ofs;
		ofs.open(ossFileName.str().c_str(), std::ios::out);
*/

		int unit = 55;
		int filenamelength = strlen(ossFileName.str().c_str());
		bplt3d_open_file_((char*)ossFileName.str().c_str(), &filenamelength, &unit);

		BlockBase* block0 = blockManager.getBlock(0);
		::Vec3i size0 = block0->getSize();
		int ix = size0.x;
		int jx = size0.y;
		int kx = size0.z;
		int nvar = 4;
		int ngrid = blockManager.getNumBlock();

/*
		ofs.write((char*)&ngrid, sizeof(int));
    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
			ofs.write((char*)&ix, sizeof(int));
			ofs.write((char*)&jx, sizeof(int));
			ofs.write((char*)&kx, sizeof(int));
			ofs.write((char*)&nvar, sizeof(int));
		}
*/

		bplt3d_write_func_header_(&ix, &jx, &kx, &nvar, &ngrid, &unit);

    float* dataP  = new float[ix*jx*kx];
    float* dataUX = new float[ix*jx*kx];
    float* dataUY = new float[ix*jx*kx];
    float* dataUZ = new float[ix*jx*kx];
    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
			::Vec3i size = block->getSize();
			Vec3r origin = block->getOrigin();
			Vec3r blockSize = block->getBlockSize();
			Vec3r cellSize = block->getCellSize();
			int level = block->getLevel();
			REAL_TYPE dx = cellSize.x;

      Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
      T* sDataP = sp->getData();
      Scalar3D<T>* sux = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UX));
      T* sDataUX = sux->getData();
      Scalar3D<T>* suy = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UY));
      T* sDataUY = suy->getData();
      Scalar3D<T>* suz = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UZ));
      T* sDataUZ = suz->getData();

			int ix0 = ix + 2*vc;
			int jx0 = jx + 2*vc;
			int kx0 = kx + 2*vc;
			for(int k=0; k<kx; k++) {
				for(int j=0; j<jx; j++) {
					for(int i=0; i<ix; i++) {
						int m = i + ix*(j + jx*k);
						int m0 = i + vc + ix0*(j + vc + jx0*(k + vc));
						dataP[m]  = sDataP[m0];
						dataUX[m] = sDataUX[m0];
						dataUY[m] = sDataUY[m0];
						dataUZ[m] = sDataUZ[m0];
					}
				}
			}
/*
			ofs.write((char*)dataP,  sizeof(float)*ix*jx*kx);
			ofs.write((char*)dataUX, sizeof(float)*ix*jx*kx);
			ofs.write((char*)dataUY, sizeof(float)*ix*jx*kx);
			ofs.write((char*)dataUZ, sizeof(float)*ix*jx*kx);
*/
			bplt3d_write_func_block_(dataP, dataUX, dataUY, dataUZ, &ix, &jx, &kx, &unit);
    }

/*
		ofs.close();
*/
		bplt3d_close_file_(&unit);

		delete [] dataP;
		delete [] dataUX;
		delete [] dataUY;
		delete [] dataUZ;
  }

  template <typename T>
  void writeData_binary(
						int dataClassID_P,
						int dataClassID_UX, int dataClassID_UY, int dataClassID_UZ,
						int vc,
						const std::string& name,
						int step,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Vec3r rootOrigin,
						REAL_TYPE rootLength) {
		ostringstream ossFileNameTime;
		ossFileNameTime << "./PLOT3D/";
		mkdir(ossFileNameTime.str().c_str(), 0755);
		ossFileNameTime.width(10);
		ossFileNameTime.setf(ios::fixed);
		ossFileNameTime.fill('0');
		ossFileNameTime << step;
		mkdir(ossFileNameTime.str().c_str(), 0755);

    const ::Vec3i& size = blockManager.getSize();
    int myrank = comm.Get_rank();

		ostringstream ossFileName;
		ossFileName << "./PLOT3D/";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << "/";
		ossFileName << "data-";
		ossFileName << name.c_str();
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << myrank;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".func";

		std::ofstream ofs;
		ofs.open(ossFileName.str().c_str(), std::ios::out);

		BlockBase* block0 = blockManager.getBlock(0);
		::Vec3i size0 = block0->getSize();
		int ix = size0.x;
		int jx = size0.y;
		int kx = size0.z;
		int nvar = 4;
		int ngrid = blockManager.getNumBlock();
		ofs.write((char*)&ngrid, sizeof(int));

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
			ofs.write((char*)&ix, sizeof(int));
			ofs.write((char*)&jx, sizeof(int));
			ofs.write((char*)&kx, sizeof(int));
			ofs.write((char*)&nvar, sizeof(int));
		}

    float* dataP  = new float[ix*jx*kx];
    float* dataUX = new float[ix*jx*kx];
    float* dataUY = new float[ix*jx*kx];
    float* dataUZ = new float[ix*jx*kx];
    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
			::Vec3i size = block->getSize();
			Vec3r origin = block->getOrigin();
			Vec3r blockSize = block->getBlockSize();
			Vec3r cellSize = block->getCellSize();
			int level = block->getLevel();
			REAL_TYPE dx = cellSize.x;

      Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
      T* sDataP = sp->getData();
      Scalar3D<T>* sux = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UX));
      T* sDataUX = sux->getData();
      Scalar3D<T>* suy = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UY));
      T* sDataUY = suy->getData();
      Scalar3D<T>* suz = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UZ));
      T* sDataUZ = suz->getData();

			int ix0 = ix + 2*vc;
			int jx0 = jx + 2*vc;
			int kx0 = kx + 2*vc;
			for(int k=0; k<kx; k++) {
				for(int j=0; j<jx; j++) {
					for(int i=0; i<ix; i++) {
						int m = i + ix*(j + jx*k);
						int m0 = i + vc + ix0*(j + vc + jx0*(k + vc));
						dataP[m]  = sDataP[m0];
						dataUX[m] = sDataUX[m0];
						dataUY[m] = sDataUY[m0];
						dataUZ[m] = sDataUZ[m0];
					}
				}
			}
			float mach[1] = {0.1};
			float alpha[1] = {0.0};
			float re[1] = {100.0};
			float time[1] = {1.0};
			ofs.write((char*)dataP,  sizeof(float)*ix*jx*kx);
			ofs.write((char*)dataUX, sizeof(float)*ix*jx*kx);
			ofs.write((char*)dataUY, sizeof(float)*ix*jx*kx);
			ofs.write((char*)dataUZ, sizeof(float)*ix*jx*kx);
    }
		ofs.close();
		delete [] dataP;
		delete [] dataUX;
		delete [] dataUY;
		delete [] dataUZ;
  }

  template <typename T>
  void writeDataQ(
						int dataClassID_P,
						int dataClassID_UX, int dataClassID_UY, int dataClassID_UZ,
						int vc,
						const std::string& name,
						int step,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Vec3r rootOrigin,
						REAL_TYPE rootLength) {
		ostringstream ossFileNameTime;
		ossFileNameTime << "./PLOT3D/";
		mkdir(ossFileNameTime.str().c_str(), 0755);
		ossFileNameTime.width(10);
		ossFileNameTime.setf(ios::fixed);
		ossFileNameTime.fill('0');
		ossFileNameTime << step;
		mkdir(ossFileNameTime.str().c_str(), 0755);

    const ::Vec3i& size = blockManager.getSize();
    int myrank = comm.Get_rank();

		ostringstream ossFileName;
		ossFileName << "./PLOT3D/";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << "/";
		ossFileName << "data-";
		ossFileName << name.c_str();
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << myrank;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".q";

		std::ofstream ofs;
		ofs.open(ossFileName.str().c_str(), std::ios::out);

		BlockBase* block0 = blockManager.getBlock(0);
		::Vec3i size0 = block0->getSize();
		int ix = size0.x;
		int jx = size0.y;
		int kx = size0.z;
		int nvar = 5;
		int ngrid = blockManager.getNumBlock();
		ofs.write((char*)&ngrid, sizeof(int));

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
			ofs.write((char*)&ix, sizeof(int));
			ofs.write((char*)&jx, sizeof(int));
			ofs.write((char*)&kx, sizeof(int));
//			ofs.write((char*)&nvar, sizeof(int));
		}

    float* dataP  = new float[ix*jx*kx];
    float* dataUX = new float[ix*jx*kx];
    float* dataUY = new float[ix*jx*kx];
    float* dataUZ = new float[ix*jx*kx];
    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
			::Vec3i size = block->getSize();
			Vec3r origin = block->getOrigin();
			Vec3r blockSize = block->getBlockSize();
			Vec3r cellSize = block->getCellSize();
			int level = block->getLevel();
			REAL_TYPE dx = cellSize.x;

      Scalar3D<T>* sp = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_P));
      T* sDataP = sp->getData();
      Scalar3D<T>* sux = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UX));
      T* sDataUX = sux->getData();
      Scalar3D<T>* suy = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UY));
      T* sDataUY = suy->getData();
      Scalar3D<T>* suz = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID_UZ));
      T* sDataUZ = suz->getData();

			int ix0 = ix + 2*vc;
			int jx0 = jx + 2*vc;
			int kx0 = kx + 2*vc;
			for(int k=0; k<kx; k++) {
				for(int j=0; j<jx; j++) {
					for(int i=0; i<ix; i++) {
						int m = i + ix*(j + jx*k);
						int m0 = i + vc + ix0*(j + vc + jx0*(k + vc));
						dataP[m]  = sDataP[m0];
						dataUX[m] = sDataUX[m0];
						dataUY[m] = sDataUY[m0];
						dataUZ[m] = sDataUZ[m0];
					}
				}
			}
			float mach[1] = {0.1};
			float alpha[1] = {0.0};
			float re[1] = {100.0};
			float time[1] = {1.0};
			ofs.write((char*)mach, sizeof(float));
			ofs.write((char*)alpha, sizeof(float));
			ofs.write((char*)re, sizeof(float));
			ofs.write((char*)time, sizeof(float));
			ofs.write((char*)dataP,  sizeof(float)*ix*jx*kx);
			ofs.write((char*)dataUX, sizeof(float)*ix*jx*kx);
			ofs.write((char*)dataUY, sizeof(float)*ix*jx*kx);
			ofs.write((char*)dataUZ, sizeof(float)*ix*jx*kx);
			ofs.write((char*)dataP,  sizeof(float)*ix*jx*kx);
    }
		ofs.close();
		delete [] dataP;
		delete [] dataUX;
		delete [] dataUY;
		delete [] dataUZ;
  }


};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // PLOT3D_WRITER_H
