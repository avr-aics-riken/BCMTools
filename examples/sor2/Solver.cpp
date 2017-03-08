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

#include "Solver.h"
#include "Scalar3DUpdater.h"
#include "Timing.h"

#include "SiloWriter.h"


Solver::Solver(const Config& conf, const std::vector<REAL_TYPE>& boundaryValue)
  : blockManager(BlockManager::getInstance()),
    comm(blockManager.getCommunicator()),
    vc(conf.vc),
    boundaryValue(boundaryValue),
    nLoopInner(conf.nLoopInner), nLoopOuter(conf.nLoopOuter),
    omega(conf.omega), separateVCUpdate(conf.separate),
    outputFileName(conf.output)
{
  size = blockManager.getSize();
  nx = size[0];
  ny = size[1];
  nz = size[2];

  // データクラス変数<f>の生成・登録
  id_f = blockManager.setDataClass<Scalar3D<REAL_TYPE>, Scalar3DUpdater<REAL_TYPE> >(vc);

  // データクラス変数<f>の仮想セル同期準備
  blockManager.prepareForVCUpdate(id_f, tag_f, separateVCUpdate);

  // 作業用データクラス変数<work>
  work = new Scalar3D<REAL_TYPE>(size, vc);   // ソース項に使用
}



Solver::~Solver()
{
  delete work;
}


void Solver::initialize()
{
  setInitialCondition();

  if (separateVCUpdate) {
    blockManager.beginUpdateVC_X(id_f);
    blockManager.beginUpdateVC_Y(id_f);
    blockManager.beginUpdateVC_Z(id_f);
    setBoundaryCondition();
    blockManager.endUpdateVC_X(id_f);
    blockManager.endUpdateVC_Y(id_f);
    blockManager.endUpdateVC_Z(id_f);
  } else {
    blockManager.beginUpdateVC(id_f);
    setBoundaryCondition();
    blockManager.endUpdateVC(id_f);
  }
}


void Solver::run()
{
  REAL_TYPE* wData = work->getData();

  for (int iOuter = 0; iOuter < nLoopOuter; iOuter++) {

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
//    BlockBase* block = blockManager.getBlock(id);
      Block* block = dynamic_cast<Block*>(blockManager.getBlock(id));

      Scalar3D<REAL_TYPE>* f = dynamic_cast<Scalar3D<REAL_TYPE>*>(
                                             block->getDataClass(id_f));
      REAL_TYPE* fData = f->getData();
      Index3DS fIndex = f->getIndex();

      const BoundaryInfo* boundaryInfo = block->getBoundaryInfo();

      const Vec3r& cellSize = block->getCellSize();
      REAL_TYPE cx = 1.0 / (cellSize.x * cellSize.x);
      REAL_TYPE cy = 1.0 / (cellSize.y * cellSize.y);
      REAL_TYPE cz = 1.0 / (cellSize.z * cellSize.z);

      REAL_TYPE c = 1.0 / (2.0 * (cx + cy + cz));

      const Vec3r& orig = block->getOrigin();
      setSource(nx, ny, nz, wData, fIndex, orig, cellSize);


      TimingStart(SOR);
      for (int iInner = 0; iInner < nLoopInner; iInner++) {
        calcSorInBlock(nx, ny, nz, fData, wData, fIndex, omega, c, cx, cy, cz);
        setBoundaryConditionInBlock(nx, ny, nz, fData, fIndex, boundaryInfo);
      }
      TimingStop(SOR);
    }

    TimingStart(VCUPDATE);
    if (separateVCUpdate) {
      blockManager.beginUpdateVC_X(id_f);
      blockManager.beginUpdateVC_Y(id_f);
      blockManager.beginUpdateVC_Z(id_f);
      TimingStart(BC);
      setBoundaryCondition();
      TimingStop(BC);
      blockManager.endUpdateVC_X(id_f);
      blockManager.endUpdateVC_Y(id_f);
      blockManager.endUpdateVC_Z(id_f);
    } else {
      blockManager.beginUpdateVC(id_f);
      TimingStart(BC);
      setBoundaryCondition();
      TimingStop(BC);
      blockManager.endUpdateVC(id_f);
    }
    TimingStop(VCUPDATE);
  }
}


void Solver::calcSorInBlock(int nx, int ny, int nz,
                            REAL_TYPE* fData, const REAL_TYPE* sData, Index3DS fIndex,
                            REAL_TYPE omega,
                            REAL_TYPE c0, REAL_TYPE cx, REAL_TYPE cy, REAL_TYPE cz)
{
  for (int c = 0; c < 2; c++) {

#pragma loop noalias
#pragma loop novrec
#pragma omp parallel for schedule(static,1)
    for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
        for (int i = (c+j+k) % 2; i < nx; i += 2) {
           REAL_TYPE fNew = c0 * (
                 + cx * (fData[fIndex(i-1,j,k)] + fData[fIndex(i+1,j,k)])
                 + cy * (fData[fIndex(i,j-1,k)] + fData[fIndex(i,j+1,k)])
                 + cz * (fData[fIndex(i,j,k-1)] + fData[fIndex(i,j,k+1)])
                 - sData[fIndex(i,j,k)] );
           fData[fIndex(i,j,k)] = omega * fNew + (1.0 - omega) * fData[fIndex(i,j,k)];
        }
      }
    }

  }
}



void Solver::setInitialCondition()
{
  for (int id = 0; id < blockManager.getNumBlock(); ++id) {
    BlockBase* block = blockManager.getBlock(id);
    Scalar3D<REAL_TYPE>* f = dynamic_cast<Scalar3D<REAL_TYPE>*>(
                              block->getDataClass(id_f));
    REAL_TYPE* fData = f->getData();
    Index3DS fIndex = f->getIndex();

    setInitialConditionInBlock(nx, ny, nz, fData, fIndex, 0.0);
  }
}

void Solver::setInitialConditionInBlock(int nx, int ny, int nz,
                                        REAL_TYPE* fData, Index3DS fIndex,
                                        REAL_TYPE value)
{
#pragma loop noalias
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        fData[fIndex(i,j,k)] = value;
      }
    }
  }
}


void Solver::setBoundaryCondition()
{
  for (int id = 0; id < blockManager.getNumBlock(); ++id) {
    Block* block = dynamic_cast<Block*>(blockManager.getBlock(id));
    const BoundaryInfo* boundaryInfo = block->getBoundaryInfo();
    Scalar3D<REAL_TYPE>* f = dynamic_cast<Scalar3D<REAL_TYPE>*>(
                              block->getDataClass(id_f));
    REAL_TYPE* fData = f->getData();
    Index3DS fIndex = f->getIndex();

    setBoundaryConditionInBlock(nx, ny, nz, fData, fIndex, boundaryInfo);
  }
}


void Solver::setBoundaryConditionInBlock(int nx, int ny, int nz,
                                         REAL_TYPE* fData, Index3DS fIndex,
                                         const BoundaryInfo* boundaryInfo)
{
  for (int i = 0; i < NUM_FACE; ++i) {
    Face face = Face(i);
    switch (boundaryInfo[face].getType()) {
      case BoundaryInfo::DIRICHLET:
        {
        REAL_TYPE v = boundaryValue[boundaryInfo[face].getID()];
        setDirichletBoundaryInBlock(nx, ny, nz, fData, fIndex, face, v);
        }
        break;
      case BoundaryInfo::NEUMANN:
      //{
      //REAL_TYPE v = boundaryValue[boundaryInfo[face].getID()];
      //setNeumannBoundaryInBlock(nx, ny, nz, fData, fIndex, face, v);
      //}
      //break;
        std::cout << "*** Neumann boundary condition is not implemented yet" << std::endl;
        Exit(EX_FAILURE);
      case BoundaryInfo::INNER:
      case BoundaryInfo::PERIODIC:
      default:
        // do nothing
        break;
    }
  }
}


void Solver::setDirichletBoundaryInBlock(int nx, int ny, int nz,
                                         REAL_TYPE* fData, Index3DS fIndex,
                                         Face face, REAL_TYPE value)
{
  switch (face) {
    case X_M:
#pragma loop noalias
#pragma loop novrec
      for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
          fData[fIndex(-1,j,k)] = 2.0 * value - fData[fIndex(0,j,k)];
        }
      }
      break;
    case X_P:
#pragma loop noalias
#pragma loop novrec
      for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
          fData[fIndex(nx,j,k)] = 2.0 * value - fData[fIndex(nx-1,j,k)];
        }
      }
      break;
    case Y_M:
#pragma loop noalias
#pragma loop novrec
      for (int k = 0; k < nz; ++k) {
        for (int i = 0; i < nx; ++i) {
          fData[fIndex(i,-1,k)] = 2.0 * value - fData[fIndex(i,0,k)];
        }
      }
      break;
    case Y_P:
#pragma loop noalias
#pragma loop novrec
      for (int k = 0; k < nz; ++k) {
        for (int i = 0; i < nx; ++i) {
          fData[fIndex(i,ny,k)] = 2.0 * value - fData[fIndex(i,ny-1,k)];
        }
      }
      break;
    case Z_M:
#pragma loop noalias
#pragma loop novrec
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          fData[fIndex(i,j,-1)] = 2.0 * value - fData[fIndex(i,j,0)];
        }
      }
      break;
    case Z_P:
#pragma loop noalias
#pragma loop novrec
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          fData[fIndex(i,j,nz)] = 2.0 * value - fData[fIndex(i,j,nz-1)];
        }
      }
      break;
    default:
      break;
  }
}


void Solver::checkResult(bool verbose)
{

  std::string outputFile(outputFileName);
  SiloWriter writer(outputFile, "mesh", false);
  writer.writeDomain("block_mesh", "domain");

  // 誤差を格納するためのデータクラス
  int id_e = blockManager.setDataClass<Scalar3D<REAL_TYPE> >(vc);

  if (comm.Get_rank() == 0) std::cout << std::endl;
  comm.Barrier();

  REAL_TYPE errorMax = 0.0;
  for (int id = 0; id < blockManager.getNumBlock(); ++id) {
    BlockBase* block = blockManager.getBlock(id);
    Vec3r origin = block->getOrigin();
    Vec3r delta = block->getCellSize();

    REAL_TYPE errorMax_inBlock = 0.0;
    Scalar3D<REAL_TYPE>* f = dynamic_cast<Scalar3D<REAL_TYPE>*>(
                            block->getDataClass(id_f));
    REAL_TYPE* fData = f->getData();
    Index3DS fIndex = f->getIndex();

    Scalar3D<REAL_TYPE>* e = dynamic_cast<Scalar3D<REAL_TYPE>*>(
                            block->getDataClass(id_e));
    assert(e);
    REAL_TYPE* eData = e->getData();
    assert(eData);


    REAL_TYPE a = 2 * M_PI;
//  REAL_TYPE a = M_PI;

#pragma loop noalias
    for (int k = 0; k < nz; ++k) {
      REAL_TYPE z = origin[2] + (k + 0.5) * delta[2];
      REAL_TYPE sz = sin(a * z);
      for (int j = 0; j < ny; ++j) {
        REAL_TYPE y = origin[1] + (j + 0.5) * delta[1];
        REAL_TYPE sy = sin(a * y);
        for (int i = 0; i < nx; ++i) {
          REAL_TYPE x = origin[0] + (i + 0.5) * delta[0];
          REAL_TYPE sx = sin(a * x);
//        std::cout << Vec3i(i,j,k) << ": " << f(i,j,k) << std::endl;
          REAL_TYPE f0 = sx * sy * sz;
          REAL_TYPE err = fabs(fData[fIndex(i,j,k)] - f0);
          eData[fIndex(i,j,k)] = err;
          errorMax_inBlock = std::max(err, errorMax_inBlock);
        }
      }
    }

    if (verbose) {
      std::cout << "Block " << comm.Get_rank() << "-" << id
                << " : errorMax = " << errorMax_inBlock
                << " (level=" <<  block->getLevel() << ", orig="
                << block->getOrigin() << ")" << std::endl;
    }

    errorMax = std::max(errorMax_inBlock, errorMax);
  }

  if (comm.Get_rank() == 0) {
    if (sizeof(REAL_TYPE) == 8) {
      comm.Reduce(MPI::IN_PLACE, &errorMax, 1, MPI::DOUBLE, MPI::MAX, 0);
    }
    else {
      comm.Reduce(MPI::IN_PLACE, &errorMax, 1, MPI::FLOAT, MPI::MAX, 0);
    }
  } else {
    if (sizeof(REAL_TYPE) == 8) {
      comm.Reduce(&errorMax, &errorMax, 1, MPI::DOUBLE, MPI::MAX, 0);
    }
    else {
      comm.Reduce(&errorMax, &errorMax, 1, MPI::FLOAT, MPI::MAX, 0);
    }
  }

  if (comm.Get_rank() == 0) {
    std::cout << "errorMax = " << errorMax << std::endl;
  }

  writer.writeScalar<REAL_TYPE>(id_f, "result");
  writer.writeScalar<REAL_TYPE>(id_e, "error");

}


void Solver::setSource(int nx, int ny, int nz,
                       REAL_TYPE* sData, Index3DS sIndex,
                       const Vec3r& orig, const Vec3r& delta)
{
  REAL_TYPE a = 2 * M_PI;
//REAL_TYPE a = M_PI;
  REAL_TYPE c = - 3.0 * a * a;

#pragma loop noalias
  for (int k = 0; k < nz; ++k) {
    REAL_TYPE z = orig[2] + (k + 0.5) * delta[2];
    REAL_TYPE sz = sin(a * z);
    for (int j = 0; j < ny; ++j) {
      REAL_TYPE y = orig[1] + (j + 0.5) * delta[1];
      REAL_TYPE sy = sin(a * y);
      for (int i = 0; i < nx; ++i) {
        REAL_TYPE x = orig[0] + (i + 0.5) * delta[0];
        REAL_TYPE sx = sin(a * x);
        sData[sIndex(i,j,k)] = c * sx * sy * sz;
      }
    }
  }
}


void Solver::dumpDataClass(int dataClassID)
{
  for (int id = 0; id < blockManager.getNumBlock(); ++id) {
    BlockBase* block = blockManager.getBlock(id);
    Scalar3D<REAL_TYPE>* f = dynamic_cast<Scalar3D<REAL_TYPE>*>(
                              block->getDataClass(id_f));
    REAL_TYPE* fData = f->getData();
    Index3DS fIndex = f->getIndex();

    std::cout << "Block " << comm.Get_rank() << "-" << id << std::endl;

    for (int k = -vc; k < nz+vc; ++k) {
      for (int j = -vc; j < ny+vc; ++j) {
        for (int i = -vc; i < nx+vc; ++i) {
          std::cout << Vec3i(i,j,k) << ": " << fData[fIndex(i,j,k)] << std::endl;
        }
      }
    }
  }
}
