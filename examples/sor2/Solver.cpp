#include "Solver.h"
#include "Scalar3DUpdater.h"
#include "Timing.h"

#include "SiloWriter.h"


Solver::Solver(const Config& conf, const std::vector<double>& boundaryValue)
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
  id_f = blockManager.setDataClass<Scalar3D<double>, Scalar3DUpdater<double> >(vc);

  // データクラス変数<f>の仮想セル同期準備
  blockManager.prepareForVCUpdate(id_f, tag_f, separateVCUpdate);

  // 作業用データクラス変数<work>
  work = new Scalar3D<double>(size, vc);   // ソース項に使用
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
  double* wData = work->getData();

  for (int iOuter = 0; iOuter < nLoopOuter; iOuter++) {

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
//    BlockBase* block = blockManager.getBlock(id);
      Block* block = dynamic_cast<Block*>(blockManager.getBlock(id));

      Scalar3D<double>* f = dynamic_cast<Scalar3D<double>*>(
                                             block->getDataClass(id_f));
      double* fData = f->getData();
      Index3DS fIndex = f->getIndex();

      const BoundaryInfo* boundaryInfo = block->getBoundaryInfo();

      const Vec3r& cellSize = block->getCellSize();
      double cx = 1.0 / (cellSize.x * cellSize.x);
      double cy = 1.0 / (cellSize.y * cellSize.y);
      double cz = 1.0 / (cellSize.z * cellSize.z);

      double c = 1.0 / (2.0 * (cx + cy + cz));

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
                            double* fData, const double* sData, Index3DS fIndex,
                            double omega,
                            double c0, double cx, double cy, double cz)
{
  for (int c = 0; c < 2; c++) {

#pragma loop noalias
#pragma loop novrec
#pragma omp parallel for schedule(static,1)
    for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
        for (int i = (c+j+k) % 2; i < nx; i += 2) {
           double fNew = c0 * (
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
    Scalar3D<double>* f = dynamic_cast<Scalar3D<double>*>(
                              block->getDataClass(id_f));
    double* fData = f->getData();
    Index3DS fIndex = f->getIndex();

    setInitialConditionInBlock(nx, ny, nz, fData, fIndex, 0.0);
  }
}

void Solver::setInitialConditionInBlock(int nx, int ny, int nz,
                                        double* fData, Index3DS fIndex,
                                        double value)
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
    Scalar3D<double>* f = dynamic_cast<Scalar3D<double>*>(
                              block->getDataClass(id_f));
    double* fData = f->getData();
    Index3DS fIndex = f->getIndex();

    setBoundaryConditionInBlock(nx, ny, nz, fData, fIndex, boundaryInfo);
  } 
}


void Solver::setBoundaryConditionInBlock(int nx, int ny, int nz,
                                         double* fData, Index3DS fIndex,
                                         const BoundaryInfo* boundaryInfo)
{
  for (int i = 0; i < NUM_FACE; ++i) {
    Face face = Face(i);
    switch (boundaryInfo[face].getType()) {
      case BoundaryInfo::DIRICHLET:
        {
        double v = boundaryValue[boundaryInfo[face].getID()];
        setDirichletBoundaryInBlock(nx, ny, nz, fData, fIndex, face, v);
        }
        break;
      case BoundaryInfo::NEUMANN:
      //{
      //double v = boundaryValue[boundaryInfo[face].getID()];
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
                                         double* fData, Index3DS fIndex,
                                         Face face, double value)
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
  int id_e = blockManager.setDataClass<Scalar3D<double> >(vc);
  
  if (comm.Get_rank() == 0) std::cout << std::endl;
  comm.Barrier();

  double errorMax = 0.0;
  for (int id = 0; id < blockManager.getNumBlock(); ++id) {
    BlockBase* block = blockManager.getBlock(id);
    Vec3r origin = block->getOrigin();
    Vec3r delta = block->getCellSize();

    double errorMax_inBlock = 0.0;
    Scalar3D<double>* f = dynamic_cast<Scalar3D<double>*>(
                            block->getDataClass(id_f));
    double* fData = f->getData();
    Index3DS fIndex = f->getIndex();

    Scalar3D<double>* e = dynamic_cast<Scalar3D<double>*>(
                            block->getDataClass(id_e));
    assert(e);
    double* eData = e->getData();
    assert(eData);


    double a = 2 * M_PI;
//  double a = M_PI;

#pragma loop noalias
    for (int k = 0; k < nz; ++k) {
      double z = origin[2] + (k + 0.5) * delta[2];
      double sz = sin(a * z);
      for (int j = 0; j < ny; ++j) {
        double y = origin[1] + (j + 0.5) * delta[1];
        double sy = sin(a * y);
        for (int i = 0; i < nx; ++i) {
          double x = origin[0] + (i + 0.5) * delta[0];
          double sx = sin(a * x);
//        std::cout << Vec3i(i,j,k) << ": " << f(i,j,k) << std::endl;
          double f0 = sx * sy * sz;
          double err = fabs(fData[fIndex(i,j,k)] - f0);
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
    comm.Reduce(MPI::IN_PLACE, &errorMax, 1, MPI::DOUBLE, MPI::MAX, 0);
  } else {
    comm.Reduce(&errorMax, &errorMax, 1, MPI::DOUBLE, MPI::MAX, 0);
  }

  if (comm.Get_rank() == 0) {
    std::cout << "errorMax = " << errorMax << std::endl;
  }

  writer.writeScalar<double>(id_f, "result");
  writer.writeScalar<double>(id_e, "error");

}


void Solver::setSource(int nx, int ny, int nz,
                       double* sData, Index3DS sIndex,
                       const Vec3r& orig, const Vec3r& delta)
{
  double a = 2 * M_PI;
//double a = M_PI;
  double c = - 3.0 * a * a;

#pragma loop noalias
  for (int k = 0; k < nz; ++k) {
    double z = orig[2] + (k + 0.5) * delta[2];
    double sz = sin(a * z);
    for (int j = 0; j < ny; ++j) {
      double y = orig[1] + (j + 0.5) * delta[1];
      double sy = sin(a * y);
      for (int i = 0; i < nx; ++i) {
        double x = orig[0] + (i + 0.5) * delta[0];
        double sx = sin(a * x);
        sData[sIndex(i,j,k)] = c * sx * sy * sz;
      }
    }
  }
}


void Solver::dumpDataClass(int dataClassID)
{
  for (int id = 0; id < blockManager.getNumBlock(); ++id) {
    BlockBase* block = blockManager.getBlock(id);
    Scalar3D<double>* f = dynamic_cast<Scalar3D<double>*>(
                              block->getDataClass(id_f));
    double* fData = f->getData();
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
