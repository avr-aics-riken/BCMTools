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

#include <cmath>
#include "Solver.h"
#include "Scalar3D.h"
#include "Vector3D.h"

#include "SiloWriter.h"

Solver::Solver(const Config& conf)
  : blockManager(BlockManager::getInstance()),
    comm(blockManager.getCommunicator()),
//  vc(conf.vc),
    outputFile(conf.output)
{
  size = blockManager.getSize();
  nx = size[0];
  ny = size[1];
  nz = size[2];

  // データクラス変数<f>の生成・登録
  id_s = blockManager.setDataClass<Scalar3D<double> >(vc);
  id_v = blockManager.setDataClass<Vector3D<double> >(vc);
}


Solver::~Solver()
{
}


void Solver::initialize()
{
  for (int id = 0; id < blockManager.getNumBlock(); ++id) {
    BlockBase* block = blockManager.getBlock(id);

    Scalar3D<double>* s
        = dynamic_cast<Scalar3D<double>*>(block->getDataClass(id_s));
    double* sData = s->getData();
    Index3DS sIndex = s->getIndex();

    Vector3D<double>* v
        = dynamic_cast<Vector3D<double>*>(block->getDataClass(id_v));
    double* vData = v->getData();
    Index3DV vIndex = v->getIndex();

    const Vec3d& orig = block->getOrigin();
    const Vec3d& pitch = block->getCellSize();

    for (int k = -vc; k < nz+vc; k++) {
      for (int j = -vc; j < ny+vc; j++) {
        for (int i = -vc; i < nx+vc; i++) {
          double x = orig.x + (i + 0.5) * pitch.x;
          double y = orig.y + (j + 0.5) * pitch.y;
          double z = orig.z + (k + 0.5) * pitch.z;
          sData[sIndex(i,j,k)] = sqrt(x*x + y*y + z*z);
          vData[vIndex(i,j,k)+0] = x;
          vData[vIndex(i,j,k)+1] = y;
          vData[vIndex(i,j,k)+2] = z;
        }
      }
    }

  }
}


void Solver::run()
{
  SiloWriter writer(outputFile, "mesh");

  writer.writeDomain("block_mesh", "domain");

  writer.writeScalar<double>(id_s, "scalar");

  writer.writeVector<double>(id_v, "vector");
}
