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

#ifndef BOUNDARY_CONDITION_SETTER_H
#define BOUNDARY_CONDITION_SETTER_H

#include "BCMTools.h"
#include "BoundaryConditionSetterBase.h"
#include "BCMOctree.h"
#include "BoundaryInfo.h"
#include "Config.h"

class BoundaryConditionSetter : public BoundaryConditionSetterBase {

  const Config* conf;

public:

  BoundaryConditionSetter(const Config* conf) : conf(conf) {}

  ~BoundaryConditionSetter() {}

  BoundaryInfo* operator() (const Node* node, const BCMOctree* tree) {
    BoundaryInfo* boundaryInfo = new BoundaryInfo[NUM_FACE];

    // BoundaryInfoの初期値は，type=INNER, id=-1

    // 境界条件の設定
    for (int i = 0; i < NUM_FACE; i++) {
      Face face = Face(i);
      if (tree->checkOnOuterBoundary(node, face)) {
        boundaryInfo[face].setType(BoundaryInfo::PERIODIC);
      }
    }

    switch (conf->type) {
      case 'x':
        if (tree->checkOnOuterBoundary(node, X_M)) {
          boundaryInfo[X_M].setType(BoundaryInfo::DIRICHLET);
          boundaryInfo[X_M].setID(0);
        }
        if (tree->checkOnOuterBoundary(node, X_P)) {
          boundaryInfo[X_P].setType(BoundaryInfo::DIRICHLET);
          boundaryInfo[X_P].setID(1);
        }
        break;
      case 'y':
        if (tree->checkOnOuterBoundary(node, Y_M)) {
          boundaryInfo[Y_M].setType(BoundaryInfo::DIRICHLET);
          boundaryInfo[Y_M].setID(0);
        }
        if (tree->checkOnOuterBoundary(node, Y_P)) {
          boundaryInfo[Y_P].setType(BoundaryInfo::DIRICHLET);
          boundaryInfo[Y_P].setID(1);
        }
        break;
      case 'z':
        if (tree->checkOnOuterBoundary(node, Z_M)) {
          boundaryInfo[Z_M].setType(BoundaryInfo::DIRICHLET);
          boundaryInfo[Z_M].setID(0);
        }
        if (tree->checkOnOuterBoundary(node, Z_P)) {
          boundaryInfo[Z_P].setType(BoundaryInfo::DIRICHLET);
          boundaryInfo[Z_P].setID(1);
        }
        break;
      default:
        assert(0);
        break;
    }

    return boundaryInfo;
  }

};


#endif // BOUNDARY_CONDITION_SETTER_H
