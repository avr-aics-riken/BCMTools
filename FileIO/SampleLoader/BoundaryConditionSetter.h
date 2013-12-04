/*
 * BCMTools
 *
 * Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef BOUNDARY_CONDITION_SETTER_H
#define BOUNDARY_CONDITION_SETTER_H

#include "BCMTools.h"
#include "BoundaryConditionSetterBase.h"
#include "BCMOctree.h"
#include "BoundaryInfo.h"

class BoundaryConditionSetter : public BoundaryConditionSetterBase {


public:

  BoundaryConditionSetter() {}

  ~BoundaryConditionSetter() {}

  BoundaryInfo* operator() (const Node* node, const BCMOctree* tree) {
    BoundaryInfo* boundaryInfo = new BoundaryInfo[NUM_FACE];

    // BoundaryInfoの初期値は，type=INNER, id=-1

    // 境界条件の設定
    for (int i = 0; i < NUM_FACE; i++) {
      Face face = Face(i);
      if (tree->checkOnOuterBoundary(node, face)) {
        boundaryInfo[face].setType(BoundaryInfo::DIRICHLET);
        boundaryInfo[face].setID(0);
      }
    }

    return boundaryInfo;
  }
  
};


#endif // BOUNDARY_CONDITION_SETTER_H
