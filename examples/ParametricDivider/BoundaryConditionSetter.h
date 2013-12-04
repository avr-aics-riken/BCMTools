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
        boundaryInfo[face].setType(BoundaryInfo::DIRICHLET);
        boundaryInfo[face].setID(0);
      }
    }

    return boundaryInfo;
  }
  
};


#endif // BOUNDARY_CONDITION_SETTER_H
