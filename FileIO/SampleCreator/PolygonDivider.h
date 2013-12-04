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

#ifndef POLYGON_DIVIDER_H
#define POLYGON_DIVIDER_H

#include <string>
#include <vector>
#include <algorithm> // for max, min
#include "BCMPolylib.h"

#include "BCMTools.h"
#include "MultiRootDivider.h"
#include "BoundingBox.h"

class PolygonDivider : public MultiRootDivider {

  int minLevel;

  const PolylibNS::BCMPolylib* pl;

  const std::vector<PolygonGroupSpec>& polygonGroupList;

  const std::vector<BoundingBoxSpec>& boundingBoxList; 

  double marginRatio;

public:

  PolygonDivider(const Vec3r& origin, double rootLength, const RootGrid* rootGrid, 
                 int minLevel, const PolylibNS::BCMPolylib* pl,
                 const std::vector<PolygonGroupSpec>& polygonGroupList,
                 const std::vector<BoundingBoxSpec>& boundingBoxList, 
                 double marginRatio = 0.0)
    : MultiRootDivider(origin, rootLength, rootGrid),
      minLevel(minLevel), pl(pl),
      polygonGroupList(polygonGroupList), boundingBoxList(boundingBoxList),
      marginRatio(marginRatio) {
  }

  ~PolygonDivider() {}

  NodeType operator() (const Pedigree& pedigree) {

    int level = pedigree.getLevel();

    if  (level < minLevel) return BRANCH;

    NodeType ret = LEAF_ACTIVE;

    for (std::vector<BoundingBoxSpec>::const_iterator it = boundingBoxList.begin();
         it != boundingBoxList.end(); ++it) {
      int maxLevel = it->level;
      if (level < maxLevel) {
        BoundingBox box = it->boundingBox;
        BoundingBox region = defineSearchRegion(pedigree, maxLevel);
        region.setMargin(marginRatio / (1 << maxLevel));
        if (box.intersects(region)) ret = BRANCH;
      }
    }


    for (std::vector<PolygonGroupSpec>::const_iterator it = polygonGroupList.begin();
         it != polygonGroupList.end(); ++it) {
      int maxLevel = it->level;
      if (level < maxLevel) {
        const std::string& polygonGroup = it->polygonGroupName;
        BoundingBox region = defineSearchRegion(pedigree, maxLevel);
        region.setMargin(marginRatio / (1 << maxLevel));
        PolylibNS::Vec3f min(region.getMin().x, region.getMin().y, region.getMin().z);
        PolylibNS::Vec3f max(region.getMax().x, region.getMax().y, region.getMax().z);
        std::vector<PolylibNS::Triangle*>* polygonList 
                        = pl->search_polygons(polygonGroup, min, max, false);
        int nPolygon = polygonList->size();
        delete polygonList;
        if (nPolygon > 0) ret = BRANCH;
      }
    }

    return ret;
  }

};

#endif // POLYGON_DIVIDER_H
