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

#ifndef BLOCK_BOUNDING_BOX_H
#define BLOCK_BOUNDING_BOX_H

#include "BoundingBox.h"
#include "BlockBase.h"
#include "BCMOctree.h"

/// ブロックからポリゴン分配用のバウンディングボックスを計算するクラス.
///
/// 外部境界条件の種類によってマージンに変化させたい場合は，
/// このクラスをカスタマイズすること.
///
class BlockBoundingBox {

  const BCMOctree* tree;
  const Vec3r& rootOrigin;
  double rootLength;
  const Vec3i& size;
  int vc;

public:

  BlockBoundingBox(const BCMOctree* tree, 
                   const Vec3r& rootOrigin, double rootLength, const Vec3i& size,
                   int vc)
    : tree(tree), rootOrigin(rootOrigin), rootLength(rootLength), size(size), vc(vc)
  {}

  BoundingBox getBoundingBox(const Node* node) const {
    Vec3r origin = tree->getOrigin(node) * rootLength + rootOrigin;
    Vec3r blockSize = node->getBlockSize() * rootLength;
    Vec3r margin((blockSize.x/size.x)*vc,
                 (blockSize.y/size.y)*vc,
                 (blockSize.z/size.z)*vc);
    return BoundingBox(origin-margin, origin+blockSize+margin);
  }

};


#endif // BLOCK_BOUNDING_BOX_H
