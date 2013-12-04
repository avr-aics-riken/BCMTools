#ifndef BLOCK_FACTORY_H
#define BLOCK_FACTORY_H

#include "Block.h"
#include "Octree.h"
#include "Partition.h"
#include "Config.h"

class Node;

class BlockFactory {

public:

  static Vec3r getOrigin(const Node* node);

  static Vec3r getBlockSize(const Node* node);

  static BoundaryInfo* makeBoundaryInfo(const Node* node, const Octree* tree,
                                        const Config* conf);

  static NeighborInfo* makeNeighborInfo(const Node* node, const Octree* tree,
                                        const BoundaryInfo* boundaryInfo,
                                        const Partition* partition);

private:

  static bool onOuterBoundary(const Pedigree& pedigree, Face face);

};

#endif // BLOCK_FACTORY_H
