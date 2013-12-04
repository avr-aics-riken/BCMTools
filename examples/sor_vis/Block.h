#ifndef BLOCK_H
#define BLOCK_H

#include "BlockBase.h"
#include "BoundaryInfo.h"


class Block : public BlockBase {

  BoundaryInfo* boundaryInfo;

public:

  Block(const Vec3i& size, const Vec3r& origin, const Vec3r& blockSize,
        int level, NeighborInfo* neighborInfo, BoundaryInfo* boundaryInfo) 
   : BlockBase(size, origin, blockSize, level, neighborInfo),
     boundaryInfo(boundaryInfo) {}

  virtual ~Block() {
    delete[] boundaryInfo;
  }

  const BoundaryInfo* getBoundaryInfo() const { return boundaryInfo; }

};

#endif // BLOCK_H
