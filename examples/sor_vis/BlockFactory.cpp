#include "BlockFactory.h"


Vec3r BlockFactory::getOrigin(const Node* node)
{
  int level = node->pedigree.level;
  int max = 1 << level;
  return Vec3r((double)node->pedigree.x/max,
               (double)node->pedigree.y/max,
               (double)node->pedigree.z/max);
}


Vec3r BlockFactory::getBlockSize(const Node* node)
{
  int level = node->pedigree.level;
  int max = 1 << level;
  return Vec3r(1.0/max, 1.0/max, 1.0/max);
}


BoundaryInfo* BlockFactory::makeBoundaryInfo(const Node* node, const Octree* tree,
                                             const Config* conf)
{
  BoundaryInfo* boundaryInfo = new BoundaryInfo[NUM_FACE];

  // BoundaryInfoの初期値は，type=INNER, id=-1

  // 境界条件の設定
  for (int i = 0; i < NUM_FACE; i++) {
    Face face = Face(i);
    if (onOuterBoundary(node->pedigree, face))
//    boundaryInfo[face].setType(BoundaryInfo::PERIODIC);
      boundaryInfo[face].setType(BoundaryInfo::DIRICHLET);
      boundaryInfo[face].setID(0);
  }

#if 0
  switch (conf->type) {
    case 'x':
      if (onOuterBoundary(node->pedigree, X_M)) {
        boundaryInfo[X_M].setType(BoundaryInfo::DIRICHLET);
        boundaryInfo[X_M].setID(0);
      }
      if (onOuterBoundary(node->pedigree, X_P)) {
        boundaryInfo[X_P].setType(BoundaryInfo::DIRICHLET);
        boundaryInfo[X_P].setID(1);
      }
      break;
    case 'y':
      if (onOuterBoundary(node->pedigree, Y_M)) {
        boundaryInfo[Y_M].setType(BoundaryInfo::DIRICHLET);
        boundaryInfo[Y_M].setID(0);
      }
      if (onOuterBoundary(node->pedigree, Y_P)) {
        boundaryInfo[Y_P].setType(BoundaryInfo::DIRICHLET);
        boundaryInfo[Y_P].setID(1);
      }
      break;
    case 'z':
      if (onOuterBoundary(node->pedigree, Z_M)) {
        boundaryInfo[Z_M].setType(BoundaryInfo::DIRICHLET);
        boundaryInfo[Z_M].setID(0);
      }
      if (onOuterBoundary(node->pedigree, Z_P)) {
        boundaryInfo[Z_P].setType(BoundaryInfo::DIRICHLET);
        boundaryInfo[Z_P].setID(1);
      }
      break;
    default:
      assert(0);
      break;
  }
#endif

  return boundaryInfo;
}


NeighborInfo* BlockFactory::makeNeighborInfo(const Node* node, const Octree* tree,
                                             const BoundaryInfo* boundaryInfo,
                                             const Partition* partition)
{
  int level = node->pedigree.level;

  NeighborInfo* neighborInfo = new NeighborInfo[NUM_FACE];

  // 隣接情報
  for (int i = 0; i < NUM_FACE; ++i) {
    Face face = Face(i);
    BoundaryInfo::Type type = boundaryInfo[face].getType();
    if (type != BoundaryInfo::INNER) neighborInfo[face].setOuterBoundary();

    Node* neighbor;
    if (type == BoundaryInfo::INNER) {
      neighbor = tree->searchNode(node, face);
    } else if (type == BoundaryInfo::PERIODIC) {
      neighbor = tree->searchNode(node, face, true);
    } else {
      continue;
    }
    assert(neighbor);
    int neighborLevel = neighbor->pedigree.level;
    int levelDiff = neighborLevel - level;

    if (levelDiff == -1) {
      neighborInfo[face].setLevelDifference(-1);
      neighborInfo[face].setID(neighbor->id);
      neighborInfo[face].setRank(partition->getRank(neighbor->id));
      Subface subface = NeighborInfo::childIdToSubface(face,
                                                       node->pedigree.getChildId(level));
      neighborInfo[face].setNeighborSubface(subface);
    }
    else if (levelDiff == 0) {
      if (neighbor->isLeafNode()) {
        neighborInfo[face].setLevelDifference(0);
        neighborInfo[face].setID(neighbor->id);
        neighborInfo[face].setRank(partition->getRank(neighbor->id));
      }
      else {
        neighborInfo[face].setLevelDifference(+1);
        for (int j = 0; j < NUM_SUBFACE; j++) {
          Subface subface = Subface(j);
          Node* child = neighbor->child[NeighborInfo::getNeighborChild(face, subface)];
          assert(child->isLeafNode());
          neighborInfo[face].setID(subface, child->id);
          neighborInfo[face].setRank(subface, partition->getRank(child->id));
        }
      }
    }

  }
  return neighborInfo;
}



bool BlockFactory::onOuterBoundary(const Pedigree& pedigree, Face face)
{
  switch (face) {
    case X_M:
      if (pedigree.x == 0) return true;
      else return false;
    case X_P:
      if (pedigree.x + 1 == 1 << pedigree.level) return true;
      else return false;
    case Y_M:
      if (pedigree.y == 0) return true;
      else return false;
    case Y_P:
      if (pedigree.y + 1 == 1 << pedigree.level) return true;
      else return false;
    case Z_M:
      if (pedigree.z == 0) return true;
      else return false;
    case Z_P:
      if (pedigree.z + 1 == 1 << pedigree.level) return true;
      else return false;
    default:
      assert(0);
      return false;
  }
}
