/*
 * BCMTools
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <algorithm>  // for random_shuffle
#include "BCMOctree.h"


/// コンストラクタ.
BCMOctree::BCMOctree(RootGrid* rootGrid, Divider* divider, Ordering ordering)
  : rootGrid(rootGrid), divider(divider), ordering(ordering)
{
  int nRoot = rootGrid->getSize();
  rootNodes = new Node*[nRoot];

  for (int id = 0; id < nRoot; id++) {
    Node* root = new Node(id);
    makeNode(root);
    if (ordering == HILBERT) {
      pickupLeafNodeHilbertOrdering(root, 0);
    } else {
      pickupLeafNodeZOrdering(root);
    }
    rootNodes[id] = root;
  }

  if (ordering == RANDOM) randomShuffle();
}


/// コンストラクタ(リーフノードのPedigreeリストから).
BCMOctree::BCMOctree(RootGrid* rootGrid, Ordering ordering,
                     int numLeafNode, const unsigned char* buf)
  : rootGrid(rootGrid), divider(0), ordering(ordering)
{
  int nRoot = rootGrid->getSize();
  rootNodes = new Node*[nRoot];

  buildTreeFromPedigreeList(numLeafNode, buf);

  for (int id = 0; id < nRoot; id++) {
    if (ordering == HILBERT) {
      pickupLeafNodeHilbertOrdering(rootNodes[id], 0);
    } else {
      pickupLeafNodeZOrdering(rootNodes[id]);
    }
  }

  if (ordering == RANDOM) randomShuffle();
}

/// コンストラクタ(ファイルロード用)
BCMOctree::BCMOctree(RootGrid* rootGrid, const std::vector<Pedigree>& pedigrees)
 : rootGrid(rootGrid), divider(0), ordering(PEDIGREELIST)
{
	using namespace std;
	int nRoot = rootGrid->getSize();
	rootNodes = new Node*[nRoot];
	for(int i = 0; i < nRoot; i++){ rootNodes[i] = new Node(i); }
	
	leafNodeArray.clear();
	leafNodeArray.reserve(pedigrees.size());

	// Build Octree form Pedigree List
	for(vector<Pedigree>::const_iterator ped = pedigrees.begin(); ped != pedigrees.end(); ++ped) {
		unsigned int rootId = ped->getRootID();
		Node* node = rootNodes[rootId];
		for(int l = 1; l <= ped->getLevel(); l++) {
			if(node->isLeafNode()) {
				node->makeChildNodes();
				for(int cid = 0; cid < 8; cid++) {
					node->getChild(cid)->setActive(false);
				}
			}
			int cid = ped->getChildId(l);
			node = node->getChild(cid);
		}
		node->setActive(true);
		node->setBlockID(leafNodeArray.size());
		leafNodeArray.push_back(node);
	}

}

/// デストラクタ.
BCMOctree::~BCMOctree()
{
  //for (int id = 0; id < rootGrid->getSize(); id++) deleteNode(rootNodes[id]);
  for (int id = 0; id < rootGrid->getSize(); id++) delete rootNodes[id];
  delete[] rootNodes;

  delete divider;
  delete rootGrid;
}


/// 再帰呼び出しによりツリー生成.
void BCMOctree::makeNode(Node* node)
{
  Divider::NodeType nodeType = (*divider)(node->getPedigree());

  if (nodeType == Divider::LEAF_ACTIVE) return;

  if (nodeType == Divider::LEAF_NO_ACTIVE) {
    node->setActive(false);
    return;
  }

  node->makeChildNodes();
  for (int i = 0; i < 8; i++) makeNode(node->getChild(i));
}


/// Octree情報を他rankにブロードキャスト.
void BCMOctree::broadcast(MPI::Intracomm& comm)
{
  assert(comm.Get_rank() == 0);
  rootGrid->broadcast(comm);

  int numLeafNode = leafNodeArray.size();
  int ibuf[2];
  ibuf[0] = numLeafNode;
  ibuf[1] = ordering;
  comm.Bcast(&ibuf, 2, MPI::INT, 0);

  size_t size = Pedigree::GetSerializeSize();
  unsigned char* buf = new unsigned char[size * numLeafNode];

  size_t ip = 0;
  for (int id = 0; id < rootGrid->getSize(); id++) {
    packPedigrees(rootNodes[id], ip, buf);
  }

  comm.Bcast(buf, size*numLeafNode, MPI::BYTE, 0);
  delete[] buf;
}


/// Pedigree情報を通信用バッファにパック.
void BCMOctree::packPedigrees(Node* node, size_t& ip, unsigned char* buf)
{
  if (node->isLeafNode()) {
    if (node->isActive()) {
      (node->getPedigree()).serialize(&buf[ip]);
      ip += Pedigree::GetSerializeSize();
    }
    return;
  }
  for (int i = 0; i < 8; i++) packPedigrees(node->getChild(i), ip, buf);
}


#if 0
void BCMOctree::unpackPedigrees(int numLeafNode, const unsigned char* buf,
                                Pedigree* pedigrees)
{
  size_t size = Pedigree::GetSerializeSize();
  size_t p = 0;

  for (int i = 0; i < numLeafNode; i++) {
    pedigrees[i].deserialize(&buf[p]);
    p += size;
  }
}
#endif


/// rank0からOctree情報を受信.
BCMOctree* BCMOctree::ReceiveFromMaster(MPI::Intracomm& comm)
{
  assert(comm.Get_rank() != 0);

  RootGrid* rootGrid = RootGrid::ReceiveFromMaster(comm);

  int ibuf[2];
  comm.Bcast(&ibuf, 2, MPI::INT, 0);
  int numLeafNode = ibuf[0];
  Ordering ordering = Ordering(ibuf[1]);

  size_t size = Pedigree::GetSerializeSize();
  unsigned char* buf = new unsigned char[size * numLeafNode];

  comm.Bcast(buf, size*numLeafNode, MPI::BYTE, 0);

  BCMOctree* tree = new BCMOctree(rootGrid, ordering, numLeafNode, buf);
  delete[] buf;

  return tree;
}


/// リーフノードのPedigreeリストからOctreeを再構築.
void BCMOctree::buildTreeFromPedigreeList(int numLeafNode, const unsigned char* buf)
{
  for (int id = 0; id < rootGrid->getSize(); id++) {
    rootNodes[id] = 0;
  }

  Pedigree pedigree;
  size_t ip = 0;
  for (int i = 0; i < numLeafNode; i++) {
    pedigree.deserialize(&buf[ip]);
    ip += Pedigree::GetSerializeSize();
    unsigned id = pedigree.getRootID();
    if (rootNodes[id] == 0) rootNodes[id] = new Node(id);
    Node* parent = rootNodes[id];
    for (int l = 1; l <= pedigree.getLevel(); l++) {
      if (parent->isLeafNode()) {
        parent->makeChildNodes();
        for (int id = 0; id < 8; id++) parent->getChild(id)->setActive(false);
      }
      int id = pedigree.getChildId(l);
      parent = parent->getChild(id);
    }
    parent->setActive(true);
  }
}


/// Zオーダリング順にリーフノードリストを作成.
void BCMOctree::pickupLeafNodeZOrdering(Node* node)
{
  if (node->isLeafNode()) {
    if (node->isActive()) {
      node->setBlockID(leafNodeArray.size());
      leafNodeArray.push_back(node);
    }
    return;
  }
  for (int i = 0; i < 8; i++) pickupLeafNodeZOrdering(node->getChild(i));
}


/// ヒルベルトオーダリング順にリーフノードリストを作成.
void BCMOctree::pickupLeafNodeHilbertOrdering(Node* node, int orientation)
{
  if (node->isLeafNode()) {
    if (node->isActive()) {
      node->setBlockID(leafNodeArray.size());
      leafNodeArray.push_back(node);
    }
    return;
  }
  for (int i = 0; i < 8; i++) {
    int childId = HilbertOrdering[orientation][i];
    int childOrientation = HilbertOrientation[orientation][i];
    pickupLeafNodeHilbertOrdering(node->getChild(childId), childOrientation);
  }
}


/// 再帰呼び出しによるツリー消去.
void BCMOctree::deleteNode(Node* node)
{
  if (node->isLeafNode()) {
    delete node;
    return;
  }
  for (int i = 0; i < 8; i++) deleteNode(node->getChild(i));
}


/// 隣接ノード探索.
Node* BCMOctree::findNeighborNode(const Node* node, Face face) const
{
  const Pedigree& pedigree = node->getPedigree();
  int max0 = pedigree.getUpperBound() - 1; // 2^level - 1

  int level = pedigree.getLevel();
  int x = pedigree.getX();
  int y = pedigree.getY();
  int z = pedigree.getZ();
  int rootID = pedigree.getRootID();

  // 隣接ノードのPedigreeを計算する．
  // 周期境界以外の外部境界面に接している場合は，隣接ノードなしと判断して，0を返す．
  switch (face) {
    case X_M:
      x--;
      if (x < 0) {
        rootID = rootGrid->getNeighborRoot(rootID, face);
        if (rootID < 0) return 0;
        x = max0;
      }
      break;
    case X_P:
      x++;
      if (x > max0) {
        rootID = rootGrid->getNeighborRoot(rootID, face);
        if (rootID < 0) return 0;
        x = 0;
      }
      break;
    case Y_M:
      y--;
      if (y < 0) {
        rootID = rootGrid->getNeighborRoot(rootID, face);
        if (rootID < 0) return 0;
        y = max0;
      }
      break;
    case Y_P:
      y++;
      if (y > max0) {
        rootID = rootGrid->getNeighborRoot(rootID, face);
        if (rootID < 0) return 0;
        y = 0;
      }
      break;
    case Z_M:
      z--;
      if (z < 0) {
        rootID = rootGrid->getNeighborRoot(rootID, face);
        if (rootID < 0) return 0;
        z = max0;
      }
      break;
    case Z_P:
      z++;
      if (z > max0) {
        rootID = rootGrid->getNeighborRoot(rootID, face);
        if (rootID < 0) return 0;
        z = 0;
      }
      break;
    default:
      return 0;
  }

  Pedigree neighborPedigree(level, x, y, z, rootID);

  // Pedigreeからノードを求める.
  Node* neighbor = rootNodes[rootID];
  for (int l = 1; l <= level; l++) {
    if (neighbor->isLeafNode()) break;
    int ijk = neighborPedigree.getX(l)
            + neighborPedigree.getY(l) * 2
            + neighborPedigree.getZ(l) * 4;
    neighbor = neighbor->getChild(ijk);
  }
  // 外部境界の除外はすんでいるので，隣接ノードは必ず存在しなくてはいけない.
  if (neighbor == 0) {
    std::cout << "***error: cannot find neighbor node." << std::endl;
    std::cout << "      node: " << pedigree << std::endl;
    std::cout << "      face: " << face << std::endl;
    std::cout << "  neigibor: " << neighborPedigree << std::endl;
    Exit(EX_FAILURE);
  }

  return neighbor;
}


/// 指定したノードの面が外部境界(周期境界も含む)かどうかチェック.
bool BCMOctree::checkOnOuterBoundary(const Node* node, Face face) const
{
  const Pedigree& pedigree = node->getPedigree();
  int rootID = pedigree.getRootID();
  int max0 = pedigree.getUpperBound() - 1;
  switch (face) {
    case X_M:
      if (pedigree.getX() == 0 && rootGrid->isOuterBoundary(rootID, face)) return true;
      else return false;
    case X_P:
      if (pedigree.getX() == max0 && rootGrid->isOuterBoundary(rootID, face)) return true;
      else return false;
    case Y_M:
      if (pedigree.getY() == 0 && rootGrid->isOuterBoundary(rootID, face)) return true;
      else return false;
    case Y_P:
      if (pedigree.getY() == max0 && rootGrid->isOuterBoundary(rootID, face)) return true;
      else return false;
    case Z_M:
      if (pedigree.getZ() == 0 && rootGrid->isOuterBoundary(rootID, face)) return true;
      else return false;
    case Z_P:
      if (pedigree.getZ() == max0 && rootGrid->isOuterBoundary(rootID, face)) return true;
      else return false;
    default:
      assert(0);
      return false;
  }
}


/// リーフノードリストをランダムシャッフル.
void BCMOctree::randomShuffle()
{
  std::random_shuffle(leafNodeArray.begin(), leafNodeArray.end());
  for (int id = 0; id < leafNodeArray.size(); id++) leafNodeArray[id]->setBlockID(id);
}


/// 指定されたノードの原点位置を取得.
Vec3d BCMOctree::getOrigin(const Node* node) const
{
  const Pedigree& pedigree = node->getPedigree();
  int rootID = pedigree.getRootID();
  int ix = rootGrid->rootID2indexX(rootID);
  int iy = rootGrid->rootID2indexY(rootID);
  int iz = rootGrid->rootID2indexZ(rootID);
  int upperBound = pedigree.getUpperBound();
  return Vec3d(ix + (double)pedigree.getX()/upperBound,
               iy + (double)pedigree.getY()/upperBound,
               iz + (double)pedigree.getZ()/upperBound);
}


/// 指定されたノードの隣接情報を計算.
NeighborInfo* BCMOctree::makeNeighborInfo(const Node* node, const Partition* partition) const
{
  int level = node->getLevel();

  NeighborInfo* neighborInfo = new NeighborInfo[NUM_FACE];

  // 隣接情報
  for (int i = 0; i < NUM_FACE; ++i) {
    Face face = Face(i);
    Node* neighbor = findNeighborNode(node, face);

    if (!neighbor) {
      continue;
    }

    int neighborLevel = neighbor->getLevel();
    int levelDiff = neighborLevel - level;

    if (levelDiff == -1) {
      if (!neighbor->isActive()) continue;
      neighborInfo[face].setLevelDifference(-1);
      neighborInfo[face].setID(neighbor->getBlockID());
      neighborInfo[face].setRank(partition->getRank(neighbor->getBlockID()));
      Subface subface = NeighborInfo::childIdToSubface(face,
                                                       node->getPedigree().getChildId(level));
      neighborInfo[face].setNeighborSubface(subface);
    }
    else if (levelDiff == 0) {
      if (neighbor->isLeafNode()) {
        if (!neighbor->isActive()) continue;
        neighborInfo[face].setLevelDifference(0);
        neighborInfo[face].setID(neighbor->getBlockID());
        neighborInfo[face].setRank(partition->getRank(neighbor->getBlockID()));
      }
      else {
        neighborInfo[face].setLevelDifference(+1);
        for (int j = 0; j < NUM_SUBFACE; j++) {
          Subface subface = Subface(j);
          Node* child = neighbor->getChild(NeighborInfo::getNeighborChildId(face, subface));
          if (!child->isLeafNode()) {
            std::cout << "***error: 2 to 1 constraint is broken." << std::endl;
            std::cout << "      node: " << node->getPedigree() << std::endl;
            std::cout << "   subface: " << subface << std::endl;
            std::cout << "  neigibor: " << neighbor->getPedigree() << std::endl;
            Exit(EX_FAILURE);
          }
          if (!child->isActive()) continue;
          neighborInfo[face].setID(subface, child->getBlockID());
          neighborInfo[face].setRank(subface, partition->getRank(child->getBlockID()));
        }
      }
    }

  }
  return neighborInfo;
}


// ヒルベルトオーダリング 子ノード選択順テーブル
const int BCMOctree::HilbertOrdering[24][8] = { 
  {0, 1, 3, 2, 6, 7, 5, 4},
  {0, 4, 6, 2, 3, 7, 5, 1},
  {0, 1, 5, 4, 6, 7, 3, 2},
  {5, 1, 0, 4, 6, 2, 3, 7},
  {3, 7, 6, 2, 0, 4, 5, 1},
  {6, 7, 3, 2, 0, 1, 5, 4},
  {5, 1, 3, 7, 6, 2, 0, 4},
  {0, 4, 5, 1, 3, 7, 6, 2},
  {5, 4, 0, 1, 3, 2, 6, 7},
  {5, 4, 6, 7, 3, 2, 0, 1},
  {0, 2, 3, 1, 5, 7, 6, 4},
  {6, 4, 0, 2, 3, 1, 5, 7},
  {5, 7, 3, 1, 0, 2, 6, 4},
  {3, 7, 5, 1, 0, 4, 6, 2},
  {6, 4, 5, 7, 3, 1, 0, 2},
  {0, 2, 6, 4, 5, 7, 3, 1},
  {6, 2, 0, 4, 5, 1, 3, 7},
  {6, 2, 3, 7, 5, 1, 0, 4},
  {3, 2, 0, 1, 5, 4, 6, 7},
  {6, 7, 5, 4, 0, 1, 3, 2},
  {5, 7, 6, 4, 0, 2, 3, 1},
  {3, 2, 6, 7, 5, 4, 0, 1},
  {3, 1, 0, 2, 6, 4, 5, 7},
  {3, 1, 5, 7, 6, 4, 0, 2},
};


// ヒルベルトオーダリング 回転テーブル
const int BCMOctree::HilbertOrientation[24][8] = { 
  {1, 2, 0, 3, 4, 0, 5, 6},
  {0, 7, 1, 8, 5, 1, 4, 9},
  {15, 0, 2, 22, 20, 2, 19, 23},
  {20, 6, 3, 23, 15, 3, 16, 22},
  {22, 13, 4, 12, 11, 4, 1, 20},
  {11, 19, 5, 20, 22, 5, 0, 12},
  {9, 3, 6, 2, 21, 6, 17, 0},
  {10, 1, 7, 11, 12, 7, 13, 14},
  {12, 9, 8, 14, 10, 8, 18, 11},
  {6, 8, 9, 7, 17, 9, 21, 1},
  {7, 15, 10, 16, 13, 10, 12, 17},
  {5, 14, 11, 9, 0, 11, 22, 8},
  {8, 20, 12, 19, 18, 12, 10, 5},
  {18, 4, 13, 5, 8, 13, 7, 19},
  {17, 11, 14, 1, 6, 14, 23, 7},
  {2, 10, 15, 18, 19, 15, 20, 21},
  {19, 17, 16, 21, 2, 16, 3, 18},
  {14, 16, 17, 15, 23, 17, 6, 10},
  {13, 21, 18, 17, 7, 18, 8, 16},
  {16, 5, 19, 4, 3, 19, 2, 13},
  {3, 12, 20, 13, 16, 20, 15, 4},
  {23, 18, 21, 10, 14, 21, 9, 15},
  {4, 23, 22, 6, 1, 22, 11, 3},
  {21, 22, 23, 0, 9, 23, 14, 2},
};
