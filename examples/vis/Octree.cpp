#include <iostream>
#include "Octree.h"

Octree::Octree(Node* root) : root(root), maxLevel(0)
{
  checkLeafNode(root, 0);

  max = 1 << maxLevel;  // 2^maxLevel

  for (int id = 0; id < leafNodeArray.size(); id++) leafNodeArray[id]->id = id;
//for (int id = 0; id < leafNodeArray.size(); id++) {
//  std::cout << id << ": pedigree=" << leafNodeArray[id]->pedigree << std::endl;
//}
}


void Octree::checkLeafNode(Node* node, int level)
{
  if (node->isLeafNode()) {
    leafNodeArray.push_back(node);
    if (level > maxLevel) maxLevel = level;
    return;
  }
  for (int i = 0; i < 8; i++) checkLeafNode(node->child[i], level+1);
}


Octree::~Octree()
{
  deleteNode(root);
}


void Octree::deleteNode(Node* node)
{
  if (node->isLeafNode()) {
    delete node;
    return;
  }
  for (int i = 0; i < 8; i++) deleteNode(node->child[i]);
}


Node* Octree::searchNode(const Node* node, Face face, bool periodic) const
{
  Pedigree pedigree = node->pedigree.findNeighbor(face, max, periodic);
//std::cout << pedigree;

  if (pedigree.level < 0) return 0;

  Node* neighbor = root;
  for (int level = 1; level <= pedigree.level; level++) {
    int ijk = pedigree.getX(level)
            + pedigree.getY(level) * 2
            + pedigree.getZ(level) * 4;
//  std::cout  << "(ijk=" << ijk << ")";
    if (neighbor->child[ijk] == 0) break;
    neighbor = neighbor->child[ijk];
  }
  return neighbor;
}
