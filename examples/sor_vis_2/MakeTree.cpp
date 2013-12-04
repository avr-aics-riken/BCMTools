#include "MakeTree.h"
#include "Octree.h"

namespace {

void makeFlatTreeNodes(Node* parent, int maxLevel)
{
  if (parent->pedigree.level >= maxLevel) return;

  for (int ijk = 0; ijk < 8; ijk++) {
    Node* child = new Node(parent, ijk);
    makeFlatTreeNodes(child, maxLevel);
    parent->child[ijk] = child;
  }
}


Octree* makeFlatTree(int maxLevel)
{
  Node* root = new Node;
  makeFlatTreeNodes(root, maxLevel);

  return new Octree(root);
}


void makeSimpleTreeNodes(Node* parent, int max, int maxLevel)
{
  if (parent->pedigree.level >= maxLevel) return;

  if (!(parent->pedigree.x == 0 ||
        parent->pedigree.y == 0 ||
        parent->pedigree.z == 0 ||
        parent->pedigree.x == max-1 ||
        parent->pedigree.y == max-1 ||
        parent->pedigree.z == max-1)) return;

  for (int ijk = 0; ijk < 8; ijk++) {
    Node* child = new Node(parent, ijk);
    makeSimpleTreeNodes(child, max*2, maxLevel);
    parent->child[ijk] = child;
  }
}

Octree* makeSimpleTree(int maxLevel)
{
  Node* root = new Node;
  int max = 1;
  makeSimpleTreeNodes(root, max, maxLevel);

  return new Octree(root);
}

void makeSimpleTreeNodes2(Node* parent, int max, int maxLevel)
{
  if (parent->pedigree.level >= maxLevel) return;

/*
  if (!(parent->pedigree.x == 0 && 
        parent->pedigree.y == 0 &&
        parent->pedigree.z == 0)) return;
*/

	if (!(parent->pedigree.x >= 0 && parent->pedigree.x <= 1 &&
				parent->pedigree.y >= 0 && parent->pedigree.y <= 1 &&
				parent->pedigree.z >= 0 && parent->pedigree.z <= 1 )) return;

  for (int ijk = 0; ijk < 8; ijk++) {
    Node* child = new Node(parent, ijk);
    makeSimpleTreeNodes2(child, max*2, maxLevel);
    parent->child[ijk] = child;
  }
}

Octree* makeSimpleTree2(int maxLevel)
{
  Node* root = new Node;
  int max = 1;
  makeSimpleTreeNodes2(root, max, maxLevel);

  return new Octree(root);
}

}


Octree* makeTree(const std::string& type, int level)
{
  if      (type == "flat") {
    return makeFlatTree(level);
  } 
  else if (type == "simple") {
    return makeSimpleTree(level);
  }
  else if (type == "simple2") {
    return makeSimpleTree2(level);
  }
  else {
    exit(EX_READ_CONFIG);
  }
}

