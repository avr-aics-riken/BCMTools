#ifndef OCTREE_H
#define OCTREE_H

#include <vector>
#include "BCMTools.h"
#include "Pedigree.h"

struct Node {
  Node* parent;
  Node* child[8];

  int id;   // リーフノード以外には-1を入れる

  Pedigree pedigree;

  Node() : parent(0), id(-1) {
    for (int i = 0; i < 8; i++) child[i] = 0;
  }

  Node(Node* parent, int ijk)
   : parent(parent), pedigree(parent->pedigree, ijk) {
    for (int i = 0; i < 8; i++) child[i] = 0;
  }

  ~Node() {}

  bool isRootNode() const { return parent == 0; }

  bool isLeafNode() const { 
    for (int i = 0; i < 8; i++) if (child[i]) return false;
    return true;
  }

};

class Octree {
  Node* root;
  int maxLevel;
  int max;

  std::vector<Node*> leafNodeArray;

public:

  Octree(Node* root);

  ~Octree();

  int getMaxLevel() const { return maxLevel; }

  Node* getRootNode() const { return root; }

  int getNumLeafNode() const { return leafNodeArray.size(); }

  std::vector<Node*>& getLeafNodeArray() { return leafNodeArray; }

  const std::vector<Node*>& getLeafNodeArray() const { return leafNodeArray; }

  Node* searchNode(const Node* node, Face face, bool periodic=false) const;

private:

  void checkLeafNode(Node* node, int level);

  void deleteNode(Node* node);

};


#endif // OCTREE_H
