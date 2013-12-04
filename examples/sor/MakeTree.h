#ifndef MAKE_TREE_H
#define MAKE_TREE_H

#include <string>

class Octree;

Octree* makeTree(const std::string& type, int level);

#endif // MAKE_TREE_H
