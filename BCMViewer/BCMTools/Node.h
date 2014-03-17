/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file Node.h
/// @brief Octree用ノードクラス
/// 

#ifndef NODE_H
#define NODE_H

#include "Pedigree.h"
#include "Vec3.h"
#include <assert.h>

using namespace Vec3class;

/// Octreeノードクラス.
class Node {

  Node* parent;      ///< 親ノードへのポインタ

  Node** childList;  ///< 子ノードリスト

  bool active;       ///< アクティブノードフラグ

  int id;   ///< ブロックID(アクティブなリーフノード以外には-1を入れる)

  Pedigree pedigree;  ///< Pedigree

public:

  /// コンストラクタ(ルートノードとして生成).
  Node(int rootID = 0) 
    : parent(0), childList(0), active(true), id(-1), pedigree(rootID) {}

  /// コンストラクタ(子ノードとして生成).
  ///
  ///  @param[in] parent 親ノード
  ///  @param[in] i  子ノード番号(0〜7)
  ///
  Node(Node* parent, int i)
   : parent(parent), childList(0), active(true), id(-1), pedigree(parent->pedigree, i) {}

  /// デストラクタ.
  ~Node() {
    if (childList) {
      for (int i = 0; i < 8; i++) delete childList[i];
      delete[] childList;
    }
  }

  /// ルートノード判定.
  ///
  ///  @return ルートノードの場合true
  ///
  bool isRootNode() const { return parent == 0; }

  /// リーフノード判定.
  ///
  ///  @return リーフノードの場合true
  ///
  bool isLeafNode() const { return childList == 0; }

  /// アクティブノード判定.
  ///
  ///  @return アクティブノードの場合true
  ///
  bool isActive() const { return active; }

  /// アクティブノードフラグの設定.
  ///
  ///  @param[in] OnOff アクティブノードフラグ値
  ///
  void setActive(bool OnOff = true) { active = OnOff; }

  /// ブロックIDを取得.
  ///
  ///  @return ブロックID
  ///
  int getBlockID() const { return id; }

  /// ブロックIDを設定.
  ///
  ///  @param[in] id ブロックID
  ///
  void setBlockID(int id) { this->id = id; }

  /// Pedigreeを取得.
  ///
  ///  @return Pedigree
  ///
  const Pedigree& getPedigree() const { return pedigree; }

  /// ツリーレベルを取得.
  ///
  ///  @return ツリーレベル
  ///
  int getLevel() const { return pedigree.getLevel(); }

  /// 8つの子ノードを生成.
  void makeChildNodes() {
    assert(childList == 0);
    childList = new Node*[8];
    for (int i = 0; i < 8; i++) childList[i] = new Node(this, i);
  }

  /// 規格化されたブロックサイズを計算.
  ///
  ///  @return ブロックサイズ
  ///
  Vec3d getBlockSize() const {
    int upperBound = pedigree.getUpperBound();
    return Vec3d(1.0/upperBound, 1.0/upperBound, 1.0/upperBound);
  }

  /// 親ノードを取得.
  ///
  ///  @return 親ノードへのポインタ
  ///
  Node* getParent() { return parent; }

  /// 子ノードを所得.
  ///
  ///  @param[in] i  子ノード番号(0〜7)
  ///  @return 子ノードへのポインタ
  ///
  Node* getChild(int i) {
 // if (!childList) return 0;
    assert(childList);
    assert(0 <= i && i < 8);
    return childList[i];
  }
  
  const Node* getChild(int i ) const {
  	return childList[i];
  }

};


#endif // NODE_H
