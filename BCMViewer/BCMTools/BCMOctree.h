/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file BCMOctree.h
/// @brief BCM用マルチルートOCtreeクラス
/// 

#ifndef BCM_OCTREE_H
#define BCM_OCTREE_H

#include <vector>
//#include "BCMTools.h"
#include "Vec3.h"
#include "RootGrid.h"
//#include "Divider.h"
#include "Pedigree.h"
#include "Node.h"
//#include "NeighborInfo.h"
#include "Partition.h"

/// フェイス番号.
enum Face { X_M, X_P, Y_M, Y_P, Z_M, Z_P, NUM_FACE };

/// BCM用マルチルートOCtreeクラス.
class BCMOctree {

public:

  /// オーダリングタイプ.
  enum Ordering {
    Z,             ///< Z(Morton)オーダリング
    HILBERT,       ///< ヒルベルトオーダリング
    RANDOM,        ///< ランダムシャッフル
	PEDIGREELIST,  ///< Pedigreeリスト順(ファイルロード用)
  };

private:

  RootGrid* rootGrid; ///< ルートノード配置情報

  //Divider* divider;   ///< ブロック分割判定クラス

  Ordering ordering;  ///< オーダリング方法

  Node** rootNodes;   ///< ルートノード配列

  std::vector<Node*> leafNodeArray;  ///< リーフノードリスト

  static const int HilbertOrdering[24][8];    ///< ヒルベルトオーダリング 子ノード選択順テーブル
  static const int HilbertOrientation[24][8]; ///< ヒルベルトオーダリング 回転テーブル

public:

  /// コンストラクタ.
  ///
  ///  @param[in] rootGrid ルートノード配置情報
  ///  @param[in] divider ブロック分割判定クラス
  ///  @param[in] ordering オーダリング方法
  ///
  ///  @note rootGridとdividerは、デストラクタにより解放される.
  ///
  //BCMOctree(RootGrid* rootGrid, Divider* divider, Ordering ordering);


  /// コンストラクタ.
  ///
  ///  @param[in] rootGrid ルートノード配置情報
  ///  @param[in] pedigrees ペディグリリスト
  ///
  ///  @note rootGridは、デストラクタにより解放される．
  ///
  BCMOctree(RootGrid* rootGrid, const std::vector<Pedigree>& pedigrees);

  /// デストラクタ.
  ~BCMOctree();

/*
  /// Octree情報を他rankにブロードキャスト.
  ///
  ///  @param[in] comm MPIコミュニケータ
  ///
  ///  @note rank0のみが呼ぶこと
  ///
  void broadcast(MPI::Intracomm& comm = MPI::COMM_WORLD);

  /// rank0からOctree情報を受信.
  ///
  ///  @param[in] comm MPIコミュニケータ
  ///  @return 新たに生成したBCMOctreeインスタンス
  ///
  ///  @note rank0からは呼ばないこと
  ///
  static BCMOctree* ReceiveFromMaster(MPI::Intracomm& comm = MPI::COMM_WORLD);
*/
  /// ルートグリッドを取得
  const RootGrid* getRootGrid() const { return rootGrid; }

  /// リーフノード総数を取得.
  int getNumLeafNode() const { return leafNodeArray.size(); }

  /// リーフノードリストを取得.
  std::vector<Node*>& getLeafNodeArray() { return leafNodeArray; }

  /// リーフノードリストを取得.
  const std::vector<Node*>& getLeafNodeArray() const { return leafNodeArray; }

  const Node* getRootNode(const int rootID) const { 
    if( rootID >= rootGrid->getSize() ){
	  return NULL;
	}
    return rootNodes[rootID];
  }

  /// 指定したノードの面が外部境界(周期境界も含む)かどうかチェック.
  ///
  ///  @param[in] node ノード
  ///  @param[in] face 面
  ///  @return 外部境界ならtrue
  ///
  bool checkOnOuterBoundary(const Node* node, Face face) const;

  /// 指定されたノードの原点位置を取得.
  ///
  ///  @param[in] node ノード
  ///  @return 原点位置
  ///
  Vec3r getOrigin(const Node* node) const;

  /// 指定されたノードの隣接情報を計算.
  ///
  ///  @param[in] node ノード
  ///  @param[in] partition 領域分割情報
  ///  @return 隣接情報クラス
  ///
  //NeighborInfo* makeNeighborInfo(const Node* node, const Partition* partition) const;

private:
  /// コンストラクタ(リーフノードのPedigreeリストから).
  ///
  ///  @param[in] rootGrid ルートノード配置情報
  ///  @param[in] ordering オーダリング方法
  ///  @param[in] numLeafNode リーフノード総数
  ///  @param[in] buf シリアライズされたPedigreeリスト
  ///
  BCMOctree(RootGrid* rootGrid, Ordering ordering,
            int numLeafNode, const unsigned char* buf);

  /// リーフノードのPedigreeリストからOctreeを再構築.
  ///
  ///  @param[in] numLeafNode リーフノード総数
  ///  @param[in] buf シリアライズされたPedigreeリスト
  ///
  void buildTreeFromPedigreeList(int numLeafNode, const unsigned char* buf);

  /// Zオーダリング順にリーフノードリストを作成.
  ///
  ///  @param[in] node ルートノード
  ///
  ///  @note 再帰呼び出しされる
  ///
  void pickupLeafNodeZOrdering(Node* node);

  /// ヒルベルトオーダリング順にリーフノードリストを作成.
  ///
  ///  @param[in] node ルートノード
  ///  @param[in] orientation 回転向き
  ///
  ///  @note 再帰呼び出しされる
  ///
  void pickupLeafNodeHilbertOrdering(Node* node, int orientation);

  /// リーフノードリストをランダムシャッフル.
  void randomShuffle();

  /// 再帰呼び出しによりツリー生成.
  ///
  ///  @param[in] node ノード
  ///
  ///  @note 内部でdividerクラスにより分割判定を行っている
  ///
  void makeNode(Node* node);

  /// 再帰呼び出しによるツリー消去.
  ///
  ///  @param[in] node ノード
  ///
  void deleteNode(Node* node);

  /// 隣接ノード探索.
  ///
  ///  @param[in] node 基準ノード
  ///  @param[in] face 面
  ///  @return 隣接ノードへのポインタ
  ///
  ///  @note 隣接ノードが非アクティブな場合でも，それを返す.
  ///        周期境界条件以外の外部境界に接している場合のみ0を返す.
  ///
  //Node* findNeighborNode(const Node* node, Face face) const;

  /// Pedigree情報を通信用バッファにパック.
  ///
  ///  @param[in] node 対象ノード
  ///  @param[in] ip バッファ内の位置インデクス
  ///  @param[in] buf バッファ配列
  ///
  void packPedigrees(Node* node, size_t& ip, unsigned char* buf);

//static void unpackPedigrees(int numLeafNode, const unsigned char* buf, Pedigree* pedigrees);




};


#endif // BCM_OCTREE_H
