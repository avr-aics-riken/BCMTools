/*
 * BCMTools
 *
 * Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file BoundaryConditionSetterBase.h
/// @brief 境界条件設定クラス(基底クラス)
/// 

#ifndef BOUNDARY_CONDITION_SETTER_BASE_H
#define BOUNDARY_CONDITION_SETTER_BASE_H

class BoundaryInfo;
class Node;
class BCMOctree;

/// 境界条件設定クラス(基底クラス).
///
/// 境界条件設定クラスのインタフェースを規定.
///
class BoundaryConditionSetterBase {

public:

  /// コンストラクタ.
  BoundaryConditionSetterBase() {}

  /// デストラクタ.
  virtual ~BoundaryConditionSetterBase() {}

  /// 境界条件設定ファンクタ.
  ///
  ///  @param[in] node ノード(ブロック)
  ///  @param[in] tree Octree
  ///
  virtual BoundaryInfo* operator() (const Node* node, const BCMOctree* tree) = 0;


};

#endif // BOUNDARY_CONDITION_SETTER_BASE_H
