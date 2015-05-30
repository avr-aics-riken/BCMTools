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

///
/// @file Divider.h
/// @brief ブロック分割判定クラス(基底クラス)
/// 

#ifndef DIVIDIER_H
#define DIVIDIER_H

#include "RootGrid.h"
#include "Pedigree.h"

/// ブロック分割判定クラス(基底クラス).
class Divider {

public:

  /// ブロック(ノード)タイプ型
  enum NodeType {
    BRANCH,          ///< 枝(分割を続ける)
    LEAF_ACTIVE,     ///< アクティブなリーフノード(分割を終了)
    LEAF_NO_ACTIVE,  ///< 非アクティブなリーフノード(分割を終了)
  };

public:

  /// コンストラクタ.
  Divider() {}

  /// デストラクタ.
  virtual ~Divider() {}

  /// ブロックを分割するかどうかを判定.
  ///
  ///   @param[in] pedigree ブロックのPedigree
  ///   @return ブロックタイプ
  ///
  virtual NodeType operator() (const Pedigree& pedigree) = 0;

};

#endif // DIVIDER_H

