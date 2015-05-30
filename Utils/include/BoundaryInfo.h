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
/// @file BoundaryInfo.h
/// @brief 境界条件情報(サンプル)クラス
/// 

#ifndef BUNDARY_INFO_H
#define BUNDARY_INFO_H

#include <iostream>
#include "BCMTools.h"
#include "BCMOctree.h"

/// 境界条件情報(サンプル)クラス.
class BoundaryInfo {

public:

  /// 境界条件タイプ
  enum Type {
    INNER,      ///< 内部境界
    DIRICHLET,  ///< Dirichlet
    NEUMANN,    ///< Neumann
    PERIODIC,   ///< 周期境界
  };

private:

  Type type;  ///< 境界条件タイプ

  int id;     ///< 境界ID(内部境界の場合は-1)

public:

  /// コンストラクタ.
  ///
  ///  @note ディフォルトは内部境界
  ///
  BoundaryInfo() : type(INNER), id(-1) {}

  /// デストラクタ.
  ~BoundaryInfo() {}

  /// 境界タイプを設定.
  void setType(Type type) { this->type = type; }

  /// 境界IDを設定.
  void setID(int id) { this->id = id; }

  /// 境界タイプを取得.
  Type getType() const { return type; }

  /// 境界IDを取得.
  int getID() const { return id; }

  /// 境界情報を出力.
  void print() const {
    switch (type) {
      case INNER:     std::cout << "type=INNER, "; break;
      case DIRICHLET: std::cout << "type=DIRICHLET, "; break;
      case NEUMANN:   std::cout << "type=NEUMANN, "; break;
      case PERIODIC:  std::cout << "type=PERIODIC, "; break;
      default: break;
    }
    std::cout << "id=" << id << std::endl;
  }

  /// 6面の境界情報を出力.
  static void print(const BoundaryInfo* bInfo) {
    for (int i = 0; i < NUM_FACE; i++) {
      Face face = Face(i);
      switch (face) {
        case X_M: std::cout << "X_M: "; break;
        case X_P: std::cout << "X_P: "; break;
        case Y_M: std::cout << "Y_M: "; break;
        case Y_P: std::cout << "Y_P: "; break;
        case Z_M: std::cout << "Z_M: "; break;
        case Z_P: std::cout << "Z_P: "; break;
        default: break;
      }
      bInfo[i].print();
    }
  }

};


#endif // BUNDARY_INFO_H
