/*
###################################################################################
#
# BCMTools
#
# Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2016 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
*/

///
/// @file Partition.h
/// @brief 1次元ブロック領域分割用ユーティリティクラス
///

#ifndef PARTITION_H
#define PARTITION_H

#include "mpi.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>

/// 1次元ブロック領域分割用ユーティリティクラス.
class Partition {

  int nProcs;  ///< プロセス数
  int nItems;  ///< 全要素数

  std::vector<int> end;  ///< 各プロセスの末尾要素番号を納めたリスト

public:

  /// コンストラクタ.
  ///  @param[in] nProcs プロセス数
  ///  @param[in] nItems 全要素数
  ///
  Partition(int nProcs, int nItems) : nProcs(nProcs), nItems(nItems), end(nProcs) {
    int m = nItems / nProcs;
    int r = nItems % nProcs;
    int end0 = 0;
    for (int i = 0; i < nProcs; i++) {
      if (i < r) {
        end[i] = end0 + m + 1;
      } else {
        end[i] = end0 + m;
      }
      end0 = end[i];
    }
    assert(end[nProcs-1] == nItems);
  }

  /// デストラクタ.
  ~Partition() {}

  /// 先頭要素番号の取得.
  ///
  ///  @param[in] rank プロセス番号
  ///  @return 先頭要素番号
  ///
  int getStart(int rank) const {
    assert(0 <= rank && rank < nProcs);
    if (rank == 0) return 0;
    return end[rank-1];
  }

  /// 末尾要素番号の取得.
  ///
  ///  @param[in] rank プロセス番号
  ///  @return 末尾要素番号+1
  ///
  int getEnd(int rank) const {
    assert(0 <= rank && rank < nProcs);
    return end[rank];
  }

  /// 担当要素数を取得.
  ///
  ///  @param[in] rank プロセス番号
  ///  @return 担当要素数
  ///
  int getNum(int rank) const {
    assert(0 <= rank && rank < nProcs);
    if (rank == 0) return end[0];
    return end[rank] - end[rank-1];
  }

  /// 担当プロセスを取得
  ///
  ///  @param[in] i 要素番号
  ///  @return 担当プロセス番号
  ///
  ///  @note 範囲外の要素番号が指定された場合MPI::PROC_NULLを返す.
  ///
  int getRank(int i) const {
    if (i < 0 || i >= nItems) return MPI::PROC_NULL;
    std::vector<int>::const_iterator it = std::upper_bound(end.begin(), end.end(), i);
    assert(it != end.end());
    return it - end.begin();
  }

  /// 分割内容を出力.
  void print() const {
    int start = 0;
    for (int rank = 0; rank < nProcs; rank++) {
      std::cout << rank << ": [" << start << ":" << end[rank]-1 << "]"
             // << " [" << getStart(rank) << ":" << getEnd(rank) << ")"
                << " #" << end[rank]-start << std::endl;
      start = end[rank];
    }
  }

};

#endif // PARTITION_H
