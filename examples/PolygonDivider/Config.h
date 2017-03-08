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

#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <string>
#include <vector>
#include "ConfigBase.h"
#include "Vec3.h"

using namespace Vec3class;

class Config : public ConfigBase {

public:

  Vec3d origin;        ///< 原点座標
  double rootLength;   ///< ルートノードボックスの辺長
  Vec3i rootN;         ///< ルートノード配置

  int minLevel;  ///< 最小分割レベル
  int maxLevel;  ///< 最大分割レベル

  std::string polylibConf;  ///< Polylib設定ファイル
  std::string polygonGroup; ///< ポリゴングループ名
  bool polygonInsideOut;    ///< ポリゴンの表裏を反転

  std::string ordering;  ///< オーダリング方法

  int size;   ///< ブロック内セル分割数
  std::string output;  ///< 結果出力ファイル

  bool verbose;          ///< 冗長メッセージフラグ

private:

  void parse() {

    origin = read<Vec3d>("origin", Vec3d(0, 0, 0));
    rootLength = read<double>("rootLength", 1.0);
    rootN = read<Vec3i>("rootGrid", Vec3i(1, 1, 1));

    minLevel = read<int>("minLevel", 0);
    maxLevel = read<int>("maxLevel", 5);

    polylibConf = read<string>("polylibConfig");
    polygonGroup = read<string>("polygonGroup");
    polygonInsideOut = read<bool>("polygonInsideOut", false);

    ordering = read<string>("ordering");

    size = read<int>("size");
    output = read<string>("output");

    verbose = read<bool>("verbose", false);
  }

  bool validate() {
    bool ret = true;
    if (!(ordering == "Z" || ordering == "Hilbert" || ordering == "random")) {
      std::cout << "error: 'ordering' must be 'Z', 'Hilbert' or 'random'."
                << std::endl;
      ret = false;
    }
    if (size < 1) {
      std::cout << "error: 'size' must be > 1" << std::endl;
      ret = false;
    }
    if (!(minLevel >= 0)) {
      std::cout << "error: 'minLevel' must be >= 0." << std::endl;
      ret = false;
    }
    if (!(minLevel <= maxLevel)) {
      std::cout << "error: 'minLevel' and 'maxLevel'  must be minLevel <= maxLevel." << std::endl;
      ret = false;
    }

    return ret;
  }

public:

  void print() const {
    std::cout.setf(std::ios::showpoint);
    std::cout << "  root grid:          " << rootN << std::endl;
    std::cout << "  root block size:    " << rootLength << std::endl;
    std::cout << "  root grid:          " << rootN << std::endl;
    std::cout << "  min level:          " << minLevel << std::endl;
    std::cout << "  max level:          " << maxLevel << std::endl;
    std::cout << "  Polylib conf. file: " << polylibConf << std::endl;
    std::cout << "  polygon group:      " << polygonGroup << std::endl;
    std::cout << "  polygon inside-out: " << (polygonInsideOut ? "yes" : "no") << std::endl;
    std::cout << "  ordering:           " << ordering << std::endl;
    std::cout << "  block size:         " << size << std::endl;
    std::cout << "  output file:        " << output << std::endl;
    std::cout << "  verbose message:    " << (verbose ? "on" : "off") << std::endl;
  }

};



#endif // CONFIG_H
