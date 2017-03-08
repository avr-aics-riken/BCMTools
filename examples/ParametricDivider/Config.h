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

  Vec3i rootN;         ///< ルートノード配置

  int minLevel;  ///< 最小分割レベル
  int maxLevel;  ///< 最大分割レベル

  std::string treeType;  ///< ツリータイプ

  std::string ordering;  ///< オーダリング方法

  REAL_TYPE cx, cy, cz;
  REAL_TYPE r;

  int size;   ///< ブロック内セル分割数
  std::string output;  ///< 結果出力ファイル
  bool verbose;          ///< 冗長メッセージフラグ


private:

  void parse() {
    Vec3r sphereCenter;
    Vec3r cylinderCenter;
    REAL_TYPE sphereRadius;
    REAL_TYPE cylinderRadius;

    rootN = read<Vec3i>("rootGrid", Vec3i(1, 1, 1));

    minLevel = read<int>("minLevel", 0);
    maxLevel = read<int>("maxLevel");
    treeType = read<string>("treeType");
    ordering = read<string>("ordering");

    size = read<int>("size");
    output = read<string>("output");

    verbose = read<bool>("verbose", false);

    sphereCenter = read<Vec3r>("sphereCenter", Vec3r(0.5, 0.5, 0.5));
    sphereRadius = read<REAL_TYPE>("sphereRadius", 0.25);

    cylinderCenter = read<Vec3r>("cylinderCenter", Vec3r(0.5, 0.5, 0.0));
    cylinderRadius = read<REAL_TYPE>("cylinderRadius", 0.4);

    cx = cy = cz = r = 0.0;
    if (treeType == "sphere") {
      cx = sphereCenter.x;
      cy = sphereCenter.y;
      cz = sphereCenter.z;
      r = sphereRadius;
    }
    if (treeType == "cylinder") {
      cx = cylinderCenter.x;
      cy = cylinderCenter.y;
      r = cylinderRadius;
    }
  }

  bool validate() {
    bool ret = true;
    if (!(treeType == "flat" ||
          treeType == "simple" ||
          treeType == "sphere" ||
          treeType == "cylinder" )) {
      std::cout << "error: 'treeType' must be 'flat', 'simple', 'sphere' or 'cylinder'." << std::endl;
      ret = false;
    }

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
    std::cout << "  min level:          " << minLevel << std::endl;
    std::cout << "  max level:          " << maxLevel << std::endl;
    std::cout << "  tree type:          " << treeType;
    if (treeType == "sphere") {
      std::cout << " [c=" << Vec3r(cx, cy, cz) << ", r=" << r << "]";
    }
    if (treeType == "cylinder") {
      std::cout << " [c=(" << cx << ", " << cy << "), r=" << r << "]";
    }
    std::cout << std::endl;
    std::cout << "  ordering:           " << ordering << std::endl;
    std::cout << "  block size:         " << size << std::endl;
    std::cout << "  output file:        " << output << std::endl;
    std::cout << "  verbose message:    " << (verbose ? "on" : "off") << std::endl;
  }

};



#endif // CONFIG_H
