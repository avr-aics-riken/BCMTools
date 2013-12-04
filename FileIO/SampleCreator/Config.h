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

#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <string>
#include <vector>
#include "ConfigBase.h"
#include "BoundingBox.h"


/// (ポリゴングループ名, 分割レベル)ペア.
struct PolygonGroupSpec {
  std::string polygonGroupName;
  int level;

  friend std::ostream& operator<<(std::ostream& os, const PolygonGroupSpec& pgs) {
    return os << pgs.polygonGroupName << " -> " << pgs.level;
  }

  friend std::istream& operator>>(std::istream& is, PolygonGroupSpec& pgs) {
    return is >> pgs.polygonGroupName >> pgs.level;
  }
};

inline std::istream& operator>>(std::istream& is,
                                std::vector<PolygonGroupSpec>& polygonGroupList) {
  while (!is.eof()) {
    PolygonGroupSpec pgs;
    is >> pgs;
    polygonGroupList.push_back(pgs);
  }
  return is;
}


/// (バウンディングボックス, 分割レベル)ペア.
struct BoundingBoxSpec {
  BoundingBox boundingBox;
  int level;

  friend std::ostream& operator<<(std::ostream& os, const BoundingBoxSpec& bbs) {
    return os << bbs.boundingBox << " -> " << bbs.level;
  }

  friend std::istream& operator>>(std::istream& is, BoundingBoxSpec& bbs) {
    return is >> bbs.boundingBox >> bbs.level;
  }
};

inline std::istream& operator>>(std::istream& is,
                                std::vector<BoundingBoxSpec>& boundingBoxList) {
  while (!is.eof()) {
    BoundingBoxSpec bbs;
    is >> bbs;
    boundingBoxList.push_back(bbs);
  }
  return is;
}


class Config : public ConfigBase {

public:

  Vec3r origin;        ///< 原点座標
  double rootLength;   ///< ルートノードボックスの辺長
  Vec3i rootN;         ///< ルートノード配置

  int minLevel;  ///< 最小分割レベル
  int maxLevel;  ///< 最大分割レベル

  std::string polylibConfig;  ///< Polylib設定ファイル

  std::vector<PolygonGroupSpec> polygonGroupList;  ///< (ポリゴングループ,分割レベル)ペアのリスト

  std::vector<BoundingBoxSpec> boundingBoxList;  ///< (バウンディングボックス,分割レベル)ペアのリスト

  std::string ordering;  ///< オーダリング方法

  int size;   ///< ブロック内セル分割数
  int vc;     ///< 仮想セル数

  std::string output;  ///< 結果出力ファイル


private:

  void parse() {
    origin = read<Vec3r>("origin", Vec3r(0, 0, 0));
    rootLength = read<double>("rootLength", 1.0);
    rootN = read<Vec3i>("rootGrid", Vec3i(1, 1, 1));

    minLevel = read<int>("minLevel", 0);
    maxLevel = read<int>("maxLevel", 5);

    polylibConfig = read<string>("polylibConfig");

    polygonGroupList = read<std::vector<PolygonGroupSpec> >(
                         "polygonGroupList",
                         std::vector<PolygonGroupSpec>());
    
    boundingBoxList = read<std::vector<BoundingBoxSpec> >(
                        "boundingBoxList",
                        std::vector<BoundingBoxSpec>());

    ordering = read<string>("ordering");

    size = read<int>("size");
    vc = read<int>("vc", 1);

    output = read<string>("output");
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
    if (!(vc > 0 && vc <= size/2)) {
      std::cout << "error: 'vc' must be 0 < vc && vc <= size/2." << std::endl;
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
    std::cout << "  ordering:           " << ordering << std::endl;
    std::cout << "  block size:         " << size << std::endl;
    std::cout << "  vc width:           " << vc << std::endl;
    std::cout << "  output file:        " << output << std::endl;
    std::cout << "  polilib config:     " << polylibConfig << std::endl;
    std::cout << "  polygonGroupList:" << std::endl;
    for (int i = 0; i < polygonGroupList.size(); i++) {
      std::cout << "       " << polygonGroupList[i] << std::endl;
    }
    std::cout << "  boundingBoxList:" << std::endl;
    for (int i = 0; i < boundingBoxList.size(); i++) {
      std::cout << "       " << boundingBoxList[i] << std::endl;
    }
  }

};



#endif // CONFIG_H
