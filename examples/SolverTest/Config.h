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
#include "Vec3.h"
#include "PolygonBBoxDivider.h"

/// <ポリゴングループ名, 分割レベル>ペアリストの読み込み
inline std::istream& operator>>(std::istream& is,
                                std::vector<PolygonGroupSpec>& polygonGroupList) {
  while (!is.eof()) {
    PolygonGroupSpec pgs;
    is >> pgs;
    polygonGroupList.push_back(pgs);
  }
  return is;
}


/// <バウンディングボックス, 分割レベル>ペアリストの読み込み.
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
	string operatorname;

  Vec3r origin;        ///< 原点座標
  double rootLength;   ///< ルートノードボックスの辺長
  Vec3i rootN;         ///< ルートノード配置
  Vec3r seed;
	bool periodicX;
	bool periodicY;
	bool periodicZ;

  int minLevel;  ///< 最小分割レベル
  int maxLevel;  ///< 最大分割レベル

	double dividerOx;
	double dividerOy;
	double dividerOz;
	double dividerR;
	double dividerDR;
	double dividerX0;
	double dividerY0;
	double dividerZ0;
	double dividerX1;
	double dividerY1;
	double dividerZ1;
	bool dividerHollow;

  std::string treeType;  ///< ツリータイプ

  std::string polylibConfig;  ///< Polylib設定ファイル
  std::vector<PolygonGroupSpec> polygonGroupList;  ///< (ポリゴングループ,分割レベル)ペアのリスト
  std::vector<BoundingBoxSpec> boundingBoxList;  ///< (バウンディングボックス,分割レベル)ペアのリスト
  std::vector<BoundingBoxSpec> sphericalBoxList;  ///< (バウンディングボックス,分割レベル)ペアのリスト

  std::string ordering;  ///< オーダリング方法

//  bool separate;         ///< 方向毎仮想セル同期フラグ
	std::string updateMethod;

  int size;   ///< ブロック内セル分割数
  int vc;     ///< 仮想セル幅

	double omegaU;
	int countMaxU;
	double epsilonU;
	int countPreConditionerU;

	double omegaP;
	int countMaxP;
	double epsilonP;
	int countPreConditionerP;

	bool bHeat;
	double omegaT;
	int countMaxT;
	double epsilonT;
	int countPreConditionerT;

	double dt;
	double rhof;
	double rhos;
	double cpf;
	double cps;
	double kf;
	double ks;
	double mu;
	double gx;
	double gy;
	double gz;
	double uxs0;
	double uys0;
	double uzs0;

	int StepStart;
	int StepEnd;

	int StepPrintTime;
	int StepPrintStats;
	int StepPrintData;
	int StepPrintBin;

	std::string advection;	

  std::string output;  ///< 結果出力ファイル
  bool verbose;          ///< 冗長メッセージフラグ

	bool cutoff;
	double cutoff_epsilon;
	bool voxelization;
	bool symmetrization;
	bool holefilling;
	bool holefilling2;
	bool masking;

	int AccelDuration;
	int STEP_ACCELDURATION;

	bool vtkWriteForce;

	bool BCMFileSave;
	bool PLOT3DSave;

	int boundaryTypeUX_X_M;
	int boundaryTypeUX_X_P;
	int boundaryTypeUX_Y_M;
	int boundaryTypeUX_Y_P;
	int boundaryTypeUX_Z_M;
	int boundaryTypeUX_Z_P;
	double boundaryValueUX_X_M;
	double boundaryValueUX_X_P;
	double boundaryValueUX_Y_M;
	double boundaryValueUX_Y_P;
	double boundaryValueUX_Z_M;
	double boundaryValueUX_Z_P;

	int boundaryTypeUY_X_M;
	int boundaryTypeUY_X_P;
	int boundaryTypeUY_Y_M;
	int boundaryTypeUY_Y_P;
	int boundaryTypeUY_Z_M;
	int boundaryTypeUY_Z_P;
	double boundaryValueUY_X_M;
	double boundaryValueUY_X_P;
	double boundaryValueUY_Y_M;
	double boundaryValueUY_Y_P;
	double boundaryValueUY_Z_M;
	double boundaryValueUY_Z_P;

	int boundaryTypeUZ_X_M;
	int boundaryTypeUZ_X_P;
	int boundaryTypeUZ_Y_M;
	int boundaryTypeUZ_Y_P;
	int boundaryTypeUZ_Z_M;
	int boundaryTypeUZ_Z_P;
	double boundaryValueUZ_X_M;
	double boundaryValueUZ_X_P;
	double boundaryValueUZ_Y_M;
	double boundaryValueUZ_Y_P;
	double boundaryValueUZ_Z_M;
	double boundaryValueUZ_Z_P;

	int boundaryTypeP_X_M;
	int boundaryTypeP_X_P;
	int boundaryTypeP_Y_M;
	int boundaryTypeP_Y_P;
	int boundaryTypeP_Z_M;
	int boundaryTypeP_Z_P;
	double boundaryValueP_X_M;
	double boundaryValueP_X_P;
	double boundaryValueP_Y_M;
	double boundaryValueP_Y_P;
	double boundaryValueP_Z_M;
	double boundaryValueP_Z_P;

	int boundaryTypeT_X_M;
	int boundaryTypeT_X_P;
	int boundaryTypeT_Y_M;
	int boundaryTypeT_Y_P;
	int boundaryTypeT_Z_M;
	int boundaryTypeT_Z_P;
	double boundaryValueT_X_M;
	double boundaryValueT_X_P;
	double boundaryValueT_Y_M;
	double boundaryValueT_Y_P;
	double boundaryValueT_Z_M;
	double boundaryValueT_Z_P;

	Vec3r boundaryValuePoiseuilleCenter;

	bool GridGenerationMode;
	bool BenchMode;

private:

  void parse() {
		operatorname = read<string>("OperatorName", "v(^_^)v");

    origin = read<Vec3r>("origin", Vec3r(0, 0, 0));
    rootLength = read<double>("rootLength", 1.0);
    rootN = read<Vec3i>("rootGrid", Vec3i(1, 1, 1));
    seed   = read<Vec3r>("seed", Vec3r(origin.x, origin.y, origin.z));

		periodicX = read<bool>("PeriodicX", false);
		periodicY = read<bool>("PeriodicY", false);
		periodicZ = read<bool>("PeriodicZ", false);

    minLevel = read<int>("minLevel", 0);
    maxLevel = read<int>("maxLevel");

		dividerOx = read<double>("dividerOx", 0.5);
		dividerOy = read<double>("dividerOy", 0.5);
		dividerOz = read<double>("dividerOz", 0.5);
		dividerR  = read<double>("dividerR" , 0.25);
		dividerDR = read<double>("dividerDR", 0.0);
		dividerX0 = read<double>("dividerX0", 0.25);
		dividerY0 = read<double>("dividerY0", 0.25);
		dividerZ0 = read<double>("dividerZ0", 0.25);
		dividerX1 = read<double>("dividerX1", 0.75);
		dividerY1 = read<double>("dividerY1", 0.75);
		dividerZ1 = read<double>("dividerZ1", 0.75);
		dividerHollow = read<bool>("dividerHollow", false);

    treeType = read<string>("treeType");

    polylibConfig = read<string>("polylibConfig");
    polygonGroupList = read<std::vector<PolygonGroupSpec> >(
                         "polygonGroupList",
                         std::vector<PolygonGroupSpec>());
    boundingBoxList = read<std::vector<BoundingBoxSpec> >(
                        "boundingBoxList",
                        std::vector<BoundingBoxSpec>());
    sphericalBoxList = read<std::vector<BoundingBoxSpec> >(
                        "sphericalBoxList",
                        std::vector<BoundingBoxSpec>());

    ordering = read<string>("ordering");

    omegaU    = read<double>("omegaU", 1.0);
		countMaxU = read<int>("countMaxU", 1000);
		epsilonU  = read<double>("epsilonU", 1.0e-5);
		countPreConditionerU
							= read<int>("countPreConditionerU", 1);

    omegaP    = read<double>("omegaP", 1.0);
		countMaxP = read<int>("countMaxP", 1000);
		epsilonP  = read<double>("epsilonP", 1.0e-5);
		countPreConditionerP
							= read<int>("countPreConditionerP", 1);

		bHeat     = read<bool>("HeatProblem", true);
    omegaT    = read<double>("omegaT", 1.0);
		countMaxT = read<int>("countMaxT", 1000);
		epsilonT  = read<double>("epsilonT", 1.0e-5);
		countPreConditionerT
							= read<int>("countPreConditionerT", 1);

//    separate = read<bool>("separateVCUpdate", false);
		updateMethod = read<string>("updateMethod", "AtOnce");

    size = read<int>("size");
    vc = read<int>("vc");

		dt		= read<double>("dt", 1.0);
		rhof	= read<double>("rhof", 1.0);
		rhos	= read<double>("rhos", 1.0);
		cpf		= read<double>("cpf", 1.0);
		cps		= read<double>("cps", 1.0);
		kf		= read<double>("kf");
		ks		= read<double>("ks");
		mu		= read<double>("mu");
		gx    = read<double>("gx", 0.0);
		gy    = read<double>("gy", 0.0);
		gz    = read<double>("gz", 0.0);
		uxs0  = read<double>("uxs0", 0.0);
		uys0  = read<double>("uys0", 0.0);
		uzs0  = read<double>("uzs0", 0.0);

		StepStart      = read<int>("StepStart");
		StepEnd        = read<int>("StepEnd");

		StepPrintTime  = read<int>("StepPrintTime");
		StepPrintStats = read<int>("StepPrintStats");
		StepPrintData  = read<int>("StepPrintData");
		StepPrintBin   = read<int>("StepPrintBin");

		advection = read<string>("AdvectionScheme", "C2");

		cutoff         = read<bool>("cutoff", false);
		cutoff_epsilon = read<double>("cutoff_epsilon", 0.0);
    voxelization   = read<bool>("voxelization", false);
		symmetrization = read<bool>("symmetrization", false);
		holefilling    = read<bool>("holefilling", true);
		holefilling2   = read<bool>("holefilling2", true);
		masking        = read<bool>("masking", true);

		AccelDuration  = read<int>("AccelDuration", 0);
		STEP_ACCELDURATION
									 = read<int>("STEP_ACCELDURATION", 1);

		vtkWriteForce  = read<bool>("VtkWriteForce", false);

		BCMFileSave    = read<bool>("BCMFileSave", false);
		PLOT3DSave     = read<bool>("PLOT3DSave", false);

		boundaryTypeUX_X_M = read<int>("boundaryTypeUX_X_M");
		boundaryTypeUX_X_P = read<int>("boundaryTypeUX_X_P");
		boundaryTypeUX_Y_M = read<int>("boundaryTypeUX_Y_M");
		boundaryTypeUX_Y_P = read<int>("boundaryTypeUX_Y_P");
		boundaryTypeUX_Z_M = read<int>("boundaryTypeUX_Z_M");
		boundaryTypeUX_Z_P = read<int>("boundaryTypeUX_Z_P");
		boundaryValueUX_X_M = read<double>("boundaryValueUX_X_M");
		boundaryValueUX_X_P = read<double>("boundaryValueUX_X_P");
		boundaryValueUX_Y_M = read<double>("boundaryValueUX_Y_M");
		boundaryValueUX_Y_P = read<double>("boundaryValueUX_Y_P");
		boundaryValueUX_Z_M = read<double>("boundaryValueUX_Z_M");
		boundaryValueUX_Z_P = read<double>("boundaryValueUX_Z_P");

		boundaryTypeUY_X_M = read<int>("boundaryTypeUY_X_M");
		boundaryTypeUY_X_P = read<int>("boundaryTypeUY_X_P");
		boundaryTypeUY_Y_M = read<int>("boundaryTypeUY_Y_M");
		boundaryTypeUY_Y_P = read<int>("boundaryTypeUY_Y_P");
		boundaryTypeUY_Z_M = read<int>("boundaryTypeUY_Z_M");
		boundaryTypeUY_Z_P = read<int>("boundaryTypeUY_Z_P");
		boundaryValueUY_X_M = read<double>("boundaryValueUY_X_M");
		boundaryValueUY_X_P = read<double>("boundaryValueUY_X_P");
		boundaryValueUY_Y_M = read<double>("boundaryValueUY_Y_M");
		boundaryValueUY_Y_P = read<double>("boundaryValueUY_Y_P");
		boundaryValueUY_Z_M = read<double>("boundaryValueUY_Z_M");
		boundaryValueUY_Z_P = read<double>("boundaryValueUY_Z_P");

		boundaryTypeUZ_X_M = read<int>("boundaryTypeUZ_X_M");
		boundaryTypeUZ_X_P = read<int>("boundaryTypeUZ_X_P");
		boundaryTypeUZ_Y_M = read<int>("boundaryTypeUZ_Y_M");
		boundaryTypeUZ_Y_P = read<int>("boundaryTypeUZ_Y_P");
		boundaryTypeUZ_Z_M = read<int>("boundaryTypeUZ_Z_M");
		boundaryTypeUZ_Z_P = read<int>("boundaryTypeUZ_Z_P");
		boundaryValueUZ_X_M = read<double>("boundaryValueUZ_X_M");
		boundaryValueUZ_X_P = read<double>("boundaryValueUZ_X_P");
		boundaryValueUZ_Y_M = read<double>("boundaryValueUZ_Y_M");
		boundaryValueUZ_Y_P = read<double>("boundaryValueUZ_Y_P");
		boundaryValueUZ_Z_M = read<double>("boundaryValueUZ_Z_M");
		boundaryValueUZ_Z_P = read<double>("boundaryValueUZ_Z_P");

		boundaryTypeP_X_M = read<int>("boundaryTypeP_X_M");
		boundaryTypeP_X_P = read<int>("boundaryTypeP_X_P");
		boundaryTypeP_Y_M = read<int>("boundaryTypeP_Y_M");
		boundaryTypeP_Y_P = read<int>("boundaryTypeP_Y_P");
		boundaryTypeP_Z_M = read<int>("boundaryTypeP_Z_M");
		boundaryTypeP_Z_P = read<int>("boundaryTypeP_Z_P");
		boundaryValueP_X_M = read<double>("boundaryValueP_X_M");
		boundaryValueP_X_P = read<double>("boundaryValueP_X_P");
		boundaryValueP_Y_M = read<double>("boundaryValueP_Y_M");
		boundaryValueP_Y_P = read<double>("boundaryValueP_Y_P");
		boundaryValueP_Z_M = read<double>("boundaryValueP_Z_M");
		boundaryValueP_Z_P = read<double>("boundaryValueP_Z_P");

		boundaryTypeT_X_M = read<int>("boundaryTypeT_X_M");
		boundaryTypeT_X_P = read<int>("boundaryTypeT_X_P");
		boundaryTypeT_Y_M = read<int>("boundaryTypeT_Y_M");
		boundaryTypeT_Y_P = read<int>("boundaryTypeT_Y_P");
		boundaryTypeT_Z_M = read<int>("boundaryTypeT_Z_M");
		boundaryTypeT_Z_P = read<int>("boundaryTypeT_Z_P");
		boundaryValueT_X_M = read<double>("boundaryValueT_X_M");
		boundaryValueT_X_P = read<double>("boundaryValueT_X_P");
		boundaryValueT_Y_M = read<double>("boundaryValueT_Y_M");
		boundaryValueT_Y_P = read<double>("boundaryValueT_Y_P");
		boundaryValueT_Z_M = read<double>("boundaryValueT_Z_M");
		boundaryValueT_Z_P = read<double>("boundaryValueT_Z_P");

		boundaryValuePoiseuilleCenter = read<Vec3r>("boundaryValuePoiseuilleCenter", Vec3r(0.0, 0.0, 0.0));

		GridGenerationMode = read<bool>("GridGenerationMode", false);
		BenchMode          = read<bool>("BenchMode", false);
  }

  bool validate() {
    bool ret = true;
/*
    if (!(treeType == "flat" || treeType == "simple" || treeType == "polygon")) {
      std::cout << "error: 'treeType' must be 'flat' or 'simple' or 'polygon'." << std::endl;
      ret = false;
    }
*/
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
    std::cout << "  min level:          " << minLevel << std::endl;
    std::cout << "  max level:          " << maxLevel << std::endl;
    std::cout << "  tree type:          " << treeType << std::endl;
    std::cout << "  ordering:           " << ordering << std::endl;
    std::cout << "  block size:         " << size << std::endl;
    std::cout << "  vc width:           " << vc << std::endl;
  }

};



#endif // CONFIG_H
