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

#ifndef SOLVER_H
#define SOLVER_H

#include "BCMTools.h"
#include "BlockManager.h"

#include "real.h"
#include "LocalScalar3D.h"
#include "ILS3D.h"
#include "Config.h"

/*
#include "SiloWriter.h"

*/
#include "VtkWriter.h"
#include "Plot3DWriter.h"
#include "BCMFileSaver.h"

class Solver {
public:
	Solver();
	~Solver();

private:
	BlockManager& blockManager;
	int myrank;
	int vc;
	int diffLevel;
	int maxLevel;
	int minLevel;

	RootGrid* rootGrid;
	BCMOctree* tree;
	Partition* partition;
	Config conf;

	bool bRestart;

private:
	real dt;

	real rhof;
	real rhos;
	real cpf;
	real cps;
	real kf;
	real ks;
	real mu;

private:
	LocalScalar3D<real> *plsUX0;
	LocalScalar3D<real> *plsUX1;
	LocalScalar3D<real> *plsUXC;
	LocalScalar3D<real> *plsUXCP;
	LocalScalar3D<real> *plsUXD;
	LocalScalar3D<real> *plsUY0;
	LocalScalar3D<real> *plsUY1;
	LocalScalar3D<real> *plsUYC;
	LocalScalar3D<real> *plsUYCP;
	LocalScalar3D<real> *plsUYD;
	LocalScalar3D<real> *plsUZ0;
	LocalScalar3D<real> *plsUZ1;
	LocalScalar3D<real> *plsUZC;
	LocalScalar3D<real> *plsUZCP;
	LocalScalar3D<real> *plsUZD;

	LocalScalar3D<real> *plsVw;
	LocalScalar3D<real> *plsVe;
	LocalScalar3D<real> *plsVs;
	LocalScalar3D<real> *plsVn;
	LocalScalar3D<real> *plsVb;
	LocalScalar3D<real> *plsVt;

	LocalScalar3D<real> *plsP0;
	LocalScalar3D<real> *plsP1;
	LocalScalar3D<real> *plsLapP;

	LocalScalar3D<real> *plsT0;
	LocalScalar3D<real> *plsT1;
	LocalScalar3D<real> *plsTC;
	LocalScalar3D<real> *plsTCP;
	LocalScalar3D<real> *plsTD;

	LocalScalar3D<int> *plsPhaseId;
	LocalScalar3D<real> *plsCut0;
	LocalScalar3D<real> *plsCut1;
	LocalScalar3D<real> *plsCut2;
	LocalScalar3D<real> *plsCut3;
	LocalScalar3D<real> *plsCut4;
	LocalScalar3D<real> *plsCut5;
	LocalScalar3D<int> *plsCutId0;
	LocalScalar3D<int> *plsCutId1;
	LocalScalar3D<int> *plsCutId2;
	LocalScalar3D<int> *plsCutId3;
	LocalScalar3D<int> *plsCutId4;
	LocalScalar3D<int> *plsCutId5;

	LocalScalar3D<real> *plsAp;
	LocalScalar3D<real> *plsAw;
	LocalScalar3D<real> *plsAe;
	LocalScalar3D<real> *plsAs;
	LocalScalar3D<real> *plsAn;
	LocalScalar3D<real> *plsAb;
	LocalScalar3D<real> *plsAt;
	LocalScalar3D<real> *plsb;

	LocalScalar3D<real> *plsr;
	LocalScalar3D<real> *plsr0;
	LocalScalar3D<real> *plsp;
	LocalScalar3D<real> *plsp_;
	LocalScalar3D<real> *plsq_;
	LocalScalar3D<real> *plss;
	LocalScalar3D<real> *plss_;
	LocalScalar3D<real> *plst_;

	LocalScalar3D<real> *plsFspx;
	LocalScalar3D<real> *plsFspy;
	LocalScalar3D<real> *plsFspz;
	LocalScalar3D<real> *plsFsvx;
	LocalScalar3D<real> *plsFsvy;
	LocalScalar3D<real> *plsFsvz;

	LocalScalar3D<real> *plsM;
	LocalScalar3D<int> *plsMaskId;

private:
	ILS3D* pils;

	real omegaU;
	int countMaxU;
	real epsilonU;
	int countPreConditionerU;
	int countUX;
	real residualUX;
	int countUY;
	real residualUY;
	int countUZ;
	real residualUZ;

	real omegaP;
	int countMaxP;
	real epsilonP;
	int countPreConditionerP;
	int countP;
	real residualP;

	real omegaT;
	int countMaxT;
	real epsilonT;
	int countPreConditionerT;
	int countT;
	real residualT;

public:
	int Init(int argc, char** argv);
	int Loop();
	int Post();

private:
	int Update(int step);
	void UpdateUX();
	void UpdateUY();
	void UpdateUZ();
	void UpdateP();
	void UpdateU();
	void UpdateT();

	int Print(int step);
	double times[32];
	void PrintTime(int step);
	void PrintILS(int step);
	void PrintStats(int step);
	void PrintForce(int step);
	void PrintData(int step);
	void PrintLog(int level, const char* format, ...);

	void WritePolygon(std::ofstream& ofs, float *pv);
	void PrintCut();

private:
	void WriteDataInVTKFormat(
						const char* dataname,
						int step,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Config& conf) {
		VtkWriter writer;

		writer.writePUT<real>(
						this->plsP0->GetID(),
						this->plsUX0->GetID(),
						this->plsUY0->GetID(),
						this->plsUZ0->GetID(),
						this->plsLapP->GetID(),
						this->vc,
						string(dataname),
						step,
						difflevel,
						rootGrid,
						tree,
						partition,
						conf.origin,
						conf.rootLength);

		writer.writeVtkOverlappingAMR_PUT<real>(
						this->plsP0->GetID(),
						this->plsUX0->GetID(),
						this->plsUY0->GetID(),
						this->plsUZ0->GetID(),
						this->plsLapP->GetID(),
						this->vc,
						string(dataname),
						step,
						difflevel,
						rootGrid,
						tree,
						partition,
						conf.origin,
						conf.rootLength);

		if( conf.vtkWriteForce ) {
			writer.writePUT<real>(
							this->plsP0->GetID(),
							this->plsFspx->GetID(),
							this->plsFspy->GetID(),
							this->plsFspz->GetID(),
							this->plsLapP->GetID(),
							this->vc,
							string("fsp"),
							step,
							difflevel,
							rootGrid,
							tree,
							partition,
							conf.origin,
							conf.rootLength);

			writer.writePUT<real>(
							this->plsP0->GetID(),
							this->plsFsvx->GetID(),
							this->plsFsvy->GetID(),
							this->plsFsvz->GetID(),
							this->plsLapP->GetID(),
							this->vc,
							string("fsv"),
							step,
							difflevel,
							rootGrid,
							tree,
							partition,
							conf.origin,
							conf.rootLength);
		}

	}

	void WriteXYZInPlot3DFormat(
						const char* dataname,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Config& conf) {
		Plot3DWriter writer;
		writer.writeXYZ(
						this->plsPhaseId->GetID(),
						this->vc,
						string(dataname),
						difflevel,
						rootGrid,
						tree,
						partition,
						conf.origin,
						conf.rootLength);
	}

	void WriteDataInPlot3DFormat(
						const char* dataname,
						int step,
						int difflevel,
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Config& conf) {
		Plot3DWriter writer;
		writer.writeData<real>(
						this->plsP0->GetID(),
						this->plsUX0->GetID(),
						this->plsUY0->GetID(),
						this->plsUZ0->GetID(),
						this->vc,
						string(dataname),
						step,
						difflevel,
						rootGrid,
						tree,
						partition,
						conf.origin,
						conf.rootLength);
	}

/*
	void WriteDataInSiloFormat(
						const char* dataname,
						int step,
						Config& conf) {
		ostringstream ossFileName;
		ossFileName << "./SILO/";
		ossFileName << "data-";
		ossFileName << dataname;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".silo";


		SiloWriter writer(ossFileName.str(), "mesh", false);

		writer.writeDomain("block_mesh", "domain");
		writer.writeScalar<REAL_TYPE>(plsP0->GetID(), "p");
		writer.writeScalar<REAL_TYPE>(plsUX0->GetID(), "ux");
		writer.writeScalar<REAL_TYPE>(plsUY0->GetID(), "uy");
		writer.writeScalar<REAL_TYPE>(plsUZ0->GetID(), "uz");
	}
*/

	void Dump(const int step);
	void Load(const int step);

	void Dump2(const int step);
	void Load2(const int step);

	void Dump3(const int step);
	void Load3(const int step);


	BCMFileIO::BCMFileSaver *psaver;

	void BCMFileSaverInit(
						RootGrid* rootGrid,
						BCMOctree* tree,
						Partition* partition,
						Config& conf) {
		Vec3r origin = conf.origin;
		Vec3r region(conf.rootLength * conf.rootN.x,
									conf.rootLength * conf.rootN.y,
									conf.rootLength * conf.rootN.z);
		psaver = new BCMFileIO::BCMFileSaver(origin, region, tree, "BCM_OUT");

		BCMFileIO::IdxUnit unit;
		unit.length = std::string("m");
		unit.L0_scale = 1.0;
		unit.velocity = std::string("Dimensional");
		unit.V0_scale = 1.0;

		BCMFileIO::IdxStep timeStep(conf.StepStart, conf.StepEnd, conf.StepPrintData);

		int id_cellId = blockManager.setDataClass< Scalar3D<unsigned char> >(conf.vc);
		for (int id = 0; id < blockManager.getNumBlock(); ++id){
			BlockBase* block = blockManager.getBlock(id);
			::Vec3i sz = block->getSize();
			Scalar3D<unsigned char>* mesh = dynamic_cast< Scalar3D<unsigned char>* >(block->getDataClass(id_cellId));
			unsigned char* data = mesh->getData();
			Index3DS idx = mesh->getIndex();
			int* pPhaseId = plsPhaseId->GetBlockData(block);
			for(int z = -conf.vc; z < sz.z + conf.vc; z++){
				for(int y = -conf.vc; y < sz.y + conf.vc; y++){
					for(int x = -conf.vc; x < sz.x + conf.vc; x++){
						data[idx(x, y, z)] = pPhaseId[idx(x, y, z)];
					}
				}
			}
		}
		psaver->RegisterCellIDInformation(
						id_cellId,
						5,
						conf.vc,
						"CellID",
						"cid",
						"lb",
						"cid");

		int id_p[1] = {0};
		int id_u[3] = {0, 0, 0};
		id_p[0] = plsP0->GetID();
		id_u[0] = plsUX0->GetID();
		id_u[1] = plsUY0->GetID();
		id_u[2] = plsUZ0->GetID();

#ifdef _REAL_IS_DOUBLE_
		psaver->RegisterDataInformation(
					id_p,
					BCMFileIO::LB_SCALAR,
					BCMFileIO::LB_FLOAT64,
					conf.vc,
					"P64",
					"P64",
					"lb",
					timeStep,
					"P",
					true);
		psaver->RegisterDataInformation(
					id_u,
					BCMFileIO::LB_VECTOR3,
					BCMFileIO::LB_FLOAT64,
					conf.vc,
					"Vel64",
					"Vel64",
					"lb",
					timeStep,
					"VEL",
					true);
#else
		psaver->RegisterDataInformation(
					id_p,
					BCMFileIO::LB_SCALAR,
					BCMFileIO::LB_FLOAT32,
					conf.vc,
					"P32",
					"P32",
					"lb",
					timeStep,
					"P",
					true);
		psaver->RegisterDataInformation(
					id_u,
					BCMFileIO::LB_VECTOR3,
					BCMFileIO::LB_FLOAT32,
					conf.vc,
					"Vel32",
					"Vel32",
					"lb",
					timeStep,
					"VEL",
					true);
#endif

		psaver->SetUnit(unit);
		psaver->Save();
		psaver->SaveLeafBlock("CellID");
	}

	void BCMFileSaverPrint(int step) {
#ifdef _REAL_IS_DOUBLE_
		psaver->SaveLeafBlock("P64", step);
		psaver->SaveLeafBlock("Vel64", step);
#else
		psaver->SaveLeafBlock("P32", step);
		psaver->SaveLeafBlock("Vel32", step);
#endif
	}

private:
	double GetTime();
	void SetValues();

private:
	string GetSolverName() {
		ostringstream ossSolverName;
		ossSolverName << "FFV-BCM(alpha)";
		return ossSolverName.str();
	}
};


#endif // SOLVER_H
