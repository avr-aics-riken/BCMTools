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
/// @file  BCMFileLoader.cpp
/// @brief BCMファイルを読み込むクラス
/// 

#include "BCMFileLoader.h"

#include "BCMOctree.h"
#include "RootGrid.h"
#include "BlockManager.h"
#include "BoundaryConditionSetterBase.h"
#include "Block.h"
#include "BlockFactory.h"
#include "PartitionMapper.h"
#include "Scalar3D.h"
#include "Vector3D.h"
#include "Scalar3DUpdater.h"
#include "Vector3DUpdater.h"

#include "Vec3.h"

#include "TextParser.h"

#include <vector>
#include <string>

#include "BCMFileCommon.h"
#include "LeafBlockLoader.h"
#include "FileSystemUtil.h"
#include "ErrorUtil.h"

#include "type.h"

#define OCTREE_LOAD_ONLY_MASTER

namespace {

	int GetUniqueTag(){
		static const int tagBase = 1000;
		static int       conter = 0;

		return tagBase + conter++;
	}

	inline 
	int CompStr( const std::string& str1, const std::string& str2, bool ignorecase=true )
	{
		std::string lstr1 = str1;
		std::string lstr2 = str2;
		if( ignorecase ){
			std::transform(lstr1.begin(), lstr1.end(), lstr1.begin(), ::tolower);
			std::transform(lstr2.begin(), lstr2.end(), lstr2.begin(), ::tolower);
		}

		return lstr1.compare(lstr2);
	}

	int ReadVec3( TextParser* tp, const std::string& label, Vec3r& v)
	{
		using namespace std;

		if(!tp){ return TP_ERROR; }

		int err = TP_NO_ERROR;
		string valStr;
		vector<string> vec_valStr;
		if( (err = tp->getValue(label, valStr)) != TP_NO_ERROR){ return err; }
		tp->splitVector(valStr, vec_valStr);

		Vec3r ret_v;
		ret_v.x = tp->convertDouble(vec_valStr[0], &err);
		ret_v.y = tp->convertDouble(vec_valStr[1], &err);
		ret_v.z = tp->convertDouble(vec_valStr[2], &err);

		v = ret_v;

		return TP_NO_ERROR;
	}

	int ReadVec3( TextParser* tp, const std::string& label, Vec3i& v)
	{
		using namespace std;

		if(!tp){ return TP_ERROR; }

		int err = TP_NO_ERROR;
		string valStr;
		vector<string> vec_valStr;

		if( (err = tp->getValue(label, valStr)) != TP_NO_ERROR){ return err; }
		tp->splitVector(valStr, vec_valStr);

		Vec3i ret_v;
		ret_v.x = tp->convertInt(vec_valStr[0], &err);
		ret_v.y = tp->convertInt(vec_valStr[1], &err);
		ret_v.z = tp->convertInt(vec_valStr[2], &err);

		v.x = ret_v.x;
		v.y = ret_v.y;
		v.z = ret_v.z;

		return TP_NO_ERROR;
	}

} // namespace

namespace BCMFileIO {

	bool GetType(const std::string& typeStr, LB_DATA_TYPE &retType){
		bool status = false;
		const char *typeList[10] = {
			"Int8", "UInt8", "Int16", "UInt16", "Int32", "UInt32", "Int64", "UInt64", "Float32", "Float64"
		};
		for(int i = 0; i < 10; i++){
			if(CompStr(typeStr, typeList[i]) == 0){
				status = true;
				retType = (LB_DATA_TYPE)(i);
				break;
			}
		}
		return status;
	}

	bool LoadOctreeHeader(FILE *fp, OctHeader& header, bool& isNeedSwap)
	{
		isNeedSwap = false;
		fread(&header, sizeof(header), 1, fp);

		if( header.identifier != OCTREE_FILE_IDENTIFIER ){
			BSwap32(&header.identifier);

			if( header.identifier != OCTREE_FILE_IDENTIFIER ){
				return false;
			}

			isNeedSwap = true;
			for(int i = 0; i < 3; i++){
				BSwap64(&header.org[i]);
				BSwap64(&header.rgn[i]);
				BSwap32(&header.rootDims[i]);
			}
			BSwap32(&header.maxLevel);
			BSwap64(&header.numLeaf);
		}

		return true;
	}

	bool LoadOctreeFile(const std::string& filename, OctHeader& header, std::vector<Pedigree>& pedigrees)
	{
		using namespace std;
		FILE *fp = NULL;
		if( (fp = fopen(filename.c_str(), "rb")) == NULL ){
			LogE("open file error(%s) [%s:%d].\n", filename.c_str(), __FILE__, __LINE__);
			return false;
		}
		
		bool isNeedSwap = false;
		if( !LoadOctreeHeader(fp, header, isNeedSwap) ){
			LogE("Load Header error(%s) [%s:%d].\n", filename.c_str(), __FILE__, __LINE__);
			fclose(fp);
			return false;
		}

		pedigrees.resize(header.numLeaf);
		fread(&pedigrees[0], sizeof(Pedigree), header.numLeaf, fp);

		fclose(fp);

		if( isNeedSwap ){
			for(vector<Pedigree>::iterator it = pedigrees.begin(); it != pedigrees.end(); ++it){
				BSwap64(&(*it));
			}
		}
		return true;
	}


	BCMFileLoader::BCMFileLoader(const std::string& idxFilename, BoundaryConditionSetterBase* bcsetter)
	 : m_blockManager(BlockManager::getInstance()),
	   m_comm(m_blockManager.getCommunicator()),
	   m_octree(NULL),
	   m_pmapper(NULL)
	{
		
		std::string dir = GetDirectory(ConvertPath(idxFilename));
		
		std::string octreeFilename;
		if( reduceError( !LoadIndex(idxFilename, m_globalOrigin, m_globalRegion, octreeFilename, 
		                            m_leafBlockSize, m_idxProcList, m_idxBlockList ) ) ){
			LogE("load index file error (%s) [%s:%d].\n", idxFilename.c_str(), __FILE__, __LINE__);
			return;
		}

		if( reduceError( !LoadOctree(std::string(dir + octreeFilename), bcsetter) ) ){
			LogE("load octree file error (%s) [%s:%d].\n", std::string(dir + octreeFilename).c_str(), __FILE__, __LINE__);
			return;
		}
	}


	BCMFileLoader::~BCMFileLoader()
	{
		if(m_pmapper != NULL) delete m_pmapper;
		if(m_octree  != NULL) delete m_octree;
	}

	bool BCMFileLoader::LoadAdditionalIndex(const std::string& filepath)
	{
		Vec3r org;
		Vec3r rgn;
		std::string octname;
		Vec3i blockSize;
		std::vector<IdxProc>  idxProcList;
		std::vector<IdxBlock> idxBlockList;

		if( reduceError( !LoadIndex(filepath, org, rgn, octname, blockSize, idxProcList, idxBlockList) ) ){
			LogE("load index file error (%s) [%s:%d].\n", filepath.c_str(), __FILE__, __LINE__);
			return false;
		}
		// Check error
		if( idxProcList.size() != m_idxProcList.size() ){
			LogE("Process Info is invalid (%s) [%s:%d].\n", filepath.c_str(), __FILE__, __LINE__);
			return false;
		}
		// Check Octree
		{
			std::string dir = GetDirectory(ConvertPath(filepath));
			std::string octFilepath = dir + octname;
			FILE *fp = NULL;
			if( (fp = fopen(octFilepath.c_str(), "rb")) == NULL ){
				LogE("Octree is invalid (%s) [%s:%d].\n", filepath.c_str(), __FILE__, __LINE__);
				return false;
			}
			OctHeader hdr;
			bool isNeedSwap = false;
			if( !LoadOctreeHeader(fp, hdr, isNeedSwap) ){
				LogE("Octree is invalid (%s) [%s:%d].\n", filepath.c_str(), __FILE__, __LINE__);
				fclose(fp);
				return false;
			}
			fclose(fp);
			if( hdr.numLeaf != m_octree->getNumLeafNode() ){
				LogE("Octree is invalid (%s) [%s:%d].\n", filepath.c_str(), __FILE__, __LINE__);
				return false;
			}
		}

		
		for(std::vector<IdxBlock>::iterator it = idxBlockList.begin(); it != idxBlockList.end(); ++it){
			m_idxBlockList.push_back(*it);
		}

		return true;
	}

	bool BCMFileLoader::LoadIndex(const std::string& filename, Vec3r& globalOrigin, Vec3r& globalRegion, std::string& octreeFilename,
	                              Vec3i& blockSize, std::vector<IdxProc>& idxProcList, std::vector<IdxBlock>& idxBlockList)
	{
		using namespace std;
		TextParser *tp = new TextParser;

		if( tp->read(filename) != TP_NO_ERROR ) { 
			LogE("[%s:%d]\n", __FILE__, __LINE__); 
			delete tp;
			return false;
		}

		tp->changeNode("/BCMTree");
		// Read Octree Filename
		if( tp->getValue("TreeFile", octreeFilename) != TP_NO_ERROR ){
			LogE("[%s:%d]\n", __FILE__, __LINE__);
			delete tp;
			return false;
		}

		// Read Proc Filename
		string procFilename;
		if( tp->getValue("ProcFile", procFilename) != TP_NO_ERROR ){
			LogE("[%s:%d]\n", __FILE__, __LINE__); 
			delete tp;
			return false;
		}

		string dir = GetDirectory(ConvertPath(filename));
		string procFilepath = dir + procFilename;
		if( !LoadIndexProc(procFilepath, idxProcList) ){
			LogE("[%s:%d]\n", __FILE__, __LINE__);
			delete tp;
			return false;
		}

		// Read Domain
		tp->changeNode("/Domain");
		{
			bool hasOrigin = false;
			bool hasRegion = false;
			vector<string> lbls;
			tp->getLabels(lbls);
			for(vector<string>::iterator it = lbls.begin(); it != lbls.end(); ++it){
				string valStr;
				tp->getValue(*it, valStr);

				if( CompStr(*it, "GlobalOrigin") == 0 ){
					if( ReadVec3(tp, *it, globalOrigin) == TP_NO_ERROR ){ hasOrigin = true; }
					continue;
				}
				if( CompStr(*it, "GlobalRegion") == 0 ){
					if( ReadVec3(tp, *it, globalRegion) == TP_NO_ERROR ){ hasRegion = true; }
					continue;
				}
			}
			if( !hasOrigin || !hasRegion ){
				LogE("[%s:%d]\n", __FILE__, __LINE__);
				delete tp;
				return false;
			}
		}
		
		tp->changeNode("/LeafBlock" );
		{
			vector<string> nodes;
			tp->getNodes(nodes, 1);
			for(vector<string>::iterator nit = nodes.begin(); nit != nodes.end(); ++nit){
				// Load Data
				if( CompStr(nit->substr(0, 4), "Data") == 0 ) {
					tp->changeNode(*nit);
					IdxBlock ib;
					if( LoadIndexData(tp, &ib) ){
						ib.rootDir = FixDirectoryPath(dir);
						idxBlockList.push_back(ib);
					}
					tp->changeNode("../");
					continue;
				}
				// Load CellID
				if( CompStr(nit->substr(0, 6), "CellID") == 0 ){
					tp->changeNode(*nit);
					IdxBlock ib;
					if( LoadIndexCellID(tp, &ib) ){
						ib.rootDir = FixDirectoryPath(dir);
						idxBlockList.push_back(ib);
					}
					tp->changeNode("../");
					continue;
				}
				// Load Unit
				if( CompStr(*nit, "Unit") == 0 ){
					tp->changeNode(*nit);
					vector<string> lbls;
					tp->getLabels(lbls);
					for(vector<string>::iterator it = lbls.begin(); it != lbls.end(); ++it){
						string valStr;
						tp->getValue(*it, valStr);

						if( CompStr(*it, "Length")   == 0 ){ m_unit.length   = valStr; continue; }
						if( CompStr(*it, "Velocity") == 0 ){ m_unit.velocity = valStr; continue; }
						if( CompStr(*it, "L0") == 0 ){ m_unit.L0_scale = atof(valStr.c_str()); continue; }
						if( CompStr(*it, "V0") == 0 ){ m_unit.V0_scale = atof(valStr.c_str()); continue; }
					}
					tp->changeNode("../");
				}
			}

			if( idxBlockList.size() == 0 ){
				LogE("[%s:%d]\n", __FILE__, __LINE__);
				delete tp;
				return false;
			}

			if( ReadVec3(tp, "size", blockSize) != TP_NO_ERROR ){
				LogE("[%s:%d]\n", __FILE__, __LINE__);
				delete tp;
				return false;
			}
		}
		delete tp;
		return true;
	}

	bool BCMFileLoader::LoadIndexProc(const std::string& filepath, std::vector<IdxProc>& procList)
	{
		using namespace std;
		TextParser *tp = new TextParser;

		if( tp->read(filepath) != TP_NO_ERROR ) {
			LogE("[%s:%d]\n", __FILE__, __LINE__);
			delete tp;
			return false;
		}

		int numProcs = 0;
		{
			tp->changeNode("/MPI");
			string valStr;
			if( tp->getValue("NumberOfRank", valStr) != TP_NO_ERROR ){
				LogE("[%s:%d]\n", __FILE__, __LINE__);
				delete tp;
				return false;
			}
			numProcs = atoi(valStr.c_str());
		}

		tp->changeNode("/process");		
		vector<string> nodes;
		tp->getNodes(nodes, 1);
		for(vector<string>::iterator nit = nodes.begin(); nit != nodes.end(); ++nit){
			if( CompStr(nit->substr(0, 4), "Rank") == 0 ){
				tp->changeNode(*nit);

				vector<string> lbls;
				tp->getLabels(lbls);

				IdxProc proc;
				for(vector<string>::iterator it = lbls.begin(); it != lbls.end(); ++it){
					string valStr;
					tp->getValue(*it, valStr);
					if( CompStr(*it, "ID") == 0 ){
						proc.rank = atoi(valStr.c_str());
						continue;
					}
					if( CompStr(*it, "Hostname") == 0 ){
						proc.hostname = valStr;
						continue;
					}

					if( CompStr(*it, "BlockRange") == 0 ){
						double range[3];
						tp->splitRange(valStr, &range[0], &range[1], &range[2]);
						proc.rangeMin = static_cast<unsigned int>(range[0]);
						proc.rangeMax = static_cast<unsigned int>(range[1]);
						continue;
					}
				}
				procList.push_back(proc);
				
				tp->changeNode("../");
			}
		}
		
		if(procList.size() != numProcs){
			LogE("[%s:%d]\n", __FILE__, __LINE__);
			delete tp;
			procList.clear();
			return false;
		}

		delete tp;
		return true;
	}
	
	bool BCMFileLoader::LoadIndexStep(TextParser *tp, IdxStep* step)
	{
		using namespace std;
		vector<string> lbls;
		tp->getLabels(lbls);
		
		bool hasBase = false;
		
		for(vector<string>::iterator it = lbls.begin(); it != lbls.end(); ++it){
			string valStr;
			tp->getValue(*it, valStr);

			if( CompStr(*it, "base") == 0 ){
				double range[3];
				if( tp->splitRange(valStr, &range[0], &range[1], &range[2]) != TP_NO_ERROR ){
					return false;
				}
				step->SetRange((unsigned int)range[0], (unsigned int)range[1], (unsigned int)range[2]);
				hasBase = true;
			}

			if( CompStr(*it, "add") == 0 ){
				vector<double> list;
				if( tp->splitList(valStr, list) != TP_NO_ERROR ){
					return false;
				}

				for(vector<double>::iterator lit = list.begin(); lit != list.end(); ++lit){
					step->AddStep((unsigned int)(*lit));
				}
			}

			if( CompStr(*it, "sub") == 0 ){
				vector<double> list;
				if( tp->splitList(valStr, list) != TP_NO_ERROR ){
					return false;
				}

				for(vector<double>::iterator lit = list.begin(); lit != list.end(); ++lit){
					step->SubStep((unsigned int)(*lit));
				}

			}
		}

		if( !hasBase ){ return false; }

		return true;
	}

	bool BCMFileLoader::LoadIndexData(TextParser *tp, IdxBlock* ib)
	{
		using namespace std;
		vector<string> lbls;
		tp->getLabels(lbls);

		bool hasName          = false;
		bool hasType          = false;
		bool hasNumComponents = false;
		bool hasVC            = false;
		bool hasPrefix        = false;
		bool hasExtension     = false;

		for(vector<string>::iterator it = lbls.begin(); it != lbls.end(); ++it){
			string valStr;
			tp->getValue(*it, valStr);

			if( CompStr(*it, "name") == 0 ){
				ib->name = valStr;
				hasName = true;
				continue;
			}

			if( CompStr(*it, "type") == 0 ){
				if( !GetType(valStr, ib->dataType) ){
					LogE("value (%s) of keyword [type] is invalid", valStr.c_str());
					return false;
				}
				hasType = true;
				continue;
			}

			if( CompStr(*it, "NumberOfComponents") == 0 ){
				int val = atoi(valStr.c_str());
				if( val != 1 && val != 3 && val != 4 && val != 6 && val != 9){
					LogE("value (%d) of keyword [NumberOfComponents] is invalid", val);
					return false;
				}
				ib->kind = static_cast<LB_KIND>(val);
				hasNumComponents = true;
				continue;
			}

			if( CompStr(*it, "VirtualCellSize") == 0 ){
				ib->vc = atoi(valStr.c_str());
				hasVC = true;
				continue;
			}

			if( CompStr(*it, "Prefix") == 0 ){
				ib->prefix = valStr;
				hasPrefix = true;
				continue;
			}

			if( CompStr(*it, "Extension") == 0 ){
				ib->extension = valStr;
				hasExtension = true;
				continue;
			}

			if( CompStr(*it, "DirectoryPath") == 0){
				ib->dataDir = FixDirectoryPath(valStr);
				continue;
			}

			if( CompStr(*it, "StepSubDirectory") == 0){
				if( CompStr(valStr, "true") == 0 ){
					ib->isStepSubDir = true;
				}else{
					ib->isStepSubDir = false;
				}
			}
		}

		if( !hasName || !hasType || !hasNumComponents || !hasVC || !hasPrefix || !hasExtension ){
			LogE("Load Index File Error [%d:%s].\n", __FILE__, __LINE__);
			return false;
		}

		ib->dataClassID.clear();
		int bitWidthTable[10] = {
			8, 8, 16, 16, 32, 32, 64, 64, 32, 64
		};
		ib->bitWidth = bitWidthTable[(int)(ib->dataType)];

		if( tp->changeNode("Step") != TP_NO_ERROR ){
			return false;
		}else{
			bool err = !LoadIndexStep(tp, &ib->step);
			tp->changeNode("../");
			if( err ){ return false; }
		}

		return true;
	}

	bool BCMFileLoader::LoadIndexCellID(TextParser *tp, IdxBlock* ib)
	{
		using namespace std;
		vector<string> lbls;
		tp->getLabels(lbls);

		bool hasName       = false;
		bool hasBitWidth   = false;
		bool hasPrefix     = false;
		bool hasExtension  = false;
		bool hasGatherMode = false;

		for(vector<string>::iterator it = lbls.begin(); it != lbls.end(); ++it){
			string valStr;
			tp->getValue(*it, valStr);

			if( CompStr(*it, "name") == 0){
				ib->name = valStr;
				hasName = true;
				continue;
			}

			if( CompStr(*it, "BitWidth") == 0 ){
				ib->bitWidth = atoi(valStr.c_str());
				ib->kind     = LB_CELLID;
				ib->dataType = LB_UINT8;
				hasBitWidth  = true;
				continue;
			}

			if( CompStr(*it, "VirtualCellSize") == 0){
				ib->vc = atoi(valStr.c_str());
				continue;
			}

			if( CompStr(*it, "Prefix") == 0 ){
				ib->prefix = valStr;
				hasPrefix  = true;
				continue;
			}

			if( CompStr(*it, "Extension") == 0 ){
				ib->extension = valStr;
				hasExtension = true;
				continue;
			}

			if( CompStr(*it, "DirectoryPath") == 0){
				ib->dataDir = FixDirectoryPath(valStr);
				continue;
			}

			if( CompStr(*it, "GatherMode") == 0 ){
				if( CompStr(valStr, "distributed") == 0){
					ib->isGather  = false;
					hasGatherMode = true;
				}else if( CompStr(valStr, "gathered") == 0){
					ib->isGather  = true;
					hasGatherMode = true;
				}else{
					LogE("value (%s) of keyword [GatherMode] is invalid.\n", valStr.c_str());
					return false;
				}
				continue;
			}
		}

		if(!hasName || !hasBitWidth || !hasPrefix || !hasExtension || !hasGatherMode ){
			LogE("Load Index File Error [%d:%s].\n", __FILE__, __LINE__);
			return false;
		}

		return true;
	}


	bool BCMFileLoader::LoadOctree(const std::string& filename, BoundaryConditionSetterBase* bcsetter)
	{
		using namespace std;

		int myRank   = m_comm.Get_rank();
		int numProcs = m_comm.Get_size();

		bool load_octree_only_master = false;
#ifdef OCTREE_LOAD_ONLY_MASTER
		load_octree_only_master = true;
#endif
		
		OctHeader header;
		vector<Pedigree> pedigrees;

		if( load_octree_only_master )
		{
			unsigned char loadError = 0;
			if( myRank == 0 ){
				if( !LoadOctreeFile(filename, header, pedigrees) ) loadError = 1; 
			}

			m_comm.Bcast(&loadError, 1, MPI::CHAR, 0);
			if( loadError == 1) return false;
			m_comm.Bcast(&header, sizeof(OctHeader), MPI::CHAR, 0);

			if( myRank != 0 ){
				pedigrees.resize(header.numLeaf);
			}
			m_comm.Bcast(&pedigrees[0], sizeof(Pedigree) * pedigrees.size(), MPI::CHAR, 0);
		}
		else
		{
			if( reduceError( !LoadOctreeFile(filename, header, pedigrees)) ) return false;
		}
		

		Vec3r rootRegion( header.rgn[0] / static_cast<double>(header.rootDims[0]),
		                  header.rgn[1] / static_cast<double>(header.rootDims[1]),
					  	  header.rgn[2] / static_cast<double>(header.rootDims[2]) );

		if( fabs(rootRegion.x - rootRegion.y) >= 1.0e-10 || fabs(rootRegion.x - rootRegion.z) >= 1.0e-10 ) {
			LogE("%lf %lf %lf\n", rootRegion.x, rootRegion.y, rootRegion.z);
			LogE("RootGrid is not regular hexahedron. [%s:%d]\n", __FILE__, __LINE__);
			return false;
		}
		
		RootGrid *rootGrid = new RootGrid(header.rootDims[0], header.rootDims[1], header.rootDims[2]);
		m_octree  = new BCMOctree(rootGrid, pedigrees);
		m_pmapper = new PartitionMapper(m_idxProcList.size(), numProcs, header.numLeaf);
		
		// Make and register Block
		Partition part(numProcs, header.numLeaf);
		BlockFactory factory(m_octree, &part, bcsetter, Vec3r(header.org), rootRegion.x, m_leafBlockSize);

		vector<Node*>& leafNodeArray = m_octree->getLeafNodeArray();
		for(int id = part.getStart(myRank); id < part.getEnd(myRank); id++){
			Node* node   = leafNodeArray[id];
			Block* block = factory.makeBlock(node);
			m_blockManager.registerBlock(block);
		}

		m_blockManager.endRegisterBlock();

		m_globalOrigin = Vec3r(header.org);
		m_globalRegion = Vec3r(header.rgn);

		return true;
	}

	bool BCMFileLoader::CreateLeafBlock(int *dataClassID, const std::string& name, const unsigned int vc, const bool separateVCUpdate)
	{
		using namespace std;
		bool err = false;

		IdxBlock* ib = findIdxBlock(m_idxBlockList, name);

		if( ib == NULL ){
			LogE("No such name as \"%s\" in loaded index.[%s:%d]\n", name.c_str(), __FILE__, __LINE__);
			err = true;
		}
		if( reduceError(err) ) { return false; }

		if(ib->kind == LB_CELLID)
		{
			if( ib->dataClassID.size() == 0 ){
				ib->dataClassID.resize(1);
				ib->dataClassID[0] = m_blockManager.setDataClass< Scalar3D<unsigned char> >(vc); // TODO Cell Updater
			}
			dataClassID[0] = ib->dataClassID[0];
		}
		else
		{
			if( ib->dataClassID.size() == 0 ){
				ib->dataClassID.resize(static_cast<size_t>(ib->kind));
				for(int i = 0; i < static_cast<int>(ib->kind); i++){
					if     (ib->dataType == LB_FLOAT32){ ib->dataClassID[i] = m_blockManager.setDataClass< Scalar3D<f32>, Scalar3DUpdater<f32> >(vc); }
					else if(ib->dataType == LB_FLOAT64){ ib->dataClassID[i] = m_blockManager.setDataClass< Scalar3D<f64>, Scalar3DUpdater<f64> >(vc); }
					#if 0 // TODO
					else if(ib->dataType == LB_INT8   ){ ib->dataClassID[i] = m_blockManager.setDataClass< Scalar3D< s8>, Scalar3DUpdater< s8> >(vc); }
					else if(ib->dataType == LB_UINT8  ){ ib->dataClassID[i] = m_blockManager.setDataClass< Scalar3D< u8>, Scalar3DUpdater< u8> >(vc); }
					else if(ib->dataType == LB_INT16  ){ ib->dataClassID[i] = m_blockManager.setDataClass< Scalar3D<s16>, Scalar3DUpdater<s16> >(vc); }
					else if(ib->dataType == LB_UINT16 ){ ib->dataClassID[i] = m_blockManager.setDataClass< Scalar3D<u16>, Scalar3DUpdater<u16> >(vc); }
					else if(ib->dataType == LB_INT32  ){ ib->dataClassID[i] = m_blockManager.setDataClass< Scalar3D<s32>, Scalar3DUpdater<s32> >(vc); }
					else if(ib->dataType == LB_UINT32 ){ ib->dataClassID[i] = m_blockManager.setDataClass< Scalar3D<u32>, Scalar3DUpdater<u32> >(vc); }
					else if(ib->dataType == LB_INT64  ){ ib->dataClassID[i] = m_blockManager.setDataClass< Scalar3D<s64>, Scalar3DUpdater<s64> >(vc); }
					else if(ib->dataType == LB_UINT64 ){ ib->dataClassID[i] = m_blockManager.setDataClass< Scalar3D<u64>, Scalar3DUpdater<u64> >(vc); }
					#endif
					m_blockManager.prepareForVCUpdate(ib->dataClassID[i], GetUniqueTag(), separateVCUpdate);
					ib->separateVCUpdate = separateVCUpdate;
				}
			}
			for(int i = 0; i < static_cast<int>(ib->kind); i++){
				dataClassID[i] = ib->dataClassID[i];
			}
		}
/*
		int dataClassID = ib->dataClassID;
		// nameに対応するブロックが作成されていない場合、作成 (初期化)
		if(dataClassID < 0){
			if(ib->kind == LB_CELLID)
			{
				// CellID向けブロック
				//dataClassID = m_blockManager.setDataClass< Scalar3D<unsigned char>, Scalar3DUpdater<unsigned char> >(vc);
				dataClassID = m_blockManager.setDataClass< Scalar3D<unsigned char> >(vc); // TODO Cell Updater
			}
			else
			{
				// Scalar向けブロック(Type別)
				if     ( ib->dataType == LB_FLOAT32) { dataClassID = m_blockManager.setDataClass< Scalar3D<f32>, Scalar3DUpdater<f32> >(vc); }
				else if( ib->dataType == LB_FLOAT64) { dataClassID = m_blockManager.setDataClass< Scalar3D<f64>, Scalar3DUpdater<f64> >(vc); }
				//else if( ib->dataType == LB_INT8   ) { dataClassID = m_blockManager.setDataClass< Scalar3D< s8>, Scalar3DUpdater< s8> >(vc); }
				//else if( ib->dataType == LB_UINT8  ) { dataClassID = m_blockManager.setDataClass< Scalar3D< u8>, Scalar3DUpdater< u8> >(vc); }
				//else if( ib->dataType == LB_INT16  ) { dataClassID = m_blockManager.setDataClass< Scalar3D<s16>, Scalar3DUpdater<s16> >(vc); }
				//else if( ib->dataType == LB_UINT16 ) { dataClassID = m_blockManager.setDataClass< Scalar3D<u16>, Scalar3DUpdater<u16> >(vc); }
				//else if( ib->dataType == LB_INT32  ) { dataClassID = m_blockManager.setDataClass< Scalar3D<s32>, Scalar3DUpdater<s32> >(vc); }
				//else if( ib->dataType == LB_UINT32 ) { dataClassID = m_blockManager.setDataClass< Scalar3D<u32>, Scalar3DUpdater<u32> >(vc); }
				//else if( ib->dataType == LB_INT64  ) { dataClassID = m_blockManager.setDataClass< Scalar3D<s64>, Scalar3DUpdater<s64> >(vc); }
				//else if( ib->dataType == LB_UINT64 ) { dataClassID = m_blockManager.setDataClass< Scalar3D<u64>, Scalar3DUpdater<u64> >(vc); }
				else{
					LogE("invalid DataType (%d)[%s:%d]\n", ib->dataType, __FILE__, __LINE__);
					err = true;
				}
				m_blockManager.prepareForVCUpdate(dataClassID, GetUniqueTag(), separateVCUpdate);
				ib->separateVCUpdate = separateVCUpdate;
			}

			ib->dataClassID = dataClassID;
		}
*/
		if( reduceError(err) ){ return false; }

		return true;
	}


	bool BCMFileLoader::LoadLeafBlock(int *dataClassID, const std::string& name, const unsigned int vc, 
	                                  const unsigned int step, const bool separateVCUpdate)
	{
		using namespace std;

		const int myRank   = m_comm.Get_rank();

		bool err = false;

		IdxBlock* ib = findIdxBlock(m_idxBlockList, name);
		
		if( ib == NULL ){
			LogE("No such name as \"%s\" in loaded index.[%s:%d]\n", name.c_str(), __FILE__, __LINE__);
			err = true;
		}
		if( reduceError(err) ){ return false; }

		if( reduceError( !CreateLeafBlock(dataClassID, name, vc, separateVCUpdate) ) ){ return false; }

		// CellIDデータロード
		if(ib->kind == LB_CELLID)
		{
			LBHeader header;
			std::vector<CellIDCapsule> cidCapsules;

			if( ib->isGather ){
				// ファイルからCellIDを読み込む (GatherMode = Gathered)
				err = !Load_LeafBlock_CellID_Gather(ib->rootDir + ib->dataDir, ib, m_comm, m_pmapper, header, cidCapsules);
			}else{
				// ファイルからCellIDを読み込む (GatherMode = Distributed)
				err = !Load_LeafBlock_CellID(ib->rootDir + ib->dataDir, ib, m_comm, m_pmapper, header, cidCapsules);
			}
			// エラーチェック
			if( reduceError(err) ){
				for(vector<CellIDCapsule>::iterator it = cidCapsules.begin(); it != cidCapsules.end(); ++it){
					delete [] it->data;
				}
				return false;
			}
			
			// MxNマッピング情報取得
			vector<PartitionMapper::FDIDList> fdidlists;
			m_pmapper->GetFDIDLists(myRank, fdidlists);

			int did = 0;
			int fid = 0;
			Vec3i bsz = Vec3i(m_leafBlockSize);
			// cidCapsulesから逐次データを展開しBlockManager配下のBlockへデータをコピー
			for(vector<PartitionMapper::FDIDList>::iterator file = fdidlists.begin(); file != fdidlists.end(); ++file){

				// rleおよびbitVoxelの圧縮を展開
				unsigned char* voxels = DecompCellIDData( header, cidCapsules[fid] );
				
				// 展開したvoxelsからブロックごとにデータをコピー
				for( vector<int>::iterator fdid = file->FDIDs.begin(); fdid != file->FDIDs.end(); ++fdid){
					Vec3i ibsz( bsz.x + vc*2,     bsz.y + vc*2,     bsz.z + vc*2    );   // 内部ブロックサイズ (仮想セル込み)
					Vec3i fbsz( bsz.x + ib->vc*2, bsz.y + ib->vc*2, bsz.z + ib->vc*2);   // ファイルブロックサイズ (仮想セル込み)
					unsigned char* block = new unsigned char[ibsz.x * ibsz.y * ibsz.z];  // データコピー用一時バッファを準備
					memset(block, 0, sizeof(unsigned char) * ibsz.x * ibsz.y * ibsz.z);  // 一時バッファの0クリア
					unsigned char *pv = &voxels[ (fbsz.x * fbsz.y * fbsz.z) * (*fdid) ];
					
					// ファイルから読み込んだCellIDを一時バッファにコピー (仮想セルサイズの不一致への対応)
					if( vc > ib->vc ){
						unsigned int vcd = vc - ib->vc;
						for(int z = 0; z < fbsz.z; z++){
							for(int y = 0; y < fbsz.y; y++){
								size_t ibloc = 0 + vcd + ( (y + vcd) + (z + vcd) * ibsz.y ) * ibsz.x;
								size_t fbloc = 0 +     + (  y        +  z        * fbsz.y ) * fbsz.x;
								memcpy(&block[ibloc], &pv[fbloc], sizeof(unsigned char) * fbsz.x );
							}
						}
					}else{
						unsigned int vcd = ib->vc - vc;
						for(int z = 0; z < ibsz.z; z++){
							for(int y = 0; y < ibsz.y; y++){
								size_t ibloc = 0 +     + (  y        +  z        * ibsz.y ) * ibsz.x;
								size_t fbloc = 0 + vcd + ( (y + vcd) + (z + vcd) * fbsz.y ) * fbsz.x;
								memcpy(&block[ibloc], &pv[fbloc], sizeof(unsigned char) * ibsz.x );
							}
						}
					}
					// ブロックマネージャ配下のブロックに値をコピー
					CopyBufferToScalar3D(m_blockManager, dataClassID[0], did, vc, block);
					did++;
					delete [] block;
				}

				delete [] voxels;
				fid++;
			}
		}
		else
		{
			// ファイルからデータを読み込み、ブロックマネージャ配下のブロックに値をコピー
			if( reduceError(!Load_LeafBlock_Data( m_comm, ib, m_blockManager, m_pmapper, vc, step)) ){
				return false;
			}
			// 現在の仮想セルサイズがファイルに記載されている仮想セルサイズよりも大きい場合、仮想セルの同期を行う
			if( vc > ib->vc ){
				for(int i = 0; i < static_cast<int>(ib->kind); i++){
					if( ib->separateVCUpdate ){
						m_blockManager.updateVC_X(dataClassID[i]);
						m_blockManager.updateVC_Y(dataClassID[i]);
						m_blockManager.updateVC_Z(dataClassID[i]);
					}
					else{
						m_blockManager.updateVC(dataClassID[i]);
					}
				}
			}
		}

		return true;
	}

	const IdxStep* BCMFileLoader::GetStep(const std::string& name ) const
	{
		const IdxBlock* ib = findIdxBlock(m_idxBlockList, name);

		if( ib == NULL ) return NULL;

		return &(ib->step);
	}

} // BCMFileIO

