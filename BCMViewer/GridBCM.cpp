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


#include "GridBCM.h"
#include "BCMFileCommon.h"
#include "BitVoxel.h"
#include "RLE.h"
#include "ErrorUtil.h"
#include "FileSystemUtil.h"
#include "LeafBlockLoader.h"
#include "RootGrid.h"
#include "PartitionMapper.h"

#include "SceneGraph.h"
#include "Render.h"

#include "TextureObject.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include "TextParser.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#endif

//#define TEXTURE_SIZE_LIMIT 8192
//#define TEXTURE_SIZE_LIMIT 1024

template<typename T> inline T MAX_(T A, T B)
{
	return ( A < B ? B : A);
}

template<typename T> inline T MIN_(T A, T B)
{
	return ( A > B ? B : A);
}

void SaveImage(const int *img, const int width, const int height)
{
	FILE *fp = fopen("img.ppm", "w");

	fprintf(fp, "P3\n%d %d\n255\n", width, height);
	for(int h = 0; h < height; h++){
		for(int w = 0; w < width; w++){
			int color = img[w + h * width];
			unsigned char r = (unsigned char)(color >> 16);
			unsigned char g = (unsigned char)(color >> 8 );
			unsigned char b = (unsigned char)(color >> 0 );
			fprintf(fp, "%d %d %d\n", r, g, b);
		}
	}

	fclose(fp);
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
	ret_v.x = (REAL_TYPE)tp->convertDouble(vec_valStr[0], &err);
	ret_v.y = (REAL_TYPE)tp->convertDouble(vec_valStr[1], &err);
	ret_v.z = (REAL_TYPE)tp->convertDouble(vec_valStr[2], &err);

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

namespace BCMFileIO {
	bool LoadIndexProc(const std::string& filepath, std::vector<IdxProc>& procList)
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
						REAL_TYPE range[3];
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

	bool LoadIndexCellID(TextParser *tp, IdxBlock* ib)
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

		ib->vc = 0;

		return true;
	}

	bool LoadIndex(const std::string& filename, const std::string& targetDir,
	               Vec3r& globalOrigin, Vec3r& globalRegion, std::string& octreeFilename,
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

		string procFilepath = targetDir + procFilename;
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
				/*
				if( CompStr(nit->substr(0, 4), "Data") == 0 ) {
					tp->changeNode(*nit);
					IdxBlock ib;
					if( LoadIndexData(tp, &ib) ){
						idxBlockList.push_back(ib);
					}
					tp->changeNode("../");
				}
				*/
				// Load CellID
				if( CompStr(nit->substr(0, 6), "CellID") == 0 ){
					tp->changeNode(*nit);
					IdxBlock ib;
					if( LoadIndexCellID(tp, &ib) ){
						ib.rootDir = targetDir;
						idxBlockList.push_back(ib);
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

			// TODO Load Unit
		}

		return true;
	}
} // namespace BCMFileIO


namespace {
	void BoxSetter(const Vec3r& org, const Vec3r& rgn,
	               SG::VertexLineFormat* vertex, SG::IndexFormat* index, size_t idxOffset)
	{
		SG::VertexLineFormat vrt[8] = {
			SG::VertexLineFormat(org.x,         org.y,         org.z        ),
			SG::VertexLineFormat(org.x + rgn.x, org.y,         org.z        ),
			SG::VertexLineFormat(org.x + rgn.x, org.y + rgn.y, org.z        ),
			SG::VertexLineFormat(org.x,         org.y + rgn.y, org.z        ),
			SG::VertexLineFormat(org.x,         org.y,         org.z + rgn.z),
			SG::VertexLineFormat(org.x + rgn.x, org.y,         org.z + rgn.z),
			SG::VertexLineFormat(org.x + rgn.x, org.y + rgn.y, org.z + rgn.z),
			SG::VertexLineFormat(org.x,         org.y + rgn.y, org.z + rgn.z)
		};

		GLuint idx[12 * 2] = {
			0, 1, 1, 2, 2, 3, 3, 0,
			4, 5, 5, 6, 6, 7, 7, 4,
			0, 4, 3, 7, 1, 5, 2, 6
		};

		for(int i = 0; i < 24; i++){
			idx[i] += idxOffset;
		}

		memcpy(vertex, vrt, sizeof(SG::VertexLineFormat)    * 8);
		memcpy(index, idx, sizeof(SG::IndexFormat) * 12* 2);
	}


	void CellSetter(const Vec3r& org, const Vec3r& rgn, const GridBCM::AXIS axis, const int divU, const int divV,
					const int texIdx_u, const int texIdx_v, const int texSize_u, int texSize_v,
	                SG::VertexFaceFormat* vertex, SG::IndexFormat* index, size_t idxOffset)
	{
		float tp[4] = {
			(float)( texIdx_u        ) / (float) texSize_u,
			(float)( texIdx_u + divU ) / (float) texSize_u,
			(float)( texIdx_v        ) / (float) texSize_v,
			(float)( texIdx_v + divV ) / (float) texSize_v
		}; // TexCoord Table (Leaft, Right, Bottom, Top)

		if( axis == GridBCM::AXIS_X ){
			const SG::VertexFaceFormat vrt[4] = {
				SG::VertexFaceFormat( org[0],          org[1],          org[2],          1.0, 0.0, 0.0, tp[0], tp[2]),
				SG::VertexFaceFormat( org[0],          org[1] + rgn[1], org[2],          1.0, 0.0, 0.0, tp[1], tp[2]),
				SG::VertexFaceFormat( org[0],          org[1] + rgn[1], org[2] + rgn[2], 1.0, 0.0, 0.0, tp[1], tp[3]),
				SG::VertexFaceFormat( org[0],          org[1],          org[2] + rgn[2], 1.0, 0.0, 0.0, tp[0], tp[3])
			};
			memcpy(vertex, vrt, sizeof(SG::VertexFaceFormat) * 4);
		}
		else if( axis == GridBCM::AXIS_Y ){
			const SG::VertexFaceFormat vrt[4] = {
				SG::VertexFaceFormat( org[0],          org[1],          org[2],          0.0, 1.0, 0.0, tp[0], tp[2]),
				SG::VertexFaceFormat( org[0] + rgn[0], org[1],          org[2],          0.0, 1.0, 0.0, tp[1], tp[2]),
				SG::VertexFaceFormat( org[0] + rgn[0], org[1],          org[2] + rgn[2], 0.0, 1.0, 0.0, tp[1], tp[3]),
				SG::VertexFaceFormat( org[0],          org[1],          org[2] + rgn[2], 0.0, 1.0, 0.0, tp[0], tp[3])
			};
			memcpy(vertex, vrt, sizeof(SG::VertexFaceFormat) * 4);
		}
		else if( axis == GridBCM::AXIS_Z ){
			const SG::VertexFaceFormat vrt[4] = {
				SG::VertexFaceFormat( org[0],          org[1],          org[2],          0.0, 0.0, 1.0, tp[0], tp[2]),
				SG::VertexFaceFormat( org[0] + rgn[0], org[1],          org[2],          0.0, 0.0, 1.0, tp[1], tp[2]),
				SG::VertexFaceFormat( org[0] + rgn[0], org[1] + rgn[1], org[2],          0.0, 0.0, 1.0, tp[1], tp[3]),
				SG::VertexFaceFormat( org[0],          org[1] + rgn[1], org[2],          0.0, 0.0, 1.0, tp[0], tp[3])
			};
			memcpy(vertex, vrt, sizeof(SG::VertexFaceFormat) * 4);
		}
		GLuint idx[6] = {0, 1, 2, 2, 3, 0};
		for( int i = 0; i < 6; i++){
			idx[i] += idxOffset;
		}

		memcpy(index, idx, sizeof(SG::IndexFormat) * 6);

	}

	void BlockSetter(const Vec3r& org, const Vec3r& sz, const GridBCM::AXIS axis,
	                 SG::VertexLineFormat* vertex, SG::IndexFormat* index, size_t idxOffset)
	{
		if( axis == GridBCM::AXIS_X ){
			const GLfloat vrt[12] = {
				org[0],         org[1],         org[2],
				org[0],         org[1] + sz[1], org[2],
				org[0],         org[1] + sz[1], org[2] + sz[2],
				org[0],         org[1],         org[2] + sz[2]
			};

			memcpy(vertex, vrt, sizeof(GLfloat) * 12);
		}
		else if( axis == GridBCM::AXIS_Y ){
			const GLfloat vrt[12] = {
				org[0],         org[1],         org[2],
				org[0] + sz[0], org[1],         org[2],
				org[0] + sz[0], org[1],         org[2] + sz[2],
				org[0],         org[1],         org[2] + sz[2]
			};

			memcpy(vertex, vrt, sizeof(GLfloat) * 12);
		}
		else if( axis == GridBCM::AXIS_Z ){
			const GLfloat vrt[12] = {
				org[0],         org[1],         org[2],
				org[0] + sz[0], org[1],         org[2],
				org[0] + sz[0], org[1] + sz[1], org[2],
				org[0],         org[1] + sz[1], org[2]
			};

			memcpy(vertex, vrt, sizeof(GLfloat) * 12);
		}

		GLuint idx[4 * 2] = {0, 1, 1, 2, 2, 3, 3, 0};
		for( int i = 0; i < 8; i++){
			idx[i] += idxOffset;
		}

		memcpy(index, idx, sizeof(GLuint) * 8);
	}


	void GridSetter(const Vec3r& org, const Vec3r& rgn, const GridBCM::AXIS axis, const int divU, const int divV,
	                     SG::VertexLineFormat* vertex, SG::IndexFormat* index, size_t idxOffset)
	{
		if( axis == GridBCM::AXIS_X) {
			for(int u = 0; u < divU-1; u++){
				float pitch = rgn.y / (float)divU;
				vertex[u * 2 + 0] = SG::VertexLineFormat(org.x, org.y + pitch * (u+1), org.z);
				vertex[u * 2 + 1] = SG::VertexLineFormat(org.x, org.y + pitch * (u+1), org.z + rgn.z);
				index[u * 2 + 0] = idxOffset + u * 2;
				index[u * 2 + 1] = idxOffset + u * 2 + 1;
			}
			for(int v = 0; v < divV-1; v++){
				float pitch = rgn.z / (float)divV;
				vertex[(divU-1) * 2 + v * 2 + 0] = SG::VertexLineFormat(org.x, org.y,         org.z + pitch * (v+1));
				vertex[(divU-1) * 2 + v * 2 + 1] = SG::VertexLineFormat(org.x, org.y + rgn.y, org.z + pitch * (v+1));
				index[ (divU-1) * 2 + v * 2 + 0] = idxOffset + (divU-1) * 2 + v * 2;
				index[ (divU-1) * 2 + v * 2 + 1] = idxOffset + (divU-1) * 2 + v * 2 + 1;
			}
			return;
		}

		else if( axis == GridBCM::AXIS_Y){
			for(int u = 0; u < divU-1; u++){
				float pitch = rgn.x / (float)divU;
				vertex[u * 2 + 0] = SG::VertexLineFormat(org.x + pitch * (u+1), org.y, org.z);
				vertex[u * 2 + 1] = SG::VertexLineFormat(org.x + pitch * (u+1), org.y, org.z + rgn.z);
				index[u * 2 + 0] = idxOffset + u * 2;
				index[u * 2 + 1] = idxOffset + u * 2 + 1;
			}
			for(int v = 0; v < divV-1; v++){
				float pitch = rgn.z / (float)divV;
				vertex[(divU-1) * 2 + v * 2 + 0] = SG::VertexLineFormat(org.x,         org.y, org.z + pitch * (v+1));
				vertex[(divU-1) * 2 + v * 2 + 1] = SG::VertexLineFormat(org.x + rgn.x, org.y, org.z + pitch * (v+1));
				index[ (divU-1) * 2 + v * 2 + 0] = idxOffset + (divU-1) * 2 + v * 2;
				index[ (divU-1) * 2 + v * 2 + 1] = idxOffset + (divU-1) * 2 + v * 2 + 1;
			}
			return;
		}

		else if( axis == GridBCM::AXIS_Z){
			for(int u = 0; u < divU-1; u++){
				float pitch = rgn.x / (float)divU;
				vertex[u * 2 + 0] = SG::VertexLineFormat(org.x + pitch * (u+1), org.y,         org.z);
				vertex[u * 2 + 1] = SG::VertexLineFormat(org.x + pitch * (u+1), org.y + rgn.y, org.z);
				index[u * 2 + 0] = idxOffset + u * 2;
				index[u * 2 + 1] = idxOffset + u * 2 + 1;
			}
			for(int v = 0; v < divV-1; v++){
				float pitch = rgn.y / (float)divV;
				vertex[(divU-1) * 2 + v * 2 + 0] = SG::VertexLineFormat(org.x,         org.y + pitch * (v+1), org.z);
				vertex[(divU-1) * 2 + v * 2 + 1] = SG::VertexLineFormat(org.x + rgn.x, org.y + pitch * (v+1), org.z);
				index[ (divU-1) * 2 + v * 2 + 0] = idxOffset + (divU-1) * 2 + v * 2;
				index[ (divU-1) * 2 + v * 2 + 1] = idxOffset + (divU-1) * 2 + v * 2 + 1;
			}
			return;
		}
		return;
	}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void GetNodesOnPlane_rec(const GridBCM::AXIS axis, const Node* node, size_t planePosition, const int maxLevel,
	                         std::vector<const Node*>& leafNodes)
	{
		//std::cout << "Pedigree : " << node->getPedigree() << std::endl;
		if(node->isLeafNode()){
			if(!node->isActive()){
				std::cout << "node (" << node->getPedigree() << ") is not Active " << __func__ << ":" << __FILE__ << ":" << __LINE__ << std::endl;
				return;
			}

			leafNodes.push_back(node);
			//std::cout << "node (" << node->getPedigree() << ") is leaf " << __func__ << ":" << __FILE__ << ":" << __LINE__ << std::endl;
			return;
		}

		int lv = node->getLevel();

		int cidTbl[2][4] = {0};
		if(axis == GridBCM::AXIS_X){
			cidTbl[0][0] = 0; cidTbl[0][1] = 2; cidTbl[0][2] = 4; cidTbl[0][3] = 6;
			cidTbl[1][0] = 1; cidTbl[1][1] = 3; cidTbl[1][2] = 5; cidTbl[1][3] = 7;
		}
		if(axis == GridBCM::AXIS_Y){
			cidTbl[0][0] = 0; cidTbl[0][1] = 1; cidTbl[0][2] = 4; cidTbl[0][3] = 5;
			cidTbl[1][0] = 2; cidTbl[1][1] = 3; cidTbl[1][2] = 6; cidTbl[1][3] = 7;
		}
		if(axis == GridBCM::AXIS_Z){
			cidTbl[0][0] = 0; cidTbl[0][1] = 1; cidTbl[0][2] = 2; cidTbl[0][3] = 3;
			cidTbl[1][0] = 4; cidTbl[1][1] = 5; cidTbl[1][2] = 6; cidTbl[1][3] = 7;
		}

		//int high = ((planePosition >> (maxLevel - (lv+1)))) > 0 ? 1 : 0;
		int high = ( planePosition >> (maxLevel - (lv+1))) & 1 == 1 ? 1 : 0;
		for(int i = 0; i < 4; i++){
			const Node* n = node->getChild(cidTbl[high][i]);
			GetNodesOnPlane_rec(axis, n, planePosition, maxLevel, leafNodes);
		}
	}

	void GetNodesOnPlane(const GridBCM::AXIS axis, const Node* root, const size_t planePosition, const int maxLevel,
	                     std::vector<const Node*>& leafNodes)
	{
		GetNodesOnPlane_rec(axis, root, planePosition, maxLevel, leafNodes);
	}

	size_t SearchMaxBlockCount(const GridBCM::AXIS axis, const BCMOctree* octree, const int maxLevel){
		int iaxis = 0;
		if(axis == GridBCM::AXIS_X){ iaxis = 0; }
		if(axis == GridBCM::AXIS_Y){ iaxis = 1; }
		if(axis == GridBCM::AXIS_Z){ iaxis = 2; }

		const RootGrid *rootGrid = octree->getRootGrid();
		size_t rootDims[3] = { rootGrid->getSizeX(), rootGrid->getSizeY(), rootGrid->getSizeZ() };
		size_t numSlice = rootDims[iaxis] * (1 << maxLevel);

		size_t maxBlocks = 0;

		//const char* axisStr[3] = { "AXIS_X", "AXIS_Y", "AXIS_Z" };
		//printf("Search MaxBlock Count [%s] : Start ", axisStr[iaxis]); fflush(stdout);
		for(size_t s = 0; s < numSlice; s++){
			size_t rPos = s / (1 << maxLevel);
			size_t nPos = s % (1 << maxLevel);

			std::vector<const Node*> leafNodes;
			// Collect leafNodes which crossed SlicePosigion from Octree
			int nru = axis == GridBCM::AXIS_X ? rootDims[1] : rootDims[0];
			int urv = axis == GridBCM::AXIS_Z ? rootDims[1] : rootDims[2];
			for(int v = 0; v < urv; v++){
				for(int u = 0; u < nru; u++){
					int x = axis == GridBCM::AXIS_X ? rPos : u;
					int y = axis == GridBCM::AXIS_Y ? rPos : axis == GridBCM::AXIS_X ? u : v;
					int z = axis == GridBCM::AXIS_Z ? rPos : v;
					const Node* root = octree->getRootNode(rootGrid->index2rootID(x, y, z));
					GetNodesOnPlane(axis, root, nPos, maxLevel, leafNodes);
				}
			}

			maxBlocks = MAX_(maxBlocks, leafNodes.size());
		}
		//printf(" Complete : MaxBlocks : %ld\n", maxBlocks); fflush(stdout);

		return maxBlocks;
	}

} // namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned char GetCellID( const BCMFileIO::bitVoxelCell* bitVoxel, const unsigned int bitWidth, const size_t size[3], const size_t position[3] )
{
	using namespace BCMFileIO;
	const unsigned char vox_per_cell = (sizeof(bitVoxelCell) * 8) / bitWidth;

	unsigned char mask = 0;
	for(int i = 0; i < bitWidth; i++) mask += (1 << i);

	size_t loc = position[0] + (position[1] + position[2] * size[1]) * size[0];

	size_t       cellIdx =  loc / vox_per_cell;
	unsigned int bitIdx  = (loc % vox_per_cell) * bitWidth;

	return (bitVoxel[cellIdx] >> bitIdx) & mask;
}

class LB
{
public:
	LB( const BCMFileIO::bitVoxelCell* bitVoxel, const size_t bitVoxelSize )
	 : m_rleBuf(NULL)
	{
		using namespace BCMFileIO;
		m_rleBuf  = rleEncode<bitVoxelCell, unsigned char>( bitVoxel, bitVoxelSize * sizeof(bitVoxelCell), &m_rleSize );
		m_orgSize = bitVoxelSize;
		//m_bitVoxel = bitVoxel;
	}

	~LB() {
		delete [] m_rleBuf;
		//delete m_bitVoxel;
	}

	const BCMFileIO::bitVoxelCell* GetBitVoxel() const {
		using namespace BCMFileIO;
		return rleDecode<bitVoxelCell, unsigned char>(m_rleBuf, m_rleSize, m_orgSize * sizeof(bitVoxelCell));
		//return m_bitVoxel;
	}

	//static Vec3i BlockSize;
	//static int   bitWidth;

private:
	//const BCMFileIO::bitVoxelCell* m_bitVoxel;
	unsigned char* m_rleBuf;
	size_t         m_rleSize;
	size_t         m_orgSize;
};

class LeafBlocks
{
public:
	LeafBlocks( const BCMFileIO::IdxBlock* ib, BCMFileIO::PartitionMapper *pmapper)
	{
		m_status = LoadLeafBlockFile( ib, pmapper );
	}


	~LeafBlocks()
	{
		for(std::vector<LB*>::iterator it = m_lb.begin(); it < m_lb.end(); it++){
			delete (*it);
		}
	}

	bool GetStatus() const { return m_status; }

	const LB* GetLeafBlock( const unsigned int NodeID ) const {
		if( NodeID >= m_lb.size() ){ return NULL; }
		return m_lb[NodeID];
	}

	const size_t* GetLeafBlockSize() const { return &m_blockSize[0]; }
	const unsigned int GetBitWidth() const { return m_bitWidth; }

private:
	bool LoadLeafBlockFile( const BCMFileIO::IdxBlock* ib, BCMFileIO::PartitionMapper *pmapper )
	{
		using namespace BCMFileIO;
		if( ib->kind != LB_CELLID ) { return false; }

		LBHeader header;
		std::vector<CellIDCapsule> cidCapsules;
		printf("    Load LeafBlock File(s) Start : ");
		if( ib->isGather ){
			if( !Load_LeafBlock_CellID_Gather(ib->rootDir + ib->dataDir, ib, pmapper, header, cidCapsules ) ){
				fprintf(stderr, "[ERR] %s : %s:%d\n", __func__, __FILE__, __LINE__);
				return false;
			}
		}else{
			if( !Load_LeafBlock_CellID(ib->rootDir + ib->dataDir, ib, pmapper, header, cidCapsules ) ){
				fprintf(stderr, "[ERR] %s : %s:%d\n", __func__, __FILE__, __LINE__);
				return false;
			}
		}
		printf("Complete\n");

		int numBlocks = pmapper->GetEnd(0) - pmapper->GetStart(0);
		m_lb.reserve(numBlocks);
		m_blockSize[0] = header.size[0];
		m_blockSize[1] = header.size[1];
		m_blockSize[2] = header.size[2];
		m_bitWidth = header.bitWidth;

		int did = 0;
		int fid = 0;
		std::vector<PartitionMapper::FDIDList> fdidlists;
		pmapper->GetFDIDLists(0, fdidlists);

		printf("    Reconstructing LeafBlocks to Internal Format [Start] "); fflush(stdout);
#ifndef _OPENMP
		for(std::vector<PartitionMapper::FDIDList>::iterator file = fdidlists.begin(); file != fdidlists.end(); ++file){
			unsigned char* voxels = DecompCellIDData( header, cidCapsules[fid] );

			for( std::vector<int>::iterator fdid = file->FDIDs.begin(); fdid != file->FDIDs.end(); ++fdid){

				// Convert File format -> internal format
				Vec3i bsz( header.size[0], header.size[1], header.size[2] );
				int   vc = header.vc;

				unsigned char* block = new unsigned char[bsz.x * bsz.y * bsz.z];
				//printf("header.vc : %d\n", header.vc);

				size_t fbsize = (bsz.x + vc*2) * (bsz.y + vc*2) * (bsz.z + vc*2);
				const unsigned char* pvox = &voxels[fbsize * (*fdid)];

				for(int z = 0; z < header.size[2]; z++){
					for(int y = 0; y < header.size[1]; y++){
						size_t bloc  = 0    + ( y     +  z     *  bsz.y         ) *  bsz.x;
						size_t fbloc = 0+vc + ((y+vc) + (z+vc) * (bsz.y+(vc*2)) ) * (bsz.x+(vc*2));
						memcpy(&block[bloc], &pvox[fbloc], sizeof(unsigned char) * bsz.x);
					}
				}

				size_t bitVoxelSize = 0;
				bitVoxelCell* bitVoxel = CompressBitVoxel(&bitVoxelSize, bsz.x * bsz.y * bsz.z, block, header.bitWidth);
				LB *lb = new LB(bitVoxel, bitVoxelSize);
				m_lb.push_back(lb);

				delete [] bitVoxel;
				delete [] block;
				did++;
			}
			fid++;
			printf("*"); fflush(stdout);
		}

		if(did != numBlocks ){ return false; }
#else
		m_lb.resize(numBlocks);

		size_t *fdidSizeTable = new size_t[fdidlists.size()];
		for(int i = 0; i < fdidlists.size(); i++){
			fdidSizeTable[i] = fdidlists[i].FDIDs.size();
		}

		#pragma omp parallel for
		for(int i = 0; i < fdidlists.size(); i++){
			size_t iLoc = 0;
			for(int j = 0; j < i; j++){
				iLoc += fdidSizeTable[j];
			}

			unsigned char* voxels = DecompCellIDData( header, cidCapsules[i] );

			for(int j = 0; j < fdidSizeTable[i]; j++){
				Vec3i bsz( header.size[0], header.size[1], header.size[2] );
				int   vc = header.vc;
				//int *p
				unsigned char* block = new unsigned char[bsz.x * bsz.y * bsz.z];
				size_t fbsize = (bsz.x + vc*2) * (bsz.y + vc*2) * (bsz.z + vc*2);
				const unsigned char* pvox = &voxels[ fbsize * j ];

				for(int z = 0; z < header.size[2]; z++){
					for(int y = 0; y < header.size[1]; y++){
						size_t bloc  = 0    + ( y     +  z     *  bsz.y         ) *  bsz.x;
						size_t fbloc = 0+vc + ((y+vc) + (z+vc) * (bsz.y+(vc*2)) ) * (bsz.x+(vc*2));
						memcpy(&block[bloc], &pvox[fbloc], sizeof(unsigned char) * bsz.x);
					}
				}

				size_t bitVoxelSize = 0;
				bitVoxelCell* bitVoxel = CompressBitVoxel(&bitVoxelSize, bsz.x * bsz.y * bsz.z, block, header.bitWidth);
				LB *lb = new LB(bitVoxel, bitVoxelSize);
				m_lb[j + iLoc] = lb;

				delete [] bitVoxel;
				delete [] block;
			}
			printf("*"); fflush(stdout);
		}
		delete [] fdidSizeTable;
#endif
		printf(" [Complete]\n");
		return true;
	}

private:
	std::vector<LB*> m_lb;

	size_t       m_blockSize[3];
	unsigned int m_bitWidth;

	bool m_status;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SliceFace
{
public:
	SliceFace( const GridBCM::AXIS axis, const size_t maxBlocks, const int numGU, const int numGV, const int numGW )
	 : m_axis(axis), m_maxBlocks(maxBlocks), m_numGU(numGU), m_numGV(numGV), m_numGW(numGW), m_img(NULL)
	{
		size_t nv = ( (m_numGU - 1) * 2 + (m_numGV - 1) * 2 ) * m_maxBlocks;
		m_blocks = new SG::Geometry<SG::VertexLineFormat>(m_maxBlocks * 4, m_maxBlocks * 8);
		m_grids  = new SG::Geometry<SG::VertexLineFormat>( nv, nv );
		m_cells  = new SG::Geometry<SG::VertexFaceFormat>( m_maxBlocks * 4, m_maxBlocks * 6);

		m_cellTex = NULL;
		InitTexture();

		SetVisible(true);
	}

	~SliceFace()
	{
		delete m_blocks;
		delete m_grids;
		delete m_cells;
		delete m_cellTex;
	}

	void CreateSlice(const BCMOctree* octree, const size_t maxLevel,
	                 const LeafBlocks* leafBlocks,
	                 const Vec3r& globalOrigin, const Vec3r& globalRegion, const size_t position)
	{
		using namespace BCMFileIO;
		int iaxis = 0;
		if(m_axis == GridBCM::AXIS_X){ iaxis = 0; }
		if(m_axis == GridBCM::AXIS_Y){ iaxis = 1; }
		if(m_axis == GridBCM::AXIS_Z){ iaxis = 2; }

		const RootGrid* rootGrid = octree->getRootGrid();
		size_t rootDims[3]  = { rootGrid->getSizeX(), rootGrid->getSizeY(), rootGrid->getSizeZ() };
		size_t maxCellCount = rootDims[iaxis] * (1 << maxLevel) * m_numGW;

		// RootGrid Position
		size_t rPos = float(rootDims[iaxis]) * float(position) / float(maxCellCount);
		// Cell Position within RootGrid
		size_t lPos = position - (rPos * (1 << maxLevel) * m_numGW);
		// Node Position within RootGrid
		size_t nPos = lPos / m_numGW;

		std::vector<const Node*> leafNodes;
		// Collect leafNodes which crossed SlicePosigion from Octree
		int nru = m_axis == GridBCM::AXIS_X ? rootDims[1] : rootDims[0];
		int urv = m_axis == GridBCM::AXIS_Z ? rootDims[1] : rootDims[2];
		for(int v = 0; v < urv; v++){
			for(int u = 0; u < nru; u++){
				int x = m_axis == GridBCM::AXIS_X ? rPos : u;
				int y = m_axis == GridBCM::AXIS_Y ? rPos : m_axis == GridBCM::AXIS_X ? u : v;
				int z = m_axis == GridBCM::AXIS_Z ? rPos : v;
				const Node* root = octree->getRootNode(rootGrid->index2rootID(x, y, z));
				GetNodesOnPlane(m_axis, root, nPos, maxLevel, leafNodes);
			}
		}

		Vec3r rootRegion( globalRegion.x / static_cast<REAL_TYPE>(rootDims[0]),
		                  globalRegion.y / static_cast<REAL_TYPE>(rootDims[1]),
						  globalRegion.z / static_cast<REAL_TYPE>(rootDims[2]) );

		float cellPitch = globalRegion[iaxis] / REAL_TYPE(maxCellCount);
		// Create All block boundary line geometry by current slice face
		SG::VertexLineFormat *bvrt = m_blocks->GetVertexBuffer();
		SG::IndexFormat      *bidx = m_blocks->GetIndexBuffer();

		// Create All block grid line geometry by current slice face
		size_t gridStride = (m_numGU -1) * 2 + (m_numGV - 1) * 2;
		SG::VertexLineFormat *gvrt = m_grids->GetVertexBuffer();
		SG::IndexFormat      *gidx = m_grids->GetIndexBuffer();

		// Create All block plane geometry by current slice face
		SG::VertexFaceFormat *cvrt = m_cells->GetVertexBuffer();
		SG::IndexFormat      *cidx = m_cells->GetIndexBuffer();


		// Create All block cell texture buffer by current slice face
		int subTexSize[2] = { m_texSize[0], (leafNodes.size() * m_numGU / m_texSize[0] + 1) * m_numGV };
		m_img = new int[subTexSize[0] * subTexSize[1]];

		// Loop by collected leaf nodes
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(int cnt = 0; cnt < leafNodes.size(); cnt++){
			size_t bIdxOffset = 0;
			size_t gIdxOffset = 0;
			size_t cIdxOffset = 0;

			// Set Cell Texture
			int NodeID = leafNodes[cnt]->getBlockID();

			// Get bitVoxel corresponds NodeID from LeafBlock Pool
			// BitVoxel is allocated by GetBitVoxel() !!! because leafBlock Pool has encoded buffer.
			const bitVoxelCell* bitVoxel = leafBlocks->GetLeafBlock(NodeID)->GetBitVoxel();
			if(!bitVoxel){
				fprintf(stderr, "[ERR] : %s [%s:%d]\n", __func__, __FILE__, __LINE__);
				//return;
			}

			const Pedigree p = leafNodes[cnt]->getPedigree();

			Vec3r rootOrg(globalOrigin.x + (rootRegion.x * rootGrid->rootID2indexX(p.getRootID())),
			              globalOrigin.y + (rootRegion.y * rootGrid->rootID2indexY(p.getRootID())),
						  globalOrigin.z + (rootRegion.z * rootGrid->rootID2indexZ(p.getRootID())));

			unsigned int lv = p.getLevel();

			Vec3r rgn( rootRegion.x / (1 << lv),     rootRegion.y / (1 << lv),     rootRegion.z / (1 << lv)     );
			Vec3r org( rootOrg.x + p.getX() * rgn.x, rootOrg.y + p.getY() * rgn.y, rootOrg.z + p.getZ() * rgn.z );
			if( m_axis == GridBCM::AXIS_X) org.x = globalOrigin.x + position * cellPitch;
			if( m_axis == GridBCM::AXIS_Y) org.y = globalOrigin.y + position * cellPitch;
			if( m_axis == GridBCM::AXIS_Z) org.z = globalOrigin.z + position * cellPitch;

			// Set Block Boundary Line
			bIdxOffset = cnt * 4;
			BlockSetter(org, rgn, m_axis, &bvrt[cnt * 4], &bidx[cnt * 8], bIdxOffset);

			// Set Block Grid Line
			gIdxOffset = cnt * gridStride;
			GridSetter(org, rgn, m_axis, m_numGU, m_numGV, &gvrt[cnt * gridStride], &gidx[cnt * gridStride], gIdxOffset);

			// Set Cell Geometry (simply plane)
			const int texIdx[2] = {
				(cnt * m_numGU) % m_texSize[0],
				(cnt * m_numGU) / m_texSize[0] * m_numGV
			};

			// block plane needs texcoord
			cIdxOffset = cnt * 4;
			CellSetter(org, rgn, m_axis, m_numGU, m_numGV, texIdx[0], texIdx[1], m_texSize[0], m_texSize[1],
			           &cvrt[cnt * 4], &cidx[cnt * 6], cIdxOffset);

			size_t oPos = m_axis == GridBCM::AXIS_X ? p.getX() : m_axis == GridBCM::AXIS_Y ? p.getY() : p.getZ();
			size_t bPos = ( lPos - (oPos << (maxLevel - lv)) * m_numGW ) >> (maxLevel - lv);
			int uidx = m_axis == GridBCM::AXIS_X ? 1 : 0;
			int vidx = m_axis == GridBCM::AXIS_Z ? 1 : 2;

			size_t pos[3] = { bPos, bPos, bPos };

			for(int v = 0; v < m_numGV; v++){
				for(int u = 0; u < m_numGU; u++){
					size_t loc = (u + texIdx[0]) + (v + texIdx[1]) * subTexSize[0];
					pos[uidx] = u;
					pos[vidx] = v;
					int id = GetCellID(bitVoxel, leafBlocks->GetBitWidth(), leafBlocks->GetLeafBlockSize(), pos);
					if( id == 0 ){
						//m_img[loc] = 0;
						m_img[loc] = (255 << 0) | (255 << 8) | (255 << 16) | (255 << 24);
					}else{
						int r  = ((id >> 0) & 1) * 127 + ((id >> 3) & 1) * 127;
						int g  = ((id >> 1) & 1) * 127 + ((id >> 4) & 1) * 127;
						int b  = ((id >> 2) & 1) * 127 + ((id >> 5) & 1) * 127;
						//m_img[loc] = (b << 0) | (g << 8) | (r << 16) | (255 << 24);
						m_img[loc] = (b << 0) | (g << 8) | (r << 16) | (255 << 24);
					}
				}
			}
			// needs bitVoxel deallocation !!!
			delete [] bitVoxel;
		}

		// Set all block boundary line geometry by current slice face
		m_blocks->UploadVertices(leafNodes.size() * 4);
		m_blocks->UploadIndices(leafNodes.size() * 4 * 2);

		// Set all block grid line geometry by current slice face
		m_grids->UploadVertices(leafNodes.size() * gridStride);
		m_grids->UploadIndices(leafNodes.size() * gridStride);

		// Set all block plane geometry by current slice face
		m_cells->UploadVertices(leafNodes.size() * 4);
		m_cells->UploadIndices(leafNodes.size() * 6);

		// Set all block cell texture buffer  by current slice face
		m_cellTex->WriteImage(m_img, subTexSize[0], subTexSize[1]);

		//SaveImage(m_img, subTexSize[0], subTexSize[1]);
	}

	void SetVisible(bool visible){
		m_visible = visible;
		m_blocks->SetVisible(visible);
		m_grids->SetVisible(visible);
		m_cells->SetVisible(visible);
	}

	bool GetVisible() const { return m_visible; }

	const SG::Geometry<SG::VertexLineFormat>* GetBlockGeometry() const { return m_blocks; }
	const SG::Geometry<SG::VertexLineFormat>* GetGridGeometry()  const { return m_grids;  }
	const SG::Geometry<SG::VertexFaceFormat>* GetCellGeometry()  const { return m_cells;  }
	unsigned int GetCellTexture() const { return m_cellTex->GetTexture(); }

private:

	void InitTexture(){
		int texture_size_limit = 0;
		glGetIntegerv( GL_MAX_TEXTURE_SIZE, &texture_size_limit );

		if( m_maxBlocks * m_numGU <= texture_size_limit ){
			m_texSize[0] = m_maxBlocks * m_numGU;
			m_texSize[1] = m_numGV;
		}else{
			bool status = false;
			for(int i = 2; i <= (texture_size_limit / m_numGV); i++){
				size_t nbu = m_maxBlocks / i + (m_maxBlocks % i == 0 ? 0 : 1);
				if( nbu * m_numGU < texture_size_limit ){
					m_texSize[0] = nbu * m_numGU;
					m_texSize[1] = i   * m_numGV;
					status = true;
					break;
				}
			}
			if( !status ) {
				fprintf(stderr, "[ERR] : cannot create texture [%s:%d]\n", __FILE__, __LINE__);
				exit(0);
			}
		}

		//printf("m_maxBlocks : %5ld m_texSize[ %5d %5d ]\n", m_maxBlocks, m_texSize[0], m_texSize[1]);

		// expand Texture size to 2^n
		for(int i = 0; i < 2; i++){
			int x = 0, j = 0;
			while(x < m_texSize[i]){ x = 1 << j; j++; }
			m_texSize[i] = x;
		}

		//printf("expand force 2^n \n");
		const char* axisStr[3] = { "AXIS_X", "AXIS_Y", "AXIS_Z" };
		int iaxis = 0;
		if(m_axis == GridBCM::AXIS_X){ iaxis = 0; }
		if(m_axis == GridBCM::AXIS_Y){ iaxis = 1; }
		if(m_axis == GridBCM::AXIS_Z){ iaxis = 2; }
		printf("    Texture Size (%s) is [%5d, %5d]\n", axisStr[iaxis], m_texSize[0], m_texSize[1]);

		if( m_cellTex ) delete m_cellTex;
		m_cellTex = new TextureObject(m_texSize[0], m_texSize[1], 4, 32, -1);
		//m_img = new int[m_texSize[0] * m_texSize[1]];
	}

private:
	const GridBCM::AXIS m_axis;
	const size_t        m_maxBlocks;
	const int           m_numGU;
	const int           m_numGV;
	const int           m_numGW;

	SG::Geometry<SG::VertexLineFormat> *m_blocks;
	SG::Geometry<SG::VertexLineFormat> *m_grids;
	SG::Geometry<SG::VertexFaceFormat> *m_cells;
	TextureObject                      *m_cellTex;

	int m_texSize[2];
	int *m_img;

	bool m_visible;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GridBCM::GridBCM( const std::string& filename )
{
	for(int i = 0; i < 3; i++){
		m_slice[i] = NULL;
	}

	m_pmapper    = NULL;
	m_octree     = NULL;
	m_leafBlocks = NULL;

	m_render = new Render::SGRender();

	m_active   = AXIS_X;
	m_showGrid = true;
	LoadFile( filename );
}


GridBCM::~GridBCM()
{
	Deinit();

	if( m_pmapper    ) delete m_pmapper;
	if( m_octree     ) delete m_octree;
	if( m_leafBlocks ) delete m_leafBlocks;
}

bool GridBCM::LoadFile(const std::string& filename)
{
	using namespace BCMFileIO;
	if( m_pmapper    ) delete m_pmapper;
	if( m_octree     ) delete m_octree;
	if( m_leafBlocks ) delete m_leafBlocks;

	std::string dir       = ConvertPath(filename);
	std::string targetDir = GetDirectory(dir);

	std::string octreeFilename;
	std::vector<IdxProc>  idxProcList;
	std::vector<IdxBlock> idxBlockList;

	Vec3i blockSize;
	if( !LoadIndex(filename, targetDir, m_globalOrigin, m_globalRegion, octreeFilename, blockSize, idxProcList, idxBlockList) ){
		fprintf(stderr, "[ERR] : %s [%s:%d]\n", __func__, __FILE__, __LINE__);
		return false;
	}

	std::string octreeFilepath = targetDir + octreeFilename;

	IdxBlock *ib = NULL;
	for(std::vector<IdxBlock>::iterator it = idxBlockList.begin(); it != idxBlockList.end(); ++it){
		if(it->kind == LB_CELLID){
			ib = &(*it);
			break;
		}
	}

	printf("Load Octree File : Start\n"); fflush(stdout);
	// Load Octree
	if( !LoadFileOctree(octreeFilepath) ){
		fprintf(stderr, "[ERR] : %s [%s:%d]\n", __func__, __FILE__, __LINE__);
		return false;
	}
	printf("Load Octree File : Complete\n"); fflush(stdout);

	printf("Load LeafBLock Files : Start\n");
	// Load Leaf Block
	m_pmapper = new PartitionMapper(idxProcList.size(), 1, m_octree->getNumLeafNode());
	m_leafBlocks = new LeafBlocks(ib, m_pmapper);
	printf("Load LeafBlock Files : Complete\n");

	if( !m_leafBlocks->GetStatus() ){
		fprintf(stderr, "[ERR] : %s [%s:%d]\n", __func__, __FILE__, __LINE__);
		return false;
	}

	const size_t *sz = m_leafBlocks->GetLeafBlockSize();
	m_blockSize[0] = sz[0];
	m_blockSize[1] = sz[1];
	m_blockSize[2] = sz[2];

	for(int i = 0; i < 3; i++){
		m_maxCellCount[i] = (1 << m_maxLevel) * m_rootDims[i] * m_blockSize[i];
		m_slicePos[i] = m_maxCellCount[i] / 2;
	}

	if( !Init() ){
		fprintf(stderr, "[ERR] : %s [%s:%d]\n", __func__, __FILE__, __LINE__);
		return false;
	}

	printf("\n");
	printf("=================== Data Information ==================\n");
	std::cout << "  GlobalOrigin  : "  << m_globalOrigin << std::endl;
	std::cout << "  GlobalRegion  : "  << m_globalRegion << std::endl;
	std::cout << "  RootGrid      : (" << m_rootDims[0]  << ", " << m_rootDims[1]  << ", " << m_rootDims[2]  << ")" << std::endl;
	std::cout << "  LeafBlockSize : (" << m_blockSize[0] << ", " << m_blockSize[1] << ", " << m_blockSize[2] << ")" << std::endl;
	std::cout << "  maxLevel      : "  << m_maxLevel << std::endl;
	std::cout << "  numLeafBlocks : "  << m_octree->getNumLeafNode() << std::endl;
	printf("=======================================================\n\n");

	return true;
}

bool GridBCM::Deinit()
{
	for(int i = 0; i < 3; i++){
		if( m_slice[i] ) delete m_slice[i];
	}
	if( m_bbox ) delete m_bbox;
}


bool GridBCM::Init()
{
	if( !m_octree || !m_pmapper ){
		fprintf(stderr, "[ERR] : %s [%s:%d]\n", __func__, __FILE__, __LINE__);
		return false;
	}

	printf("Initialize Graphics : Start\n");
	static bool isInit = false;
	if( isInit ) Deinit();

	// BBox Line
	{
		m_bbox = new SG::Geometry<SG::VertexLineFormat>(8, 24);
		SG::VertexLineFormat vrt[8];
		SG::IndexFormat     idx[24];
		size_t idxOffset = 0;
		BoxSetter(m_globalOrigin, m_globalRegion, vrt, idx, idxOffset);
		m_bbox->SetVertices(vrt, 8);
		m_bbox->SetIndices(idx, 24);
	}

	for(int i = 0; i < 3; i++){
		AXIS axis = (AXIS)( 1 << i );
		const int numRootUVW[3] = {
			axis == AXIS_X ? m_rootDims[1] : m_rootDims[0],
			axis == AXIS_Z ? m_rootDims[1] : m_rootDims[2],
			axis == AXIS_X ? m_rootDims[0] : axis == AXIS_Y ? m_rootDims[1] : m_rootDims[2]
		};
		const int numBlockUVW[3] = {
			axis == AXIS_X ? m_blockSize[1] : m_blockSize[0],
			axis == AXIS_Z ? m_blockSize[1] : m_blockSize[2],
			axis == AXIS_X ? m_blockSize[0] : axis == AXIS_Y ? m_blockSize[1] : m_blockSize[2]
		};

		//size_t maxBlocks = numRootUVW[0] * numRootUVW[1] * (1 << m_maxLevel ) * (1 << m_maxLevel);
		size_t maxBlocks = SearchMaxBlockCount( axis, m_octree, m_maxLevel );

		m_slice[i] = new SliceFace( axis, maxBlocks, numBlockUVW[0], numBlockUVW[1], numBlockUVW[2]);
		m_slice[i]->SetVisible(false);
	}
	printf("Initialize Graphics : End\n");

	isInit = true;
	return true;
}

bool GridBCM::LoadFileOctree(const std::string& filename)
{
	using namespace std;
	using namespace BCMFileIO;

	OctHeader header;
	vector<Pedigree> pedigrees;

	FILE *fp = NULL;

	if( (fp = fopen(filename.c_str(), "rb")) == NULL ){
		fprintf(stderr, "open file error(%s) [%s:%d].\n", filename.c_str(), __FILE__, __LINE__);
		return false;
	}

	bool isNeedSwap = false;

	fread(&header, sizeof(header), 1, fp);

	if( header.identifier != OCTREE_FILE_IDENTIFIER ){
		BSwap32(&header.identifier);

		if( header.identifier != OCTREE_FILE_IDENTIFIER ){
			fprintf(stderr, "%s is not Octree file [%s:%d].\n", filename.c_str(), __FILE__, __LINE__);
			return false;
		}

		isNeedSwap = true;

		for(int i = 0; i < 3; i++){
			BSwap64(&header.org[i]);
			BSwap64(&header.rgn[i]);
			BSwap64(&header.rootDims[i]);
		}
		BSwap32(&header.maxLevel);
		BSwap64(&header.numLeaf);
	}

	pedigrees.resize(header.numLeaf);
	fread(&pedigrees[0], sizeof(Pedigree), header.numLeaf, fp);
	fclose(fp);

	if( isNeedSwap ){
		for(vector<Pedigree>::iterator it = pedigrees.begin(); it != pedigrees.end(); ++it){
			BSwap64(&(*it));
		}
	}

	Vec3r rootRegion( header.rgn[0] / static_cast<REAL_TYPE>(header.rootDims[0]),
	                  header.rgn[1] / static_cast<REAL_TYPE>(header.rootDims[1]),
			          header.rgn[2] / static_cast<REAL_TYPE>(header.rootDims[2]));

	if( fabs(rootRegion.x - rootRegion.y) >= 1.0e-10 || fabs(rootRegion.x - rootRegion.z) >= 1.0e-10 ) {
		fprintf(stderr, "[ERR] %s [%s:%d]\n", __func__, __FILE__, __LINE__);
		return false;
	}

	RootGrid *rootGrid = new RootGrid(header.rootDims[0], header.rootDims[1], header.rootDims[2]);
	m_octree  = new BCMOctree(rootGrid, pedigrees);


	m_globalOrigin = Vec3r(header.org);
	m_globalRegion = Vec3r(header.rgn);
	m_maxLevel = header.maxLevel;
	m_rootDims[0] = header.rootDims[0];
	m_rootDims[1] = header.rootDims[1];
	m_rootDims[2] = header.rootDims[2];

	return true;
}


size_t GridBCM::SetSlicePosition(size_t position)
{
	int i = 0;
	if(m_active == GridBCM::AXIS_X) i = 0;
	if(m_active == GridBCM::AXIS_Y) i = 1;
	if(m_active == GridBCM::AXIS_Z) i = 2;

	size_t pos = MAX_(MIN_(position, (m_maxCellCount[i]-1)), (size_t)0);
	m_slicePos[i] = pos;

	printf("slice position : (%5ld, %5ld, %5ld)\n", m_slicePos[0], m_slicePos[1], m_slicePos[2]);

	UpdateGeometry(m_active);
	return pos;
}


bool GridBCM::UpdateGeometry()
{
	for(int i = 0; i < 3; i++) {
		AXIS axis = (AXIS)(1 << i);
		if( !UpdateGeometry(axis) ){ return false; }
	}

	return true;
}

bool GridBCM::UpdateGeometry(const AXIS axis )
{
	for(int i = 0; i < 3; i++){
		if( (axis >> i) & 1 ) m_slice[i]->CreateSlice(m_octree, m_maxLevel, m_leafBlocks, m_globalOrigin, m_globalRegion, m_slicePos[i]);
	}
	return true;
}

void GridBCM::AddSlicePlane(const AXIS axis)
{
	for(int i = 0; i < 3; i++){
		if( (axis >> i) & 1 ) m_slice[i]->SetVisible(true);
	}

	UpdateGeometry(axis);
}

void GridBCM::DeleteSlicePlane(const AXIS axis)
{
	for(int i = 0; i < 3; i++){
		if( (axis >> i) & 1 ) m_slice[i]->SetVisible(false);
	}

	UpdateGeometry(axis);
}


void GridBCM::Render()
{
	glEnableClientState(GL_VERTEX_ARRAY);

	if( m_bbox ){
		if(m_bbox->GetVisible()){
			glColor4d(0.0, 0.0, 0.0, 1.0);
			glLineWidth(1.0);
			m_render->RenderLine( m_bbox );
		}
	}

	glLineWidth(1.0);
	for( int i = 0; i < 3; i++){
		if( !m_slice[i] ) continue;
		if( m_slice[i]->GetVisible()) {
			glColor4d(0.0, 0.0, 0.0, 1.0);
			m_render->RenderLine(m_slice[i]->GetBlockGeometry());
			if( m_showGrid) m_render->RenderLine(m_slice[i]->GetGridGeometry());
		}
	}

	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glEnable(GL_ALPHA_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_TEXTURE_2D);

	for(int i = 0; i < 3; i++){
		if( !m_slice[i] ) continue;
		if( m_slice[i]->GetVisible() ){
			glBindTexture(GL_TEXTURE_2D, m_slice[i]->GetCellTexture());
			m_render->RenderFace(m_slice[i]->GetCellGeometry());
		}
	}

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glDisable(GL_ALPHA_TEST);

	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glDisableClientState(GL_VERTEX_ARRAY);
}
