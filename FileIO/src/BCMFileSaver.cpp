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
/// @file  BCMFileSaver.cpp
/// @brief BCMファイルを出力するクラス
/// 
#include "BCMOctree.h"
#include "RootGrid.h"
#include "BlockManager.h"
#include "Partition.h"

#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "ErrorUtil.h"
#include "FileSystemUtil.h"

#include "BCMFileCommon.h"
#include "BCMFileSaver.h"
#include "LeafBlockSaver.h"

#include "Scalar3D.h"
#include "Vector3D.h"

#define ENABLE_RLE_ENCODE

namespace BCMFileIO {

	unsigned char* GetCellIDBlock(const IdxBlock* ib, BlockManager& blockManager)
	{
		int vc = ib->vc;

		const Vec3i size = blockManager.getSize();
		const size_t numBlock = blockManager.getNumBlock();
		const size_t tsz = (size.x + vc*2) * (size.y + vc*2) * (size.z + vc*2) * numBlock;

		unsigned char* grids = new unsigned char[tsz];

		for(int id = 0; id < numBlock; ++id){
			BlockBase* block = blockManager.getBlock(id);

			Scalar3D<unsigned char>* mesh = dynamic_cast< Scalar3D<unsigned char>* >(block->getDataClass(ib->dataClassID[0]));
			unsigned char* data = mesh->getData();
			Index3DS idx = mesh->getIndex();

			Vec3i sz( size.x + vc*2, size.y + vc*2, size.z + vc*2);
			for(int z = 0; z < sz.z; z++){
				for(int y = 0; y < sz.y; y++){
					for(int x = 0; x < sz.x; x++){
						size_t loc = x + (y + (z + id * sz.z) * sz.y ) * sz.x;
						grids[loc] = data[ idx( x-vc, y-vc, z-vc ) ];
					}
				}
			}
		}
		return grids;
	}


	BCMFileSaver::BCMFileSaver( const Vec3r& globalOrigin, const Vec3r& globalRegion, const BCMOctree* octree, const std::string dir )
	 : m_blockManager(BlockManager::getInstance()), m_comm(m_blockManager.getCommunicator()),
	   m_octree(octree), m_globalOrigin(globalOrigin), m_globalRegion(globalRegion)
	{
		m_targetDir = FixDirectoryPath(dir);

		m_unit.length   = std::string("NonDimensional");
		m_unit.L0_scale = 1.0;
		m_unit.velocity = std::string("NonDimensional");
		m_unit.V0_scale = 1.0;
	}

	BCMFileSaver::~BCMFileSaver()
	{

	}
	
	bool BCMFileSaver::RegisterCellIDInformation( const int          dataClassID, 
                                                  const unsigned int bitWidth, 
												  const short        vc,
                                                  const std::string& name,
                                                  const std::string& prefix,
								                  const std::string& extension,
												  const std::string& dataDir,
								                  const bool         gather
												  )
	{
		if(dataClassID < 0){ return false; }

		if( findIdxBlock(m_idxBlockList, dataClassID) != NULL ){
			LogE("dataClassID(%s) is already registerd [%s:%d]\n", dataClassID, __FILE__, __LINE__);
			return false;
		}

		IdxBlock ib;
		ib.dataClassID.resize(1);
		ib.dataClassID[0] = dataClassID;
		ib.rootDir     = m_targetDir;
		ib.dataDir     = FixDirectoryPath(dataDir);
		ib.kind        = LB_CELLID;
		ib.dataType    = LB_UINT8;
		ib.bitWidth    = bitWidth;
		ib.vc          = vc;
		ib.name        = name;
		ib.prefix      = prefix;
		ib.extension   = extension;
		ib.isGather    = gather;

		m_idxBlockList.push_back(ib);

		return true;
	}

	bool BCMFileSaver::RegisterDataInformation( const int          *dataClassID,
	                                            const LB_KIND       kind,
										        const LB_DATA_TYPE  dataType,
												const short         vc,
										        const std::string&  name,
										        const std::string&  prefix,
										        const std::string&  extension,
										        const IdxStep&      step,
												const std::string&  dataDir,
												const bool          stepSubDir )
	{
		if(!dataClassID){ return false; }
		for(int i = 0; i < static_cast<int>(kind); i++){
			if(dataClassID[i] < 0){
				LogE("dataClassID(%d) is invalid) [%s:%d]\n", dataClassID[i], __FILE__, __LINE__);
				return false;
			}

			if( findIdxBlock(m_idxBlockList, dataClassID[i]) != NULL ){
				LogE("dataClassID(%d) is already registerd [%s:%d]\n", dataClassID[i], __FILE__, __LINE__);
				return false;
			}
		}


		unsigned int bitWidthTable[10] = {
			8, 8, 16, 16, 32, 32, 64, 64, 32, 64
		};


		IdxBlock ib;

		ib.dataClassID.resize(static_cast<int>(kind));
		for(int i = 0; i < static_cast<int>(kind); i++){
			ib.dataClassID[i]  = dataClassID[i];
		}

		ib.rootDir      = m_targetDir;
		ib.dataDir      = FixDirectoryPath(dataDir);
		ib.kind         = kind;
		ib.dataType     = dataType;
		ib.bitWidth     = bitWidthTable[(int)dataType];
		ib.vc           = vc;
		ib.name         = name;
		ib.prefix       = prefix;
		ib.extension    = extension;
		ib.isStepSubDir = stepSubDir;
		ib.step         = step;

		m_idxBlockList.push_back(ib);
		
		return true;
	}

	bool BCMFileSaver::SetUnit( const IdxUnit& unit )
	{
		m_unit = unit;
		return true;
	}


	bool BCMFileSaver::Save()
	{
		using namespace std;

		string octFilename("tree.oct");
		string octFilepath = m_targetDir + octFilename;
		//if( m_comm.Get_rank() == 0) LogI("Output Files are : IDX[%s], OCT[%s]\n", idxFilepath.c_str(), octFilepath.c_str());

		bool err = false;

		if( m_comm.Get_rank() == 0 ){
			err = !CreateDirectory(m_targetDir, m_targetDir.find("/") == 0 ? true : false);
		}else{
			err = false;
		}
		if( reduceError(err) ){
			return false;
		}

		err = reduceError( !SaveIndex(octFilename, m_octree->getNumLeafNode()) );
		if( err ){
			LogE("faild to save index file. [%s:%d]\n", __FILE__, __LINE__);
			return false;
		}

		err = reduceError( !SaveOctree(octFilepath, m_octree) );
		if( err ){
			LogE("faild to save octree file. [%s:%d]\n", __FILE__, __LINE__);
			return false;
		}

		return true;
	}


	bool BCMFileSaver::SaveLeafBlock(const char* name, unsigned int step )
	{
		using namespace std;

		bool err = false;

		IdxBlock *ib = findIdxBlock(m_idxBlockList, name);

		if( ib == NULL ){
			LogE("%s is not registerd. [%s:%d]\n", name, __FILE__, __LINE__);
			err = true;
		}

		if( reduceError(err) ){ return false; }
		err = false;

		if( ib->kind == LB_CELLID )
		{
			string lbdir = ib->rootDir + ib->dataDir;
			if( ib->isGather ){
				if( m_comm.Get_rank() == 0 ) err = !CreateDirectory(lbdir, lbdir.find("/") == 0 ? true : false);
			}else{
				err = !CreateDirectory(lbdir, lbdir.find("/") == 0 ? true : false);
			}

			if( reduceError(err) ){ 
				LogE("Cannot Create Output Directory (%s) [%s:%d]\n", lbdir.c_str(), __FILE__, __LINE__);
				return false;
			}

			unsigned char *data = GetCellIDBlock(ib, m_blockManager);
			if( data == NULL ){
				err = true;
			} else {
#ifdef ENABLE_RLE_ENCODE
				err = !Save_LeafBlock_CellID(m_comm, ib, m_blockManager.getSize(), m_blockManager.getNumBlock(), data, true);
#else
				err = !Save_LeafBlock_CellID(m_comm, ib, m_blockManager.getSize(), m_blockManager.getNumBlock(), data, false);
#endif // ENABLE_RLE_ENCODE
				delete [] data;
			}

			if( reduceError(err) ){
				LogE("Save Leaf Block (CellID) [%s:%d]\n", __FILE__, __LINE__);
				return false;
			}
		}
		else
		{
			err = !Save_LeafBlock_Data(m_comm, ib, m_blockManager, step);

			if( reduceError(err) ){
				LogE("Save Leaf Block (Scalar) [%s:%d]\n", __FILE__, __LINE__);
				return false;
			}
		}
		
		return true;

	}

	bool BCMFileSaver::SaveIndexProc(const std::string& filepath, const Partition& part)
	{
		using namespace std;

		// Gather Hostnames
		char hostname[MPI_MAX_PROCESSOR_NAME + 1] = {0};
		int nameLen;
		MPI_Get_processor_name(hostname, &nameLen);
		hostname[nameLen] = '\0';
		nameLen++;

		int* nameLenTable = NULL;
		if(m_comm.Get_rank() == 0 ){
			nameLenTable = new int[m_comm.Get_size()];
		}
		m_comm.Gather(&nameLen, 1, MPI::INT, nameLenTable, 1, MPI::INT, 0);


		char** hostnameList = NULL;
		if(m_comm.Get_rank() == 0){
			hostnameList = new char*[m_comm.Get_size()];
			hostnameList[0] = new char[nameLenTable[0]];
			memset(hostnameList[0], 0, sizeof(char) * nameLenTable[0]);
			memcpy(hostnameList[0], hostname, sizeof(char) * nameLen);

			for(int i = 1; i < m_comm.Get_size(); i++){
				hostnameList[i] = new char[nameLenTable[i]];
				memset(hostnameList[i], 0, sizeof(char) * nameLenTable[i]);
				m_comm.Recv(hostnameList[i], nameLenTable[i], MPI::CHAR, i, i);
			}
			delete [] nameLenTable;
		}else{
			m_comm.Send(hostname, nameLen, MPI::CHAR, 0, m_comm.Get_rank());
		}
		
		if( m_comm.Get_rank() != 0){
			return true;
		}
		
		stringstream os;
		os.setf(ios::scientific);
		os.precision(6);

		int  NumberOfRank  = 0;
		int  NumberOfGroup = 1;
		int  RankID        = 0;
		int  GroupID       = 0;

		MPI_Comm_size(MPI_COMM_WORLD, &NumberOfRank);
		MPI_Comm_rank(MPI_COMM_WORLD, &RankID);

		os << "MPI {" << endl;
		os << "  NumberOfRank   = " << NumberOfRank   << endl;
		os << "  NumberOfGroup  = " << NumberOfGroup  << endl;
		os << "  RankID         = " << RankID         << endl;
		os << "  GroupID        = " << GroupID        << endl;
		os << "}" << endl;
		os << endl;
		os << "Process { " << endl;
		for(int proc = 0; proc < m_comm.Get_size(); proc++){
			os << "  Rank[@] { " << endl;
			os << "    ID         = "   << proc << endl;
			os << "    HostName   = \"" << hostnameList[proc] << "\"" << endl;
			os << "    BlockRange = @range(" << part.getStart(proc) <<  "," << part.getEnd(proc) - 1 << ")" << endl;
			os << "  }" << endl << endl;
		}
		os << "}" << endl;

		ofstream ofs(filepath.c_str());
		if( !ofs ){
			LogE("failed to open file (%s) .[%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
			for(int i = 0; i < m_comm.Get_size(); i++) delete [] hostnameList[i];
			delete hostnameList;
			return false;
		}

		ofs << os.str(); 

		for(int i = 0; i < m_comm.Get_size(); i++) delete [] hostnameList[i];
		delete hostnameList;

		return true;
	}

	bool BCMFileSaver::SaveIndexCellID(const std::string& procName, const std::string& octName)
	{
		using namespace std;
		if( m_comm.Get_rank() != 0 ){ return true; }

		// search CellID IdxBlock
		IdxBlock *ib = NULL;
		for(vector<IdxBlock>::iterator it = m_idxBlockList.begin(); it != m_idxBlockList.end(); ++it){
			if( it->kind == LB_CELLID ){
				ib = &(*it);
				break;
			}
		}
		if( !ib ) { 
			return true;
		}

		Vec3i leafSize = m_blockManager.getSize();

		ostringstream os;
		os.setf(ios::scientific);
		os.precision(6);

		// write CellID Information
		os << "Domain { " << endl;
		os << "  GlobalOrigin = " << m_globalOrigin << endl;
		os << "  GlobalRegion = " << m_globalRegion << endl;
		os << "}" << endl;
		os << endl;
		os << "BCMTree {" << endl;
		os << "  TreeFile  = \"" << octName  << "\"" << endl;
		os << "  ProcFile  = \"" << procName << "\"" << endl;
		os << "}" << endl;
		os << endl;
		os << "LeafBlock {" << endl;
		os << "  size = " << leafSize << endl << endl;
		os << "  Unit {" << endl;
		os << "    Length = \"" << m_unit.length << "\"" << endl;
		os << "    L0     = " << m_unit.L0_scale << endl;
		os << "  }" << endl << endl;
		os << "  CellID { " << endl;
		os << "    name            = \"" << ib->name      << "\"" << endl;
		os << "    BitWidth        = "   << ib->bitWidth  << endl;
		os << "    VirtualCellSize = "   << ib->vc        << endl;
		os << "    DirectoryPath   = \"" << ib->dataDir   << "\"" << endl;
		os << "    Prefix          = \"" << ib->prefix    << "\"" << endl;
		os << "    Extension       = \"" << ib->extension << "\"" << endl;
		os << "    GatherMode      = \"" << (ib->isGather ? string("gathered") : string("distributed")) << "\"" << endl;
		os << "  }" << endl;
		os << "}" << endl;
		
		std::string filepath = m_targetDir + std::string("cellid.bcm");
		ofstream ofs( filepath.c_str() );
		if( !ofs ){
			LogE("failed to open file (%s) .[%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
			return false;
		}
		
		ofs << os.str();
		return true;
	}

	bool BCMFileSaver::SaveIndexData(const std::string& procName, const std::string& octName)
	{
		using namespace std;
		if( m_comm.Get_rank() != 0 ){ return true; }
		// search Data IdxBlocks
		vector<IdxBlock*> ibs;
		for(vector<IdxBlock>::iterator it = m_idxBlockList.begin(); it != m_idxBlockList.end(); ++it){
			if( it->kind != LB_CELLID ){
				ibs.push_back(&(*it));
			}
		}
		if( ibs.size() == 0) {
			return true;
		}

		Vec3i leafSize = m_blockManager.getSize();

		ostringstream os;
		os.setf(ios::scientific);
		os.precision(6);

		const char *typeStr[10] = {
			"Int8", "UInt8", "Int16", "UInt16", "Int32", "UInt32", "Int64", "UInt64", "Float32", "Float64"
		};

		// write Data Information
		os << "Domain { " << endl;
		os << "  GlobalOrigin = " << m_globalOrigin << endl;
		os << "  GlobalRegion = " << m_globalRegion << endl;
		os << "}" << endl;
		os << endl;
		os << "BCMTree { " << endl;
		os << "  TreeFile  = \"" << octName  << "\"" << endl;
		os << "  ProcFile  = \"" << procName << "\"" << endl;
		os << "}" << endl;
		os << endl;
		os << "LeafBlock {" << endl;
		os << "  size = " << leafSize << endl << endl;
		os << "  Unit { " << endl;
		os << "    Length   = \"" << m_unit.length << "\"" << endl;
		os << "    L0       = "   << m_unit.L0_scale << endl;
		os << "    Velocity = \"" << m_unit.velocity << "\"" << endl;
		os << "    V0       = "  << m_unit.V0_scale << endl;
		os << "  }" << endl << endl;
		for(vector<IdxBlock*>::iterator it = ibs.begin(); it != ibs.end(); ++it){
			os << "  Data[@] { " << endl;
			os << "    name               = \"" << (*it)->name                        << "\"" << endl;
			os << "    NumberOfComponents = "   << ((*it)->kind == LB_SCALAR ? 1 : 3) << endl;
			os << "    Type               = \"" << typeStr[(int)((*it)->dataType)]    << "\"" << endl;
			os << "    VirtualCellSize    = "   << (*it)->vc                          << endl;
			os << "    DirectoryPath      = \"" << (*it)->dataDir                     << "\"" << endl;
			os << "    Prefix             = \"" << (*it)->prefix                      << "\"" << endl;
			os << "    Extension          = \"" << (*it)->extension                   << "\"" << endl;
			os << "    StepSubDirectory   = \"" << ((*it)->isStepSubDir ? string("true") : string("false")) << "\"" << endl;

			os << endl;
			unsigned int stepRange[3] = { (*it)->step.GetRangeMin(), (*it)->step.GetRangeMax(), (*it)->step.GetRangeInterval() };
			const vector<unsigned int>& stepAdds = (*it)->step.GetAddStepList();
			const vector<unsigned int>& stepSubs = (*it)->step.GetSubStepList();

			os << "    Step {" << endl;
			os << "      base = @range(" << stepRange[0] << "," << stepRange[1] << "," << stepRange[2] << ")" << endl;
			if( stepAdds.size() != 0){
				os << "      Add  = @list(" << stepAdds[0];
				//for(vector<unsigned int>::const_iterator it = ++(stepAdds.begin()); it != stepAdds.end(); ++it){ os << "," << (*it); }
				for(size_t i = 1; i < stepAdds.size(); i++){ os << "," << stepAdds[i]; }
				os << ")" << endl;
			}
			if( stepSubs.size() != 0){
				os << "      Sub  = @list(" << stepSubs[0];
				//for(vector<unsigned int>::const_iterator it = ++(stepSubs.begin()); it != stepSubs.end(); ++it){ os << "," << (*it); }
				for(size_t i = 1; i < stepSubs.size(); i++){ os << "," << stepSubs[i]; }
				os << ")" << endl;
			}
			os << endl;
			os << "      Time   = " << (*it)->step.GetInitialTime() << endl;
			os << "      DeltaT = " << (*it)->step.GetDeltaT() << endl;
			os << "    }" <<endl;

			os << "  }" << endl << endl;
		}
		os << "}" << endl;

		std::string filepath = m_targetDir + std::string("data.bcm");
		ofstream ofs( filepath.c_str() );
		if( !ofs ){
			LogE("failed to open file (%s) .[%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
			return false;
		}
		
		ofs << os.str();
		return true;
	}


	bool BCMFileSaver::SaveIndex(const std::string& octName, const int numLeaf)
	{
		using namespace std;

		if( numLeaf < m_comm.Get_size() ){
			LogE("less number of leafs than number of procs. [%s:%d]\n", __FILE__, __LINE__);
			return false;
		}
		
		ostringstream os;
		os.setf(ios::scientific);
		os.precision(6);

		Partition part(m_comm.Get_size(), numLeaf);

		std::string procName("proc.bcm");
		std::string procPath = m_targetDir + procName;
		if( reduceError(!SaveIndexProc(procPath, part)) )      { LogE("%s:%s:%d\n", __func__, __FILE__, __LINE__); return false; }
		if( reduceError(!SaveIndexCellID(procName, octName)) ) { LogE("%s:%s:%d\n", __func__, __FILE__, __LINE__); return false; }
		if( reduceError(!SaveIndexData(procName, octName)) )   { LogE("%s:%s:%d\n", __func__, __FILE__, __LINE__); return false; }

		return true;
	}
	
	bool BCMFileSaver::SaveOctree( const std::string& filepath, const BCMOctree* octree)
	{
		using namespace std;

		if( m_comm.Get_rank() != 0 ){
			return true;
		}
		
		FILE *fp = NULL;
		if( (fp = fopen(filepath.c_str(), "wb")) == NULL ){
			LogE("faild open file (%s) .[%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
			return false;
		}

		OctHeader header;
		header.identifier = OCTREE_FILE_IDENTIFIER;

		const RootGrid* rootGrid = octree->getRootGrid();

		header.org[0]      = m_globalOrigin.x;
		header.org[1]      = m_globalOrigin.y;
		header.org[2]      = m_globalOrigin.z;
		header.rgn[0]      = m_globalRegion.x;
		header.rgn[1]      = m_globalRegion.y;
		header.rgn[2]      = m_globalRegion.z;
		header.rootDims[0] = rootGrid->getSizeX();
		header.rootDims[1] = rootGrid->getSizeY();
		header.rootDims[2] = rootGrid->getSizeZ();

		header.numLeaf = octree->getNumLeafNode();

		const vector<Node*> nodes  = octree->getLeafNodeArray();
		vector<Pedigree>    pedigs;
		pedigs.reserve(octree->getNumLeafNode());
		
		int maxLevel = 0;
		for(vector<Node*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it){
			const Node* n = *it;
			maxLevel = n->getLevel() > maxLevel ? n->getLevel() : maxLevel;
			pedigs.push_back(n->getPedigree());
		}
		
		header.maxLevel = maxLevel;

		fwrite(&header, sizeof(header), 1, fp);
		fwrite(&pedigs[0], sizeof(Pedigree), pedigs.size(), fp);

		fclose(fp);

		return true;
	}


} // namespace BCMFileIO

