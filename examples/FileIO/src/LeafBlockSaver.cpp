///
/// @file  LeafBlockSaver.cpp
/// @brief LeafBlockファイルを出力する関数群
///

#include "LeafBlockSaver.h"
#include "BCMFileCommon.h"
#include "BitVoxel.h"
#include "RLE.h"
#include "ErrorUtil.h"

#include "FileSystemUtil.h"

#include "BlockManager.h"
#include "Scalar3D.h"

#include "type.h"

namespace BCMFileIO {

    bool Save_LeafBlock_CellID( const MPI::Intracomm& comm,
	                            const IdxBlock*       ib,
							    const Vec3i&          size,
							    const size_t          numBlock,
							    const unsigned char*  datas,
							    bool                  rle )
	{
		using namespace std;

		int rank = comm.Get_rank();

		// リーフブロックヘッダを準備
		LBHeader header;
		header.identifier = LEAFBLOCK_FILE_IDENTIFIER;
		header.kind       = static_cast<unsigned char>(ib->kind);
		header.dataType   = static_cast<unsigned char>(ib->dataType);
		header.bitWidth   = static_cast<unsigned short>(ib->bitWidth);
		header.vc         = ib->vc;
		header.size[0]    = size.x;
		header.size[1]    = size.y;
		header.size[2]    = size.z;
		header.numBlock   = 0; // temporary

		int vc = ib->vc;
		// 自プロセスの担当ブロックの保存用一時バッファサイズを計算
		const size_t tsz = (size.x + vc*2) * (size.y + vc*2) * (size.z + vc*2) * numBlock;

		// 自プロセスの担当ブロックをBitVoxel化
		size_t bitVoxelSize = 0;
		bitVoxelCell* bitVoxel = CompressBitVoxel(&bitVoxelSize, tsz, datas, ib->bitWidth);

		// リーフブロックのCellIDヘッダを準備
		LBCellIDHeader ch;
		ch.numBlock = numBlock;

		// BitVoxelバッファのサイズ(Byte単位)を計算
		size_t bvs = bitVoxelSize * sizeof(bitVoxelCell);
		unsigned char* rleBuf = NULL;
		unsigned char* dp     = NULL;

		if( rle ){ // RLE圧縮の場合
			size_t rleSize = 0;
			// BitVoxelをRLE圧縮
			rleBuf = rleEncode<bitVoxelCell, unsigned char>(bitVoxel, bvs, &rleSize);
			// リーフブロックのCellIDヘッダに圧縮サイズを記載
			ch.compSize = rleSize;
			dp = rleBuf;
		}else{ // RLEなしの場合
			// リーフブロックのCellIDヘッダに圧縮サイズ(=0)を記載
			ch.compSize = 0;
			dp = reinterpret_cast<unsigned char*>(bitVoxel);
		}
		
		if( ib->isGather ){ // GatherMode = "Gathered"
			int *numBlockTable      = NULL;
			int *leafBlockSizeTable = NULL;
			int bSz = ch.compSize == 0 ? bvs : static_cast<int>(ch.compSize);
			int nb  = static_cast<int>(ch.numBlock);
			if(rank == 0 ){
				numBlockTable      = new int[comm.Get_size()]; // 各ランクのブロック数取得バッファ
				leafBlockSizeTable = new int[comm.Get_size()]; // 各ランクのデータサイズ取得バッファ
			}

			// 各ランクのブロック数を集約
			comm.Gather(&bSz, 1, MPI::INT, leafBlockSizeTable, 1, MPI::INT, 0);
			// 各ランクのデータサイズを集約
			comm.Gather(&nb,  1, MPI::INT, numBlockTable,      1, MPI::INT, 0);
			
			unsigned char* rcvBuf = NULL;
			int *displs = NULL;
			
			uint64_t allSz = 0;
			if( rank == 0 ){
				displs = new int[comm.Get_size()];
				for(int i = 0; i < comm.Get_size(); i++){
					displs[i] = allSz;
					allSz += leafBlockSizeTable[i];
				}
				rcvBuf = new unsigned char[allSz];
			}
			// 各ランクの圧縮済みCellIDバッファを集約
			comm.Gatherv(dp, bSz, MPI::UNSIGNED_CHAR, rcvBuf, leafBlockSizeTable, displs, MPI::UNSIGNED_CHAR, 0);

			if( rank == 0 ){
				char filename[128];
				sprintf(filename, "%s.%s", ib->prefix.c_str(), ib->extension.c_str());
				string filepath = ib->rootDir + ib->dataDir + string(filename);
				FILE *fp = NULL;
				if( (fp = fopen(filepath.c_str(), "wb")) == NULL) {
					LogE("fileopen error <%s>. [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
					return false;
				}
				
				for(int i = 0; i < comm.Get_size(); i++){
					header.numBlock += numBlockTable[i];
				}
				// ヘッダを出力
				fwrite(&header, sizeof(header), 1, fp);

				uint64_t number_of_procs = comm.Get_size();
				// プロセス数を出力
				fwrite(&number_of_procs, sizeof(uint64_t), 1, fp);

				for(int i = 0; i < comm.Get_size(); i++){
					LBCellIDHeader ch;
					ch.numBlock = numBlockTable[i];
					ch.compSize = rle ? leafBlockSizeTable[i] : 0;
					// 圧縮サイズを出力
					fwrite(&ch, sizeof(LBCellIDHeader), 1, fp);
				}
				// 集約したCellIDバッファを出力
				fwrite(rcvBuf, sizeof(unsigned char), allSz, fp);

				fclose(fp);

			}
			// 各メモリの解放
			if( rank == 0){
				delete [] numBlockTable;
				delete [] leafBlockSizeTable;
				delete [] displs;
				delete [] rcvBuf;
			}

		}else{ // GatherMode = "Distributed"

			// write file
			char filename[128];
			sprintf(filename, "%s_%06d.%s", ib->prefix.c_str(), rank, ib->extension.c_str());
			string filepath = ib->rootDir + ib->dataDir + string(filename);
			FILE *fp = NULL;
			if( (fp = fopen(filepath.c_str(), "wb")) == NULL) {
				LogE("fileopen error <%s>. [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				return false;
			}
	
			header.numBlock = numBlock;
	
			fwrite(&header, sizeof(LBHeader),     1, fp);
			fwrite(&ch,     sizeof(LBCellIDHeader), 1, fp);

			size_t dsz = ch.compSize == 0 ? bvs : ch.compSize;
			
			fwrite(dp, sizeof(unsigned char), dsz, fp);

			fclose(fp);
		}

		delete [] bitVoxel;
		if( rleBuf != NULL ) delete [] rleBuf;
	
		return true;
	}


	/////////////////////////////////////////////////////////////////////////////////////////////
	template<typename T>
	bool CopyScalar3DToBuffer(BlockManager& blockManager, const int dataClassID, const int dataID, const int vc, T* buf)
	{
		Vec3i size = blockManager.getSize();
		
		BlockBase* block = blockManager.getBlock(dataID);
		Scalar3D<T>* mesh = dynamic_cast< Scalar3D<T>* >(block->getDataClass(dataClassID));
		T* data      = mesh->getData();
		Index3DS idx = mesh->getIndex();

		for(int z = -vc; z < size.z + vc; z++){
			for(int y = -vc; y < size.y + vc; y++){
				for(int x = -vc; x < size.x + vc; x++){
					size_t loc = ( (x+vc) + ((y+vc) + (z+vc) * (size.y + (vc*2))) * (size.x + (vc*2)) );
					buf[loc] = data[idx(x, y, z)];
				}
			}
		}
		return true;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	template<typename T>
	bool _Save_LeafBlock_Data(const MPI::Intracomm& comm,
							  const IdxBlock*       ib, 
							  BlockManager&         blockManager, 
							  const unsigned int    step)
	{
		using namespace std;
		int rank = comm.Get_rank();

		Vec3i size = blockManager.getSize();

		LBHeader header;
		header.identifier = LEAFBLOCK_FILE_IDENTIFIER;
		header.kind       = static_cast<unsigned char>(ib->kind);
		header.dataType   = static_cast<unsigned char>(ib->dataType);
		header.bitWidth   = static_cast<unsigned short>(ib->bitWidth);
		header.vc         = ib->vc;
		header.size[0]    = size.x;
		header.size[1]    = size.y;
		header.size[2]    = size.z;
		header.numBlock   = blockManager.getNumBlock();

		int vc = ib->vc;
		const size_t sz = (size.x + vc*2) * (size.y + vc*2) * (size.z + vc*2);
		
		string outputDir = ib->rootDir + ib->dataDir;
		if(ib->isStepSubDir){
			char stepDirName[128];
			sprintf(stepDirName, "%010d/", step);
			outputDir += std::string(stepDirName);
		}
		
		CreateDirectory(outputDir, outputDir.find("/") == 0 ? true : false);

		FILE *fp = NULL;
		char filename[128];
		sprintf(filename, "%s_%010d_%06d.%s", ib->prefix.c_str(), step, rank, ib->extension.c_str());

		string filepath = outputDir + string(filename);
		
		if( (fp = fopen(filepath.c_str(), "wb")) == NULL){
			LogE("fileopen err <%s>. [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
			return false;
		}

		fwrite(&header, sizeof(header), 1, fp);
		
		for(int id = 0; id < blockManager.getNumBlock(); ++id){
			for(int comp = 0; comp < static_cast<int>(ib->dataClassID.size()); comp++){
				T* buf = new T[sz];
				CopyScalar3DToBuffer(blockManager, ib->dataClassID[comp], id, vc, buf);
				fwrite(buf, sizeof(T), sz, fp);
				delete buf;
			}
		}

		fclose(fp);

		return true;
	}

	bool Save_LeafBlock_Data(const MPI::Intracomm& comm,
							 const IdxBlock*       ib, 
							 BlockManager&         blockManager, 
							 const unsigned int    step)
	{
		bool status = false;
		if     ( ib->dataType == LB_INT8   ) { status = _Save_LeafBlock_Data< s8>(comm, ib, blockManager, step); }
		else if( ib->dataType == LB_UINT8  ) { status = _Save_LeafBlock_Data< u8>(comm, ib, blockManager, step); }
		else if( ib->dataType == LB_INT16  ) { status = _Save_LeafBlock_Data<s16>(comm, ib, blockManager, step); }
		else if( ib->dataType == LB_UINT16 ) { status = _Save_LeafBlock_Data<u16>(comm, ib, blockManager, step); }
		else if( ib->dataType == LB_INT32  ) { status = _Save_LeafBlock_Data<s32>(comm, ib, blockManager, step); }
		else if( ib->dataType == LB_UINT32 ) { status = _Save_LeafBlock_Data<u32>(comm, ib, blockManager, step); }
		else if( ib->dataType == LB_INT64  ) { status = _Save_LeafBlock_Data<s64>(comm, ib, blockManager, step); }
		else if( ib->dataType == LB_UINT64 ) { status = _Save_LeafBlock_Data<u64>(comm, ib, blockManager, step); }
		else if( ib->dataType == LB_FLOAT32) { status = _Save_LeafBlock_Data<f32>(comm, ib, blockManager, step); }
		else if( ib->dataType == LB_FLOAT64) { status = _Save_LeafBlock_Data<f64>(comm, ib, blockManager, step); }
		else{
			LogE("invalid DataType (%d)[%s:%d]\n", ib->dataType, __FILE__, __LINE__);
			return false;
		}

		return status;
	}

} // BCMFIleIO
