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
/// @file  LeafBlockLoader.cpp
/// @brief LeafBlockファイルを読み込むクラス
///

#include <vector>
#include <string>

#include "LeafBlockLoader.h"
#include "BitVoxel.h"
#include "RLE.h"
#include "ErrorUtil.h"

#include "type.h"

namespace BCMFileIO {
	inline void DUMMY(void*){}

	inline size_t GetBitVoxelSize( const LBHeader& hdr, size_t numBlocks ){
		size_t blockSize = (hdr.size[0] + hdr.vc * 2) * (hdr.size[1] + hdr.vc * 2) * (hdr.size[2] + hdr.vc * 2);
		return GetBitVoxelSize(blockSize * numBlocks, hdr.bitWidth);
	}

	inline bool Load_LeafBlock_Header( FILE *fp, LBHeader& hdr, bool& isNeedSwap)
	{
		fread(&hdr, sizeof(LBHeader), 1, fp);

		if( hdr.identifier != LEAFBLOCK_FILE_IDENTIFIER ){
			BSwap32(&hdr.identifier);

			if( hdr.identifier != LEAFBLOCK_FILE_IDENTIFIER ){
				return false;
			}

			isNeedSwap = true;

			BSwap16(&hdr.bitWidth);
			BSwap32(&hdr.vc);
			BSwap32(&hdr.size[0]);
			BSwap32(&hdr.size[1]);
			BSwap32(&hdr.size[2]);
			BSwap64(&hdr.numBlock);
		}

		//LogI("Header [bw : %d, vc : %d, sz : (%3d, %3d, %3d), nb : %d\n", 
		//       hdr.bitWidth, hdr.vc, hdr.size[0], hdr.size[1], hdr.size[2], hdr.numBlock);

		return true;
	}

	inline bool Load_LeafBlock_CellIDHeader( FILE *fp, LBCellIDHeader& chdr, const bool isNeedSwap )
	{
		fread(&chdr, sizeof(LBCellIDHeader), 1, fp);
		if( isNeedSwap ){
			BSwap64(&chdr.numBlock);
			BSwap64(&chdr.compSize);
		}
		if(chdr.compSize != 0 && (chdr.compSize % sizeof(GridRleCode)) != 0){
			LogE("compress size is invalid\n");
			return false;
		}
		return true;
	}
	
	inline bool Load_LeafBlock_CellIDData( FILE *fp, unsigned char** data, const LBHeader& hdr, const LBCellIDHeader& chdr, const bool isNeedSwap)
	{
		size_t sz = 0;
		if( chdr.compSize == 0 ){
			sz = GetBitVoxelSize(hdr, chdr.numBlock) * sizeof(bitVoxelCell);
		}else{
			sz = chdr.compSize;
		}

		*data = new unsigned char[sz];

		fread(*data, sizeof(unsigned char), sz, fp);

		if( isNeedSwap ){
			if( chdr.compSize == 0 ){
				size_t bitVoxelSize = GetBitVoxelSize(hdr, chdr.numBlock);
				bitVoxelCell* bitVoxel = reinterpret_cast<bitVoxelCell*>(*data);
				for(int i = 0; i < bitVoxelSize; i++){
					BSwap32(&bitVoxel[i]);
				}
			}else{
				GridRleCode* p = reinterpret_cast<GridRleCode*>(*data);
				for(int i = 0; i < sz / sizeof(GridRleCode); i++){
					BSwap32(&(p[i].c));
				}
			}
		}
		
		return true;
	}

	bool Load_LeafBlock_CellID( const std::string&          dir, 
	                            const IdxBlock*             ib, 
							    const MPI::Intracomm&       comm, 
							    PartitionMapper*            pmapper,
							    LBHeader&                   header, 
							    std::vector<CellIDCapsule>& cidCapsules )
	{
		using namespace std;

		cidCapsules.clear();
		
		vector<PartitionMapper::FDIDList> fdidlists;
		pmapper->GetFDIDLists(comm.Get_rank(), fdidlists);

		cidCapsules.reserve(fdidlists.size());
		
		bool err = false;

		for(vector<PartitionMapper::FDIDList>::const_iterator file = fdidlists.begin(); file != fdidlists.end(); ++file){
			char filename[128];
			sprintf(filename, "%s_%06d.%s", ib->prefix.c_str(), file->FID, ib->extension.c_str());
			string filepath = dir + string(filename);

			FILE *fp = NULL;
			if( (fp = fopen(filepath.c_str(), "rb")) == NULL ){
				LogE("file open error \"%s\" [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				err = true;
				break;
			}

			bool isNeedSwap = false;
			
			LBHeader hdr;

			if( !Load_LeafBlock_Header(fp, hdr, isNeedSwap) ){
				LogE("%s is not leafBlock file [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				fclose(fp); err = true; break;
			}

			if(hdr.kind != static_cast<unsigned char>(LB_CELLID)){
				LogE("%s is not CellID file [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				fclose(fp); err = true; break;
			}
			if(hdr.bitWidth < 1 && header.bitWidth > 5) {
				LogE("%s is not CEllID file [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				fclose(fp); err = true; break;
			}

			CellIDCapsule cc;
			if( !Load_LeafBlock_CellIDHeader(fp, cc.header, isNeedSwap) ){
				fclose(fp); err = true; break;
			}
			
			Load_LeafBlock_CellIDData(fp, &cc.data, hdr, cc.header, isNeedSwap);

			cidCapsules.push_back(cc);
			header = hdr;
			fclose(fp);
		}

		if( err ){
			LogE("Clear cidCapsules [%s:%d]\n", __FILE__, __LINE__);
			// データロードに失敗したためロード済みのデータを破棄
			for(vector<CellIDCapsule>::iterator it = cidCapsules.begin(); it != cidCapsules.end(); ++it){
				if( it->data != NULL) delete [] it->data;
			}
			cidCapsules.clear();
			return false;
		}

		return true;
	}

	bool Load_LeafBlock_CellID_Gather( const std::string&        dir, 
	                                   const IdxBlock*           ib, 
									   const MPI::Intracomm&     comm, 
									   PartitionMapper*          pmapper,
							           LBHeader&                 header, 
									   std::vector<CellIDCapsule>& cidCapsules )
	{
		using namespace std;

		cidCapsules.clear();

		int rank = comm.Get_rank();

		// rank 0でファイルをロード
		if( rank == 0 ){
			// ファイルロード
			char filename[128];
			sprintf(filename, "%s.%s", ib->prefix.c_str(), ib->extension.c_str());
			string filepath = dir + string(filename);

			FILE *fp = NULL;
			if( (fp = fopen(filepath.c_str(), "rb")) == NULL ){
				printf("err : file open error \"%s\" [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				
				// ファイルロードエラーを全プロセスに通知
				unsigned char loadError = 1; comm.Bcast(&loadError, 1, MPI::CHAR, 0);
				return false;
			}

			bool isNeedSwap = false;
			
			LBHeader hdr;

			if( !Load_LeafBlock_Header(fp, hdr, isNeedSwap) ){
				LogE("%s is not leafBlock file [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				// ファイルロードエラーを全プロセスに通知
				unsigned char loadError = 1; comm.Bcast(&loadError, 1, MPI::CHAR, 0);
				fclose(fp); return false;
			}

			if(hdr.kind != static_cast<unsigned char>(LB_CELLID)){
				LogE("%s is not Grid file [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				// ファイルロードエラーを全プロセスに通知
				unsigned char loadError = 1; comm.Bcast(&loadError, 1, MPI::CHAR, 0);
				fclose(fp); return false;
			}

			if(hdr.bitWidth < 1 && header.bitWidth > 5) {
				LogE("%s is not Grid file [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				// ファイルロードエラーを全プロセスに通知
				unsigned char loadError = 1; comm.Bcast(&loadError, 1, MPI::CHAR, 0);
				fclose(fp); return false;
			}

			header = hdr;

			uint64_t wnp = 0;
			fread(&wnp, sizeof(uint64_t), 1, fp);
			if( isNeedSwap ){
				BSwap64(&wnp);
			}

			if( wnp != pmapper->GetWriteProcs() ){
				printf("err : %s's write procs is invalid  [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				// ファイルロードエラーを全プロセスに通知
				unsigned char loadError = 1; comm.Bcast(&loadError, 1, MPI::CHAR, 0);
				fclose(fp); return false;
			}

			vector<LBCellIDHeader> chs(wnp);
			for(int i = 0; i < wnp; i++){
				if( !Load_LeafBlock_CellIDHeader(fp, chs[i], isNeedSwap) ){
					// ファイルロードエラーを全プロセスに通知
					unsigned char loadError = 1; comm.Bcast(&loadError, 1, MPI::CHAR, 0);
					fclose(fp); return false;
				}
			}

			vector<unsigned char*> contents(wnp);

			for(int i = 0; i < wnp; i++){
				contents[i] = NULL;
				Load_LeafBlock_CellIDData( fp, &contents[i], hdr, chs[i], isNeedSwap );
			}

			fclose(fp);
			
			// ファイルロード完了を全プロセスに通知	
			unsigned char loadError = 0; comm.Bcast(&loadError, 1, MPI::CHAR, 0);

			// ヘッダ情報をブロードキャスト	
			comm.Bcast(&header, sizeof(LBHeader), MPI::CHAR, 0);
			// Gridヘッダ情報をブロードキャスト
			comm.Bcast(&chs[0], wnp * sizeof(LBCellIDHeader), MPI::CHAR, 0);
			
			// 各計算ノードにデータを送信
			for(int i = 1; i < comm.Get_size(); i++){
				vector<PartitionMapper::FDIDList> fdidlists;
				pmapper->GetFDIDLists(i, fdidlists);
				for(vector<PartitionMapper::FDIDList>::iterator file = fdidlists.begin(); file != fdidlists.end(); ++file){
					size_t sz = 0;
					if( chs[file->FID].compSize == 0){
						sz = GetBitVoxelSize(header, chs[file->FID].numBlock) * sizeof(bitVoxelCell);
					}else{
						sz = chs[file->FID].compSize;
					}

					// 送信
					int tag = file->FID;
					comm.Send(contents[file->FID], sz, MPI::CHAR, i, tag);
				}
			}
			
			// Rank 0用のデータコピー

			bool freeMask[wnp];
			for( int i = 0; i < wnp; i++){ freeMask[i] = true; }

			vector<PartitionMapper::FDIDList> fdidlists;
			pmapper->GetFDIDLists(rank, fdidlists);
			for(vector<PartitionMapper::FDIDList>::iterator file = fdidlists.begin(); file != fdidlists.end(); ++file){
				CellIDCapsule cc;
				cc.header = chs[file->FID];
				cc.data   = contents[file->FID];
				cidCapsules.push_back(cc);
				freeMask[file->FID] = false;
			}
			
			// 不要なデータを解放
			for(int i = 0; i < wnp; i++){
				if(freeMask[i]) delete [] contents[i];
			}
		}
		// rank 0 以外
		else
		{
			// rank 0からの通知を受け取る
			unsigned char loadError = 0; comm.Bcast(&loadError, 1, MPI::CHAR, 0);
			if(loadError == 1){
				return false;
			}

			// ヘッダ情報を取得
			comm.Bcast(&header, sizeof(LBHeader), MPI::CHAR, 0);

			// Gridヘッダ情報を取得
			int wnp = pmapper->GetWriteProcs();
			vector<LBCellIDHeader> chs(wnp);
			comm.Bcast(&chs[0], wnp * sizeof(LBCellIDHeader), MPI::CHAR, 0);
			
			vector<PartitionMapper::FDIDList> fdidlists;
			pmapper->GetFDIDLists(rank, fdidlists);
			for(vector<PartitionMapper::FDIDList>::iterator file = fdidlists.begin(); file != fdidlists.end(); ++file){
				CellIDCapsule cc;
				cc.header = chs[file->FID];
				size_t sz = 0;
				if( chs[file->FID].compSize == 0){
					sz = GetBitVoxelSize(header, chs[file->FID].numBlock) * sizeof(bitVoxelCell);
				}else{
					sz = chs[file->FID].compSize;
				}
				cc.data = new unsigned char[sz];
				
				// 受信
				int tag = file->FID;
				comm.Recv(cc.data, sz, MPI::CHAR, 0, tag);

				cidCapsules.push_back(cc);
			}
		}
		return true;
	}

	// CellIDCapsuleのdataはここで解放されます
	unsigned char* DecompCellIDData( const LBHeader &header,  const CellIDCapsule& cc)
	{
		unsigned char* ret = NULL;

		size_t blockSize = (header.size[0] + header.vc*2) * (header.size[1] + header.vc*2) * (header.size[2] + header.vc*2);
		size_t dataSize  = blockSize * cc.header.numBlock;

		bitVoxelCell* bitVoxel = NULL;
		if( cc.header.compSize != 0){
			size_t bitVoxelSize = GetBitVoxelSize(dataSize, header.bitWidth);
			size_t dsize = bitVoxelSize * sizeof(bitVoxelCell);
			bitVoxel = rleDecode<bitVoxelCell, unsigned char>(cc.data, cc.header.compSize, dsize);
			delete [] cc.data;
		}else{
			bitVoxel = reinterpret_cast<bitVoxelCell*>(cc.data);
		}

		ret = DecompressBitVoxel(dataSize, bitVoxel, header.bitWidth);

		delete [] bitVoxel;

		return ret;
	}

	////////////////////////////////////////////////////////////////////////
	
	unsigned char* Load_BlockContents(FILE *fp, const LBHeader& hdr, const Vec3i& bsz, const int vc, const bool isNeedSwap)
	{
		const static size_t typeByteTable[10] = { 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };
		static void (*BSwap[10])(void*) = { DUMMY, DUMMY, BSwap16, BSwap16, BSwap32, BSwap32, BSwap64, BSwap64, BSwap32, BSwap64 };

		size_t typeByte = typeByteTable[hdr.dataType];
		
		Vec3i ibsz( bsz.x + vc*2,     bsz.y + vc*2,     bsz.z + vc*2);
		Vec3i fbsz( bsz.x + hdr.vc*2, bsz.y + hdr.vc*2, bsz.z + hdr.vc*2);

		unsigned char* block = new unsigned char[ibsz.x * ibsz.y * ibsz.z * typeByte];
		unsigned char* buf   = new unsigned char[fbsz.x * fbsz.y * fbsz.z * typeByte];

		memset(block, 0, sizeof(unsigned char) * ibsz.x * ibsz.y * ibsz.z * typeByte);

		fread(buf, sizeof(unsigned char), fbsz.x * fbsz.y * fbsz.z * typeByte, fp);

		if( isNeedSwap ){
			for(int i = 0; i < fbsz.x * fbsz.y * fbsz.z; i++){
				BSwap[hdr.dataType](&buf[i * typeByte]);
			}
		}

		if( vc > hdr.vc ){
			unsigned int vcd = vc - hdr.vc;
			for(int z = 0; z < fbsz.z; z++){
				for(int y = 0; y < fbsz.y; y++){
					size_t ibloc = 0 + vcd + ( (y + vcd) + (z + vcd) * ibsz.y ) * ibsz.x;
					size_t fbloc = 0 +     + (  y        +  z        * fbsz.y ) * fbsz.x;
					memcpy(&block[ibloc * typeByte], &buf[fbloc * typeByte], typeByte * fbsz.x );
				}
			}
		}else{
			unsigned int vcd = hdr.vc - vc;
			for(int z = 0; z < ibsz.z; z++){
				for(int y = 0; y < ibsz.y; y++){
					size_t ibloc = 0 +     + (  y        +  z        * ibsz.y ) * ibsz.x;
					size_t fbloc = 0 + vcd + ( (y + vcd) + (z + vcd) * fbsz.y ) * fbsz.x;
					memcpy(&block[ibloc * typeByte], &buf[fbloc * typeByte], typeByte * ibsz.x );
				}
			}
		}
		
		delete [] buf;
		return block;
	}

	
	bool Load_LeafBlock_Data(const MPI::Intracomm& comm,
							 const IdxBlock*       ib,
							 BlockManager&         blockManager,
							 PartitionMapper*      pmapper,
							 const int             vc,
							 const unsigned int    step )
	{
		using namespace std;
		vector<PartitionMapper::FDIDList> fdidlists;
		pmapper->GetFDIDLists(comm.Get_rank(), fdidlists);

		Vec3i bsz = blockManager.getSize();
		
		int did = 0;
		for(vector<PartitionMapper::FDIDList>::iterator file = fdidlists.begin(); file != fdidlists.end(); ++file){
			char filename[128];
			sprintf(filename, "%s_%010d_%06d.%s", ib->prefix.c_str(), step, file->FID, ib->extension.c_str());

			string dirpath = ib->rootDir + ib->dataDir;
			if(ib->isStepSubDir){
				char stepDirName[128];
				sprintf(stepDirName, "%010d/", step);
				dirpath += string(stepDirName);
			}

			string filepath = dirpath + string(filename);
			
			FILE *fp = NULL;
			if( (fp = fopen(filepath.c_str(), "rb")) == NULL ) {
				LogE("Cannnot open file (%s) [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				return false;
			}
			
			bool isNeedSwap = false;
			LBHeader hdr;

			if( !Load_LeafBlock_Header(fp, hdr, isNeedSwap) ){
				LogE("%s is not leafBlock file [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
				return false;
			}

			if(hdr.kind != static_cast<unsigned char>(ib->kind) ){
				LogE("%s's kind(%d) is not corresponds IndexFile(%d) [%s:%d]\n", filepath.c_str(), hdr.kind, ib->kind, __FILE__, __LINE__);
				return false;
			}

			if(hdr.dataType != static_cast<unsigned char>(ib->dataType)){
				LogE("%s's Type(%d) is not corresponds IndexFile(%d). [%s:%d]\n", filepath.c_str(), hdr.dataType, ib->dataType, __FILE__, __LINE__);
				return false;
			}

			if(hdr.vc != ib->vc ){
				LogE("%s's vc(%d) is not corresponds IndexFile(%d). [%s:%d]\n", filepath.c_str(), hdr.vc, ib->vc, __FILE__, __LINE__);
				return false;
			}

			if(hdr.size[0] != bsz.x || hdr.size[1] != bsz.y || hdr.size[2] != bsz.z){
				LogE("%s's size(%3d, %3d, %3d) is not corresponds IndexFile(%3d, %3d, %3d). [%s:%d]\n",
				     filepath.c_str(), hdr.size[0], hdr.size[1], hdr.size[2], bsz.x, bsz.y, bsz.z, __FILE__, __LINE__);
				return false;
			}
			
			const static size_t typeByteTable[10] = { 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };
			size_t typeByte = typeByteTable[hdr.dataType];
			Vec3i fbsz( bsz.x + hdr.vc*2, bsz.y + hdr.vc*2, bsz.z + hdr.vc*2);

			fseek(fp, file->FDIDs[0] * typeByte * (fbsz.x * fbsz.y * fbsz.z) * static_cast<size_t>(ib->kind), SEEK_CUR);

			for(vector<int>::iterator fdid = file->FDIDs.begin(); fdid != file->FDIDs.end(); ++fdid){
				for(int i = 0; i < static_cast<int>(ib->kind); i++){
					int dcid = ib->dataClassID[i];
					unsigned char* block = Load_BlockContents(fp, hdr, bsz, vc, isNeedSwap);

					if     (ib->dataType == LB_INT8   ){ CopyBufferToScalar3D(blockManager, dcid, did, vc, reinterpret_cast< s8*>(block)); }
					else if(ib->dataType == LB_UINT8  ){ CopyBufferToScalar3D(blockManager, dcid, did, vc, reinterpret_cast< u8*>(block)); }
					else if(ib->dataType == LB_INT16  ){ CopyBufferToScalar3D(blockManager, dcid, did, vc, reinterpret_cast<s16*>(block)); }
					else if(ib->dataType == LB_UINT16 ){ CopyBufferToScalar3D(blockManager, dcid, did, vc, reinterpret_cast<u16*>(block)); }
					else if(ib->dataType == LB_INT32  ){ CopyBufferToScalar3D(blockManager, dcid, did, vc, reinterpret_cast<s32*>(block)); }
					else if(ib->dataType == LB_UINT32 ){ CopyBufferToScalar3D(blockManager, dcid, did, vc, reinterpret_cast<u32*>(block)); }
					else if(ib->dataType == LB_INT64  ){ CopyBufferToScalar3D(blockManager, dcid, did, vc, reinterpret_cast<s64*>(block)); }
					else if(ib->dataType == LB_UINT64 ){ CopyBufferToScalar3D(blockManager, dcid, did, vc, reinterpret_cast<u64*>(block)); }
					else if(ib->dataType == LB_FLOAT32){ CopyBufferToScalar3D(blockManager, dcid, did, vc, reinterpret_cast<f32*>(block)); }
					else if(ib->dataType == LB_FLOAT64){ CopyBufferToScalar3D(blockManager, dcid, did, vc, reinterpret_cast<f64*>(block)); }

					delete [] block;
				}

				did++;
			}

			fclose(fp);
		}

		return true;
	}


} // BCMFileIO

