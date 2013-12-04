///
/// @file  LeafBlockLoader.h
/// @brief LeafBlockファイルを読み込むクラス
///

#ifndef __BCMTOOLS_LEAF_BLOCK_LOADER_H__
#define __BCMTOOLS_LEAF_BLOCK_LOADER_H__

#include <mpi.h>

#include "BCMFileCommon.h"
#include "PartitionMapper.h"
#include "BlockManager.h"
#include "Scalar3D.h"

namespace BCMFileIO {

	/// グリッドヘッダとデータを一括りにした構造体
	struct CellIDCapsule
	{
		LBCellIDHeader   header; ///< リーフブロックのグリッドヘッダ
		unsigned char*   data;   ///< リーフブロックデータ
		CellIDCapsule() : data(NULL){}
	};
	
	/// LeafBlockファイル(CellID)の読み込み (Gatherなし)
	/// 
	/// @param[in]  dir         入力元ディレクトリ
	/// @param[in]  ib          ブロック情報
	/// @param[in]  comm        MPIコミュニケータ
	/// @param[in]  pmapper     MxNデータマッパ
	/// @param[out] header      LeafBlockファイルヘッダ
	/// @param[out] cidCapsules 自プロセスが担当するデータのリスト
	/// @return 成功した場合true, 失敗した場合false
	/// 
	/// @note cidCapsulesに入る値はファイルから読み込んだ生の状態．
	///       圧縮符号やbitVoxelの展開，真に必要なデータの選択は
	///       ここでは実施しない
	/// 
	bool Load_LeafBlock_CellID( const std::string&          dir, 
	                            const IdxBlock*             ib, 
							    const MPI::Intracomm&       comm, 
							    PartitionMapper*            pmapper,
							    LBHeader&                   header, 
							    std::vector<CellIDCapsule>& cidCapsules );


	/// LeafBlockファイル(CellID)の読み込み (Gatherあり)
	/// 
	/// @param[in]  dir         入力元ディレクトリ
	/// @param[in]  ib          ブロック情報
	/// @param[in]  comm        MPIコミュニケータ
	/// @param[in]  pmapper     MxNデータマッパ
	/// @param[out] header      LeafBlockファイルヘッダ
	/// @param[out] cidCapsules 自プロセスが担当するデータのリスト
	/// @return 成功した場合true, 失敗した場合false
	/// 
	/// @note cidCapsulesに入る値はファイルから読み込んだ生の状態．
	///       圧縮符号やbitVoxelの展開，真に必要なデータの選択は
	///       ここでは実施しない
	/// 
	bool Load_LeafBlock_CellID_Gather( const std::string&          dir, 
	                                   const IdxBlock*             ib, 
									   const MPI::Intracomm&       comm, 
									   PartitionMapper*            pmapper,
							           LBHeader&                   header, 
									   std::vector<CellIDCapsule>& cidCapsules );
	/// 圧縮データの展開
	/// 
	/// @param[in] header     LeafBlockファイルヘッダ
	/// @param[in] cidCapsule CellIDカプセル
	/// @return 展開後のデータ (失敗した場合NULLを返す．)
	/// 
	/// @note cidCapsuleのdataはこの関数内で解放される．
	/// 
	unsigned char* DecompCellIDData( const LBHeader &header,  const CellIDCapsule& cidCapsule);


	/// LeafBlockファイル(Scalar)の読み込み
	/// 
	/// @param[in] comm           MPIコミュニケータ
	/// @param[in] ib             ブロック情報
	/// @param[in] blockManager   ブロックマネージャ
	/// @param[in] pmapper        MxNデータマッパ
	/// @param[in] vc             内部構造の仮想セルサイズ
	/// @param[in] step           読み込むデータのタイムステップインデックス番号
	///
	/// @return 成功した場合true, 失敗した場合false
	///
	/// @note このメソッドの実行により、LeafBlockファイルの中身がBlockManager配下のBlockに格納されます。
	///       仮想セルサイズは、ファイルに記載されている仮想セルサイズと関係なく指定できます。
	///       ブロック間の仮想セルの同期は行っていないため、ファイルの仮想セルサイズよりも大きい値を入れた
	///       場合、読めない値は0で埋まります。
	///
	bool Load_LeafBlock_Data(const MPI::Intracomm& comm,
							 const IdxBlock*       ib,
							 BlockManager&         blockManager,
							 PartitionMapper*      pmapper,
							 const int             vc,
							 const unsigned int    step );
	

	/// データバッファからBlockManager配下のBlockにデータをコピー
	///
	/// @param[in] blockManager ブロックマネージャ
	/// @param[in] dataClassID  データクラスID
	/// @param[in] blockID      ブロックID (プロセス内部でのブロック番号)
	/// @param[in] vc           仮想セルサイズ
	/// @param[in] buf          コピー元データバッファ
	///
	/// @return 成功した場合true, 失敗した場合false
	///
	template<typename T>
	bool CopyBufferToScalar3D(BlockManager& blockManager, const int dataClassID, const int blockID, const int vc, const T*  buf)
	{
		Vec3i size = blockManager.getSize();

		BlockBase* block = blockManager.getBlock(blockID);
		Scalar3D<T>* mesh = dynamic_cast< Scalar3D<T>* >(block->getDataClass(dataClassID));
		T* data      = mesh->getData();
		Index3DS idx = mesh->getIndex();

		for(int z = -vc; z < size.z + vc; z++){
			for(int y = -vc; y < size.y + vc; y++){
				for(int x = -vc; x < size.x + vc; x++){
					size_t loc = ( (x+vc) + ((y+vc) + (z+vc) * (size.y + (vc*2))) * (size.x + (vc*2)) );
					data[idx(x, y, z)] = buf[loc];
				}
			}
		}
		return true;
	}


} // namespace BCMFileIO

#endif// __BCMTOOLS_LEAF_BLOCK_LOADER_H__

