/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file  LeafBlockLoader.h
/// @brief LeafBlockファイルを読み込むクラス
///

#ifndef __BCMTOOLS_LEAF_BLOCK_LOADER_H__
#define __BCMTOOLS_LEAF_BLOCK_LOADER_H__

//#include <mpi.h>

#include "BCMFileCommon.h"
#include "PartitionMapper.h"

namespace BCMFileIO {

	/// グリッドヘッダとデータを一括りにした構造体
	struct CellIDCapsule
	{
		LBCellIDHeader   header; ///< リーフブロックのグリッドヘッダ
		unsigned char*   data;   ///< リーフブロックデータ
		CellIDCapsule() : data(NULL){}
	};
	
	/// LeafBlockファイルの読み込み (Gatherなし)
	/// 
	/// @param[in]  dir         入力元ディレクトリ
	/// @param[in]  ib          ブロック情報
	/// @param[in]  comm        MPIコミュニケータ
	/// @param[in]  pmapper     MxNデータマッパ
	/// @param[out] header      LeafBlockファイルヘッダ
	/// @param[out] cidCapsules 自プロセスが担当するデータのリスト
	/// @return 成功した場合true, 失敗した場合false
	/// 
	/// @note gridCapsulesに入る値はファイルから読み込んだ生の状態．
	///       圧縮符号やbitVoxelの展開，真に必要なデータの選択は
	///       ここでは実施しない
	/// 
	bool Load_LeafBlock_CellID( const std::string&          dir, 
	                            const IdxBlock*             ib, 
							    //const MPI::Intracomm&       comm, 
							    PartitionMapper*            pmapper,
							    LBHeader&                   header, 
							    std::vector<CellIDCapsule>& cidCapsules );


	/// LeafBlockファイルの読み込み (Gatherあり)
	/// 
	/// @param[in]  dir         入力元ディレクトリ
	/// @param[in]  ib          ブロック情報
	/// @param[in]  comm        MPIコミュニケータ
	/// @param[in]  pmapper     MxNデータマッパ
	/// @param[out] header      LeafBlockファイルヘッダ
	/// @param[out] cidCapsules 自プロセスが担当するデータのリスト
	/// @return 成功した場合true, 失敗した場合false
	/// 
	/// @note gridCapsulesに入る値はファイルから読み込んだ生の状態．
	///       圧縮符号やbitVoxelの展開，真に必要なデータの選択は
	///       ここでは実施しない
	/// 
	bool Load_LeafBlock_CellID_Gather( const std::string&          dir, 
	                                   const IdxBlock*             ib, 
									   //const MPI::Intracomm&       comm, 
									   PartitionMapper*            pmapper,
							           LBHeader&                   header, 
									   std::vector<CellIDCapsule>& cidCapsules );
	/// 圧縮データの展開
	/// 
	/// @param[in] header     LeafBlockファイルヘッダ
	/// @param[in] cidCapsule CellIDカプセル
	/// @return 展開後のデータ (失敗した場合NULLを返す．)
	/// 
	/// @note gridCapsuleのdataはこの関数内で解放される．
	/// 
	unsigned char* DecompCellIDData( const LBHeader &header,  const CellIDCapsule& cidCapsule);

} // namespace BCMFileIO

#endif// __BCMTOOLS_LEAF_BLOCK_LOADER_H__

