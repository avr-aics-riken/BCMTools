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
/// @file  LeafBlockSaver.h
/// @brief LeafBlockファイルを出力する関数群
///

#ifndef __BCMTOOLS_LEAFBLOCK_SAVER_H__
#define __BCMTOOLS_LEAFBLOCK_SAVER_H__

#include <mpi.h>

#include "BCMFileCommon.h"
#include "Vec3.h"

class BlockManager;

namespace BCMFileIO {
	
	/// LeafBlockファイル(CellID)の出力
	/// 
	/// @param[in] comm      MPIコミュニケータ
	/// @param[in] ib        ブロック情報
	/// @param[in] size      リーフブロックサイズ
	/// @param[in] numBlock  総ブロック数
	/// @param[in] datas     CellIDが格納されたデータバッファ
	/// @param[in] rle       RLE圧縮フラグ (trueの場合RLE圧縮を行う)
	///
	/// @return 成功した場合true, 失敗した場合false
	/// 
    bool Save_LeafBlock_CellID( const MPI::Intracomm& comm,
	                            const IdxBlock*       ib,
							    const Vec3i&          size,
							    const size_t          numBlock,
							    const unsigned char*  datas,
							    bool                  rle);


	/// LeafBlockファイル(物理量)の出力
	/// 
	/// @param[in] comm         MPIコミュニケータ
	/// @param[in] ib           ブロック情報
	/// @param[in] blockManager ブロックマネージャ
	/// @param[in] step         出力タイムステップのインデックス番号
	/// 
	/// @return 成功した場合true, 失敗した場合false
	///
	bool Save_LeafBlock_Data(const MPI::Intracomm& comm,
							 const IdxBlock*       ib, 
							 BlockManager&         blockManager, 
							 const unsigned int    step);
	
} // BCMFileIO

#endif // __BCMTOOLS_LEAFBLOCK_SAVER_H__

