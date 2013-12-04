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
/// @file  BCMFileSaver.h
/// @brief BCMファイルを出力するクラス
/// 


#ifndef __BCMTOOLS_BCM_FILE_SAVER_H__
#define __BCMTOOLS_BCM_FILE_SAVER_H__

#include <mpi.h>
#include <string>
#include <vector>

#include <Vec3.h>

#include "BCMFileCommon.h"

class BCMOctree;
class RootGrid;
class BlockManager;

namespace BCMFileIO {

	/// BCMファイルを出力するクラス
	class BCMFileSaver
	{
	public:

		/// コンストラクタ
		/// 
		/// @param[in] globalOrigin 計算空間全体の起点座標
		/// @param[in] globalRegion 計算空間全体の領域サイズ
		/// @param[in] octree 出力Octree
		/// @param[in] dir ファイル出力先ディレクトリ(省略した場合、カレントディレクトリ)
		///
		BCMFileSaver( const Vec3r& globalOrigin, const Vec3r& globalRegion, const BCMOctree* octree, const std::string dir = std::string("") );

		/// デストラクタ
		~BCMFileSaver();
	
		/// 出力対象のリーフブロック情報を登録 (CellID用)
		/// 
		/// @param[in] dataClassID データクラスID
		/// @param[in] bitWidth    CellIDの表現ビット幅
		/// @param[in] vc          仮想セルサイズ
		/// @param[in] name        系の名称
		/// @param[in] prefix      リーフブロックファイルのPrefix
		/// @param[in] extension   リーフブロックファイルの拡張子
		/// @param[in] dataDir     リーフブロックファイルの出力ディレクトリを指定 (コンストラクタで指定した出力ディレクトリからの相対パス)
		/// @param[in] gatherMode  集約モード (trueの場合、Rank 0に集約)
		/// @return 成功した場合true, 失敗した場合false
		///
		bool RegisterCellIDInformation( const int          dataClassID, 
	                                    const unsigned int bitWidth, 
										const short        vc,
	                                    const std::string& name,
	                                    const std::string& prefix,
									    const std::string& extension,
										const std::string& dataDir = std::string("./"),
									    const bool         gatherMode = true );

		/// 出力対象のリーフブロック情報を登録 (Data用)
		/// 
		/// @param[in] dataClassID  データクラスID (配列)
		/// @param[in] kind         データの種類
		/// @param[in] dataType     データの型
		/// @param[in] vc           仮想セルサイズ
		/// @param[in] name         系の名称
		/// @param[in] prefix       リーフブロックファイルのPrefix
		/// @param[in] extension    リーフブロックファイルの拡張子
		/// @param[in] step         タイムステップ情報
		/// @param[in] dataDir      リーフブロックファイルの出力ディレクトリを指定 (コンストラクタで指定した出力ディレクトリからの相対パス)
		/// @param[in] stepSubDir   タイムステップごとの出力ディレクトリフラグ (trueの場合、タイムステップごとのディレクトリを作成)
		/// @return 成功した場合true, 失敗した場合false
		/// 
		bool RegisterDataInformation( const int          *dataClassID,
		                              const LB_KIND       kind,
									  const LB_DATA_TYPE  dataType,
									  const short         vc,
									  const std::string&  name,
									  const std::string&  prefix,
									  const std::string&  extension,
									  const IdxStep&      step,
									  const std::string&  dataDir = std::string("./"),
									  const bool          stepSubDir = false );

		/// 出力対象データの単位系設定
		/// 
		/// @param[in] unit 単位系
		/// @return 成功した場合true, 失敗した場合false
		/// @note 設定を省略した場合、デフォルト値で動作
		///
		bool SetUnit( const IdxUnit& unit );
		
		/// ファイル出力を実行
		/// 
		/// @return 成功した場合true, 失敗した場合false
		///
		/// @note ここで出力される情報は、インデックスファイルとOctree。
		///       RegisterCellIDInformation()でCellIDが登録されていない場合、cellid.bcmは出力されない
		///       RegisterDataInformation()でDataが登録されていない場合、data.bcmは出力されない
		///       RegisterDataInformation()で複数のDataを登録している場合、data.bcmにまとめて記載
		/// 
		bool Save();

		/// リーフブロックファイルを出力
		/// 
		/// @param[in] name 系の名称 (Registerした際に設定した名前)
		/// @param[in] step タイムステップ
		/// @return 成功した場合true, 失敗した場合false
		/// 
		/// @note ファイル出力する対象がCellIDの場合、stepは無視される．
		/// 
		bool SaveLeafBlock(const char* name, unsigned int step = 0);
	
	private:
		/// インデックスファイルを出力
		/// 
		/// @param[in] octName  Octreeファイル名
		/// @param[in] numLeaf  総リーフブロック数
		/// @return 成功した場合true, 失敗した場合false
		/// 
		bool SaveIndex(const std::string& octName, const int numLeaf);

		/// プロセス情報ファイルを出力
		///
		/// @param[in] filepath プロセス情報ファイルの名前(相対パス)
		/// @param[in] part     リーフブロックの1次元分割情報
		/// @return 成功した場合true, 失敗した場合false
		/// 
		bool SaveIndexProc(const std::string& filepath, const Partition& part);
		
		/// CellID情報ファイルを出力
		///
		/// @param[in] procName プロセス情報ファイルの名前
		/// @param[in] octName  Octreeファイルの名前
		/// @return 成功した場合true, 失敗した場合false
		/// 
		bool SaveIndexCellID(const std::string& procName, const std::string& octName);

		/// Data情報ファイルを出力
		/// 
		/// @param[in] procName プロセス情報ファイルの名前
		/// @param[in] octName  Octreeファイルの名前
		/// @return 成功した場合true, 失敗した場合false
		/// 
		bool SaveIndexData(const std::string& procName, const std::string& octName);

		/// Octreeファイルを出力
		/// 
		/// @param[in] filepath Octreeファイル名 (相対パス)
		/// @param[in] octree   Octree
		/// 
		/// @return 成功した場合true, 失敗した場合false
		/// 
		bool SaveOctree(const std::string& filepath, const BCMOctree* octree);
	
	private:
		BlockManager&          m_blockManager;   ///< ブロックマネージャ
		const MPI::Intracomm&  m_comm;           ///< MPIコミュニケータ
		const BCMOctree*       m_octree;         ///< BCMOctree
		const Vec3r            m_globalOrigin;   ///< 計算空間の起点座標
		const Vec3r            m_globalRegion;   ///< 計算空間全体の領域サイズ
		IdxUnit                m_unit;           ///< 単位系
		std::string            m_targetDir;      ///< ファイル出力ターゲットディレクトリ名
		std::vector<IdxBlock>  m_idxBlockList;   ///< 登録されたブロック情報リスト
	};

} // namespace BCMFileIO


#endif // __BCMTOOLS_BCM_FILE_SAVER_H__

