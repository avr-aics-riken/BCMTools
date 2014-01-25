/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file  BCMFileCommon.h
/// @brief BCMファイルIO用共通クラス群
///


#ifndef __BCMTOOLS_BCM_FILE_HEADER_H__
#define __BCMTOOLS_BCM_FILE_HEADER_H__

#include <limits.h>
#include <vector>
#include <string>
#include <list>

#include "BitVoxel.h"

#if defined(_WIN32)
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#else
#include <stdint.h>
#endif

/// Octreeファイルのエンディアン識別子 (OC01)
#define OCTREE_FILE_IDENTIFIER    (('O' | ('C' << 8) | ('0' << 16) | ('1' << 24)))

/// LeafBlockファイルのエンディアン識別子 (LB01)
#define LEAFBLOCK_FILE_IDENTIFIER (('L' | ('B' << 8) | ('0' << 16) | ('1' << 24)))


namespace BCMFileIO {

#ifdef __GNUC__
#pragma pack(push, 1)
#define ALIGNMENT __attribute__((packed))
#else
#pragma pack(1)
#define ALIGNMENT
#endif // __GNUC__
	
	/// Octreeファイルヘッダ構造体
	struct OctHeader
	{
		unsigned int identifier;   ///< エンディアン識別子
		double       org[3];       ///< 原点座標
		double       rgn[3];       ///< 領域サイズ
		unsigned int rootDims[3];  ///< ルート分割数
		unsigned int maxLevel;     ///< Octree最大分割レベル
		uint64_t     numLeaf;      ///< リーフノード数
		uint64_t     padding;      ///< 16バイトアライメント用パディング

		OctHeader() : padding(0) {}
	
	} ALIGNMENT;

	/// LeafBlockファイルヘッダ構造体
	struct LBHeader
	{
		unsigned int   identifier; ///< エンディアン識別子
		unsigned char  kind;       ///< ブロックファイル種類
		unsigned char  dataType;   ///< 1セルあたりのサイズ
		unsigned short bitWidth;   ///< 1セルあたりのビット幅
		unsigned int   vc;         ///< 仮想セルサイズ
		unsigned int   size[3];    ///< ブロックサイズ
		uint64_t       numBlock;   ///< ファイルに記載されている総ブロック数

	} ALIGNMENT;
	
	/// LeafBlockのCellIDヘッダ構造体
	struct LBCellIDHeader
	{
		uint64_t numBlock; ///< ブロック数
		uint64_t compSize; ///< 圧縮符号サイズ (バイト単位)

	} ALIGNMENT;

	/// RLE圧縮符号の走査用構造体
	struct GridRleCode
	{
		bitVoxelCell    c; ///< データ
		unsigned char len; ///< ラン長

	} ALIGNMENT;


#ifdef __GNUC__
#pragma pack(pop)
#else  // __GNUC__
#pragma pack()
#endif // __GNUC__

	/// リーフブロックデータタイプ
	enum LB_KIND
	{
		LB_CELLID = 0, ///< グリッド
		LB_SCALAR = 1, ///< スカラ
		LB_VECTOR = 3, ///< ベクター
		LB_TENSOR = 9, ///< テンソル
	};
	
	/// リーフセルのデータ識別子
	enum LB_DATA_TYPE
	{
		LB_INT8    =  0,
		LB_UINT8   =  1,
		LB_INT16   =  2,
		LB_UINT16  =  3,
		LB_INT32   =  4,
		LB_UINT32  =  5,
		LB_INT64   =  6,
		LB_UINT64  =  7,
		LB_FLOAT32 =  8,
		LB_FLOAT64 =  9
	};

	/// インデックスファイル用単位系情報
	struct IdxUnit
	{
		std::string length;    ///< 長さ単位 (NonDimensional, m, cm, mm)
		double      L0_scale;  ///< 規格化に用いたスケール (単位:指定単位)
		std::string velocity;  ///< 時間単位 (NonDimensional, Dimensional)
		double      V0_scale;  ///< 規格化に用いた時間スケール (単位:Dimensionalの場合m/s)
	};

	/// インデックスファイル用プロセス情報
	struct IdxProc
	{
		std::string  hostname; ///< ホスト名
		unsigned int rank;     ///< ランク番号
		unsigned int rangeMin; ///< ブロックIDのレンジ最小値
		unsigned int rangeMax; ///< ブロックIDのレンジ最大値
	};

	/// インデックスファイル用タイムステップ情報
	class IdxStep
	{
	public:
		IdxStep() : m_rangeMin(0), m_rangeMax(0), m_rangeInterval(0) { }
		~IdxStep(){ }

		bool SetRange(const unsigned int rangeMin, const unsigned int rangeMax, const unsigned int rangeInterval=1)
		{
			if(rangeMax <= rangeMin){ return false; }
			m_rangeMin      = rangeMin;
			m_rangeMax      = rangeMax;
			m_rangeInterval = rangeInterval;
		}

		void AddStep(const unsigned int step ){
			m_adds.push_back(step);
		}

		void SubStep(const unsigned int step ){
			m_subs.push_back(step);
		}

		// TODO !TEST
		/// @note リストはメソッド内で確保するため、リスト取得後、不要になったら解放してください。
		std::list<unsigned int>* GetStepLists() const
		{
			std::list<unsigned int>*steps = new std::list<unsigned int>;

			for(unsigned int i = m_rangeMin; i < m_rangeMax; i+= m_rangeInterval){
				steps->push_back(i);
			}
			
			// 追加リストからのステップ追加
			for(std::vector<unsigned int>::const_iterator it = m_adds.begin(); it != m_adds.end(); ++it){
				steps->push_back((*it));
			}

			// 削除リストからステップ削除
			for(std::vector<unsigned int>::const_iterator it = m_subs.begin(); it != m_subs.end(); ++it){
				steps->remove((*it));
			}
			
			steps->sort();

			return steps;
		}

		unsigned int GetRangeMin()      const { return m_rangeMin; }
		unsigned int GetRangeMax()      const { return m_rangeMax; }
		unsigned int GetRangeInterval() const { return m_rangeInterval; }

		const std::vector<unsigned int>& GetAddStepList() const { return m_adds; }
		const std::vector<unsigned int>& GetSubStepList() const { return m_subs; }

	private:
		unsigned int m_rangeMin;          ///< タイムステップレンジ (Min)
		unsigned int m_rangeMax;          ///< タイムステップレンジ (Max)
		unsigned int m_rangeInterval;     ///< タイムステップレンジ (Interval)
		std::vector<unsigned int> m_adds; ///< 追加タイムステップリスト
		std::vector<unsigned int> m_subs; ///< 削除タイムステップリスト
	};
	
	/// インデックスファイル用ブロック情報
	struct IdxBlock
	{
		std::string  rootDir;      ///< インデックスファイルのディレクトリ
		std::string  dataDir;      ///< データディレクトリ

		int          dataClassID;  ///< データクラスID(ファイルには記載しない)
		LB_DATA_TYPE dataType;     ///< セルのデータ識別子
		std::string  name;         ///< 系の名称
		LB_KIND      kind;         ///< リーフブロックタイプ
		unsigned int bitWidth;     ///< セルあたりのビット幅
		unsigned int vc;           ///< 仮想セルサイズ
		std::string  prefix;       ///< ファイル名Prefix
		std::string  extension;    ///< ファイル拡張子
		bool         isGather;     ///< Gatherフラグ
		bool         isSubStepDir; ///< ステップごとのサブディレクトリフラグ

		IdxStep*     stepInfo;

		IdxBlock() :
			rootDir(std::string("")),
			dataDir(std::string("")),
			dataClassID(-1),
			isGather(false),
			isSubStepDir(false),
			stepInfo(NULL)
		{}
	};

	/// データクラスIDからブロック情報を取得するユーティリティ関数
	/// 
	/// @param[in] idxBlockList ブロック情報リスト
	/// @param[in] dataClassID  データクラスID
	/// @return ブロック情報のポインタ
	///
	inline IdxBlock* findIdxBlock(std::vector<IdxBlock>& idxBlockList, const int dataClassID ){
		for(std::vector<IdxBlock>::iterator it = idxBlockList.begin(); it != idxBlockList.end(); ++it){
			if( it->dataClassID == dataClassID){
				return &(*it);
			}
		}
		return NULL;
	}

	/// 系の名称からブロック情報を取得するユーティリティ関数
	/// 
	/// @param[in] idxBlockList ブロック情報リスト
	/// @param[in] name         系の名称
	/// @return ブロック情報のポインタ
	///
	inline IdxBlock* findIdxBlock(std::vector<IdxBlock>& idxBlockList, const std::string& name ){
		for(std::vector<IdxBlock>::iterator it = idxBlockList.begin(); it != idxBlockList.end(); ++it){
			if( it->name == name){
				return &(*it);
			}
		}
		return NULL;
	}

	
	/// 2byte用エンディアンスワップ
	static inline void BSwap16(void* a){
		unsigned short* x = (unsigned short*)a;
		*x = (unsigned short)( ((((*x) & 0xff00) >> 8 ) | (((*x) & 0x00ff) << 8)) );
	}

	/// 4byte用エンディアンスワップ
	static inline void BSwap32(void* a){
		unsigned int* x = (unsigned int*)a;
		*x = ( (((*x) & 0xff000000) >> 24) | (((*x) & 0x00ff0000) >> 8 ) |
		       (((*x) & 0x0000ff00) <<  8) | (((*x) & 0x000000ff) << 24) );
	}

	/// 8byte用エンディアンスワップ
	static inline void BSwap64(void* a){
		uint64_t* x = (uint64_t*)a;
		*x=( (((*x) & 0xff00000000000000ull) >> 56) | (((*x) & 0x00ff000000000000ull) >> 40) |
		     (((*x) & 0x0000ff0000000000ull) >> 24) | (((*x) & 0x000000ff00000000ull) >>  8) |
		     (((*x) & 0x00000000ff000000ull) <<  8) | (((*x) & 0x0000000000ff0000ull) << 24) |
		     (((*x) & 0x000000000000ff00ull) << 40) | (((*x) & 0x00000000000000ffull) << 56) );
	}


} // namespace BCMFileIO

#endif // __BCMTOOLS_BCM_FILE_HEADER_H__

