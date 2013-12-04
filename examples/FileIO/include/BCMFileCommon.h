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
		LB_CELLID  = 0, ///< グリッド
		LB_SCALAR  = 1, ///< スカラ
		LB_VECTOR3 = 3, ///< ベクトル (3要素)
		LB_VECTOR4 = 4, ///< ベクトル (4要素)
		LB_VECTOR6 = 6, ///< ベクトル (6要素)
		LB_TENSOR  = 9, ///< テンソル (9要素)
	};
	
	/// リーフセルのデータ識別子
	enum LB_DATA_TYPE
	{
		LB_INT8    =  0, ///< 符号付き 8bit整数型
		LB_UINT8   =  1, ///< 符号なし 8bit整数型
		LB_INT16   =  2, ///< 符号付き16bit整数型
		LB_UINT16  =  3, ///< 符号なし16bit整数型
		LB_INT32   =  4, ///< 符号付き32bit整数型
		LB_UINT32  =  5, ///< 符号なし32bit整数型
		LB_INT64   =  6, ///< 符号付き64bit整数型
		LB_UINT64  =  7, ///< 符号なし64bit整数型
		LB_FLOAT32 =  8, ///< 32bit浮動小数点 (単精度浮動小数点)
		LB_FLOAT64 =  9  ///< 64bit浮動小数点 (倍精度浮動小数点)
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
		/// コンストラクタ
		IdxStep() : m_rangeMin(0), m_rangeMax(0), m_rangeInterval(0), m_time(.0f), m_deltaT(.1f) { }
		
		/// コンストラクタ
		///
		/// @param[in] rangeMin      タイムステップレンジの開始インデックス
		/// @param[in] rangeMax      タイムステップレンジの終了インデックス
		/// @param[in] rangeInterval ステップ間隔
		///
		IdxStep(const unsigned int rangeMin, const unsigned int rangeMax, const unsigned int rangeInterval=1)
		  : m_rangeMin(0), m_rangeMax(0), m_rangeInterval(0), m_time(.0f), m_deltaT(.1f)
		{
			SetRange(rangeMin, rangeMax, rangeInterval);
		}
		
		/// デストラクタ
		~IdxStep(){ }

		/// タイムステップレンジ設定
		///
		/// @param[in] rangeMin      タイムステップレンジの開始インデックス
		/// @param[in] rangeMax      タイムステップレンジの終了インデックス
		/// @param[in] rangeInterval ステップ間隔
		/// @return 成功した場合true, 失敗した場合false
		///
		bool SetRange(const unsigned int rangeMin, const unsigned int rangeMax, const unsigned int rangeInterval=1)
		{
			if(rangeMax <= rangeMin){ return false; }
			m_rangeMin      = rangeMin;
			m_rangeMax      = rangeMax;
			m_rangeInterval = rangeInterval;

			return true;
		}

		/// 追加ステップの設定
		/// 
		/// @param[in] step 追加ステップ
		///
		void AddStep(const unsigned int step ){
			m_adds.push_back(step);
		}

		/// 削除ステップの設定
		///
		/// @param[in] step 削除ステップ
		///
		void SubStep(const unsigned int step ){
			m_subs.push_back(step);
		}

		/// ステップが設定したリストに含まれるかを判定
		///
		/// @param[in] step 判定するステップ番号
		/// @return stepがリストに含まれる場合true, 含まれない場合false
		///
		bool IsCorrect(const unsigned int step) const
		{
			for(std::vector<unsigned int>::const_iterator it = m_adds.begin(); it != m_adds.end(); ++it){
				if( *it == step ){ return true; }
			}
			for(std::vector<unsigned int>::const_iterator it = m_subs.begin(); it != m_subs.end(); ++it){
				if( *it == step ){ return false; }
			}

			if( m_rangeMin <= step && m_rangeMax >= step ){
				if( ((step - m_rangeMin) % m_rangeInterval) == 0 ){
					return true;
				}else{
					return false;
				}
			}
			return false;
		}

		/// Step = 0における時刻を設定
		///
		/// @param[in] time Step = 0における時刻
		///
		void SetInitalTime(float time){ m_time   = time;   }
		
		/// Step間の時刻幅を設定
		///
		/// @param[in] deltaT Step間の時刻幅
		///
		void SetDeltaT(float deltaT)  { m_deltaT = deltaT; }

		/// 設定したStepのリストを取得
		///
		/// @return ステップのリスト
		/// @note リストはメソッド内で確保するため、リスト取得後、不要になったら解放してください。
		///
		std::list<unsigned int>* GetStepList() const
		{
			if( m_rangeInterval == 0 ){ return NULL; }

			std::list<unsigned int>*steps = new std::list<unsigned int>;

			for(unsigned int i = m_rangeMin; i <= m_rangeMax; i+= m_rangeInterval){
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
		
		/// ステップの開始インデックスを取得
		///
		/// @return 開始インデックス
		///
		unsigned int GetRangeMin()      const { return m_rangeMin; }
		
		/// ステップの終了インデックスを取得
		///
		/// @return 終了インデックス
		///
		unsigned int GetRangeMax()      const { return m_rangeMax; }

		/// ステップ間隔を取得
		///
		/// @return ステップ間隔
		///
		unsigned int GetRangeInterval() const { return m_rangeInterval; }

		/// 追加ステップリストを取得
		///
		/// @return 追加ステップリスト
		///
		const std::vector<unsigned int>& GetAddStepList() const { return m_adds; }
		
		/// 削除ステップリストを取得
		/// 
		/// @return 削除ステップリスト
		///
		const std::vector<unsigned int>& GetSubStepList() const { return m_subs; }
		
		/// Step = 0における時刻を取得
		/// 
		/// @return Step = 0における時刻
		///
		float GetInitialTime() const { return m_time;   }
		
		/// Step間の時間幅を取得
		///
		/// @return Step間の時間幅
		///
		float GetDeltaT()      const { return m_deltaT; }

	private:
		unsigned int m_rangeMin;          ///< タイムステップレンジ (Min)
		unsigned int m_rangeMax;          ///< タイムステップレンジ (Max)
		unsigned int m_rangeInterval;     ///< タイムステップレンジ (Interval)
		std::vector<unsigned int> m_adds; ///< 追加タイムステップリスト
		std::vector<unsigned int> m_subs; ///< 削除タイムステップリスト
		float        m_time;              ///< Step = 0 における時刻
		float        m_deltaT;            ///< Step間の時間幅
	};
	
	/// インデックスファイル用ブロック情報
	struct IdxBlock
	{
		std::string      rootDir;      ///< インデックスファイルのディレクトリ
		std::string      dataDir;      ///< データディレクトリ

		std::vector<int> dataClassID;  ///< データクラスID(マルチコンポーネント対応のため配列)
		LB_DATA_TYPE     dataType;     ///< セルのデータ識別子
		std::string      name;         ///< 系の名称
		LB_KIND          kind;         ///< リーフブロックタイプ
		unsigned int     bitWidth;     ///< セルあたりのビット幅
		unsigned int     vc;           ///< 仮想セルサイズ
		std::string      prefix;       ///< ファイル名Prefix
		std::string      extension;    ///< ファイル拡張子
		bool             isGather;     ///< Gatherフラグ
		bool             isStepSubDir; ///< ステップごとのサブディレクトリフラグ
		IdxStep          step;         ///< タイムステップ情報

		bool         separateVCUpdate;

		IdxBlock() :
			rootDir(std::string("")),
			dataDir(std::string("")),
			vc(0),
			isGather(false),
			isStepSubDir(false),
			separateVCUpdate(false)
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
			for(size_t i = 0; i < it->dataClassID.size(); i++){
				if( it->dataClassID[i] == dataClassID){
					return &(*it);
				}
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

	/// データクラスIDからブロック情報を取得するユーティリティ関数 (const用)
	/// 
	/// @param[in] idxBlockList ブロック情報リスト
	/// @param[in] dataClassID  データクラスID
	/// @return ブロック情報のポインタ
	///
	inline const IdxBlock* findIdxBlock(const std::vector<IdxBlock>& idxBlockList, const int dataClassID ){
		for(std::vector<IdxBlock>::const_iterator it = idxBlockList.begin(); it != idxBlockList.end(); ++it){
			for(size_t i = 0; i < it->dataClassID.size(); i++){
				if( it->dataClassID[i] == dataClassID){
					return &(*it);
				}
			}
		}
		return NULL;
	}

	/// 系の名称からブロック情報を取得するユーティリティ関数 (const用)
	/// 
	/// @param[in] idxBlockList ブロック情報リスト
	/// @param[in] name         系の名称
	/// @return ブロック情報のポインタ
	///
	inline const IdxBlock* findIdxBlock(const std::vector<IdxBlock>& idxBlockList, const std::string& name ){
		for(std::vector<IdxBlock>::const_iterator it = idxBlockList.begin(); it != idxBlockList.end(); ++it){
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

