///
/// @file  PartitionMapper.h
/// @brief MxNデータロードのためのマッピングクラス
/// 
/// 用語リスト
/// - did   : グローバルなデータID
/// - FID   : ファイルID
/// - FDID  : ファイル内の相対データID

#ifndef __BCMTOOLS_PARTITION_MAPPER_H__
#define __BCMTOOLS_PARTITION_MAPPER_H__

#include "Partition.h"

namespace BCMFileIO {

	class PartitionMapper
	{
	public:
		
		/// ファイルIDとファイル内のデータIDリスト構造体
		struct FDIDList
		{
			int              FID;   ///< FID
			std::vector<int> FDIDs; ///< FDIDリスト
		};
	
		/// コンストラクタ
		/// 
		/// @param[in] writeProcs ファイル出力時の並列数
		/// @param[in] readProcs  ファイル読込時の並列数
		/// @param[in] numLeaf    総リーフブロック数
		/// 
		PartitionMapper(const int writeProcs, const int readProcs, const int numLeaf)
			: partW(writeProcs, numLeaf),
			  partR(readProcs, numLeaf),
			  writeProcs(writeProcs),
			  readProcs(readProcs)
		{
	
		}

		/// デストラクタ
		~PartitionMapper()
		{
	
		}

		/// ファイル読込時の並列数を取得
		int GetReadProcs() const { return readProcs; }
		
		/// ファイル出力時の並列数を取得
		int GetWriteProcs() const { return writeProcs; }
	
		/// ファイル読込時の先頭didを取得
		/// 
		/// @param[in] rank プロセス番号
		/// @return 先頭did
		///
		int GetStart(const int rank) const { return partR.getStart(rank); }

		/// ファイル読込時の末尾didを取得
		///
		/// @param[in] rank プロセス番号
		/// @return 末尾did
		///
		int GetEnd(const int rank) const   { return partR.getEnd(rank);   }
	
		/// グローバルデータID(did)が保存されているファイルのID(FID)を取得
		/// 
		/// @param[in] did did
		/// @return didが保存されているFID
		///
		int GetFID(int did ){ return partW.getRank(did); }
	
		/// グローバルデータID(did)のファイル相対データID(FDID)を取得
		/// 
		/// @param[in] did did
		/// @return didのFDID
		/// 
		int GetFDID(int did )   { return did - partW.getStart( this->GetFID(did) ); }
		
		/// FDIDリストのリストを取得
		/// 
		/// @param[in]  myRank    プロセス番号
		/// @param[out] fdidlists FDIDリストのリスト
		/// @return 成功した場合true, 失敗した場合false
		/// 
		bool GetFDIDLists(const int myRank, std::vector<FDIDList>& fdidlists )
		{
			// 自ノードが担当するブロックのレンジを取得
			int didRange[2] = { GetStart(myRank), GetEnd(myRank) };
	
			// 自ノードが担当するDIDが保存されているFIDをリストアップ
			std::vector<int> fids;
			fids.reserve(didRange[1] - didRange[0]);
			for(int did = didRange[0]; did < didRange[1]; did++){
				fids.push_back(GetFID(did));
			}
	
			// FIDの重複削除
			sort(fids.begin(), fids.end());
			fids.erase( unique(fids.begin(), fids.end()), fids.end());
	
	
			/// FDIDListのリストを作成
			fdidlists.resize(fids.size());
	
			for(size_t i = 0; i < fids.size(); i++){
				fdidlists[i].FID = fids[i];
			}
	
			for(int did = didRange[0]; did < didRange[1]; did++){
				fdidlists[ (GetFID(did) - fids[0]) ].FDIDs.push_back(GetFDID(did));
			}
	
			return true;
		}
	
	
	private:
		Partition partW;      ///< ファイル出力時のPartition
		Partition partR;      ///< ファイル読込時のPartition
		const int writeProcs; ///< ファイル出力時の並列数
		const int readProcs;  ///< ファイル読込時の並列数
	};

} // namespace BCMFileIO



#endif // __BCMTOOLS_PARTITION_MAPPER_H__

