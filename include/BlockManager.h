///
/// @file BlockManager.h
/// @brief ブロックマネージャクラス
///

#ifndef BLOCK_MANAGER_H
#define BLOCK_MANAGER_H

#include <vector>
#include <map>
#include <algorithm>
#include "mpi.h"
#include "BlockBase.h"
#include "CommBuffer.h"
#include "BCMTools.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

/// ブロックマネージャ (シングルトン).
/// 並列計算ノード内の全ブロックを管理する．
///
///  @todo 計算セルを一つも含まないブロックをスキップできる仕組みを追加．
///        例えば，BlockBaseクラスににフラグを設けるとか．
///
class BlockManager {

  /// ブロック面を識別するための構造体.
  struct FaceID {
    int id;           ///< ブロックID(ローカル)
    Face face;        ///< フェイス番号
    Subface subface;  ///< サブフェイス番号
    FaceID(int id, Face face, Subface subface = Subface(0))
     : id(id), face(face), subface(subface) {}
  };


  /// 通信バッファクラステーブル.
  struct CommBufferTable {
    int tag;                  ///< 通信タグ番号
    SendBuffer* sendBuffer;   ///< 送信用バッファクラス配列
    RecvBuffer* recvBuffer;   ///< 受信用バッファクラス配列

    CommBufferTable(int tag, bool separate = false) : tag(tag) {
      sendBuffer = new SendBuffer(tag);
      recvBuffer = new RecvBuffer(tag);
    }

    ~CommBufferTable() {
      delete sendBuffer;
      delete recvBuffer;
    }
  };


  typedef std::vector<BlockBase*> BlockList;
  typedef std::vector<FaceID> FaceList;
  typedef std::map<int, FaceList> FaceListMap;
  typedef std::map<int, CommBufferTable*> CommBufferTableMap;


  int startID;     ///< 先頭ブロックID
  int numBlock;    ///< 担当ブロック数

  MPI::Intracomm comm;   ///< MPIコミュニケータ

  BlockList blockList;   ///< ブロックのリスト

  Vec3i size;  ///< 各ブロック内のセル数

  bool faceListPrepared;
  bool separateFaceListPrepared;

  FaceList localFaceList;     ///< 同一ノード内のブロックと接する面のリスト
  FaceListMap sendFaceList;   ///< 他ノードへ送信する面のリスト
  FaceListMap recvFaceList;   ///< 他ノードから受信する面のリスト

  FaceList localSeparateFaceList[3];    ///< 同一ノード内のブロックと接する面のリスト(方向別)
  FaceListMap sendSeparateFaceList[3];  ///< 他ノードへ送信する面のリスト(方向別)
  FaceListMap recvSeparateFaceList[3];  ///< 他ノードから受信する面のリスト(方向別)


  CommBufferTableMap commBufferTableMap;  ///< 「データクラスID→通信バッファ」のマップ

  CommBufferTableMap commBufferSeparateTableMap[3]; ///< 「データクラスID→通信バッファ」のマップ(方向別)

  /// recvFaceListソート用比較ファンクタ.
  /// これを用いて送信側sendFaceListのFaceIDオーダに合わせてソートする
  class RecvFaceComp {
    const BlockList& bList;
  public:
    RecvFaceComp(const BlockList& bList) : bList(bList) {};
    bool operator() (const FaceID& a, const FaceID& b) {
      return getOrder(a) < getOrder(b);
    }
  private:
    int getOrder(const FaceID& c) const {
      const NeighborInfo* neighborInfo = bList[c.id]->getNeighborInfo();
      return neighborInfo[c.face].getID(c.subface) * (NUM_FACE*NUM_SUBFACE)
             + NeighborInfo::reverseFace(c.face) * NUM_SUBFACE
             + neighborInfo[c.face].getNeighborSubface();
    }
  };

public:

  /// ブロックマネージャ インスタンスの取得.
  static BlockManager& getInstance() {
    static BlockManager instance;
    return instance;
  }

  /// MPIコミュニケータ設定.
  void setCommunicator(const MPI::Intracomm& comm) { this->comm = comm; }

  /// MPIコミュニケータ取得.
  const MPI::Intracomm& getCommunicator() { return comm; }

  /// ブロックを登録.
  void registerBlock(BlockBase* block);

  /// ブロック登録の完了.
  /// 全担当ブロックをregisterBlockメソッドにより登録後，このメソッドを呼ぶ必要がある
  void endRegisterBlock();

  /// ブロックサイズ(全ブロックで等しい)取得.
  const Vec3i& getSize() const { assert(size != 0); return size; }

  /// 先頭ブロックIDの取得.
  int getStartID() const { return startID; }

  /// 担当ブロック数の取得.
  int getNumBlock() const { return numBlock;; }

//const BlockList& getBlockList() const { return blockList; }

  /// ブロックの取得.
  ///
  ///  @param[in] localID ローカルブロックID
  ///
  BlockBase* getBlock(int localID) {
    if (localID < 0 || localID >= numBlock) {
      std::cout << "***error: illegal local block ID" << std::endl;
      Exit(EX_FAILURE);
    }
    return blockList[localID];
  }

  /// ブロック配置情報を出力.
  void printBlockLayoutInfo();

  /// データクラスDの生成・登録.
  ///
  ///  @param[in] vc 仮想セル幅
  ///
  template <typename D>
  int setDataClass(int vc) {
    int id = -1;
    BlockList::const_iterator it = blockList.begin();
    for (; it != blockList.end(); ++it) {
      BlockBase* block = *it;
      Vec3i size = block->getSize();
      int id0 = block->setDataClass(new D(size, vc)); 
      if (it == blockList.begin()) {
        id = id0;
      } else {
        if (id != id0) {
          std::cout << "error: DataClass register inconsistency" << std::endl;
          Exit(EX_FAILURE);
        }
      }
    }
    return id;
  }

  /// 仮想セルアップデータUを持つデータクラスDの生成・登録.
  ///
  ///  @param[in] vc 仮想セル幅
  ///
  template <typename D, typename U>
  int setDataClass(int vc) {
    int id = -1;
    BlockList::const_iterator it = blockList.begin();
    for (; it != blockList.end(); ++it) {
      BlockBase* block = *it;
      Vec3i size = block->getSize();
      D* dataClass = new D(size, vc);
      U* updater = new U(block->getNeighborInfo());
      dataClass->setUpdater(updater);
      int id0 = block->setDataClass(dataClass); 
      if (it == blockList.begin()) {
        id = id0;
      } else {
        if (id != id0) {
          std::cout << "error: DataClass register inconsistency" << std::endl;
          Exit(EX_FAILURE);
        }
      }
    }
    return id;
  }

  /// 仮想セル同期の準備.
  ///
  ///  @param[in] dataClassID 対象データクラスのID
  ///  @param[in] tag 通信タグ番号
  ///  @param[in] separate true=3方向同時/false=3方向別々通信
  ///
  void prepareForVCUpdate(int dataClassID, int tag, bool separate = false);

  /// 仮想セル同期(3方向同時).
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void updateVC(int dataClassID) {
    assert(commBufferTableMap.find(dataClassID) != commBufferTableMap.end());
    beginUpdateVC(dataClassID);
    endUpdateVC(dataClassID);
  }

  /// 仮想セル同期(x方向).
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void updateVC_X(int dataClassID) { updateVC(dataClassID, 0); }

  /// 仮想セル同期(y方向).
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void updateVC_Y(int dataClassID) { updateVC(dataClassID, 1); }

  /// 仮想セル同期(z方向).
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void updateVC_Z(int dataClassID) { updateVC(dataClassID, 2); }

  /// 仮想セル同期(個別方向).
  ///
  ///  @param[in] dataClassID データクラスID
  ///  @param[in] xyz 方向(0:x, 1:y, 2:z)
  ///
  void updateVC(int dataClassID, int xyz) {
    assert(commBufferSeparateTableMap[xyz].find(dataClassID)
           != commBufferSeparateTableMap[xyz].end());
    beginUpdateVC(dataClassID, xyz);
    endUpdateVC(dataClassID, xyz);
  }

  /// 仮想セル同期(3方向同時)通信開始.
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void beginUpdateVC(int dataClassID) {
    assert(commBufferTableMap.find(dataClassID) != commBufferTableMap.end());
    recvVCBegin(dataClassID);
    copyVCToSendBuffer(dataClassID);
    sendVCBegin(dataClassID);
    copyVCFromNeighbor(dataClassID);
  }

  /// 仮想セル同期(x方向)通信開始.
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void beginUpdateVC_X(int dataClassID) { beginUpdateVC(dataClassID, 0); }

  /// 仮想セル同期(y方向)通信開始.
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void beginUpdateVC_Y(int dataClassID) { beginUpdateVC(dataClassID, 1); }

  /// 仮想セル同期(z方向)通信開始.
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void beginUpdateVC_Z(int dataClassID) { beginUpdateVC(dataClassID, 2); }

  /// 仮想セル同期(個別方向)通信開始.
  ///
  ///  @param[in] dataClassID データクラスID
  ///  @param[in] xyz 方向(0:x, 1:y, 2:z)
  ///
  void beginUpdateVC(int dataClassID, int xyz) {
    assert(commBufferSeparateTableMap[xyz].find(dataClassID) 
           != commBufferSeparateTableMap[xyz].end());
    recvVCBegin(dataClassID, xyz);
    copyVCToSendBuffer(dataClassID, xyz);
    sendVCBegin(dataClassID, xyz);
    copyVCFromNeighbor(dataClassID, xyz);
  }

  /// 仮想セル同期(3方向同時)通信完了待ち.
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void endUpdateVC(int dataClassID) {
    assert(commBufferTableMap.find(dataClassID) != commBufferTableMap.end());
    recvVCEnd(dataClassID);
    copyVCFromRecvBuffer(dataClassID);
    sendVCEnd(dataClassID);
  }

  /// 仮想セル同期(x方向)通信完了待ち.
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void endUpdateVC_X(int dataClassID) { endUpdateVC(dataClassID, 0); }

  /// 仮想セル同期(y方向)通信完了待ち.
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void endUpdateVC_Y(int dataClassID) { endUpdateVC(dataClassID, 1); }

  /// 仮想セル同期(z方向)通信完了待ち.
  ///
  ///  @param[in] dataClassID データクラスID
  ///
  void endUpdateVC_Z(int dataClassID) { endUpdateVC(dataClassID, 2); }

  /// 仮想セル同期(個別方向)通信完了待ち.
  ///
  ///  @param[in] dataClassID データクラスID
  ///  @param[in] xyz 方向(0:x, 1:y, 2:z)
  ///
  void endUpdateVC(int dataClassID, int xyz) {
    assert(commBufferSeparateTableMap[xyz].find(dataClassID)
           != commBufferSeparateTableMap[xyz].end());
    recvVCEnd(dataClassID, xyz);
    copyVCFromRecvBuffer(dataClassID, xyz);
    sendVCEnd(dataClassID, xyz);
  }


#if 0
  void printFaceList(const FaceList& faceList) {
    for (int i = 0; i < faceList.size(); ++i) {
      std::cout << "(" << faceList[i].id +startID << "," << faceList[i].face << ","
                       << faceList[i].subface << ")";
    }
    std::cout << std::endl;
  }

  void printFaceListMap(const FaceListMap& faceListMap) {
    FaceListMap::const_iterator it = faceListMap.begin();
    for (; it != faceListMap.end(); ++it) {
      std::cout << "  " << it->first << ": ";
      printFaceList(it->second);
    }
  }

  void printFaceListPeer(const FaceList& faceList) {
    for (int i = 0; i < faceList.size(); ++i) {
      int id = faceList[i].id;
      Face face = faceList[i].face;
      Subface subface = faceList[i].subface;
      const NeighborInfo* neighborInfo = blockList[id]->getNeighborInfo();
      std::cout << "(" << neighborInfo[face].getID(subface)
                << "," << NeighborInfo::reverseFace(face)
                << "," << neighborInfo[face].getNeighborSubface() << ")";
    }
    std::cout << std::endl;
  }

  void printFaceListMapPeer(const FaceListMap& faceListMap) {
    FaceListMap::const_iterator it = faceListMap.begin();
    for (; it != faceListMap.end(); ++it) {
      std::cout << "  " << it->first << ": ";
      printFaceListPeer(it->second);
    }
  }
#endif

private:

  /// コピーコンストラクタ(コピー禁止).
  BlockManager(const BlockManager& rhs);

  /// 代入演算子(コピー禁止).
  BlockManager& operator=(const BlockManager& rhs);

  /// フェイスリストの設定.
  void setFaceLists();

  /// フェイスリスト(3方向別)の設定.
  void setSeparateFaceLists();

  /// 各ブロックで，指定したデータクラスの隣接データクラス情報を設定.
  void setDataClassNeighbor(int dataClassID);

  /// 各ブロックで，指定したデータクラスの隣接データクラス情報(3方向別)を設定.
  void setDataClassNeighbor(int dataClassID, int xyz);

  /// 各ブロックで，指定したデータクラスの送信バッファポインタを設定.
  void setSendBufferPointers(SendBuffer* sendBuffer, int dataClassID);

  /// 各ブロックで，指定したデータクラスの送信バッファポインタ(3方向別)を設定.
  void setSendBufferPointers(SendBuffer* sendBuffer, int dataClassID, int xyz);

  /// 各ブロックで，指定したデータクラスの受信バッファポインタを設定.
  void setRecvBufferPointers(RecvBuffer* recvBuffer, int dataClassID);

  /// 各ブロックで，指定したデータクラスの受信バッファポインタ(3方向別)を設定.
  void setRecvBufferPointers(RecvBuffer* recvBuffer, int dataClassID, int xyz);

  /// 各ブロックで，指定したデータクラスの仮想セルデータを隣接ブロックからコピー.
  void copyVCFromNeighbor(int dataClassID);

  /// 各ブロックで，指定したデータクラスの仮想セルデータを隣接ブロックからコピー(3方向別).
  void copyVCFromNeighbor(int dataClassID, int xyz);

  /// 各ブロックで，指定したデータクラスの仮想セルデータを送信バッファにコピー.
  void copyVCToSendBuffer(int dataClassID);

  /// 各ブロックで，指定したデータクラスの仮想セルデータを送信バッファにコピー(3方向別).
  void copyVCToSendBuffer(int dataClassID, int xyz);

  /// 各ブロックで，指定したデータクラスの仮想セルデータを受信バッファからコピー.
  void copyVCFromRecvBuffer(int dataClassID);

  /// 各ブロックで，指定したデータクラスの仮想セルデータを受信バッファからコピー(3方向別).
  void copyVCFromRecvBuffer(int dataClassID, int xyz);

  /// 仮想セル同期データ送信.
  void sendVC(int dataClassID) {
    commBufferTableMap[dataClassID]->sendBuffer->send();
  }

  /// 仮想セル同期データ送信(3方向別).
  void sendVC(int dataClassID, int xyz) {
    commBufferSeparateTableMap[xyz][dataClassID]->sendBuffer->send();
  }

  /// 仮想セル同期データ送信開始.
  void sendVCBegin(int dataClassID) {
    commBufferTableMap[dataClassID]->sendBuffer->sendBegin();
  }

  /// 仮想セル同期データ送信開始(3方向別).
  void sendVCBegin(int dataClassID, int xyz) {
    commBufferSeparateTableMap[xyz][dataClassID]->sendBuffer->sendBegin();
  }

  /// 仮想セル同期データ送信完了待ち.
  void sendVCEnd(int dataClassID) {
    commBufferTableMap[dataClassID]->sendBuffer->sendEnd();
  }

  /// 仮想セル同期データ送信完了待ち(3方向別).
  void sendVCEnd(int dataClassID, int xyz) {
    commBufferSeparateTableMap[xyz][dataClassID]->sendBuffer->sendEnd();
  }

  /// 仮想セル同期データ受信.
  void recvVC(int dataClassID) {
    commBufferTableMap[dataClassID]->recvBuffer->recv();
  }

  /// 仮想セル同期データ受信(3方向別).
  void recvVC(int dataClassID, int xyz) {
    commBufferSeparateTableMap[xyz][dataClassID]->recvBuffer->recv();
  }

  /// 仮想セル同期データ受信開始.
  void recvVCBegin(int dataClassID) {
    commBufferTableMap[dataClassID]->recvBuffer->recvBegin();
  }

  /// 仮想セル同期データ受信開始(3方向別).
  void recvVCBegin(int dataClassID, int xyz) {
    commBufferSeparateTableMap[xyz][dataClassID]->recvBuffer->recvBegin();
  }

  /// 仮想セル同期データ受信完了待ち.
  void recvVCEnd(int dataClassID) {
    commBufferTableMap[dataClassID]->recvBuffer->recvEnd();
  }

  /// 仮想セル同期データ受信完了待ち(3方向別).
  void recvVCEnd(int dataClassID, int xyz) {
    commBufferSeparateTableMap[xyz][dataClassID]->recvBuffer->recvEnd();
  }

  int faceToXYZ(Face face) {
    switch (face) {
      case X_M: return 0;
      case X_P: return 0;
      case Y_M: return 1;
      case Y_P: return 1;
      case Z_M: return 2;
      case Z_P: return 2;
      default: Exit(EX_FAILURE);
    }
    /* NOTREACHED */
  }

protected:

  /// コンストラクタ.
  ///
  ///  @note シングルトンとして継承可能なようにprivateではなくprotectedにしておく
  ///
  BlockManager()
   : numBlock(0), faceListPrepared(false), separateFaceListPrepared(false) {
    comm = MPI::COMM_WORLD;
  }

  /// デストラクタ.
  /// 登録されているブロックのメモリ解放も行う．
  ///
  ///  @note シングルトンとして継承可能なようにprivateではなくprotectedにしておく
  ///
  ~BlockManager() {
    CommBufferTableMap::iterator it = commBufferTableMap.begin();
    for (; it != commBufferTableMap.end(); it++) delete it->second;
  }

};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // BLOCK_MANAGER_H
