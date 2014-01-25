/*
 * BCMTools
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file CommBuffer.h
/// @brief 通信バッファクラス
///

#ifndef COMM_BUFFER_H
#define COMM_BUFFER_H

#include <map>
#include <vector>
#include "mpi.h"
#include "PointerSetter.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

///
/// 通信バッファクラス(通信機能なし).
///
class CommBuffer {

private:

  /// Byte型.
//typedef uint8_t Byte;
  typedef unsigned char Byte;

protected:

  /// バッファ管理構造体.
  struct Buffer {
    size_t size;    ///< バッファサイズ(byte)
    Byte* data;     ///< データ領域
    bool dirty;     ///< 要再アロケートフラグ
    std::vector<PointerSetterBase*> pointerTable;   ///< PointerSetterテーブル
    std::vector<size_t> offsetTable;    ///< データ領域先頭から各ポインタへのByte数

    /// コンストラクタ.
    Buffer() : size(0), data(0), dirty(true) {}

    /// デストラクタ.
    ~Buffer() { 
      delete[] data;
      for (int i = 0; i < pointerTable.size(); ++i) {
        pointerTable[i]->setPointer(0);
        delete pointerTable[i];
      }
    }
  };

  /// バッファテーブル型(「ランク番号→バッファ管理構造体」のマップ).
  typedef std::map<int, Buffer*> BufferTable;

  BufferTable bufferTable;  ///< バッファテーブル

  int nPeer;   ///< 通信相手ノード数
  MPI::Status* status;     ///< MPIステータス配列
  MPI::Request* request;   ///< MPIリクエスト配列

//public:
protected: 

  /// コンストラクタ.
  CommBuffer();

  /// デストラクタ.
  ~CommBuffer();

public: 

  /// バッファにT型データ領域を登録.
  ///
  ///  @param[in] rank ランク番号
  ///  @param[in] size 領域サイズ(T型の配列長)
  ///  @param[in] pointer 「データ領域先頭へのポインタ」へのポインタ
  ///
  ///  @note バッファ領域確保後にデータ領域先頭アドレスが*pointerに格納される.
  ///
  template <typename T>
  void setData(int rank, size_t size, T** pointer);

  /// バッファにデータ領域を登録(PointerSetter使用).
  ///
  ///  @param[in] rank ランク番号
  ///  @param[in] size 領域サイズ(byte)
  ///  @param[in] pointerSetter 指定された型のPointerSetter
  ///
  ///  @note pointerSetterはCommBufferのデストラクタ内で解放される.
  ///
  void setData(int rank, size_t size, PointerSetterBase* pointerSetter);

  /// バッファデータ領域の確保.
  ///
  ///  @return 通信相手ノード数
  ///
  int allocateBuffer();

  /// ランクを指定してバッファを消去.
  ///
  ///  @param[in] rank ランク番号
  ///
  void deleteBuffer(int rank);

  /// 全バッファを消去.
  void deleteBuffer();

private:

  /// コピーコンストラクタ(コピー禁止).
  CommBuffer(const CommBuffer& rhs);

  /// 代入演算子(コピー禁止).
  CommBuffer& operator=(const CommBuffer& rhs);

};


///
/// 送信用通信バッファクラス.
///
class SendBuffer : public CommBuffer {

  int tag;                ///< タグ値
  const MPI::Comm& comm;  ///< コミュニケータ

public:

  /// コンストラクタ.
  ///
  ///  @param[in] tag タグ値
  ///  @param[in] comm コミュニケータ(ディフォルトMPI_COMM_WORLD)
  ///
  SendBuffer(int tag, const MPI::Comm& comm = MPI::COMM_WORLD);

  /// デストラクタ.
  ~SendBuffer();

  /// データ送信.
  void send();

  /// ノンブロッキング送信開始.
  void sendBegin();

  /// ノンブロッキング送信終了待ち.
  void sendEnd();

};


///
/// 受信用通信バッファクラス.
///
class RecvBuffer : public CommBuffer {

  int tag;                ///< タグ値
  const MPI::Comm& comm;  ///< コミュニケータ

public:

  /// コンストラクタ.
  ///
  ///  @param[in] tag タグ値
  ///  @param[in] comm コミュニケータ(ディフォルトMPI_COMM_WORLD)
  ///
  RecvBuffer(int tag, const MPI::Comm& comm = MPI::COMM_WORLD);

  /// デストラクタ.
  ~RecvBuffer();

  /// データ受信.
  void recv();

  /// ノンブロッキング受信開始.
  void recvBegin();

  /// ノンブロッキング受信終了待ち.
  void recvEnd();

};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // COM_BUFFER_H

