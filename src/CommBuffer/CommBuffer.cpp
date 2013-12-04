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

#include "CommBuffer.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif


// CommBuffer -- 通信バッファクラス(通信機能なし). //////////

/// コンストラクタ.
CommBuffer::CommBuffer() : nPeer(0), status(0), request(0)
{
}


/// デストラクタ.
CommBuffer::~CommBuffer()
{
  BufferTable::iterator it = bufferTable.begin();
  for (; it != bufferTable.end(); ++it) {
    delete it->second;
  }
  delete[] status;
  delete[] request;
}


/// バッファにT型データ領域を登録.
template <typename T>
void CommBuffer::setData(int rank, size_t size, T** pointer)
{
  size = size * sizeof(T); // sizeをバイト単位に変換
  addBufferSize(rank, size, new PointerSetter<T>(pointer));
}


/// バッファにデータ領域を登録(PointerSetter使用).
void CommBuffer::setData(int rank, size_t size, PointerSetterBase* pointerSetter)
{
  // sizeはバイト単位
  if (size % 8 != 0) size = size - (size % 8) + 8; // 64ビット境界にパディング
  BufferTable::iterator it = bufferTable.find(rank);
  if (it == bufferTable.end()) {
    bufferTable[rank] = new Buffer;
  }
  bufferTable[rank]->pointerTable.push_back(pointerSetter);
  bufferTable[rank]->offsetTable.push_back(bufferTable[rank]->size);
  bufferTable[rank]->size += size;
  bufferTable[rank]->dirty = true;
}


/// バッファデータ領域の確保.
int CommBuffer::allocateBuffer()
{
  BufferTable::iterator it = bufferTable.begin();
  for (; it != bufferTable.end(); ++it) {
    Buffer* buffer = it->second;
    if (buffer->dirty) {
      delete[] buffer->data;
      buffer->data = new Byte[buffer->size];
      for (int i = 0; i < buffer->pointerTable.size(); ++i) {
        buffer->pointerTable[i]->setPointer(buffer->data + buffer->offsetTable[i]);
      }
      buffer->dirty = false;
    }
  }

  nPeer = bufferTable.size();
  if (nPeer > 0) {
    delete[] status;
    delete[] request;
    status = new MPI::Status[nPeer];
    request = new MPI::Request[nPeer];
  }
  return nPeer;
}


/// ランクを指定してバッファを消去.
void CommBuffer::deleteBuffer(int rank)
{
  BufferTable::iterator it = bufferTable.find(rank);
  if (it != bufferTable.end()) {
    delete it->second;
    bufferTable.erase(it);
  } 
  nPeer = bufferTable.size();
}


/// 全バッファを消去.
void CommBuffer::deleteBuffer()
{
  BufferTable::iterator it = bufferTable.begin();
  for (; it != bufferTable.end(); ++it) {
    delete it->second;
  }
  bufferTable.clear();
  nPeer = 0;
}


// SendBuffer -- 送信用通信バッファクラス. //////////

/// コンストラクタ.
SendBuffer::SendBuffer(int tag, const MPI::Comm& comm) : tag(tag), comm(comm)
{
}


/// デストラクタ.
SendBuffer::~SendBuffer()
{
}


/// データ送信.
void SendBuffer::send()
{
  BufferTable::const_iterator it = bufferTable.begin();
  for (; it != bufferTable.end(); ++it) {
    comm.Send(it->second->data, it->second->size, MPI::BYTE, it->first, tag);
  }
}


/// ノンブロッキング送信開始.
void SendBuffer::sendBegin()
{
  BufferTable::const_iterator it = bufferTable.begin();
  for (int i = 0; it != bufferTable.end(); ++it, ++i) {
    request[i] = comm.Isend(it->second->data, it->second->size,
                            MPI::BYTE, it->first, tag);
  }
}


/// ノンブロッキング送信終了待ち.
void SendBuffer::sendEnd()
{
  MPI::Request::Waitall(nPeer, request, status);
}


// RecvBuffer -- 受信用通信バッファクラス. //////////

/// コンストラクタ.
RecvBuffer::RecvBuffer(int tag, const MPI::Comm& comm) : tag(tag), comm(comm)
{
}


/// デストラクタ.
RecvBuffer::~RecvBuffer()
{
}


/// データ受信.
void RecvBuffer::recv()
{
  BufferTable::iterator it = bufferTable.begin();
  for (int i = 0; it != bufferTable.end(); ++it, ++i) {
    comm.Recv(it->second->data, it->second->size, MPI::BYTE, it->first, tag,
              status[i]);
  }
}


/// ノンブロッキング受信開始.
void RecvBuffer::recvBegin()
{
  BufferTable::iterator it = bufferTable.begin();
  for (int i = 0; it != bufferTable.end(); ++it, ++i) {
    request[i] = comm.Irecv(it->second->data, it->second->size,
                            MPI::BYTE, it->first, tag);
  }
}


/// ノンブロッキング受信終了待ち.
void RecvBuffer::recvEnd()
{
  MPI::Request::Waitall(nPeer, request, status);
}


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif
