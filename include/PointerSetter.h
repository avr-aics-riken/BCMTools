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
/// @file PointerSetter.h
/// @brief 通信バッファクラス用のポインタキャストクラス
///

#ifndef POINTER_SETTER_H
#define POINTER_SETTER_H

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

///
/// PointerSetterテンプレートクラスのポリモーフィズムのための基底クラス.
///
class PointerSetterBase {

public:

  /// デストラクタ.
  virtual ~PointerSetterBase() {}

  /// ポインタのキャスト.
  ///
  ///   @param[in] voidPointer void型ポインタ
  ///
  virtual void setPointer(void* voidPointer) = 0;

};


///
/// 任意の基本型Tへのポインタにキャストするためのテンプレートクラス.
///
template <typename T>
class PointerSetter : public PointerSetterBase {

  T** pointer;  ///< 「アドレス値設定対象となるポインタ」へのポインタ

public:

  /// コンストラクタ.
  ///
  ///  @param[in] pointer 「アドレス値設定対象となるポインタ」へのポインタ
  ///
  PointerSetter(T** pointer) : pointer(pointer) {}

  /// ポインタ値のT*型へのキャスト.
  ///
  ///   @param[in] voidPointer void型ポインタ
  ///
  void setPointer(void* voidPointer) {
    *pointer = static_cast<T*>(voidPointer);
  }

};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif
