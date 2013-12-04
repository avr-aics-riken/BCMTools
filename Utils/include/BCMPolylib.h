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
/// @file BCMPolylib.h
/// @brief BCM用に拡張し，非シングルトン化したMPIPolylibクラス
/// 

#ifndef BCM_POLYLIB_H
#define BCM_POLYLIB_H

#include "mpi.h"
#include "MPIPolylib.h"
#include "BoundingBox.h"

namespace PolylibNS {

/// BCM用に拡張し，非シングルトン化したMPIPolylibクラス.
///
/// 以下のような使い方を想定:
///   -# コンストラクタによりPolylibインスタンスを生成(全rank)
///   -# 設定ファイルおよびポリゴンデータの読み込み(rank0のみ)
///   -# (rank0におけるツリー構築，BCMブロック配置決定，領域分割)
///   -# rank0から他rankへ，設定ファイル内容と担当領域内ポリゴンデータを送信
///   -# 不要になったら，デストラクタによりPolylibインスタンスとポリゴンデータの解放
///
class BCMPolylib : public MPIPolylib {

  /// Polylib設定ファイル内容の文字列
  std::string m_config_contents;

public:

  /// コンストラクタ.
  ///
  ///   @param[in] comm MPIコミュニケータ
  ///
  BCMPolylib(MPI::Comm& comm = MPI::COMM_WORLD);

  /// デストラクタ.
  ~BCMPolylib();

  /// 設定ファイルおよびポリゴンデータの読み込み.
  ///
  /// rank0の計算ノードのみから呼ぶこと.
  ///
  ///   @param[in] config_filename Polylib設定XMLファイル
  ///
  POLYLIB_STAT load(std::string config_filename, float scale=1.0);

  /// rank0から設定ファイル内容と担当ポリゴンデータを受信.
  ///
  /// rank0以外の計算ノードから呼ぶこと.
  ///
  POLYLIB_STAT load_from_rank0();

  /// 各rankの担当領域のバウンディングボックスを拡大.
  ///
  /// rank0の計算ノードのみから呼ぶこと.
  ///
  ///   @param[in] rank rank番号
  ///   @param[in] min 追加バウンディングボックス下端点座標
  ///   @param[in] max 追加バウンディングボックス上端点座標
  ///
  POLYLIB_STAT set_bounding_box(int rank, const Vec3f& min, const Vec3f& max);

  /// 各rankの担当領域のバウンディングボックスを拡大.
  ///
  /// rank0の計算ノードのみから呼ぶこと.
  ///
  ///   @param[in] rank rank番号
  ///   @param[in] box 追加バウンディングボックス
  ///
  POLYLIB_STAT set_bounding_box(int rank, const BoundingBox& box) {
    Vec3f min(box.getMin().x, box.getMin().y, box.getMin().z);
    Vec3f max(box.getMax().x, box.getMax().y, box.getMax().z);
    return set_bounding_box(rank, min, max);
  }

  /// 各rankに設定ファイル内容と担当ポリゴンデータを送信.
  ///
  /// rank0の計算ノードのみから呼ぶこと.
  /// データ送信後，rank0の保持するポリゴンデータは，
  /// 自分の担当分を残して解放する．
  ///
  POLYLIB_STAT send_to_all();

protected:


private:

  /// 呼び出し禁止にしたMPIPolylibの公開メソッド.
  static MPIPolylib* get_instance();

  /// 呼び出し禁止にしたMPIPolylibの公開メソッド.
  POLYLIB_STAT init_parallel_info(MPI_Comm comm,
                      float bpos[3], unsigned int bbsize[3],
                      unsigned int gcsize[3], float dx[3]);

  /// 呼び出し禁止にしたMPIPolylibの公開メソッド.
  POLYLIB_STAT load_rank0(std::string config_filename = "");

  /// 呼び出し禁止にしたMPIPolylibの公開メソッド.
  POLYLIB_STAT load_parallel(std::string config_filename = "",
                       ID_FORMAT id_format = ID_BIN);

  /// 呼び出し禁止にしたMPIPolylibの公開メソッド.
  POLYLIB_STAT move(PolylibMoveParams &params);

  /// 呼び出し禁止にしたMPIPolylibの公開メソッド.
  POLYLIB_STAT migrate();
  
};

}

#endif // BCM_POLYLIB_H
