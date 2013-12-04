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
/// @file ConfigBase.h
/// @brief 設定パラメータクラス(基底クラス)
/// 

#ifndef CONFIG_BASE_H
#define CONFIG_BASE_H

#include "ConfigFile.h"
#include "mpi.h"


/// 設定パラメータクラス(基底クラス).
///
/// このクラスを次ぎのように継承してカスタマイズすること
///  @li 設定パラメータをpublicメンバに追加
///  @li parseメソッド内で、readメソッドを用いて各パラメータ値を入力
///  @li 必要なら、validateメソッド内でパラメータ値の検証
///
class ConfigBase {

  MPI::Comm& comm;   ///< MPIコミュニケータ

  ConfigFile* configFile;  ///< ConfigFileオブジェクト

public:
  /// コンストラクタ.
  ///
  ///  @param[in] comm MPIコミュニケータ
  ///
  ConfigBase(MPI::Comm& comm = MPI::COMM_WORLD);

  /// デストラクタ.
  ~ConfigBase();

  /// 設定ファイル読み込み.
  ///
  ///  @param[in] file 設定ファイルパス
  ///
  ///  @note rank0のみがファイルから読み込み，他rankに転送する.
  ///
  void load(const char* file);


private:

  /// パラメータのパース(派生クラスでカスタマイズ).
  virtual void parse() = 0;

  /// パラメータ値のチェック(派生クラスでカスタマイズ).
  virtual bool validate() { return true; }

  /// ConfigFileオブジェクトの内容をrank0からブロードキャスト.
  void broadcastConfigFile(const ConfigFile* configFile);

  /// ConfigFileオブジェクトの内容をrank0から受信.
  void receiveConfigFile(ConfigFile* configFile);

protected:

  /// keyに対応したパラメータの読み込み.
  template<class T> T read(const std::string& key) const {
    return configFile->read<T>(key);
  }

  /// keyに対応したパラメータの読み込み(ディフォルト値あり).
  template<class T> T read(const std::string& key, const T& value) const {
    return configFile->read<T>(key, value);
  }

  /// エラー終了.
  void errorExit(const char* message, int code = 1); 

};


#endif // CONFIG_BASE_H
