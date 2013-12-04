#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include "mpi.h"


/// 設定パラメータクラス.
class Config {

public:

  int level;             ///< ツリーレベル
  std::string treeType;  ///< ツリータイプ

  char type;             ///< 固定境界の向き
  int size;              ///< ブロック内セル数
  int vc;
  double b0, b1;         ///< 境界値
  std::string output;    ///< 結果出力ファイル(未使用)

  int nLoopInner;
  int nLoopOuter;

  double omega;

  bool randomShuffle;    ///< ランダムに領域分割

  bool separate;         ///< 方向毎仮想セル同期フラグ

  bool verbose;          ///< 冗長メッセージフラグ

public:
  Config();

  /// コンストラクタ.
  /// 設定ファイルのパース
  /// @param file 設定ファイル
  Config(const char* file);

  void load(const char* file);

  /// 設定パラメータをコンソール出力.
  void print() const;

  void bcast(const MPI::Comm& comm = MPI::COMM_WORLD);

};

#endif // CONFIG_H
