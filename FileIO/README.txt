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
       2013/02/13



1. はじめに
　BCMTools FileIOはBCMTools向けのファイル入出力機能追加です．


2. BCMToolsに対する変更箇所
　FileIO向けに以下のファイルに変更を加えました．

　変更ファイル
　・examples/util/include/BCMOctree.h
　・examples/util/src/BCMOctree.cpp
　・src/Block/BlockManager.cpp

　変更内容
　[BCMOctree.[h|cpp]]
　・OctreeファイルのペディグリリストからOctreeを復元するルーチン(コンストラクタ)の追加．
　・ファイル出力用にRootGridポインタの取得．
　・メモリ解放ルーチン (メモリリークの修正)
　
　[BlockManager.cpp]
　・Reduce/Allreduceで使用しているMPI data TypeをMPI::INTEGERからMPI::INTに変更
　
3. FileIO用ファイル構成
　ファイル入出力関連のソースコードは，すべて以下に配置されています．
　examples/FileIO/[inlcude|src]
　
　各ファイルの概要を説明します．
　・RLE.h
　  -> ランレングス圧縮/展開ライブラリ
　・BitVoxel.[h|cpp]
　  -> ビットボクセル圧縮/展開ライブラリ
　・ErrorUtil.[h|cpp]
　  -> エラー処理関連関数群
　・FileSystemUtil.[h|cpp]
　  -> ファイル操作関連関数群
　・PartitionMapper.h
　  -> MxNデータロード用データマッピングユーティリティ
　・type.h
　  -> データ型宣言群
　・BCMFileCommon.h
　  -> ファイル入出力共通関数群
　・BCMFileSaver.[h|cpp]
　  -> BCMファイル出力クラス
　・LeafBlockSaver.[h|cpp]
　  -> LeafBlockファイル出力関数群
　・BCMFileLoader.[h|cpp]
　  -> BCMファイル読込クラス
　・LeafBlockLoader.[h|cpp]
　  -> LeafBlockファイル読込関数群


4. ファイル保存機能使用方法
　ファイル保存機能は以下のように使用します．
　
　　1. BCMFileIO::BCMFileSaverのインスタンス生成
　　2. 出力するデータ(ブロック)の情報を登録
　　3. インデックスファイル/Octreeファイル出力
　　4. リーフブロックファイル出力

　=== 例 ===
　// BCMFileSaverのインスタンス生成
　BCMFileIO::BCMFileSaver saver(globalOrigin, // 計算空間全体の起点座標
　                              globalRegion, // 計算空間全体の領域サイズ
　                              octree,       // BCMOctreeのポインタ
　                              "out"         // ファイル出力先ディレクトリ (省略した場合"./")
　                              );
　
　// 出力するタイムステップ情報を設定
　BCMFileIO::IdxStep step(0, 100, 5);  // 0から100まで5刻みのリストを設定
　
　// 出力するデータの情報を登録 (CellID)
　saver.RegisterCellIDInformation( dcid_cid, // データクラスID
　                                 bitWidth, // セルのビット幅
　                                 "CellID", // 系の名称
　                                 "cid",    // ファイル名のPrefix
　                                 "lb",     // 拡張子
　                                 "cid",    // 出力先ディレクトリ ※1
　                                 true      // Gatherフラグ (省略した場合true)
　                                );
　
　// 出力するデータの情報を登録 (Data) ※2
　saver.RegisterDataInformation( &dcid_tmp,              // データクラスID (ポインタ) ※3
　                                BCMFileIO::LB_SCALAR   // データの種別
　                                BCMFioeIO::LB_FLOAT64, // データの型
　                                virtualCell,           // 仮想セルサイズ
　                                "Temperature",         // 系の名称
　                                "tmp",                 // ファイル名のPrefix
　                                "lb",                  // 拡張子
　                                step,                  // タイムステップリスト
　                                "tmp",                 // 出力先ディレクトリ ※1
　                                true                   // タイムステップごとのサブディレクトリフラグ ※4
　                               );
　
　※1 出力先ディレクトリは、コンストラクタで指定したディレクトリからの相対パスとなります．
　
　※2 複数のDataを登録する場合、同様にRegisterDataInformationを複数呼び出します．
　
　※3 複数コンポーネント(ベクトル・テンソル等)の場合、各要素をそれぞれScalar3D<T>として保持するため、 データを登録する
　    際には、要素毎のデータクラスIDを配列にし、先頭アドレスを渡してください．その場合、配列のサイズは、データの種別で
　    指定する要素数と同じにしてください．
　
　※4 タイムステップごとのサブディレクトリフラグをtrueにした場合、出力先として指定したディレクトリの下に、タイムステッ
　    プ毎のディレクトリが作成され、リーフブロックファイルは群はその下に出力されます．
　
　// インデックスファイルとOctreeファイルを出力
　// 登録されたデータ(CellID/Data)の情報とOctreeをファイルに保存
　saver.Save();
　             
　// リーフブロックファイルを出力
　saver.SaveLeafBlock("CellID");        // 出力するデータの系の名称
　saver.SaveLeafBlock("Temprature", 0); // 出力するデータの系の名称, タイムステップ番号
　============
　
5. ファイル読込機能使用方法
　ファイル読込機能は以下のように使用します．
　
　　1. BCMFileIO::BCMFileLoaderのインスタンス生成とメタ情報読込
　　2. 複数のインデックスファイルを読む場合、追加のインデックスファイルを読み込み
　　3. リーフブロックファイル読込
　
　※ メタ情報
　　- インデックスファイルの情報
　　- Octree
　
　
　=== 例 ===
　// BCMFileLoaderのインスタンス生成とメタ情報
　BCMFileIO::BCMFileLoader loader("cellId.bcm", // インデックスファイルのパス
　                                bcsetter      // 境界条件設定ファンクタ
　                               );
　
　// 追加のインデックスファイルの読み込み ※1
　loader.LoadAdditionalIndex("data.bcm");
　
　// リーフブロックファイル読込 (CellID)
　loader.LoadLeafBlock(&cid_dcid, // データクラスID ※2
　                     "CellID",  // 系の名称
　                     vc         // 仮想セルサイズ
　                     );         // returnはデータクラスID
　
　// リーフブロックファイル読込 (Data)
　loader.LoadLeafBlock(&tmp_dcid,     // データクラスID ※2
　                     "Temperature", // 系の名称
　                     vc,            // 仮想セルサイズ
　                     0              // タイムステップ
　                     );             // returnはデータクラスID
　
　※2 LoadLeafBlockを実行すると、系の名称に対応するブロックが生成されます．生成したブロックのデータクラスIDは、
　    LoadLeafBlockの第一引数として出力されます．また、すでに作成されている場合は、そのデータクラスIDを出力します．
　    系の名称に対応するデータが複数コンポーネントの場合、各要素に対応したデータクラスIDを出力するため、第一引数に渡す
　    バッファは配列として、事前に確保してください．
　    また、データを読み込まずブロックの生成だけを実行したい場合、loader.CreateLeafBlock(&dcid, name, vc)を呼び出してく
　    ださい．
　
　==========
　
　
6. サンプルプログラム
　BCMTools FileIOを使ったサンプルを2種類作成しました。
　
　・SampleCreator
　  -> ParallelMeshGenerationのメッシュ生成機能ルーチンを使い、メッシュ生成を行い、
　     結果(CellID/Scalar/Vecotr)をファイルとして出力
　
　・SampleLoader
　  -> SampleCreatorで作成したBCMファイル群を読み込み、結果をファイルに出力
　
