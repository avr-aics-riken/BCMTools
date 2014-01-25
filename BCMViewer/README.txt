/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */


1. はじめに
BCMViewerはBCMTools FileIOにより出力されたCellIDデータを読み込み、
スライス断面表示による可視化を行う簡易ツールです．

2. 動作プラットホーム
OpenGL 2.1 以降が動作するLinuxとMac OS X

3. ビルド方法

3.1 Makefileの編集
 13行目の"TEXTPARSER = "のパスをTextParser-1.2のインストールされて
 いるディレクトリに変更

3.2 ビルド
 - ターミナルを起動し、BCMViewerのディレクトリへ移動
 - makeを実行
 - 実行ファイル"main"が出力されればビルド成功

4. 起動方法
以下の要領で起動
> cd ${BCMViewer}
> ./main <cellid.bcmファイルのファイルパス>


5. 操作方法

- 視点操作
　マウス左ボタンドラッグ        : 回転
　Shft + マウス右ボタンドラッグ : 平行移動
　Ctrl + マウス左ボタンドラッグ : 拡大/縮小

- スライス断面操作
　x/y/z キー  : 操作対象の面を変更
　              (xキーを入力した場合、X断面が操作対象)
　nキー       : 操作対象の面を+1スライド
　Nキー       : 操作対象の面を+10スライド
　pキー       : 操作対象の面を-1スライド
　Pキー       : 操作対象の面を-10スライド
　Dキー       : 操作対象の面を非表示に設定
　gキー       : リーフブロックのグリッド表示をON/OFF(トグル)

- その他
　ESCキー     : 終了


