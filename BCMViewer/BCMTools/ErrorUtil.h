/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file  ErrorUtil.h
/// @brief エラー処理関連のユーティリティ
///

#ifndef __BCMTOOLS_ERROR_UTIL_H__
#define __BCMTOOLS_ERROR_UTIL_H__

//#include <mpi.h>

namespace BCMFileIO {

	/// ログ出力 (エラー)
	void LogE( const char *format, ...);
	
	/// ログ出力 (警告)
	void LogW( const char *format, ...);
	
	/// ログ出力 (情報)
	void LogI( const char *format, ...);
	
	/// ログ出力 (デバッグ情報)
	void LogD( const char *format, ...);

	/// 全プロセスに対しエラー情報を配信
	/// 
	/// @param[in] err  エラーがある場合trueを入力
	/// @param[in] comm MPIコミュニケータ
	/// @return 1プロセスでもエラーがある場合trueを返す．全プロセスでエラーが無い場合false
	/// 
	//bool reduceError( const bool err, MPI::Intracomm& comm = MPI::COMM_WORLD );
	
} // namespace BCMFileIO

#endif // __BCMTOOLS_ERROR_UTIL_H__

