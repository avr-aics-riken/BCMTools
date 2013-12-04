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
/// @file FileSystemUtil.h
/// @brief ファイル操作関連ユーティリティ
///

#ifndef __BCMTOOLS_FILESYSTEM_UTIL_H__
#define __BCMTOOLS_FILESYSTEM_UTIL_H__

#include <string>

namespace BCMFileIO {

	/// ファイルパスの "\\" を "/" に変換
	/// return : "/"に変換されたファイルパス
	inline std::string ConvertPath(const std::string& path)
	{
		std::string r = path;
		size_t p = r.find("\\");
		while (p != std::string::npos)
		{
			*(r.begin() + p) = '/';
			p = r.find("\\");
		}
		return r;
	}

	/// pathからディレクトリ名を抜き出す
	/// return : ディレクトリ名 (最後の"/"は残す)
	inline std::string GetDirectory(const std::string& path)
	{
		std::string cpath = ConvertPath(path);

		std::string dir;
		size_t p = cpath.rfind("/");
		if(p != std::string::npos)
		{
			dir = cpath.substr(0, p+1);
		}else{
			dir = std::string("./");
		}
		return dir;
	}
	
	/// pathから拡張子以前のファイル名を抜き出す
	inline std::string GetFilePrefix(const std::string& path)
	{
		std::string cpath = ConvertPath(path);
		std::string filename = cpath.substr(cpath.rfind("/")+1);
		
		return filename.substr(0, filename.rfind("."));
	}

	/// dirをディレクトリ名として整形
	/// 空の場合 : "./"
	/// 文字列が入っている場合最後に "/"を追加
	inline std::string FixDirectoryPath(const std::string& dir)
	{
		if( dir == std::string("") ){
			return std::string("./");
		}else if( dir.rfind("/") != dir.length()-1 ){
			return dir + std::string("/");
		}else{
			return dir;
		}
	}

	/// pathで指定したディレクトリを作成 (mkdir -p 相当)
	/// absolutePath = trueの場合、ルートディレクトリ (/) からディレクトリを作成
	/// ディレクトリ作成に失敗した場合、falseを返す
	bool CreateDirectory(const std::string& path, bool absolutePath = false);


} // namespace BCMFileIO


#endif // __BCMTOOLS_FILESYSTEM_UTIL_H__

