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
/// @file  ErrorUtil.cpp
/// @brief エラー処理関連のユーティリティ
///

#include "ErrorUtil.h"

#include <stdio.h>
#include <stdarg.h>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

namespace BCMFileIO {

	enum LOG_LEVEL
	{
		LOG_ERROR = 0,
		LOG_WARN,
		LOG_INFO,
		LOG_DEBUG
	};
	
	
	void Log( enum LOG_LEVEL level, const std::string& msg )
	{
		static const char *log_header[4] = {
			"[ERR]", "[WRN]", "[INF]", "[DBG]"
		};
	
		MPI::Intracomm& comm = MPI::COMM_WORLD;
		char rankstr[128];
		sprintf(rankstr, "[RANK:%6d]", comm.Get_rank());
		
		std::string logstr(rankstr);
		logstr += std::string(" ") + std::string(log_header[level]) + std::string(" ") + msg;
	
	#ifdef _WIN32
		OutputDebugStringA(logstr.c_str());
	#else
		printf("%s", logstr.c_str());	
	#endif
	}
	
	void LogE (const char *format, ...)
	{
		char buf[256];
		va_list arg;
		va_start(arg, format);
		vsprintf(buf, format, arg);
		va_end(arg);
		
		Log(LOG_ERROR, buf);
	}
	
	void LogW (const char *format, ...)
	{
		char buf[256];
		va_list arg;
		va_start(arg, format);
		vsprintf(buf, format, arg);
		va_end(arg);
		
		Log(LOG_WARN, buf);
	}
	
	void LogI (const char *format, ...)
	{
		char buf[256];
		va_list arg;
		va_start(arg, format);
		vsprintf(buf, format, arg);
		va_end(arg);
		
		Log(LOG_INFO, buf);
	}
	
	void LogD (const char *format, ...)
	{
		char buf[256];
		va_list arg;
		va_start(arg, format);
		vsprintf(buf, format, arg);
		va_end(arg);
		
		Log(LOG_DEBUG, buf);
	}


	bool reduceError(const bool err, MPI::Intracomm& comm) {
		int ierr = err ? 1 : 0;
		comm.Allreduce(&ierr, &ierr, 1, MPI::INT, MPI::BOR);
		return ierr == 0 ? false : true;
	}


} // namespace BCMFileIO
	
