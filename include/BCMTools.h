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
/// @file BCMTools.h
/// @brief BCM Tools 共通ヘッダ
///

#ifndef BCM_TOOLS_H
#define BCM_TOOLS_H

#include<cstdio>
#include<cstdlib>

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

/// 呼び出し箇所が分かるexit関数マクロ(assertの代わりに使用).
#define Exit(x) \
  ((void)printf("exit at %s:%u\n", __FILE__, __LINE__), exit((x)))

/// DEBUGマクロの定義時のみassertマクロを有効に.
#ifndef DEBUG
#define NDEBUG
#endif
#include <cassert>

/// フェイス番号.
enum Face { X_M, X_P, Y_M, Y_P, Z_M, Z_P, NUM_FACE };

/// サブフェイス番号.
enum Subface { SF_00, SF_01, SF_10, SF_11, NUM_SUBFACE };

/// リターンコード.
enum ExitStatus {
    EX_SUCCESS     = 0,    ///< successful termination
    EX_USAGE       = 16,   ///< incorrect command arguments
    EX_MEMORY      = 17,   ///< memory allocation failure
    EX_OPEN_FILE   = 18,   ///< open file error
    EX_READ_CONFIG = 19,   ///< read configuration file error
    EX_READ_DATA   = 20,   ///< read data error
    EX_WRITE_DATA  = 21,   ///< write data error
    EX_FAILURE     = 1,    ///< other failure
};

#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // BCM_TOOLS_H
