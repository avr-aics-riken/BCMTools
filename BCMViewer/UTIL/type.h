/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __OZ_TYPE_H__
#define __OZ_TYPE_H__


#ifdef _WIN32

typedef short          int16_t;
typedef unsigned short uint16_t;
typedef int            int32_t;
typedef unsigned int   uint32_t;
typedef long           int64_t;
typedef unsigned long  uint64_t;

#else  // _WIN32

#include <stdint.h>

#endif // _WIN32

typedef bool          b8;
typedef char          s8;
typedef unsigned char u8;
typedef int16_t       s16;  
typedef uint16_t      u16;
typedef int32_t       s32;
typedef uint32_t      u32;
typedef int64_t       s64;
typedef uint64_t      u64;
typedef float         f32;
typedef double        f64;

#endif // __OZ_TYPE_H__
