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

#ifndef PM_H
#define PM_H

#ifdef __K_FAPP_L0
#include <fj_tool/fapp.h>
#endif

#ifdef __K_FAPP_L1
#include <fj_tool/fapp.h>
#endif

#ifdef __K_FAPP_L2
#include <fj_tool/fapp.h>
#endif

#ifdef __K_FAPP_L3
#include <fj_tool/fapp.h>
#endif

#ifdef __K_FAPP_L4
#include <fj_tool/fapp.h>
#endif

#ifdef __K_FAPP_L5
#include <fj_tool/fapp.h>
#endif


void PM_Start_L0(const char* name, int id);
void PM_Stop_L0(const char* name, int id);
void PM_Start_L1(const char* name, int id);
void PM_Stop_L1(const char* name, int id);
void PM_Start_L2(const char* name, int id);
void PM_Stop_L2(const char* name, int id);
void PM_Start_L3(const char* name, int id);
void PM_Stop_L3(const char* name, int id);
void PM_Start_L4(const char* name, int id);
void PM_Stop_L4(const char* name, int id);

enum timing_key {
	tm_Init_LoadSTL,
	tm_Init_DivideDomain,
	tm_Init_OrderBlock,
	tm_Update,
	tm_UpdateT,
	tm_UpdateUX,
	tm_UpdateUY,
	tm_UpdateUZ,
	tm_UpdateP,
	tm_UpdateU,
	tm_DOT,
	tm_DOT_Calc,
	tm_DOT_Comm,
	tm_JacobiSmoother,
	tm_JacobiSmoother_Calc,
	tm_JacobiSmoother_Swap,
	tm_JacobiSmoother_Comm,
	tm_END,
};

void PM_Start(unsigned key, int level=0, int id=0, bool bBarrier=false);
void PM_Stop(unsigned key, int level=0, int id=0, double flopPerTask=0.0, unsigned iterationCount=1);

#endif

