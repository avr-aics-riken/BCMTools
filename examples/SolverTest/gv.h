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

#include "Config.h"
#include <PerfMonitor.h>

#ifdef GLOBAL_VALUE_DEFINE
#define GLOBAL_VARIABLE
#define GLOBAL_VALUE(v) = (v)
#else
#define GLOBAL_VARIABLE extern
#define GLOBAL_VALUE(v)
#endif

GLOBAL_VARIABLE Config* g_pconf;
GLOBAL_VARIABLE pm_lib::PerfMonitor* g_pPM;

#undef GLOBAL_VARIABLE
#undef GLOBAL_VALUE


void PrintLog(int myrank, int level, const char* format, ...);
double RandomUniform();
double RandomNormal(double mu, double sigma);

