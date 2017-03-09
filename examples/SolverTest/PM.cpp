/*
###################################################################################
#
# BCMTools
#
# Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2016 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
*/

#include "PM.h"
#include <PerfMonitor.h>
#include "gv.h"

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

void PM_Start_L0(const char* name, int id) {
#ifdef __K_FAPP_L0
	fapp_start(name, id, 0);
#endif
}

void PM_Stop_L0(const char* name, int id) {
#ifdef __K_FAPP_L0
	fapp_stop(name, id, 0);
#endif
}

void PM_Start_L1(const char* name, int id) {
#ifdef __K_FAPP_L1
	fapp_start(name, id, 1);
#endif
}

void PM_Stop_L1(const char* name, int id) {
#ifdef __K_FAPP_L1
	fapp_stop(name, id, 1);
#endif
}

void PM_Start_L2(const char* name, int id) {
#ifdef __K_FAPP_L2
	fapp_start(name, id, 2);
#endif
}

void PM_Stop_L2(const char* name, int id) {
#ifdef __K_FAPP_L2
	fapp_stop(name, id, 2);
#endif
}

void PM_Start_L3(const char* name, int id) {
#ifdef __K_FAPP_L3
	fapp_start(name, id, 3);
#endif
}

void PM_Stop_L3(const char* name, int id) {
#ifdef __K_FAPP_L3
	fapp_stop(name, id, 3);
#endif
}

void PM_Start_L4(const char* name, int id) {
#ifdef __K_FAPP_L4
	fapp_start(name, id, 4);
#endif
}

void PM_Stop_L4(const char* name, int id) {
#ifdef __K_FAPP_L4
	fapp_stop(name, id, 4);
#endif
}

void PM_Start_L5(const char* name, int id) {
#ifdef __K_FAPP_L5
	fapp_start(name, id, 5);
#endif
}

void PM_Stop_L5(const char* name, int id) {
#ifdef __K_FAPP_L5
	fapp_stop(name, id, 5);
#endif
}

void PM_Start(unsigned key, int level, int id, bool bBarrier) {
	if( bBarrier == true && g_pconf->BenchMode == true ) {
		MPI_Barrier(MPI_COMM_WORLD);
	}
	g_pPM->start(key);
}

void PM_Stop(unsigned key, int level, int id, double flopPerTask, unsigned iterationCount) {
	g_pPM->stop(key, flopPerTask, iterationCount);
}
