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

#include <iostream>
#include <cstdlib>
#include <sys/stat.h>

#include "BCMTools.h"

#include "Solver.h"

int main(int argc, char** argv)
{
	Solver *pSolver = new Solver();

	int nResultInit = pSolver->Init(argc, argv);
	switch( nResultInit ) {
		case EX_SUCCESS : {
			break;
		}
		case EX_FAILURE : {
			delete pSolver;
			return EX_SUCCESS;
			break;
		}
		default : {
			break;
		}
	}

	int nResultLoop = pSolver->Loop();
	switch( nResultLoop ) {
		case EX_SUCCESS : {
			break;
		}
		case EX_FAILURE : {
			delete pSolver;
			return EX_SUCCESS;
			break;
		}
		default : {
			break;
		}
	}

	int nResultPost = pSolver->Post();
	switch( nResultPost ) {
		case EX_SUCCESS : {
			break;
		}
		case EX_FAILURE : {
			delete pSolver;
			return EX_SUCCESS;
			break;
		}
		default : {
			break;
		}
	}

	delete pSolver;

  return EX_SUCCESS;
}

