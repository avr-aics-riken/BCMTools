#include <iostream>
#include <cstdlib>
#include <sys/stat.h>

#include "BCMTools.h"

#include "Solver.h"

int main(int argc, char** argv)
{
	mkdir("./VTI", 0755);
	mkdir("./BIN", 0755);

	Solver *pSolver = new Solver();

	bool bResult = false;

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
		case 1 : {
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
		case 1 : {
			break;
		}
		default : {
			break;
		}
	}

	delete pSolver;

  return EX_SUCCESS;
}

