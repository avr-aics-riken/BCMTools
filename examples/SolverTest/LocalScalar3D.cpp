#include "LocalScalar3D.h"

#include "bsf3d.h"

template <>
void LocalScalar3D<real>::CalcStats(BlockManager& blockManager) {
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		BlockBase* block = blockManager.getBlock(id);
		::Vec3i size = block->getSize();
		::Vec3r origin = block->getOrigin();
		::Vec3r blockSize = block->getBlockSize();
		::Vec3r cellSize = block->getCellSize();

		int sz[3] = {size.x, size.y, size.z};
		int g[1] = {vc};
		real dx = cellSize.x;
	
		real* pData = GetBlockData(block);
		if( id==0 ) {
			int m0 = vc + (size.x+2*vc)*(vc + (size.y+2*vc)*vc);
			sum = 0.0;
			max = pData[m0];
			min = pData[m0];
			absmax = abs(pData[m0]);
			absmin = abs(pData[m0]);
		}
		sf3d_calc_stats_(
					&sum,
					&max,
					&min,
					&absmax,
					&absmin,
					pData,
					sz, g);
	}

	real sum_tmp = sum;
#ifdef _REAL_IS_DOUBLE_		
	MPI_Allreduce(&sum_tmp, &sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&sum_tmp, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif

	real max_tmp = max;
#ifdef _REAL_IS_DOUBLE_		
	MPI_Allreduce(&max_tmp, &max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&max_tmp, &max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
#endif

	real min_tmp = min;
#ifdef _REAL_IS_DOUBLE_		
	MPI_Allreduce(&min_tmp, &min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&min_tmp, &min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif

	real absmax_tmp = absmax;
#ifdef _REAL_IS_DOUBLE_		
	MPI_Allreduce(&absmax_tmp, &absmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&absmax_tmp, &absmax, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
#endif

	real absmin_tmp = absmin;
#ifdef _REAL_IS_DOUBLE_		
	MPI_Allreduce(&absmin_tmp, &absmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&absmin_tmp, &absmin, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif

/*
	real n_grids = (real)(m_iNX*m_iNY*m_iNZ);
	real n_grids_tmp = n_grids;
#ifdef _REAL_IS_DOUBLE_
	MPI_Allreduce(&n_grids_tmp, &n_grids, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
#else
	MPI_Allreduce(&n_grids_tmp, &n_grids, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif
*/

/*
	m_rAve = sum/n_grids;
	m_rMax = max;
	m_rMin = min;
	m_rAbsMax = absmax;
	m_rAbsMin = absmin;
*/
}

template <>
void LocalScalar3D<real>::Dump(BlockManager& blockManager, const int step, const char* label) {
	ImposeBoundaryCondition(blockManager);
	MPI::Intracomm comm = blockManager.getCommunicator();

#ifdef _LARGE_BLOCK_
#else
#endif
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		BlockBase* block = blockManager.getBlock(id);

		::Vec3i size = block->getSize();
		Vec3r origin = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize = block->getCellSize();
		int level = block->getLevel();

		ostringstream ossFileName;
		ossFileName << "./BIN/";
		ossFileName << "dump-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << comm.Get_rank();
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << id;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".bin";

		int cx = size.x + 2*vc;
		int cy = size.y + 2*vc;
		int cz = size.z + 2*vc;
		int iNE = 1;

		real* pData = GetBlockData(block);

		ofstream ofs;
		ofs.open(ossFileName.str().c_str(), ios::out | ios::binary);
		ofs.write((char*)&size.x, sizeof(int));
		ofs.write((char*)&size.y, sizeof(int));
		ofs.write((char*)&size.z, sizeof(int));
		ofs.write((char*)&vc    , sizeof(int));
		ofs.write((char*)&iNE   , sizeof(int));
		ofs.write((char*)pData  , sizeof(real)*cx*cy*cz);
		ofs.close();
	}
}

template <>
void LocalScalar3D<real>::Load(BlockManager& blockManager, const int step, const char* label) {
	MPI::Intracomm comm = blockManager.getCommunicator();

#ifdef _LARGE_BLOCK_
#else
#endif
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		BlockBase* block = blockManager.getBlock(id);

		::Vec3i size = block->getSize();
		Vec3r origin = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize = block->getCellSize();
		int level = block->getLevel();

		ostringstream ossFileName;
		ossFileName << "./BIN/";
		ossFileName << "dump-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << comm.Get_rank();
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << id;
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".bin";

		int nx = 0;
		int ny = 0;
		int nz = 0;
		int nv = 0;
		int ne = 0;

		real* pData = GetBlockData(block);

		ifstream ifs;
		ifs.open(ossFileName.str().c_str(), ios::in | ios::binary);
		ifs.read((char*)&nx, sizeof(int));
		ifs.read((char*)&ny, sizeof(int));
		ifs.read((char*)&nz, sizeof(int));
		ifs.read((char*)&nv, sizeof(int));
		ifs.read((char*)&ne, sizeof(int));
		if( nx == size.x && ny == size.y && nz == size.z && nv == vc && ne == 1 ) {
			int cx = nx + 2*nv;
			int cy = ny + 2*nv;
			int cz = nz + 2*nv;
			ifs.read((char*)pData, sizeof(real)*cx*cy*cz);
		} else if( 2*nx == size.x && 2*ny == size.y && 2*nz == size.z && nv == vc && ne == 1 ) {
			int cx = nx + 2*nv;
			int cy = ny + 2*nv;
			int cz = nz + 2*nv;
			real *pDataS = new real [cx*cy*cz];
			ifs.read((char*)pDataS, sizeof(real)*cx*cy*cz);

			int sz[3] = {2*nx, 2*ny, 2*nz};
			sf3d_copy_x2_(
					(real*)pData,
					(real*)pDataS,
					(int*)sz,
					(int*)&vc);

			delete [] pDataS;
		} else {
			Exit(0);
		}
		ifs.close();
	}

	ImposeBoundaryCondition(blockManager);
}

template <>
void LocalScalar3D<real>::Dump2(BlockManager& blockManager, const int step, const char* label) {
	ImposeBoundaryCondition(blockManager);
	MPI::Intracomm comm = blockManager.getCommunicator();

	ostringstream ossFileName;
	ossFileName << "./BIN/";
	ossFileName << "dump-";
	ossFileName << label;
	ossFileName << "-";
	ossFileName.width(5);
	ossFileName.setf(ios::fixed);
	ossFileName.fill('0');
	ossFileName << comm.Get_rank();
	ossFileName << "-";
	ossFileName.width(10);
	ossFileName.setf(ios::fixed);
	ossFileName.fill('0');
	ossFileName << step;
	ossFileName << ".bin";

	BlockBase* block = blockManager.getBlock(0);
	::Vec3i size = block->getSize();
	int cx = size.x + 2*vc;
	int cy = size.y + 2*vc;
	int cz = size.z + 2*vc;
	int iNE = 1;
	int iNB = blockManager.getNumBlock();

	ofstream ofs;
	ofs.open(ossFileName.str().c_str(), ios::out | ios::binary);
	ofs.write((char*)&size.x, sizeof(int));
	ofs.write((char*)&size.y, sizeof(int));
	ofs.write((char*)&size.z, sizeof(int));
	ofs.write((char*)&vc    , sizeof(int));
	ofs.write((char*)&iNE   , sizeof(int));
	ofs.write((char*)&iNB   , sizeof(int));

#ifdef _LARGE_BLOCK_
#else
#endif
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		block = blockManager.getBlock(id);

		real* pData = GetBlockData(block);

		ofs.write((char*)pData  , sizeof(real)*cx*cy*cz);
	}

	ofs.close();
}

template <>
void LocalScalar3D<real>::Load2(BlockManager& blockManager, const int step, const char* label) {
	MPI::Intracomm comm = blockManager.getCommunicator();

	ostringstream ossFileName;
	ossFileName << "./BIN/";
	ossFileName << "dump-";
	ossFileName << label;
	ossFileName << "-";
	ossFileName.width(5);
	ossFileName.setf(ios::fixed);
	ossFileName.fill('0');
	ossFileName << comm.Get_rank();
	ossFileName << "-";
	ossFileName.width(10);
	ossFileName.setf(ios::fixed);
	ossFileName.fill('0');
	ossFileName << step;
	ossFileName << ".bin";

	int nx = 0;
	int ny = 0;
	int nz = 0;
	int nv = 0;
	int ne = 0;
	int nb = 0;

	ifstream ifs;
	ifs.open(ossFileName.str().c_str(), ios::in | ios::binary);
	ifs.read((char*)&nx, sizeof(int));
	ifs.read((char*)&ny, sizeof(int));
	ifs.read((char*)&nz, sizeof(int));
	ifs.read((char*)&nv, sizeof(int));
	ifs.read((char*)&ne, sizeof(int));
	ifs.read((char*)&nb, sizeof(int));

	int cx = nx + 2*nv;
	int cy = ny + 2*nv;
	int cz = nz + 2*nv;

	BlockBase* block = blockManager.getBlock(0);
	::Vec3i size = block->getSize();

	if( nx == size.x && ny == size.y && nz == size.z && nv == vc && ne == 1 && nb == blockManager.getNumBlock() ) {
		for (int id = 0; id < blockManager.getNumBlock(); ++id) {
			block = blockManager.getBlock(id);

			real* pData = GetBlockData(block);

			ifs.read((char*)pData, sizeof(real)*cx*cy*cz);
		}
	} else if( 2*nx == size.x && 2*ny == size.y && 2*nz == size.z && nv == vc && ne == 1 && nb == blockManager.getNumBlock() ) {
		real *pDataS = new real [cx*cy*cz];
		for (int id = 0; id < blockManager.getNumBlock(); ++id) {
			block = blockManager.getBlock(id);

			real* pData = GetBlockData(block);

			ifs.read((char*)pDataS, sizeof(real)*cx*cy*cz);

			int sz[3] = {2*nx, 2*ny, 2*nz};
			sf3d_copy_x2_(
					(real*)pData,
					(real*)pDataS,
					(int*)sz,
					(int*)&vc);
		}
		delete [] pDataS;
	} else {
		Exit(0);
	}

	ifs.close();

	ImposeBoundaryCondition(blockManager);
}

template <>
void LocalScalar3D<real>::Dump3(BlockManager& blockManager, const int step, const char* label, Partition* partition, int myrank) {
	ImposeBoundaryCondition(blockManager);
	MPI::Intracomm comm = blockManager.getCommunicator();

#ifdef _LARGE_BLOCK_
#else
#endif
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		BlockBase* block = blockManager.getBlock(id);

		::Vec3i size = block->getSize();
		Vec3r origin = block->getOrigin();
		Vec3r blockSize = block->getBlockSize();
		Vec3r cellSize = block->getCellSize();
		int level = block->getLevel();

		ostringstream ossFileName;
		ossFileName << "./BIN/";
		ossFileName << "dump-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << id + partition->getStart(myrank);
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".bin";

		int cx = size.x + 2*vc;
		int cy = size.y + 2*vc;
		int cz = size.z + 2*vc;
		int iNE = 1;

		real* pData = GetBlockData(block);

		ofstream ofs;
		ofs.open(ossFileName.str().c_str(), ios::out | ios::binary);
		ofs.write((char*)&size.x, sizeof(int));
		ofs.write((char*)&size.y, sizeof(int));
		ofs.write((char*)&size.z, sizeof(int));
		ofs.write((char*)&vc    , sizeof(int));
		ofs.write((char*)&iNE   , sizeof(int));
		ofs.write((char*)pData  , sizeof(real)*cx*cy*cz);
		ofs.close();
	}
}

template <>
void LocalScalar3D<real>::Load3(BlockManager& blockManager, const int step, const char* label, Partition* partition, int myrank) {
	MPI::Intracomm comm = blockManager.getCommunicator();

#ifdef _LARGE_BLOCK_
#else
#endif
	for (int id = 0; id < blockManager.getNumBlock(); ++id) {
		ostringstream ossFileName;
		ossFileName << "./BIN/";
		ossFileName << "dump-";
		ossFileName << label;
		ossFileName << "-";
		ossFileName.width(5);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << id + partition->getStart(myrank);
		ossFileName << "-";
		ossFileName.width(10);
		ossFileName.setf(ios::fixed);
		ossFileName.fill('0');
		ossFileName << step;
		ossFileName << ".bin";

		int nx = 0;
		int ny = 0;
		int nz = 0;
		int nv = 0;
		int ne = 0;
		int nb = 0;

		ifstream ifs;
		ifs.open(ossFileName.str().c_str(), ios::in | ios::binary);
		ifs.read((char*)&nx, sizeof(int));
		ifs.read((char*)&ny, sizeof(int));
		ifs.read((char*)&nz, sizeof(int));
		ifs.read((char*)&nv, sizeof(int));
		ifs.read((char*)&ne, sizeof(int));
		ifs.read((char*)&nb, sizeof(int));

		int cx = nx + 2*nv;
		int cy = ny + 2*nv;
		int cz = nz + 2*nv;

		BlockBase* block = blockManager.getBlock(id);
		::Vec3i size = block->getSize();
		real* pData = GetBlockData(block);
		real *pDataS = new real [cx*cy*cz];
		if( nx == size.x && ny == size.y && nz == size.z && nv == vc && ne == 1 && nb == blockManager.getNumBlock() ) {
			ifs.read((char*)pData, sizeof(real)*cx*cy*cz);
		} else if( 2*nx == size.x && 2*ny == size.y && 2*nz == size.z && nv == vc && ne == 1 && nb == blockManager.getNumBlock() ) {
			ifs.read((char*)pDataS, sizeof(real)*cx*cy*cz);
			int sz[3] = {2*nx, 2*ny, 2*nz};
			sf3d_copy_x2_(
					(real*)pData,
					(real*)pDataS,
					(int*)sz, (int*)&vc);
		} else {
			Exit(0);
		}
		ifs.close();
		delete [] pDataS;
	}
	ImposeBoundaryCondition(blockManager);
}

