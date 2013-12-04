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

#ifndef BCM_POLYLIB_H
#define BCM_POLYLIB_H

#include "mpi.h"
#include "MPIPolylib.h"
#include "BoundingBox.h"

namespace PolylibNS {

class BCMPolylib : public MPIPolylib {

  std::string m_config_contents;

public:

  BCMPolylib(MPI::Comm& comm = MPI::COMM_WORLD);

  ~BCMPolylib();

  POLYLIB_STAT load(std::string config_filename);

  POLYLIB_STAT load_from_rank0();

  POLYLIB_STAT set_bounding_box(int rank, const Vec3f& min, const Vec3f& max);

  POLYLIB_STAT set_bounding_box(int rank, const BoundingBox& box) {
    Vec3f min(box.getMin().x, box.getMin().y, box.getMin().z);
    Vec3f max(box.getMax().x, box.getMax().y, box.getMax().z);
    return set_bounding_box(rank, min, max);
  }

  POLYLIB_STAT send_to_all();

protected:


private:

  // 以下は、呼び出し禁止にするMPIPolylibの公開メソッド

	static MPIPolylib* get_instance();

  POLYLIB_STAT init_parallel_info(MPI_Comm comm,
                      float bpos[3], unsigned int bbsize[3],
                      unsigned int gcsize[3], float dx[3]);

	POLYLIB_STAT load_rank0(std::string config_filename = "");

	POLYLIB_STAT load_parallel(std::string config_filename = "",
                       ID_FORMAT id_format = ID_BIN);

	POLYLIB_STAT move(PolylibMoveParams &params);

	POLYLIB_STAT migrate();
  
};

}

#endif // BCM_POLYLIB_H
