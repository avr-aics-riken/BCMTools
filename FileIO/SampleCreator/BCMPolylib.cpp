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

#include "BCMPolylib.h"

using namespace PolylibNS;

BCMPolylib::BCMPolylib(MPI::Comm& comm)
{
  m_mycomm = comm;  // ちゃんとMPI_Commにキャストされる(see openmpi/ompi/mpi/cxx/comm.h)
  MPI_Comm_rank(m_mycomm, &m_myrank);
  MPI_Comm_size(m_mycomm, &m_numproc);
}


BCMPolylib::~BCMPolylib()
{
  // MPIPolilibのprotectedだったデストラクタが呼ばれる
}


POLYLIB_STAT BCMPolylib::load(std::string config_filename)
{
  // rank0のみPolylibクラスのloadメソッドを呼ぶことを許可する
  if (m_myrank != 0) return PLSTAT_NG;

  POLYLIB_STAT ret;

  // 設定ファイル読み込み。
  if ((ret = Polylib::load_config_file(&m_config_contents, config_filename)) != PLSTAT_OK ) {
    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():Polylib::load_config() faild. returns:"
              << PolylibStat2::String(ret) << std::endl;
    return ret;
  }

  // グループ階層構造構築。
  if ((ret = Polylib::make_group_tree(m_config_contents)) != PLSTAT_OK ) {
    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():Polylib::make_group_tree() faild. returns:"
              << PolylibStat2::String(ret) << std::endl;
    return ret;
  }

  // ポリゴン情報を構築 (三角形IDファイルは不要なので、第二引数はダミー)
  if ((ret = load_polygons(false, ID_BIN)) != PLSTAT_OK) {
    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():load_polygons() faild. returns:"
              << PolylibStat2::String(ret) << std::endl;
    return ret;
  }

  return PLSTAT_OK;
}


POLYLIB_STAT BCMPolylib::load_from_rank0()
{
  if (m_myrank == 0) return PLSTAT_NG;

  POLYLIB_STAT ret;

  // 設定ファイルの内容をrank0から受信する
  if ((ret = broadcast_config_from_rank0()) != PLSTAT_OK ) {
    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():broadcast_config_from_rank0() faild. returns:"
              << PolylibStat2::String(ret) << std::endl;
    return ret;
  }

  // ポリゴン情報をrank0から受信する。
  if ((ret = receive_polygons_from_rank0()) != PLSTAT_OK ) {
    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():receive_polygons_from_rank0() faild. returns:"
              << PolylibStat2::String(ret) << std::endl;
    return ret;
  }

  return PLSTAT_OK;
}


POLYLIB_STAT BCMPolylib::set_bounding_box(int rank,
                                          const Vec3f& min, const Vec3f& max)
{
  if (m_myrank != 0) return PLSTAT_NG;

  // ParallelInfoの必要最低限のメンバのみを流用
  if (rank == 0) {
    m_myproc.m_comm = m_mycomm;
    m_myproc.m_rank = rank;
    m_myproc.m_area.m_gcell_min = min;
    m_myproc.m_area.m_gcell_max = max;
    m_myproc.m_area.m_gcell_bbox.add(min);
    m_myproc.m_area.m_gcell_bbox.add(max);
  }
  else {
    ParallelInfo* proc = new ParallelInfo;
    proc->m_comm = m_mycomm;
    proc->m_rank = rank;
    proc->m_area.m_gcell_min = min;
    proc->m_area.m_gcell_max = max;
    proc->m_area.m_gcell_bbox.add(min);
    proc->m_area.m_gcell_bbox.add(max);
    m_other_procs.push_back(proc);
  }

  return PLSTAT_OK;
}


POLYLIB_STAT BCMPolylib::send_to_all()
{
  if (m_myrank != 0) return PLSTAT_NG;

  POLYLIB_STAT ret;

  // 設定ファイルの内容を他PEへブロードキャストする
  if ((ret = broadcast_config(m_config_contents )) != PLSTAT_OK ) {
    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():broadcast_config() faild. returns:"
              << PolylibStat2::String(ret) << std::endl;
    return ret;
  }

  // ポリゴン情報を他PEへ配信する
  if ((ret = send_polygons_to_all()) != PLSTAT_OK ) {
    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():send_polygons_to_all() faild. returns:"
              << PolylibStat2::String(ret) << std::endl;
      return ret;
   }

  // 他PE領域ポリゴン情報を削除して自領域分のみでデータ構造再構築
  if ((ret = erase_outbounded_polygons()) != PLSTAT_OK ) {
    PL_ERROSH << "[ERROR]MPIPolylib::load_rank0():erase_outbounded_polygons() faild. returns:"
              << PolylibStat2::String(ret) << std::endl;
    return ret;
  }

  return PLSTAT_OK;
}
