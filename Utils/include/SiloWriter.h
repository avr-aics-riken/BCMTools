/*
 * BCMTools
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file SiloWriter.h
/// @brief SiloフォーマットによるBCMデータ出力クラス
/// 
/// @note HDF5ライブラリをリンクしない場合は，
///       コンパイルオプションにUSE_PDB_DRIVERを追加すること．
///

#ifndef SILO_WRITER_H
#define SILO_WRITER_H

#include <string>
#include <sstream>
#include <cstring>
#include "mpi.h"
#include "silo.h"
#include "BCMTools.h"
#include "BlockManager.h"
#include "Scalar3D.h"
#include "Vector3D.h"

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif


/// SiloフォーマットによるBCMデータ出力クラス.
class SiloWriter {

  BlockManager& blockManager;   ///< ブロックマネージャ

  const MPI::Intracomm& comm;   ///< MPIコミュニケータ

  DBfile* file;         ///< 分散出力ファイル
  DBfile* masterFile;   ///< マスターファイル

  const std::string& fileName;  ///< ファイル名(マスターファイルの)

  bool withGhostZones;  ///< ゴーストセル使用フラグ

  int* blockIdTable;  ///< 各ランクでの開始ブロック番号(ランク0のみ)

  const std::string meshName;  ///< メッシュ名

public:

  /// コンストラクタ.
  /// 各ノードで分散出力ファイルをオープンしてメッシュデータを出力．
  /// ノード0は，マスターファイルもオープン．
  ///
  ///  @param[in] fileName (マスターファイルの)ファイル名
  ///  @param[in] meshName メッシュ名
  ///  @param[in] withGhostZones ゴーストセル使用フラグ
  ///
  ///  @note マスタファイル名xxx.extに対して，分散ファイル名はxxx_ランク番号.ext
  ///
  SiloWriter(const std::string& fileName, const std::string& meshName, bool withGhostZones = true)
   : blockManager(BlockManager::getInstance()),
     comm(blockManager.getCommunicator()),
     fileName(fileName), withGhostZones(withGhostZones), meshName(meshName) {

    int myrank = comm.Get_rank();

#ifndef USE_PDB_DRIVER
    DBSetCompression("METHOD=GZIP");  // 圧縮はHDF5ドライバのみ対応
#endif

    file = DBCreate(makeLocalFileName(fileName, myrank).c_str(),
#ifdef USE_PDB_DRIVER
                    DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
#else
                    DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);
#endif
    assert(file);

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);

      // 各ブロックに対応して，block_nというディレクトリを作成(nはグローバルID)
      std::ostringstream dirName;
      dirName << "block" << id + blockManager.getStartID();

      DBMkDir(file, dirName.str().c_str());
      DBSetDir(file, dirName.str().c_str());

      // メッシュデータ出力
      writeMesh(meshName, block);

      // ルートディレクトリに戻る
      DBSetDir(file, "/");
    }


    if (myrank == 0) {
      // マスターファイル
#ifdef USE_PDB_DRIVER
      masterFile = DBCreate(fileName.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
#else
      masterFile = DBCreate(fileName.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);
#endif
      assert(masterFile);
      blockIdTable = new int[comm.Get_size() + 1];
      blockIdTable[0] = 0;
    } else {
      masterFile = 0;
      blockIdTable = 0;
    }

    // blockIdTableには，各ランクでの先頭ブロックIDを格納．
    // ただし，末尾の余分な要素には全ブロック数を入れる．
    int endId = blockManager.getStartID() + blockManager.getNumBlock();
    comm.Gather(&endId, 1, MPI::INT, &blockIdTable[1], 1, MPI::INT, 0);

    if (myrank == 0) writeExternalArrays();

    if (myrank == 0) writeMultiMesh(meshName);
    
  }

  /// デストラクタ.
  ~SiloWriter() {
    DBClose(file);
    if (masterFile) DBClose(masterFile);
    delete[] blockIdTable;
  }


  /// 領域分割情報の出力.
  /// ブロック単位で，担当ランク番号(ノード番号)を出力．
  ///
  ///  @param[in] blockMeshName ブロック境界を定めるメッシュデータのメッシュ名
  ///  @param[in] name ランク番号データ名
  ///
  void writeDomain(const std::string& blockMeshName, const std::string& name) {
    int myRank = comm.Get_rank();
    int dims[3] = { 1, 1, 1 };

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);

      std::ostringstream dirName;
      dirName << "block" << id + blockManager.getStartID();
      DBSetDir(file, dirName.str().c_str());
      writeBlockBoundary(blockMeshName, block);
      DBPutQuadvar1(file, name.c_str(), blockMeshName.c_str(), &myRank, dims, 3,
                    NULL, 0, DB_INT, DB_ZONECENT, NULL);
      DBSetDir(file, "/");
    }

    if (myRank == 0) {
      writeMultiVar(name);
      writeMultiMesh(blockMeshName);
    }
  }


  /// スカラーデータの出力.
  ///
  ///  @param[in] dataClassID データクラスID
  ///  @param[in] name データ名
  ///
  ///  @note セル中心で定義されたScalar3Dのみ対応．整数型もfloatに変換して出力.
  ///
  template <typename T>
  void writeScalar(int dataClassID, const std::string& name) {
    const Vec3i& size = blockManager.getSize();
    float* data = new float[(size.x+2) * (size.y+2) * (size.z+2)];

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
      int gz_m[3], gz_p[3];
      setGhostZones(block, gz_m, gz_p);
      Vec3i vc_m(gz_m);
      Vec3i vc_p(gz_p);
      Vec3i size0 = size + vc_m + vc_p;
      int dims[3] = { size0.x, size0.y, size0.z };

      Scalar3D<T>* s = dynamic_cast<Scalar3D<T>*>(block->getDataClass(dataClassID));
      assert(s);
      T* sData = s->getData();
      Index3DS sIndex = s->getIndex();
      int ijk = 0;
      for (int k = -vc_m.z; k < size.z+vc_p.z; k++) {
        for (int j = -vc_m.y; j < size.y+vc_p.y; j++) {
          for (int i = -vc_m.x; i < size.x+vc_p.x; i++) {
            data[ijk++] = sData[sIndex(i,j,k)];
          }
        }
      }

      std::ostringstream dirName;
      dirName << "block" << id + blockManager.getStartID();
      DBSetDir(file, dirName.str().c_str());
      DBPutQuadvar1(file, name.c_str(), meshName.c_str(), data, dims, 3,
                    NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);
      DBSetDir(file, "/");
    }

    delete[] data;

    if (comm.Get_rank() == 0) writeMultiVar(name);
  }


  /// ベクトルデータの出力.
  ///
  ///  @param[in] dataClassID データクラスID
  ///  @param[in] name データ名
  ///
  ///  @note セル中心で定義されたVector3Dのみ対応．整数型もfloatに変換して出力.
  ///
  template <typename T>
  void writeVector(int dataClassID, const std::string& name) {
    const Vec3i& size = blockManager.getSize();
    float* vx = new float[(size.x+2) * (size.y+2) * (size.z+2)];
    float* vy = new float[(size.x+2) * (size.y+2) * (size.z+2)];
    float* vz = new float[(size.x+2) * (size.y+2) * (size.z+2)];
    float* data[3] = { vx, vy, vz };

    std::string name_x = name + std::string("_x");
    std::string name_y = name + std::string("_y");
    std::string name_z = name + std::string("_z");
    char* coordNames[3];
    coordNames[0] = strdup(name_x.c_str());
    coordNames[1] = strdup(name_y.c_str());
    coordNames[2] = strdup(name_z.c_str());

    for (int id = 0; id < blockManager.getNumBlock(); ++id) {
      BlockBase* block = blockManager.getBlock(id);
      int gz_m[3], gz_p[3];
      setGhostZones(block, gz_m, gz_p);
      Vec3i vc_m(gz_m);
      Vec3i vc_p(gz_p);
      Vec3i size0 = size + vc_m + vc_p;
      int dims[3] = { size0.x, size0.y, size0.z };

      Vector3D<T>* v = dynamic_cast<Vector3D<T>*>(block->getDataClass(dataClassID));
      assert(v);
      T* vData = v->getData();
      Index3DV vIndex = v->getIndex();
      int ijk = 0;
      for (int k = -vc_m.z; k < size.z+vc_p.z; k++) {
        for (int j = -vc_m.y; j < size.y+vc_p.y; j++) {
          for (int i = -vc_m.x; i < size.x+vc_p.x; i++) {
            vx[ijk] = vData[vIndex(i,j,k)+0];
            vy[ijk] = vData[vIndex(i,j,k)+1];
            vz[ijk] = vData[vIndex(i,j,k)+2];
            ijk++;
          }
        }
      }

      std::ostringstream dirName;
      dirName << "block" << id + blockManager.getStartID();
      DBSetDir(file, dirName.str().c_str());
      DBPutQuadvar(file, name.c_str(), meshName.c_str(), 3, coordNames, data, dims, 3,
                   NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);
      DBSetDir(file, "/");
    }

    delete[] vx;
    delete[] vy;
    delete[] vz;
    free(coordNames[0]);
    free(coordNames[1]);
    free(coordNames[2]);

    if (comm.Get_rank() == 0) writeMultiVar(name);
  }



private:

  ///  Nameschemeから参照されるexternal arrayを書き出す.
  ///
  ///  @note External arrayの名前ははVisItのソースコード
  ///        src/database/Silo/avtSiloMBObjectCache.Cによる.
  ///
  void writeExternalArrays() {
    int nRank = comm.Get_size();
    int nBlock = blockIdTable[nRank];
    int ndims = 1;
    int dims;

    int* procs = new int[nBlock];
    for (int iRank = 0; iRank < nRank; iRank++) {
      for (int id = blockIdTable[iRank]; id < blockIdTable[iRank+1]; id++) {
        procs[id] = iRank;
      }
    }
    dims = nBlock;
    DBWrite(masterFile, "procs", procs, &dims, ndims, DB_INT);
    delete[] procs;

    dims = 1;
    int dummy = 0;
    DBWrite(masterFile, "DomainFiles", &dummy, &dims, ndims, DB_INT);
  }

  /// ローカルファイル(分散ファイル)のファイル名を決定する.
  ///  @param[in] fileName マスターファイル名
  ///  @param[in] rank ランク番号
  ///  @return 分散ファイル名
  ///
  ///  @note マスタファイル名xxx.extに対して，分散ファイル名はxxx_ランク番号.ext
  ///
  const std::string makeLocalFileName(const std::string& fileName, int rank) {
    std::string localFileName(fileName);
    std::ostringstream rankStr;
    rankStr << "_" << rank;
    std::string::size_type pos = localFileName.rfind(".");
    localFileName.insert(pos, rankStr.str());
    return localFileName;
  }
  

  /// ブロック単位でメッシュデータを出力．
  ///
  ///  @param[in] name メッシュ名
  ///  @param[in] block ブロック
  ///
  void writeMesh(const std::string& name, const BlockBase* block) {
    const Vec3i& size = block->getSize();
    const Vec3d& orig = block->getOrigin();
    const Vec3d& pitch = block->getCellSize();

//  DBoptlist* optList = DBMakeOptlist(2);
    DBoptlist* optList = DBMakeOptlist(5);
    int gz_m[3], gz_p[3];
    setGhostZones(block, gz_m, gz_p);
    DBAddOption(optList, DBOPT_LO_OFFSET, gz_m);
    DBAddOption(optList, DBOPT_HI_OFFSET, gz_p);

    DBAddOption(optList, DBOPT_XLABEL, (void*)" ");
    DBAddOption(optList, DBOPT_YLABEL, (void*)" ");
    DBAddOption(optList, DBOPT_ZLABEL, (void*)" ");

    Vec3i vc_m(gz_m);
    Vec3i vc_p(gz_p);

    Vec3i size0 = size + vc_m + vc_p;
    int dims[3] = { size0.x+1, size0.y+1, size0.z+1 };
 // double x[size0.x+1];
 // double y[size0.y+1];
 // double z[size0.z+1];
    double* x = new double[size0.x+1];
    double* y = new double[size0.y+1];
    double* z = new double[size0.z+1];
    double *coords[3] = {x, y, z };
    for (int i = 0; i < size0.x+1; i++) { x[i] = orig.x + (i - vc_m.x) * pitch.x; }
    for (int j = 0; j < size0.y+1; j++) { y[j] = orig.y + (j - vc_m.y) * pitch.y; }
    for (int k = 0; k < size0.z+1; k++) { z[k] = orig.z + (k - vc_m.z) * pitch.z; }

    DBPutQuadmesh(file, name.c_str(), NULL, coords, dims, 3,
                  DB_DOUBLE, DB_COLLINEAR, optList);

    DBFreeOptlist(optList);
    delete[] x;
    delete[] y;
    delete[] z;
  }


  /// ブロック単位でブロック境界メッシュデータを出力．
  ///
  ///  @param[in] name メッシュ名
  ///  @param[in] block ブロック
  ///
  void writeBlockBoundary(const std::string& name, const BlockBase* block) {
    const Vec3d& orig = block->getOrigin();
    const Vec3d& blockSize = block->getBlockSize();
    double x[2] = { orig.x, orig.x + blockSize.x };
    double y[2] = { orig.y, orig.y + blockSize.y };
    double z[2] = { orig.z, orig.z + blockSize.z };
    int dims[3] = { 2, 2, 2 };
    double *coords[3] = {x, y, z };
    DBoptlist* optList = DBMakeOptlist(3);
    DBAddOption(optList, DBOPT_XLABEL, (void*)" ");
    DBAddOption(optList, DBOPT_YLABEL, (void*)" ");
    DBAddOption(optList, DBOPT_ZLABEL, (void*)" ");
    DBPutQuadmesh(file, name.c_str(), NULL, coords, dims, 3,
                  DB_DOUBLE, DB_COLLINEAR, optList);
    DBFreeOptlist(optList);
  }


  /// メッシュのマスターデータを出力.
  ///
  ///  @param[in] name メッシュ名
  ///
  void writeMultiMesh(const std::string& name) {
    int nRank = comm.Get_size();
    int nBlock = blockIdTable[nRank];

#if 1
    DBoptlist* optList = DBMakeOptlist(3);

    int meshType = DB_QUAD_RECT;
    DBAddOption(optList, DBOPT_MB_BLOCK_TYPE, &meshType);

    const std::string fileNsStr = makeFileNameschemeString(fileName);
    DBAddOption(optList, DBOPT_MB_FILE_NS, (void*)fileNsStr.c_str());

    const std::string blockNsStr = makeBlockNameschemeString(name);
    DBAddOption(optList, DBOPT_MB_BLOCK_NS, (void*)blockNsStr.c_str());

    DBPutMultimesh(masterFile, name.c_str(), nBlock, NULL, NULL, optList);

    DBFreeOptlist(optList);
#else

    char** meshNames = new char*[nBlock];
    int* meshTypes = new int[nBlock];

    for (int iRank = 0; iRank < nRank; iRank++) {
      const std::string localFileName = makeLocalFileName(fileName, iRank);
      for (int id = blockIdTable[iRank]; id < blockIdTable[iRank+1]; id++) {
        std::ostringstream localName;
        localName << localFileName << ":" << "block" << id << "/" << name;
      //std::cout << "localName: " << localName.str() << std::endl;

        meshNames[id] = strdup(localName.str().c_str());
        meshTypes[id] = DB_QUAD_RECT;
      }
    }
    DBPutMultimesh(masterFile, name.c_str(), nBlock, meshNames, meshTypes, NULL);
    for (int id = 0; id < nBlock; id++) delete[] meshNames[id];
    delete[] meshNames;
    delete[] meshTypes;

#endif
  }

  /// セルデータのマスターデータを出力.
  ///
  ///  @param[in] name データ名
  ///
  void writeMultiVar(const std::string&  name) {
    int nRank = comm.Get_size();
    int nBlock = blockIdTable[nRank];

#if 1
    DBoptlist* optList = DBMakeOptlist(3);

    int meshType = DB_QUADVAR;
    DBAddOption(optList, DBOPT_MB_BLOCK_TYPE, &meshType);

    const std::string fileNsStr = makeFileNameschemeString(fileName);
    DBAddOption(optList, DBOPT_MB_FILE_NS, (void*)fileNsStr.c_str());

    const std::string blockNsStr = makeBlockNameschemeString(name);
    DBAddOption(optList, DBOPT_MB_BLOCK_NS, (void*)blockNsStr.c_str());

    DBPutMultivar(masterFile, name.c_str(), nBlock, NULL, NULL, optList);

    DBFreeOptlist(optList);

#else
    char** varNames = new char*[nBlock];
    int* varTypes = new int[nBlock];

    for (int iRank = 0; iRank < nRank; iRank++) {
      const std::string localFileName = makeLocalFileName(fileName, iRank);
      for (int id = blockIdTable[iRank]; id < blockIdTable[iRank+1]; id++) {
        std::ostringstream localName;
        localName << localFileName << ":" << "block" << id << "/" << name;
        varNames[id] = strdup(localName.str().c_str());
        varTypes[id] = DB_QUADVAR;
      }
    }
    DBPutMultivar(masterFile, name.c_str(), nBlock, varNames, varTypes, NULL);
    for (int id = 0; id < nBlock; id++) free(varNames[id]);
    delete[] varNames;
    delete[] varTypes;
#endif
  }

  /// ブロック毎に，ゴーストゾーンの幅を設定.
  ///
  ///  @param[in] block ブロック
  ///  @param[out] gz_m[] マイナス側ゴーストゾーン幅
  ///  @param[out] gz_p[] プラス側ゴーストゾーン幅
  ///
  void setGhostZones(const BlockBase* block, int gz_m[], int gz_p[]) {
    if (withGhostZones) {
      const NeighborInfo* neighborInfo = block->getNeighborInfo();
      // 隣接ブロックが存在したら、ゴーストゾーンの幅を1にする.
      gz_m[0] = neighborInfo[X_M].exists() ? 1 : 0;
      gz_p[0] = neighborInfo[X_P].exists() ? 1 : 0;
      gz_m[1] = neighborInfo[Y_M].exists() ? 1 : 0;
      gz_p[1] = neighborInfo[Y_P].exists() ? 1 : 0;
      gz_m[2] = neighborInfo[Z_M].exists() ? 1 : 0;
      gz_p[2] = neighborInfo[Z_P].exists() ? 1 : 0;
    } else {
      gz_m[0] = gz_m[1] = gz_m[2] = 0;
      gz_p[0] = gz_p[1] = gz_p[2] = 0;
    }
  }


  ///  in:  fileName = "<file>" + "." + "<ext>"
  ///  out: "@<file>_%-d.<ext>@#P[n]"
  std::string makeFileNameschemeString(const std::string& fileName) {
//  const char* fileNs = "@test_%-d.silo@#P[n]";
    std::string fileNameStr(fileName);
    std::string::size_type pos = fileNameStr.rfind(".");
    fileNameStr.insert(pos, "_%-d");
    return std::string("@") + fileNameStr + std::string("@#P[n]");
  }

  /// "@bclock%-d/<name>@n"
  std::string makeBlockNameschemeString(const std::string& name) {
//  const char* blockNs = "@block%-d/mesh@n";
    return std::string("@block%-d/") + name + std::string("@n");
  }

};


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // SILO_WRITER_H
