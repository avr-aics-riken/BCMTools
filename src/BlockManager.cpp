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

#include "BlockManager.h"
#include "limits.h"   // for INT_MAX

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif

#include <fstream>

/// ブロックを登録.
void BlockManager::registerBlock(BlockBase* block)
{
  blockList.push_back(block);
  numBlock++;
}


/// ブロック登録の完了.
void BlockManager::endRegisterBlock()
{
  if (numBlock == 0) {
    std::cout << "***error: number of blocks must be > 0" << std::endl;
    Exit(EX_FAILURE);
  }

  for (int i = 0; i < numBlock; i++) {
    BlockBase* block = blockList[i];
    if (i == 0) {
      size = block->getSize();
    } else {
      // 各ブロックでのセル数が等しいことを確認
      if (size != block->getSize()) {
        std::cout << "***error: block size inconsistent" << std::endl;
        Exit(EX_FAILURE);
      }
    }
  }

  int numBlockScanned;
  comm.Scan(&numBlock, &numBlockScanned, 1, MPI::INT, MPI::SUM);
  startID = numBlockScanned - numBlock;
}


/// 仮想セル同期の準備.
void BlockManager::prepareForVCUpdate(int dataClassID, int tag, VCUpdateMethod::Type method)
{
  switch (method) {
    case VCUpdateMethod::AtOnce:
      if (!faceListPrepared) setFaceLists();
      setDataClassNeighbor(dataClassID, localFaceList);
      assert(commBufferTableMap.find(dataClassID) == commBufferTableMap.end());
      commBufferTableMap[dataClassID] = new CommBufferTable(tag);
      setSendBufferPointers(commBufferTableMap[dataClassID]->sendBuffer,
                            dataClassID, sendFaceList);
      setRecvBufferPointers(commBufferTableMap[dataClassID]->recvBuffer,
                            dataClassID, recvFaceList);
      commBufferTableMap[dataClassID]->sendBuffer->allocateBuffer();
      commBufferTableMap[dataClassID]->recvBuffer->allocateBuffer();
      break;
    case VCUpdateMethod::SeparateXYZ:
      if (!separateFaceListPrepared) setSeparateFaceLists();
      for (int xyz = 0; xyz < 3; xyz++) {
        setDataClassNeighbor(dataClassID, localSeparateFaceList[xyz]);
        assert(commBufferSeparateTableMap[xyz].find(dataClassID)
               == commBufferSeparateTableMap[xyz].end());
        commBufferSeparateTableMap[xyz][dataClassID] = new CommBufferTable(tag);
        setSendBufferPointers(commBufferSeparateTableMap[xyz][dataClassID]->sendBuffer,
                              dataClassID, sendSeparateFaceList[xyz]);
        setRecvBufferPointers(commBufferSeparateTableMap[xyz][dataClassID]->recvBuffer,
                              dataClassID, recvSeparateFaceList[xyz]);
        commBufferSeparateTableMap[xyz][dataClassID]->sendBuffer->allocateBuffer();
        commBufferSeparateTableMap[xyz][dataClassID]->recvBuffer->allocateBuffer();
      }
      break;
    case VCUpdateMethod::SeparateXYZ_SeparateLevelDiff:
      if (!separateLevelDiffFaceListPrepared) setSeparateLevelDiffFaceLists();
      for (int xyz = 0; xyz < 3; xyz++) {
        for (int ld = 0; ld < 3; ld++) {
          setDataClassNeighbor(dataClassID, localSeparateLevelDiffFaceList[xyz][ld]);
          assert(commBufferSeparateLevelDiffTableMap[xyz][ld].find(dataClassID)
                 == commBufferSeparateLevelDiffTableMap[xyz][ld].end());
          commBufferSeparateLevelDiffTableMap[xyz][ld][dataClassID] = new CommBufferTable(tag);
          setSendBufferPointers(commBufferSeparateLevelDiffTableMap[xyz][ld][dataClassID]->sendBuffer,
                                dataClassID, sendSeparateLevelDiffFaceList[xyz][ld]);
          setRecvBufferPointers(commBufferSeparateLevelDiffTableMap[xyz][ld][dataClassID]->recvBuffer,
                                dataClassID, recvSeparateLevelDiffFaceList[xyz][ld]);
          commBufferSeparateLevelDiffTableMap[xyz][ld][dataClassID]->sendBuffer->allocateBuffer();
          commBufferSeparateLevelDiffTableMap[xyz][ld][dataClassID]->recvBuffer->allocateBuffer();
        }
      }
      break;
    default:
      assert(0);
      break;
  }
}


/// ブロック配置情報を出力.
void BlockManager::printBlockLayoutInfo()
{
  int myrank = comm.Get_rank();
  int nprocs = comm.Get_size();

  if (myrank == 0) {
    std::cout << "Block Layout Information" << std::endl;
  }

  // check level of blocks

  int levelMin = INT_MAX;
  int levelMax = 0;
  for (int localID = 0; localID < numBlock; localID++) {
    int level = blockList[localID]->getLevel();
    levelMin = std::min(level, levelMin);
    levelMax = std::max(level, levelMax);
  }

//std::cout << myrank << ": Lmin=" << levelMin << " Lmax=" << levelMax << std::endl;

  comm.Allreduce(MPI_IN_PLACE, &levelMin, 1, MPI::INT, MPI::MIN);
  comm.Allreduce(MPI_IN_PLACE, &levelMax, 1, MPI::INT, MPI::MAX);

  if (myrank == 0) {
    std::cout << "  Min level: " << levelMin << std::endl;
    std::cout << "  Max level: " << levelMax << std::endl;
  }


  // count blocks

  int nLevel = levelMax - levelMin + 1;
  int* nBlock = new int[nLevel];
  for (int i = 0; i < nLevel; i++) nBlock[i] = 0;

  for (int localID = 0; localID < numBlock; localID++) {
    BlockBase* block = blockList[localID];
    int level = block->getLevel();
    nBlock[level-levelMin]++;
  }

  int numBlockSum;
  comm.Reduce(&numBlock, &numBlockSum, 1, MPI::INT, MPI::SUM, 0);
  if (myrank == 0) {
    comm.Reduce(MPI_IN_PLACE, nBlock, nLevel, MPI::INT, MPI::SUM, 0);
  } else {
    comm.Reduce(nBlock, nBlock, nLevel, MPI::INT, MPI::SUM, 0);
  }

  if (myrank == 0) {
    std::cout << "  Number of blocks" << std::endl;
    for (int level = levelMin; level <= levelMax; level++) {
      std::cout << "    L=" << level << ": " << nBlock[level-levelMin] << std::endl;
    }
    std::cout << "    total: " << numBlockSum << std::endl;
  }

  delete[] nBlock;

  int numBlockMin, numBlockMax;
  comm.Reduce(&numBlock, &numBlockMin, 1, MPI::INT, MPI::MIN, 0);
  comm.Reduce(&numBlock, &numBlockMax, 1, MPI::INT, MPI::MAX, 0);
  int numBlock2Sum = numBlock * numBlock;
  if (myrank == 0) {
    comm.Reduce(MPI_IN_PLACE, &numBlock2Sum, 1, MPI::INT, MPI::SUM, 0);
  } else {
    comm.Reduce(&numBlock2Sum, &numBlock2Sum, 1, MPI::INT, MPI::SUM, 0);
  }

  if (myrank == 0) {
    std::cout << "  Number of total blocks / node" << std::endl;
    REAL_TYPE ave = (REAL_TYPE)numBlockSum / nprocs;
    std::cout << "    ave: " << ave << std::endl;
    std::cout << "    min: " << numBlockMin << std::endl;
    std::cout << "    max: " << numBlockMax << std::endl;
    std::cout << "     sd: " << sqrt((REAL_TYPE)numBlock2Sum/nprocs - ave*ave)<< std::endl;
  }


  // count faces

  int nFaceInter[3] = { 0, 0, 0 };
  int nFaceIntra[3] = { 0, 0, 0 };

  for (int localID = 0; localID < numBlock; localID++) {
    BlockBase* block = blockList[localID];
    const NeighborInfo* neighborInfo = block->getNeighborInfo();
    for (int i = 0; i < NUM_FACE; i++) {
      Face face = Face(i);
      if (neighborInfo[face].exists()) {
        int levelDiff = neighborInfo[face].getLevelDifference();
        if (levelDiff == 0) {
          if (neighborInfo[face].getRank() == myrank) {
            nFaceIntra[LD_0]++;
          } else {
            nFaceInter[LD_0]++;
          }
        }
        else if (levelDiff == 1) {
          for (int j = 0; j < NUM_SUBFACE; j++) {
            Subface subface = Subface(j);
            if (neighborInfo[face].getRank(subface) == myrank) {
              nFaceIntra[LD_P1]++;
            } else {
              nFaceInter[LD_P1]++;
            }
          }
        }
        else if (levelDiff == -1) {
          if (neighborInfo[face].getRank() == myrank) {
            nFaceIntra[LD_M1]++;
          } else {
            nFaceInter[LD_M1]++;
          }
        }
        else {
          Exit(EX_FAILURE);
        }
      }
    }
  }

  int nFaceInterSum[3] = { 0, 0, 0 };
  int nFaceIntraSum[3] = { 0, 0, 0 };

  comm.Reduce(nFaceInter, nFaceInterSum, 3, MPI::INT, MPI::SUM, 0);
  comm.Reduce(nFaceIntra, nFaceIntraSum, 3, MPI::INT, MPI::SUM, 0);

  if (myrank == 0) {
    std::cout << "  Number of faces" << std::endl;
    std::cout << "    intra-node " << std::endl;
    std::cout << "        dLevel=-1: " << nFaceIntraSum[LD_M1] << std::endl;
    std::cout << "        dLevel= 0: " << nFaceIntraSum[LD_0] << std::endl;
    std::cout << "        dLevel=+1: " << nFaceIntraSum[LD_P1] << std::endl;
    std::cout << "            total: "
           << nFaceIntraSum[LD_M1] + nFaceIntraSum[LD_0] + nFaceIntraSum[LD_P1]
           << std::endl;
    std::cout << "    inter-node " << std::endl;
    std::cout << "        dLevel=-1: " << nFaceInterSum[LD_M1] << std::endl;
    std::cout << "        dLevel= 0: " << nFaceInterSum[LD_0] << std::endl;
    std::cout << "        dLevel=+1: " << nFaceInterSum[LD_P1] << std::endl;
    std::cout << "            total: "
           << nFaceInterSum[LD_M1] + nFaceInterSum[LD_0] + nFaceInterSum[LD_P1]
           << std::endl;
  }

  int nFaceInterTotal = nFaceInter[LD_M1] + nFaceInter[LD_0] + nFaceInter[LD_P1];

  int nFaceInterTotalMin, nFaceInterTotalMax;
  comm.Reduce(&nFaceInterTotal, &nFaceInterTotalMin, 1, MPI::INT, MPI::MIN, 0);
  comm.Reduce(&nFaceInterTotal, &nFaceInterTotalMax, 1, MPI::INT, MPI::MAX, 0);

  int nFaceInterTotal2Sum = nFaceInterTotal * nFaceInterTotal;
  if (myrank == 0) {
    comm.Reduce(MPI_IN_PLACE, &nFaceInterTotal2Sum, 1, MPI::INT, MPI::SUM, 0);
  } else {
    comm.Reduce(&nFaceInterTotal2Sum, &nFaceInterTotal2Sum, 1, MPI::INT, MPI::SUM, 0);
  }

  if (myrank == 0) {
    std::cout << "  Number of total inter-node faces / node" << std::endl;
    REAL_TYPE ave
      = (REAL_TYPE)(nFaceInterSum[LD_M1] + nFaceInterSum[LD_0] + nFaceInterSum[LD_P1]) / nprocs;
    std::cout << "    ave: " << ave << std::endl;
    std::cout << "    min: " << nFaceInterTotalMin << std::endl;
    std::cout << "    max: " << nFaceInterTotalMax << std::endl;
    std::cout << "     sd: " << sqrt((REAL_TYPE)nFaceInterTotal2Sum/nprocs - ave*ave) << std::endl;
  }
}

/// ブロック配置情報を出力(ファイル).
void BlockManager::printBlockLayoutInfo(const char* filename)
{
  int myrank = comm.Get_rank();
  int nprocs = comm.Get_size();

	std::ofstream ofs;
  if (myrank == 0) {
		ofs.open(filename, std::ios::out);
	}

  if (myrank == 0) {
    ofs << "Block Layout Information" << std::endl;
  }

  // check level of blocks

  int levelMin = INT_MAX;
  int levelMax = 0;
  for (int localID = 0; localID < numBlock; localID++) {
    int level = blockList[localID]->getLevel();
    levelMin = std::min(level, levelMin);
    levelMax = std::max(level, levelMax);
  }

//ofs << myrank << ": Lmin=" << levelMin << " Lmax=" << levelMax << std::endl;

  comm.Allreduce(MPI_IN_PLACE, &levelMin, 1, MPI::INT, MPI::MIN);
  comm.Allreduce(MPI_IN_PLACE, &levelMax, 1, MPI::INT, MPI::MAX);

  if (myrank == 0) {
    ofs << "  Min level: " << levelMin << std::endl;
    ofs << "  Max level: " << levelMax << std::endl;
  }


  // count blocks

  int nLevel = levelMax - levelMin + 1;
  int* nBlock = new int[nLevel];
  for (int i = 0; i < nLevel; i++) nBlock[i] = 0;

  for (int localID = 0; localID < numBlock; localID++) {
    BlockBase* block = blockList[localID];
    int level = block->getLevel();
    nBlock[level-levelMin]++;
  }

  int numBlockSum;
  comm.Reduce(&numBlock, &numBlockSum, 1, MPI::INT, MPI::SUM, 0);
  if (myrank == 0) {
    comm.Reduce(MPI_IN_PLACE, nBlock, nLevel, MPI::INT, MPI::SUM, 0);
  } else {
    comm.Reduce(nBlock, nBlock, nLevel, MPI::INT, MPI::SUM, 0);
  }

  if (myrank == 0) {
    ofs << "  Number of blocks" << std::endl;
    for (int level = levelMin; level <= levelMax; level++) {
      ofs << "    L=" << level << ": " << nBlock[level-levelMin] << std::endl;
    }
    ofs << "    total: " << numBlockSum << std::endl;
  }

  delete[] nBlock;

  int numBlockMin, numBlockMax;
  comm.Reduce(&numBlock, &numBlockMin, 1, MPI::INT, MPI::MIN, 0);
  comm.Reduce(&numBlock, &numBlockMax, 1, MPI::INT, MPI::MAX, 0);
  int numBlock2Sum = numBlock * numBlock;
  if (myrank == 0) {
    comm.Reduce(MPI_IN_PLACE, &numBlock2Sum, 1, MPI::INT, MPI::SUM, 0);
  } else {
    comm.Reduce(&numBlock2Sum, &numBlock2Sum, 1, MPI::INT, MPI::SUM, 0);
  }

  if (myrank == 0) {
    ofs << "  Number of total blocks / node" << std::endl;
    REAL_TYPE ave = (REAL_TYPE)numBlockSum / nprocs;
    ofs << "    ave: " << ave << std::endl;
    ofs << "    min: " << numBlockMin << std::endl;
    ofs << "    max: " << numBlockMax << std::endl;
    ofs << "     sd: " << sqrt((REAL_TYPE)numBlock2Sum/nprocs - ave*ave)<< std::endl;
  }


  // count faces

  int nFaceInter[3] = { 0, 0, 0 };
  int nFaceIntra[3] = { 0, 0, 0 };

  for (int localID = 0; localID < numBlock; localID++) {
    BlockBase* block = blockList[localID];
    const NeighborInfo* neighborInfo = block->getNeighborInfo();
    for (int i = 0; i < NUM_FACE; i++) {
      Face face = Face(i);
      if (neighborInfo[face].exists()) {
        int levelDiff = neighborInfo[face].getLevelDifference();
        if (levelDiff == 0) {
          if (neighborInfo[face].getRank() == myrank) {
            nFaceIntra[LD_0]++;
          } else {
            nFaceInter[LD_0]++;
          }
        }
        else if (levelDiff == 1) {
          for (int j = 0; j < NUM_SUBFACE; j++) {
            Subface subface = Subface(j);
            if (neighborInfo[face].getRank(subface) == myrank) {
              nFaceIntra[LD_P1]++;
            } else {
              nFaceInter[LD_P1]++;
            }
          }
        }
        else if (levelDiff == -1) {
          if (neighborInfo[face].getRank() == myrank) {
            nFaceIntra[LD_M1]++;
          } else {
            nFaceInter[LD_M1]++;
          }
        }
        else {
          Exit(EX_FAILURE);
        }
      }
    }
  }

  int nFaceInterSum[3] = { 0, 0, 0 };
  int nFaceIntraSum[3] = { 0, 0, 0 };

  comm.Reduce(nFaceInter, nFaceInterSum, 3, MPI::INT, MPI::SUM, 0);
  comm.Reduce(nFaceIntra, nFaceIntraSum, 3, MPI::INT, MPI::SUM, 0);

  if (myrank == 0) {
    ofs << "  Number of faces" << std::endl;
    ofs << "    intra-node " << std::endl;
    ofs << "        dLevel=-1: " << nFaceIntraSum[LD_M1] << std::endl;
    ofs << "        dLevel= 0: " << nFaceIntraSum[LD_0] << std::endl;
    ofs << "        dLevel=+1: " << nFaceIntraSum[LD_P1] << std::endl;
    ofs << "            total: "
           << nFaceIntraSum[LD_M1] + nFaceIntraSum[LD_0] + nFaceIntraSum[LD_P1]
           << std::endl;
    ofs << "    inter-node " << std::endl;
    ofs << "        dLevel=-1: " << nFaceInterSum[LD_M1] << std::endl;
    ofs << "        dLevel= 0: " << nFaceInterSum[LD_0] << std::endl;
    ofs << "        dLevel=+1: " << nFaceInterSum[LD_P1] << std::endl;
    ofs << "            total: "
           << nFaceInterSum[LD_M1] + nFaceInterSum[LD_0] + nFaceInterSum[LD_P1]
           << std::endl;
  }

  int nFaceInterTotal = nFaceInter[LD_M1] + nFaceInter[LD_0] + nFaceInter[LD_P1];

  int nFaceInterTotalMin, nFaceInterTotalMax;
  comm.Reduce(&nFaceInterTotal, &nFaceInterTotalMin, 1, MPI::INT, MPI::MIN, 0);
  comm.Reduce(&nFaceInterTotal, &nFaceInterTotalMax, 1, MPI::INT, MPI::MAX, 0);

  int nFaceInterTotal2Sum = nFaceInterTotal * nFaceInterTotal;
  if (myrank == 0) {
    comm.Reduce(MPI_IN_PLACE, &nFaceInterTotal2Sum, 1, MPI::INT, MPI::SUM, 0);
  } else {
    comm.Reduce(&nFaceInterTotal2Sum, &nFaceInterTotal2Sum, 1, MPI::INT, MPI::SUM, 0);
  }

  if (myrank == 0) {
    ofs << "  Number of total inter-node faces / node" << std::endl;
    REAL_TYPE ave
      = (REAL_TYPE)(nFaceInterSum[LD_M1] + nFaceInterSum[LD_0] + nFaceInterSum[LD_P1]) / nprocs;
    ofs << "    ave: " << ave << std::endl;
    ofs << "    min: " << nFaceInterTotalMin << std::endl;
    ofs << "    max: " << nFaceInterTotalMax << std::endl;
    ofs << "     sd: " << sqrt((REAL_TYPE)nFaceInterTotal2Sum/nprocs - ave*ave) << std::endl;
  }

	if( myrank == 0 ) {
		ofs.close();
	}
}

/// フェイスリストの設定.
void BlockManager::setFaceLists()
{
  if (faceListPrepared) return;

  int myrank = comm.Get_rank();

  for (int localID = 0; localID < numBlock; localID++) {
    BlockBase* block = blockList[localID];
    const NeighborInfo* neighborInfo = block->getNeighborInfo();

    for (int i = 0; i < NUM_FACE; i++) {
      Face face = Face(i);

      if (neighborInfo[face].getLevelDifference() == 0 ||
          neighborInfo[face].getLevelDifference() == -1) {
        int rank = neighborInfo[face].getRank();
        if (rank == myrank) {
          localFaceList.push_back(FaceID(localID, face));
        }
        else if (rank != MPI::PROC_NULL) {
          sendFaceList[rank].push_back(FaceID(localID, face));
          recvFaceList[rank].push_back(FaceID(localID, face));
        }
      }

      if (neighborInfo[face].getLevelDifference() == 1) {
        for (int j = 0; j < NUM_SUBFACE; j++) {
          Subface subface = Subface(j);
          int rank = neighborInfo[face].getRank(subface);
          if (rank == myrank) {
            localFaceList.push_back(FaceID(localID, face, subface));
          }
          else if (rank != MPI::PROC_NULL) {
            sendFaceList[rank].push_back(FaceID(localID, face, subface));
            recvFaceList[rank].push_back(FaceID(localID, face, subface));
          }
        }
      }
    }
  }

  // 送信側のFaceIDオーダーに合わせるため，RecvFaceListをソートする
  FaceListMap::iterator it = recvFaceList.begin();
  for (; it != recvFaceList.end(); ++it) {
    std::sort(it->second.begin(), it->second.end(), RecvFaceComp(blockList));
  }

#if 0
  std::cout << std::endl << "SendFaceList" << std::endl;
  printFaceListMap(sendFaceList);
  std::cout << std::endl << "RecvFaceList" << std::endl;
  printFaceListMap(recvFaceList);
  std::cout << std::endl << "SendFaceListPeer" << std::endl;
  printFaceListMapPeer(sendFaceList);
  std::cout << std::endl << "RecvFaceListPeer" << std::endl;
  printFaceListMapPeer(recvFaceList);
#endif

  faceListPrepared = true;
}


/// フェイスリスト(3方向別)の設定.
void BlockManager::setSeparateFaceLists() {
  if (separateFaceListPrepared) return;

  int myrank = comm.Get_rank();

  for (int localID = 0; localID < numBlock; localID++) {
    BlockBase* block = blockList[localID];
    const NeighborInfo* neighborInfo = block->getNeighborInfo();

    for (int i = 0; i < NUM_FACE; i++) {
      Face face = Face(i);
      int xyz = faceToXYZ(face);

      if (neighborInfo[face].getLevelDifference() == 0 ||
          neighborInfo[face].getLevelDifference() == -1) {
        int rank = neighborInfo[face].getRank();
        if (rank == myrank) {
          localSeparateFaceList[xyz].push_back(FaceID(localID, face));
        }
        else if (rank != MPI::PROC_NULL) {
          sendSeparateFaceList[xyz][rank].push_back(FaceID(localID, face));
          recvSeparateFaceList[xyz][rank].push_back(FaceID(localID, face));
        }
      }

      if (neighborInfo[face].getLevelDifference() == 1) {
        for (int j = 0; j < NUM_SUBFACE; j++) {
          Subface subface = Subface(j);
          int rank = neighborInfo[face].getRank(subface);
          if (rank == myrank) {
            localSeparateFaceList[xyz].push_back(FaceID(localID, face, subface));
          }
          else if (rank != MPI::PROC_NULL) {
            sendSeparateFaceList[xyz][rank].push_back(FaceID(localID, face, subface));
            recvSeparateFaceList[xyz][rank].push_back(FaceID(localID, face, subface));
          }
        }
      }
    }
  }

  // 送信側のFaceIDオーダーに合わせるため，RecvFaceListをソートする
  for (int xyz = 0; xyz < 3; xyz++) {
    FaceListMap::iterator it = recvSeparateFaceList[xyz].begin();
    for (; it != recvSeparateFaceList[xyz].end(); ++it) {
      std::sort(it->second.begin(), it->second.end(), RecvFaceComp(blockList));
    }
  }

  separateFaceListPrepared = true;
}


/// フェイスリスト(3方向別・レベル差別)の設定.
void BlockManager::setSeparateLevelDiffFaceLists() {
  if (separateLevelDiffFaceListPrepared) return;

  int myrank = comm.Get_rank();

  for (int localID = 0; localID < numBlock; localID++) {
    BlockBase* block = blockList[localID];
    const NeighborInfo* neighborInfo = block->getNeighborInfo();

    for (int i = 0; i < NUM_FACE; i++) {
      Face face = Face(i);
      int xyz = faceToXYZ(face);

      if (neighborInfo[face].getLevelDifference() == 0) {
        int rank = neighborInfo[face].getRank();
        if (rank == myrank) {
          localSeparateLevelDiffFaceList[xyz][LD_0].push_back(FaceID(localID, face));
        }
        else if (rank != MPI::PROC_NULL) {
          sendSeparateLevelDiffFaceList[xyz][LD_0][rank].push_back(FaceID(localID, face));
          recvSeparateLevelDiffFaceList[xyz][LD_0][rank].push_back(FaceID(localID, face));
        }
      }
      else if (neighborInfo[face].getLevelDifference() == -1) {
        int rank = neighborInfo[face].getRank();
        if (rank == myrank) {
          localSeparateLevelDiffFaceList[xyz][LD_M1].push_back(FaceID(localID, face));
        }
        else if (rank != MPI::PROC_NULL) {
          sendSeparateLevelDiffFaceList[xyz][LD_M1][rank].push_back(FaceID(localID, face));
          recvSeparateLevelDiffFaceList[xyz][LD_M1][rank].push_back(FaceID(localID, face));
        }
      }
      else if (neighborInfo[face].getLevelDifference() == 1) {
        for (int j = 0; j < NUM_SUBFACE; j++) {
          Subface subface = Subface(j);
          int rank = neighborInfo[face].getRank(subface);
          if (rank == myrank) {
            localSeparateLevelDiffFaceList[xyz][LD_P1].push_back(FaceID(localID, face, subface));
          }
          else if (rank != MPI::PROC_NULL) {
            sendSeparateLevelDiffFaceList[xyz][LD_P1][rank].push_back(FaceID(localID, face, subface));
            recvSeparateLevelDiffFaceList[xyz][LD_P1][rank].push_back(FaceID(localID, face, subface));
          }
        }
      }
    }
  }

  // 送信側のFaceIDオーダーに合わせるため，RecvFaceListをソートする
  for (int xyz = 0; xyz < 3; xyz++) {
    for (int ld = 0; ld < 3; ld++) {
      FaceListMap::iterator it = recvSeparateLevelDiffFaceList[xyz][ld].begin();
      for (; it != recvSeparateLevelDiffFaceList[xyz][ld].end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), RecvFaceComp(blockList));
      }
    }
  }

  separateLevelDiffFaceListPrepared = true;
}


/// 各ブロックで，指定したデータクラスの隣接データクラス情報を設定.
void BlockManager::setDataClassNeighbor(int dataClassID, const FaceList& faceList) {
  FaceList::const_iterator it = faceList.begin();
  for (; it != faceList.end(); ++it) {
    int id = it->id;
    Face face = it->face;
    Subface subface = it->subface;
    BlockBase* block = blockList[id];
    UpdatableDataClass* dataClass
      = dynamic_cast<UpdatableDataClass*>(block->getDataClass(dataClassID));
    assert(dataClass);
    int neighborBlockID = (block->getNeighborInfo())[face].getID(subface)
                            - startID;
    dataClass->setNeighbor(face, subface,
                           blockList[neighborBlockID]->getDataClass(dataClassID));
  }
}


/// 各ブロックで，指定したデータクラスの送信バッファポインタを設定.
void BlockManager::setSendBufferPointers(SendBuffer* sendBuffer, int dataClassID, const FaceListMap& faceListMap)
{
  FaceListMap::const_iterator it_map = faceListMap.begin();
  for (; it_map != faceListMap.end(); ++it_map) {
    int rank = it_map->first;
    FaceList::const_iterator it = it_map->second.begin();
    for (; it != it_map->second.end(); ++it) {
      int id = it->id;
      Face face = it->face;
      Subface subface = it->subface;
      UpdatableDataClass* dataClass
        = dynamic_cast<UpdatableDataClass*>(
                          blockList[id]->getDataClass(dataClassID));
      sendBuffer->setData(rank,
                          dataClass->getSendBufferByteSize(face, subface),
                          dataClass->getSendBufferPointerSetter(face, subface));
    }
  }
}


/// 各ブロックで，指定したデータクラスの受信バッファポインタを設定.
void BlockManager::setRecvBufferPointers(RecvBuffer* recvBuffer, int dataClassID, const FaceListMap& faceListMap)
{
  FaceListMap::const_iterator it_map = faceListMap.begin();
  for (; it_map != faceListMap.end(); ++it_map) {
    int rank = it_map->first;
    FaceList::const_iterator it = it_map->second.begin();
    for (; it != it_map->second.end(); ++it) {
      int id = it->id;
      Face face = it->face;
      Subface subface = it->subface;
      UpdatableDataClass* dataClass
        = dynamic_cast<UpdatableDataClass*>(
                          blockList[id]->getDataClass(dataClassID));
      recvBuffer->setData(rank,
                          dataClass->getRecvBufferByteSize(face, subface),
                          dataClass->getRecvBufferPointerSetter(face, subface));
    }
  }
}


/// 各ブロックで，指定したデータクラスの仮想セルデータを隣接ブロックからコピー.
void BlockManager::copyVCFromNeighbor(int dataClassID, const FaceList& faceList)
{
	int nf = faceList.size();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
  for (int n=0; n<nf; n++) {
    int id = faceList[n].id;
    Face face = faceList[n].face;
    Subface subface = faceList[n].subface;
    UpdatableDataClass* dataClass
      = dynamic_cast<UpdatableDataClass*>(
                          blockList[id]->getDataClass(dataClassID));
    dataClass->copyFromNeighbor(face, subface);
	}
/*
  FaceList::const_iterator it = faceList.begin();
  for (; it != faceList.end(); ++it) {
    int id = it->id;
    Face face = it->face;
    Subface subface = it->subface;
    UpdatableDataClass* dataClass
      = dynamic_cast<UpdatableDataClass*>(
                          blockList[id]->getDataClass(dataClassID));
    dataClass->copyFromNeighbor(face, subface);
  }
*/
}


/// 各ブロックで，指定したデータクラスの仮想セルデータを送信バッファにコピー.
void BlockManager::copyVCToSendBuffer(int dataClassID, const FaceListMap& faceListMap)
{
  FaceListMap::const_iterator it_map = faceListMap.begin();
  for (; it_map != faceListMap.end(); ++it_map) {
		FaceList faceList = it_map->second;
		int nf = faceList.size();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for(int n=0; n<nf; n++) {
			int id = faceList[n].id;
			Face face = faceList[n].face;
			Subface subface = faceList[n].subface;
			UpdatableDataClass* dataClass
				= dynamic_cast<UpdatableDataClass*>(
														blockList[id]->getDataClass(dataClassID));
      dataClass->copyToCommBuffer(face, subface);
		}
	}

/*
  FaceListMap::const_iterator it_map = faceListMap.begin();
  for (; it_map != faceListMap.end(); ++it_map) {
    FaceList::const_iterator it = it_map->second.begin();
    for (; it != it_map->second.end(); ++it) {
      int id = it->id;
      Face face = it->face;
      Subface subface = it->subface;
      UpdatableDataClass* dataClass
        = dynamic_cast<UpdatableDataClass*>(
                          blockList[id]->getDataClass(dataClassID));
      dataClass->copyToCommBuffer(face, subface);
    }
  }
*/
}


/// 各ブロックで，指定したデータクラスの仮想セルデータを受信バッファからコピー.
void BlockManager::copyVCFromRecvBuffer(int dataClassID, const FaceListMap& faceListMap)
{
  FaceListMap::const_iterator it_map = faceListMap.begin();
  for (; it_map != faceListMap.end(); ++it_map) {
		FaceList faceList = it_map->second;
		int nf = faceList.size();
#ifdef _BLOCK_IS_LARGE_
#else
#pragma omp parallel for
#endif
		for(int n=0; n<nf; n++) {
			int id = faceList[n].id;
			Face face = faceList[n].face;
			Subface subface = faceList[n].subface;
			UpdatableDataClass* dataClass
				= dynamic_cast<UpdatableDataClass*>(
														blockList[id]->getDataClass(dataClassID));
      dataClass->copyFromCommBuffer(face, subface);
		}
	}

/*
  FaceListMap::const_iterator it_map = faceListMap.begin();
  for (; it_map != faceListMap.end(); ++it_map) {
    FaceList::const_iterator it = it_map->second.begin();
    for (; it != it_map->second.end(); ++it) {
      int id = it->id;
      Face face = it->face;
      Subface subface = it->subface;
      UpdatableDataClass* dataClass
        = dynamic_cast<UpdatableDataClass*>(
                          blockList[id]->getDataClass(dataClassID));
      dataClass->copyFromCommBuffer(face, subface);
    }
  }
*/
}


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif
