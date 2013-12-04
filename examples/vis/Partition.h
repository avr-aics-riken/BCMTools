#ifndef PARTITION_H
#define PARTITION_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>
#include "mpi.h"

class Partition {

  int nProcs;
  int nItems;

  std::vector<int> end;

public:

  Partition(int nProcs, int nItems) : nProcs(nProcs), nItems(nItems), end(nProcs) {
    int m = nItems / nProcs;
    int r = nItems % nProcs;
    int end0 = 0;
    for (int i = 0; i < nProcs; i++) {
      if (i < r) {
        end[i] = end0 + m + 1;
      } else {
        end[i] = end0 + m;
      }
      end0 = end[i];
    }
    assert(end[nProcs-1] == nItems);
  }

  ~Partition() {}

  int getStart(int rank) const {
    assert(0 <= rank && rank < nProcs);
    if (rank == 0) return 0;
    return end[rank-1];
  }

  int getEnd(int rank) const {
    assert(0 <= rank && rank < nProcs);
    return end[rank];
  }

  int getNum(int rank) const {
    assert(0 <= rank && rank < nProcs);
    if (rank == 0) return end[0];
    return end[rank] - end[rank-1];
  }

  int getRank(int i) const {
    if (i < 0 || i >= nItems) return MPI::PROC_NULL;
    std::vector<int>::const_iterator it = std::upper_bound(end.begin(), end.end(), i);
    assert(it != end.end());
    return it - end.begin();
  }

  void print() const {
    int start = 0;
    for (int rank = 0; rank < nProcs; rank++) {
      std::cout << rank << ": [" << start << ":" << end[rank]-1 << "]"
             // << " [" << getStart(rank) << ":" << getEnd(rank) << ")"
                << " #" << end[rank]-start << std::endl;
      start = end[rank];
    }
  }

};

#endif // PARTITION_H
