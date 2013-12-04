#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <iostream>
#include "BCMTools.h"

struct Pedigree {
  int level;
  int x;
  int y;
  int z;

  Pedigree() : level(0), x(0), y(0), z(0) {}

  Pedigree(int level, int x, int y, int z) : level(level), x(x), y(y), z(z) {
    if (level >= 0) {
      assert(x < (1 << level));
      assert(y < (1 << level));
      assert(z < (1 << level));
    }
  }

  Pedigree(const Pedigree& parent, int ijk) {
    assert(0 <= ijk && ijk < 8);
    int i = ijk & 0x01;
    int j = (ijk >> 1) & 0x01;
    int k = (ijk >> 2) & 0x01;
    x = parent.x * 2 + i;
    y = parent.y * 2 + j;
    z = parent.z * 2 + k;
    level = parent.level + 1;
  }

  ~Pedigree() {}

  int getX(int level) const { return (x >> (this->level - level)) & 0x01; }
  int getY(int level) const { return (y >> (this->level - level)) & 0x01; }
  int getZ(int level) const { return (z >> (this->level - level)) & 0x01; }

  int getChildId(int level) const {
    return getX(level) + getY(level) * 2 + getZ(level) * 4;
  }

  Pedigree findNeighbor(Face face, int max, bool periodic=false) const {
    if (level == 0) {
      if (periodic) return Pedigree(level, x, y, z);
      else return Pedigree(-1, x, y, z);
    }
    int xx = x;
    int yy = y;
    int zz = z;
    switch (face) {
      case X_M:
        --xx;
        if (xx < 0) {
          if (periodic) xx += max;
          else return Pedigree(-1, xx, yy, zz);
        }
        break;
      case X_P:
        ++xx;
        if (xx >= max) {
          if (periodic) xx -= max;
          else return Pedigree(-1, xx, yy, zz);
        }
        break;
      case Y_M:
        --yy;
        if (yy < 0) {
          if (periodic) yy += max;
          else return Pedigree(-1, xx, yy, zz);
        }
        break;
      case Y_P:
        ++yy;
        if (yy >= max) {
          if (periodic) yy -= max;
          else return Pedigree(-1, xx, yy, zz);
        }
        break;
      case Z_M:
        --zz;
        if (zz < 0) {
          if (periodic) zz += max;
          else return Pedigree(-1, xx, yy, zz);
        }
        break;
      case Z_P:
        ++zz;
        if (zz >= max) {
          if (periodic) zz -= max;
          else return Pedigree(-1, xx, yy, zz);
        }
        break;
      default:
        return Pedigree(-1, xx, yy, zz);
    }
    return Pedigree(level, xx, yy, zz);
  }
  
};


inline std::ostream& operator<<(std::ostream& os, const Pedigree& p) {
  return os << "(" << p.level << "(" << p.x << "," << p.y << "," << p.z << "))";
}

#endif // PEDIGREE_H
