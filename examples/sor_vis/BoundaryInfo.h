#ifndef BUNDARY_INFO_H
#define BUNDARY_INFO_H

#include <iostream>
#include "BCMTools.h"

class BoundaryInfo {

public:

  enum Type { INNER, DIRICHLET, NEUMANN, PERIODIC };

private:

  Type type;

  int id;

public:

  BoundaryInfo() : type(INNER), id(-1) {}

  ~BoundaryInfo() {}

  void setType(Type type) { this->type = type; }

  void setID(int id) { this->id = id; }

  Type getType() const { return type; }

  int getID() const { return id; }

  void print() const {
    switch (type) {
      case INNER:     std::cout << "type=INNER, "; break;
      case DIRICHLET: std::cout << "type=DIRICHLET, "; break;
      case NEUMANN:   std::cout << "type=NEUMANN, "; break;
      case PERIODIC:  std::cout << "type=PERIODIC, "; break;
      default: break;
    }
    std::cout << "id=" << id << std::endl;
  }

  static void print(const BoundaryInfo* bInfo) {
    for (int i = 0; i < NUM_FACE; i++) {
      Face face = Face(i);
      switch (face) {
        case X_M: std::cout << "X_M: "; break;
        case X_P: std::cout << "X_P: "; break;
        case Y_M: std::cout << "Y_M: "; break;
        case Y_P: std::cout << "Y_P: "; break;
        case Z_M: std::cout << "Z_M: "; break;
        case Z_P: std::cout << "Z_P: "; break;
        default: break;
      }
      bInfo[i].print();
    }
  }

};


#endif // BUNDARY_INFO_H
