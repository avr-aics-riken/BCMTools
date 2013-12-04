#include <iostream>
#include <cstdlib>
#include "Partition.h"

int main(int argc, char** argv)
{
  if (argc != 3) {
    std::cout << "usage: " << argv[0] << " nProcs nItems" << std::endl;
    exit(1);
  }

  int nProcs = atoi(argv[1]);
  int nItems = atoi(argv[2]);

  std::cout << "nProcs = " << nProcs << std::endl;
  std::cout << "nItems = " << nItems << std::endl;

  Partition p(nProcs, nItems);

  p.print();

  for (;;) {
    int i;
    std::cin >> i ;
    std::cout << i << "->" << p.getRank(i) << std::endl;
  }

  return 0;
}
