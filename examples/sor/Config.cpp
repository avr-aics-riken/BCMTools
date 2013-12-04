#include <iostream>
#include <cstdlib>
#include "Config.h"
#include "ConfigFile.h"
#include "BCMTools.h"

Config::Config()
{
}

Config::Config(const char* file)
{
  load(file);
}

void Config::load(const char* file)
{
//  std::cout << std::endl << "Configuration file: " << file << std::endl;

    try {
      ConfigFile configFile(file);

      level = configFile.read<int>("level");
      treeType = configFile.read<string>("treeType");

      type = configFile.read<char>("type");
      size = configFile.read<int>("size");
      vc = configFile.read<int>("vc");
      b0 = configFile.read<double>("b0", 0.0);
      b1 = configFile.read<double>("b1", 1.0);
      output = configFile.read<string>("output");

      nLoopInner = configFile.read<int>("nLoopInner");
      nLoopOuter = configFile.read<int>("nLoopOuter");

      omega = configFile.read<double>("omega", 1.0);

      randomShuffle = configFile.read<bool>("randomShuffle", false);

      separate = configFile.read<bool>("separateVCUpdate", false);

      verbose = configFile.read<bool>("verbose", false);

    }
    catch(ConfigFile::file_not_found& e) {
        std::cout << "error: cannot open configfile: " << e.filename << std::endl;
        exit(EX_OPEN_FILE);
    }
    catch(ConfigFile::key_not_found& e) {
        std::cout << "error: cannot find key: " << e.key << std::endl;
        exit(EX_READ_CONFIG);
    }

    if (!(treeType == "flat" || treeType == "simple")) {
      std::cout << "error: 'treeType' must be 'flat' or 'simple'." << std::endl;
      exit(EX_READ_CONFIG);
    }

    if (!(vc > 0 && vc <= size/2)) {
      std::cout << "error: 'vc' must be 0 < vc && vc <= size/2." << std::endl;
      exit(EX_READ_CONFIG);
    }

    switch (type) {
      case 'X':
        type = 'x';
        break;
      case 'Y':
        type = 'y';
        break;
      case 'Z':
        type = 'z';
        break;
      case 'x':
      case 'y':
      case 'z':
        break;
      default:
        std::cout << "error: 'type' must be 'x', 'y' or 'z'." << std::endl;
        exit(EX_READ_CONFIG);
        break;
    }
}


void Config::print() const
{
    std::cout.setf(std::ios::showpoint);
    std::cout << "  level:              " << level << std::endl;
    std::cout << "  tree type:          " << treeType << std::endl;
    std::cout << "  type:               " << type << std::endl;
    std::cout << "  block size:         " << size << std::endl;
    std::cout << "  vc width:           " << vc << std::endl;
    std::cout << "  boundary values:    " << b0 << ", " << b1 << std::endl;
    std::cout << "  SOR oemga:          " << omega << std::endl;
    std::cout << "  inner loop length:  " << nLoopInner << std::endl;
    std::cout << "  outer loop length:  " << nLoopOuter << std::endl;
    std::cout << "  random shuffle:     " << (randomShuffle ? "on" : "off") << std::endl;
    std::cout << "  separate vc-update: " << (separate ? "on" : "off") << std::endl;
    std::cout << "  output file:        " << output << std::endl;
    std::cout << "  verbose message:    " << (verbose ? "on" : "off") << std::endl;
}


void Config::bcast(const MPI::Comm& comm)
{
  int bufferSize;
 
  if (comm.Get_rank() == 0) {
    bufferSize = sizeof(b0) + sizeof(b1) + sizeof(omega)
               + sizeof(level)
               + sizeof(size) + sizeof(vc) + sizeof(nLoopInner) + sizeof(nLoopOuter)
               + sizeof(randomShuffle) + sizeof(separate) + sizeof(verbose)
               + sizeof(type) + (treeType.size() + 1) + (output.size() + 1);
  } 

  comm.Bcast(&bufferSize, 1, MPI::INT, 0);

  char* buffer = new char[bufferSize];

  if (comm.Get_rank() == 0) {
    int position = 0;
    MPI::DOUBLE.Pack(&b0, 1, buffer,  bufferSize, position, comm);
    MPI::DOUBLE.Pack(&b1, 1, buffer,  bufferSize, position, comm);
    MPI::DOUBLE.Pack(&omega, 1, buffer,  bufferSize, position, comm);
    MPI::INT.Pack(&level, 1, buffer, bufferSize, position, comm);
    MPI::INT.Pack(&size, 1, buffer, bufferSize, position, comm);
    MPI::INT.Pack(&vc, 1, buffer, bufferSize, position, comm);
    MPI::INT.Pack(&nLoopInner, 1, buffer, bufferSize, position, comm);
    MPI::INT.Pack(&nLoopOuter, 1, buffer, bufferSize, position, comm);
    MPI::BOOL.Pack(&randomShuffle, 1, buffer, bufferSize, position, comm);
    MPI::BOOL.Pack(&separate, 1, buffer, bufferSize, position, comm);
    MPI::BOOL.Pack(&verbose, 1, buffer, bufferSize, position, comm);
    MPI::CHAR.Pack(&type, 1, buffer, bufferSize, position, comm);
    MPI::CHAR.Pack(treeType.c_str(), treeType.size()+1, buffer, bufferSize, position, comm);
    MPI::CHAR.Pack(output.c_str(), output.size()+1, buffer, bufferSize, position, comm);
  }

  comm.Bcast(buffer, bufferSize, MPI::CHAR, 0);

  if (comm.Get_rank() != 0) {
    int position = 0;
    MPI::DOUBLE.Unpack(buffer, bufferSize, &b0, 1, position, comm);
    MPI::DOUBLE.Unpack(buffer, bufferSize, &b1, 1, position, comm);
    MPI::DOUBLE.Unpack(buffer, bufferSize, &omega, 1, position, comm);
    MPI::INT.Unpack(buffer, bufferSize, &level, 1, position, comm);
    MPI::INT.Unpack(buffer, bufferSize, &size, 1, position, comm);
    MPI::INT.Unpack(buffer, bufferSize, &vc, 1, position, comm);
    MPI::INT.Unpack(buffer, bufferSize, &nLoopInner, 1, position, comm);
    MPI::INT.Unpack(buffer, bufferSize, &nLoopOuter, 1, position, comm);
    MPI::BOOL.Unpack(buffer, bufferSize, &randomShuffle, 1, position, comm);
    MPI::BOOL.Unpack(buffer, bufferSize, &separate, 1, position, comm);
    MPI::BOOL.Unpack(buffer, bufferSize, &verbose, 1, position, comm);
    MPI::CHAR.Unpack(buffer, bufferSize, &type, 1, position, comm);
    treeType = std::string(&buffer[position]);
    position += (treeType.size() + 1);
    output = std::string(&buffer[position]);
  }

  delete[] buffer;
}
