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

      size = configFile.read<int>("size");
      output = configFile.read<string>("output");

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

    if (size < 1) {
      std::cout << "error: 'size' must be > 1" << std::endl;
      exit(EX_READ_CONFIG);
    }

}


void Config::print() const
{
    std::cout.setf(std::ios::showpoint);
    std::cout << "  level:              " << level << std::endl;
    std::cout << "  tree type:          " << treeType << std::endl;
    std::cout << "  block size:         " << size << std::endl;
    std::cout << "  output file:        " << output << std::endl;
    std::cout << "  verbose message:    " << (verbose ? "on" : "off") << std::endl;
}


void Config::bcast(const MPI::Comm& comm)
{
  int bufferSize;
 
  if (comm.Get_rank() == 0) {
    bufferSize = 
               + sizeof(level)
               + sizeof(size)
               + sizeof(verbose)
               + (treeType.size() + 1) + (output.size() + 1);
  } 

  comm.Bcast(&bufferSize, 1, MPI::INT, 0);

  char* buffer = new char[bufferSize];

  if (comm.Get_rank() == 0) {
    int position = 0;
    MPI::INT.Pack(&level, 1, buffer, bufferSize, position, comm);
    MPI::INT.Pack(&size, 1, buffer, bufferSize, position, comm);
    MPI::BOOL.Pack(&verbose, 1, buffer, bufferSize, position, comm);
    MPI::CHAR.Pack(treeType.c_str(), treeType.size()+1, buffer, bufferSize, position, comm);
    MPI::CHAR.Pack(output.c_str(), output.size()+1, buffer, bufferSize, position, comm);
  }

  comm.Bcast(buffer, bufferSize, MPI::CHAR, 0);

  if (comm.Get_rank() != 0) {
    int position = 0;
    MPI::INT.Unpack(buffer, bufferSize, &level, 1, position, comm);
    MPI::INT.Unpack(buffer, bufferSize, &size, 1, position, comm);
    MPI::BOOL.Unpack(buffer, bufferSize, &verbose, 1, position, comm);
    treeType = std::string(&buffer[position]);
    position += (treeType.size() + 1);
    output = std::string(&buffer[position]);
  }

  delete[] buffer;
}
