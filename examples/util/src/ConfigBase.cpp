#include "ConfigBase.h"
#include <iostream>
#include <sstream>


/// コンストラクタ.
ConfigBase::ConfigBase(MPI::Comm& comm) : comm(comm)
{
}


/// デストラクタ.
ConfigBase::~ConfigBase()
{
  delete configFile;
}


/// 設定ファイル読み込み.
void ConfigBase::load(const char* file)
{
  if (comm.Get_rank() == 0) {
    try {
      configFile = new ConfigFile(file);
      parse();
    }
    catch(ConfigFile::file_not_found& e) {
      std::cout << "error: cannot open configfile: " << e.filename << std::endl;
      errorExit("input config file.");
    }
    catch(ConfigFile::key_not_found& e) {
      std::cout << "error: cannot find key: " << e.key << std::endl;
      errorExit("input config file.");
    }
    if (!validate()) errorExit("input config file.");
    broadcastConfigFile(configFile);
  }
  else {
    configFile = new ConfigFile;
    receiveConfigFile(configFile);
    parse();
  }

}


/// ConfigFileオブジェクトの内容をrank0からブロードキャスト.
void ConfigBase::broadcastConfigFile(const ConfigFile* configFile)
{
  std::ostringstream outStr;
  outStr << *configFile;

  int size = outStr.str().size() + 1;
  comm.Bcast(&size, 1, MPI::INT, 0);

  comm.Bcast((void*)outStr.str().c_str(), size, MPI::CHAR, 0);
  // constを消すためにキャストが必要
}


/// ConfigFileオブジェクトの内容をrank0から受信.
void ConfigBase::receiveConfigFile(ConfigFile* configFile)
{
  int size;
  comm.Bcast(&size, 1, MPI::INT, 0);

  char* buffer = new char[size];
  comm.Bcast(buffer, size, MPI::CHAR, 0);
  
  std::istringstream inStr(buffer);
  inStr >> *configFile;

  delete[] buffer;
}


/// エラー終了.
void ConfigBase::errorExit(const char* message, int code) 
{
  std::cout << "error: " << message << std::endl;
  comm.Abort(code);
}
