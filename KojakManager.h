#ifndef _KOJAKMANAGER_H
#define _KOJAKMANAGER_H

#include "KLog.h"
#include "KParams.h"

#define VERSION "2.0.0 alpha 8"
#define BDATE "August 25 2021"

class KojakManager {
public:
  KojakManager();

  void clearFiles();

  int setFile(const char* fn);
  int setFile(std::string& s);
  void setParam(const char* p);
  void setParam(std::string& s);
  bool setParams(const char* fn);
  bool setParams(std::string& s);

  bool getBaseFileName(std::string& base, const char* fName, std::string& extP);
  int run();


private:
  std::vector<kFile> files;
  KLog log;
  std::string paramFile;
  KParams param_obj;
  kParams params;

};

#endif
