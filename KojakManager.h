#ifndef _KOJAKMANAGER_H
#define _KOJAKMANAGER_H

#include "KAnalysis.h"
#include "KData.h"
#include "KDB.h"
#include "KIons.h"
#include "KParams.h"

#define VERSION "2.0.0-dev"
#define BDATE "October 16 2019"

class KojakManager {
public:
  KojakManager();

  void clearFiles();

  int setFile(const char* fn);
  int setFile(string& s);
  void setParam(const char* p);
  void setParam(string& s);
  bool setParams(const char* fn);
  bool setParams(string& s);

  bool getBaseFileName(string& base, const char* fName, string& extP);
  int run();


private:
  vector<kFile> files;
  KParams param_obj;
  kParams params;
  
};

#endif