#ifndef _KPARAMS_H
#define _KPARAMS_H

#include "KStructs.h"

class KParams {
public:
  KParams();
  KParams(kParams* p);
  ~KParams();

  bool parseConfig(char* fname);

private:

  kParams* params;

  void parse(char* cmd);
  void warn(char* c, int i);

};

#endif
