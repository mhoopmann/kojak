/*
Copyright 2017, Michael R. Hoopmann, Institute for Systems Biology

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef _KTOPPEPS_H
#define _KTOPPEPS_H

#include "KStructs.h"
#include <list>

class KTopPeps{
public:

  KTopPeps();
  KTopPeps(const KTopPeps& c);
  ~KTopPeps();

  KTopPeps& operator=(const KTopPeps& c);

  int singletBins;
  int singletCount;
  int singletMax;
  
  kSingletScoreCard*    singletFirst;   //pointer to start of linked list
  kSingletScoreCard*    singletLast;    //pointer to end of linked list
  std::list<kSingletScoreCard*>**  singletList;

  void  checkSingletScore(kSingletScoreCard& s);
  void  resetSingletList(double mass);

};

#endif
