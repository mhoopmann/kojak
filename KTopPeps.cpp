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

#include "KTopPeps.h"

using namespace std;

KTopPeps::KTopPeps(){
  singletMax = 0;
}

void KTopPeps::checkSingletScore(kSingletScoreCard& s){

  singletList.push_back(s);
  list<kSingletScoreCard>::reverse_iterator i=singletList.rbegin();
  while(true){
    list<kSingletScoreCard>::reverse_iterator i2=i;
    ++i2;
    if(i2==singletList.rend()) break;
    if(i->simpleScore>i2->simpleScore){
      swap(*i,*i2);
    } else break;
    ++i;
  }
  if(singletList.size()>singletMax) singletList.pop_back();

}

