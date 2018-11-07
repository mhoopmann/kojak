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
  singletCount = 0;
  singletFirst = NULL;
  singletLast = NULL;
  singletMax = 0;

  singletList = NULL;
  singletBins = 0;
}

KTopPeps::KTopPeps(const KTopPeps& c){
  singletCount = c.singletCount;
  singletMax = c.singletMax;
  singletFirst = NULL;
  singletLast = NULL;
  kSingletScoreCard* sc = NULL;
  kSingletScoreCard* tmp = c.singletFirst;
  if (tmp != NULL) {
    singletFirst = new kSingletScoreCard(*tmp);
    sc = singletFirst;
    tmp = tmp->next;
    while (tmp != NULL){
      sc->next = new kSingletScoreCard(*tmp);
      sc->next->prev = sc;
      sc = sc->next;
      tmp = tmp->next;
    }
    singletLast = sc;
  }

  singletBins = c.singletBins;
  if (singletBins == 0) singletList = NULL;
  else {
    singletList = new list<kSingletScoreCard*>*[singletBins];
    for (size_t j = 0; j<singletBins; j++){
      if (c.singletList[j] == NULL) singletList[j] = NULL;
      else {
        singletList[j] = new list<kSingletScoreCard*>;
        list<kSingletScoreCard*>::iterator it = c.singletList[j]->begin();
        while (it != c.singletList[j]->end()){
          singletList[j]->emplace_back(*it);
          it++;
        }
      }
    }
  }
}

KTopPeps::~KTopPeps(){
  while (singletFirst != NULL){
    kSingletScoreCard* tmp = singletFirst;
    singletFirst = singletFirst->next;
    delete tmp;
  }
  singletLast = NULL;

  if (singletList != NULL){
    for (size_t j = 0; j<singletBins; j++){
      if (singletList[j] != NULL) delete singletList[j];
    }
    delete[] singletList;
  }
}

KTopPeps& KTopPeps::operator=(const KTopPeps& c){
  if(this!=&c){
    size_t j;

    singletCount = c.singletCount;
    singletMax = c.singletMax;
    
    while (singletFirst != NULL){
      kSingletScoreCard* tmp = singletFirst;
      singletFirst = singletFirst->next;
      delete tmp;
    }
    singletLast = NULL;

    kSingletScoreCard* sc = NULL;
    kSingletScoreCard* tmp = c.singletFirst;
    if (tmp != NULL) {
      singletFirst = new kSingletScoreCard(*tmp);
      sc = singletFirst;
      tmp = tmp->next;
      while (tmp != NULL){
        sc->next = new kSingletScoreCard(*tmp);
        sc->next->prev = sc;
        sc = sc->next;
        tmp = tmp->next;
      }
      singletLast = sc;
    }

    if (singletList != NULL){
      for (j = 0; j<singletBins; j++){
        if (singletList[j] != NULL) delete singletList[j];
      }
      delete[] singletList;
    }
    singletBins = c.singletBins;
    if (singletBins == 0) singletList = NULL;
    else {
      singletList = new list<kSingletScoreCard*>*[singletBins];
      for (size_t j = 0; j<singletBins; j++){
        if (c.singletList[j] == NULL) singletList[j] = NULL;
        else {
          singletList[j] = new list<kSingletScoreCard*>;
          list<kSingletScoreCard*>::iterator it = c.singletList[j]->begin();
          while (it != c.singletList[j]->end()){
            singletList[j]->emplace_back(*it);
            it++;
          }
        }
      }
    }
  }
  return *this;
}

void KTopPeps::checkSingletScore(kSingletScoreCard& s){

  kSingletScoreCard* sc;
  kSingletScoreCard* cur;
  size_t ind;

  //If list is empty, add the score card
  if (singletCount == 0){
    singletFirst = new kSingletScoreCard(s);
    singletLast = singletFirst;
    singletCount++;

    //add to singlet list
    ind = (int)(s.mass / 10);
    if (singletList[ind] == NULL) singletList[ind] = new list<kSingletScoreCard*>;
    singletList[ind]->emplace_back(singletFirst);

    return;
  }

  //check if we can just add to the end
  if (s.simpleScore <= singletLast->simpleScore){
    //check if we need to store the singlet
    //if(singletCount>=singletMax) return;

    singletLast->next = new kSingletScoreCard(s);
    singletLast->next->prev = singletLast;
    singletLast = singletLast->next;
    singletCount++;

    //add to singlet list
    ind = (int)(s.mass / 10);
    if (singletList[ind] == NULL) singletList[ind] = new list<kSingletScoreCard*>;
    singletList[ind]->emplace_back(singletLast);

    return;
  }

  //check if it goes in the front
  if (s.simpleScore >= singletFirst->simpleScore){
    singletFirst->prev = new kSingletScoreCard(s);
    singletFirst->prev->next = singletFirst;
    singletFirst = singletFirst->prev;

    //add to singlet list
    ind = (int)(s.mass / 10);
    if (singletList[ind] == NULL) singletList[ind] = new list<kSingletScoreCard*>;
    singletList[ind]->emplace_back(singletFirst);
    singletCount++;

    if (singletCount>singletMax){
      int i = singletCount;
      cur = singletLast;
      while (i>singletMax){ //step to singletMax position
        cur = cur->prev;
        i--;
      }
      while (cur->next != NULL && cur->next->simpleScore == cur->simpleScore){ //step to first instance of score lower than singletMax
        cur = cur->next;
      }

      //delete everything hereafter
      while (cur->next != NULL){
        sc = cur->next;
        cur->next = sc->next;
        //if(sc->next!=NULL) sc->next->prev=cur; //is this necessary if they're all going to go?

        ind = (int)(sc->mass / 10);
        if (singletList[ind]->size() == 1) {
          delete singletList[ind];
          singletList[ind] = NULL;
        } else {
          list<kSingletScoreCard*>::iterator it = singletList[ind]->begin();
          while (*it != sc) it++;
          singletList[ind]->erase(it);
        }
        delete sc;
        singletCount--;
      }

      singletLast = cur;
    }
    return;
  }


  //scan to find insertion point
  cur = singletFirst->next;
  int i = 1;
  while (s.simpleScore < cur->simpleScore){
    i++;
    cur = cur->next;
  }

  sc = new kSingletScoreCard(s);
  sc->prev = cur->prev;
  sc->next = cur;
  cur->prev->next = sc;
  cur->prev = sc;
  if (sc->prev == NULL) singletFirst = sc;

  //add to singlet list
  ind = (int)(s.mass / 10);
  if (singletList[ind] == NULL) singletList[ind] = new list<kSingletScoreCard*>;
  singletList[ind]->emplace_back(sc);
  singletCount++;

  if (singletCount>singletMax){
    int i = singletCount;
    cur = singletLast;
    while (i>singletMax){ //step to singletMax position
      cur = cur->prev;
      i--;
    }
    while (cur->next != NULL && cur->next->simpleScore == cur->simpleScore){ //step to first instance of score lower than singletMax
      cur = cur->next;
    }

    //delete everything hereafter
    while (cur->next != NULL){
      sc = cur->next;
      cur->next = sc->next;
      //if(sc->next!=NULL) sc->next->prev=cur; //is this necessary if they're all going to go?

      ind = (int)(sc->mass / 10);
      if (singletList[ind]->size() == 1) {
        delete singletList[ind];
        singletList[ind] = NULL;
      } else {
        list<kSingletScoreCard*>::iterator it = singletList[ind]->begin();
        while (*it != sc) it++;
        singletList[ind]->erase(it);
      }
      delete sc;
      singletCount--;
    }

    singletLast = cur;
  }

}

void KTopPeps::resetSingletList(double mass){
  size_t j;
  if (singletList != NULL){
    for (j = 0; j<singletBins; j++){
      if (singletList[j] != NULL) delete singletList[j];
    }
    delete[] singletList;
  }
  singletBins = (int)(mass / 10 + 1);
  singletList = new list<kSingletScoreCard*>*[singletBins];
  for (j = 0; j<singletBins; j++) singletList[j] = NULL;
}
