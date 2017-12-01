#include "KTopPeps.h"

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
  if (s.simpleScore<singletLast->simpleScore){
    //check if we need to store the singlet
    if (singletCount == singletMax) return;

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

    if (singletCount<singletMax) {
      singletCount++;
    } else {
      cur = singletLast;
      singletLast = singletLast->prev;
      singletLast->next = NULL;

      //delete expired singlet from list
      ind = (int)(cur->mass / 10);
      if (singletList[ind]->size() == 1) {
        delete singletList[ind];
        singletList[ind] = NULL;
      } else {
        list<kSingletScoreCard*>::iterator it = singletList[ind]->begin();
        while (*it != cur) it++;
        singletList[ind]->erase(it);
      }

      delete cur;
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

  if (singletCount<singletMax) {
    singletCount++;
  } else {
    cur = singletLast;
    singletLast = singletLast->prev;
    singletLast->next = NULL;

    //delete expired singlet from list
    ind = (int)(cur->mass / 10);
    if (singletList[ind]->size() == 1) {
      delete singletList[ind];
      singletList[ind] = NULL;
    } else {
      list<kSingletScoreCard*>::iterator it = singletList[ind]->begin();
      while (*it != cur) it++;
      singletList[ind]->erase(it);
    }

    delete cur;
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
