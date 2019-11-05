/*
Copyright 2014, Michael R. Hoopmann, Institute for Systems Biology

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

#include "KSpectrum.h"
#include <iostream>

using namespace std;

/*============================
  Constructors & Destructors
============================*/
KSpectrum::KSpectrum(const int& i, const double& bs, const double& os){
  binOffset=os;
  binSize=bs;
  instrumentPrecursor=false;
  invBinSize=1.0/binSize;
  charge = 0; 
  maxIntensity=0;
  mz = 0;
  precursor = new vector<kPrecursor>;
  singlets = new vector<KTopPeps>;
  spec = new vector<kSpecPoint>;
  scanNumber = 0;
  rTime = 0;
  xCorrArraySize=0;
  xCorrSparseArraySize=0;
  xCorrSparseArray=NULL;
  
  singletCount=0;
  singletFirst=NULL;
  singletLast=NULL;
  singletMax=i;

  kojakSparseArray=NULL;
  kojakBins=0;

  singletList=NULL;
  singletBins=0;

  lowScore=0;

  int j;
  for (j = 0; j<HISTOSZ; j++) histogram[j] = 0;
  histogramCount = 0;
  histoMaxIndex = 0;

  for (j = 0; j<HISTOSZ; j++) histogramSinglet[j] = 0;
  histogramSingletCount = 0;

  //diagnostics - probably temporary
  tmpIntercept=0;
  tmpSlope=0;
  tmpIStartCorr=0;
  tmpINextCorr=0;
  tmpIMaxCorr=0;
  tmpRSquare=0;
  tmpSingletIntercept = 0;
  tmpSingletSlope = 0;
  tmpSingletIStartCorr = 0;
  tmpSingletINextCorr = 0;
  tmpSingletIMaxCorr = 0;
  tmpSingletRSquare = 0;
  tmpSingCount=0;
  tmpHistCount=0;

  cc=0;
  sc=0;
}

KSpectrum::KSpectrum(const KSpectrum& p){
  unsigned int i;
  int j;
  spec = new vector<kSpecPoint>;
  for(i=0;i<p.spec->size();i++) spec->push_back(p.spec->at(i));
  for(i=0;i<20;i++) topHit[i]=p.topHit[i];
  precursor = new vector<kPrecursor>;
  for(i=0;i<p.precursor->size();i++) precursor->push_back(p.precursor->at(i));
  singlets = new vector<KTopPeps>;
  for (i = 0; i<p.singlets->size(); i++) singlets->push_back(p.singlets->at(i));

  binOffset = p.binOffset;
  binSize = p.binSize;
  instrumentPrecursor = p.instrumentPrecursor;
  invBinSize = p.invBinSize;
  charge= p.charge;
  maxIntensity = p.maxIntensity;
  mz = p.mz;
  scanNumber = p.scanNumber;
  rTime = p.rTime;
  xCorrArraySize = p.xCorrArraySize;
  xCorrSparseArraySize = p.xCorrSparseArraySize;
  lowScore=p.lowScore;

  for (i = 0; i<HISTOSZ; i++) histogram[i] = p.histogram[i];
  histogramCount = p.histogramCount;
  histoMaxIndex = p.histoMaxIndex;

  for (i = 0; i<HISTOSZ; i++) histogramSinglet[i] = p.histogramSinglet[i];
  histogramSingletCount = p.histogramSingletCount;

  //diagnostics - probably temporary
  tmpIntercept = p.tmpIntercept;
  tmpSlope = p.tmpSlope;
  tmpIStartCorr = p.tmpIStartCorr;
  tmpINextCorr = p.tmpINextCorr;
  tmpIMaxCorr = p.tmpIMaxCorr;
  tmpRSquare = p.tmpRSquare;
  tmpSingletIntercept = p.tmpSingletIntercept;
  tmpSingletSlope = p.tmpSingletSlope;
  tmpSingletIStartCorr = p.tmpSingletIStartCorr;
  tmpSingletINextCorr = p.tmpSingletINextCorr;
  tmpSingletIMaxCorr = p.tmpSingletIMaxCorr;
  tmpSingletRSquare = p.tmpSingletRSquare;
  tmpSingCount = p.tmpSingCount;
  tmpHistCount = p.tmpHistCount;

  cc=p.cc;
  sc=p.sc;

  singletCount=p.singletCount;
  singletMax=p.singletMax;
  singletFirst=NULL;
  singletLast=NULL;
  kSingletScoreCard* sc=NULL;
  kSingletScoreCard* tmp=p.singletFirst;
  if(tmp!=NULL) {
    singletFirst=new kSingletScoreCard(*tmp);
    sc=singletFirst;
    tmp=tmp->next;
    while(tmp!=NULL){
      sc->next=new kSingletScoreCard(*tmp);
      sc->next->prev=sc;
      sc=sc->next;
      tmp=tmp->next;
    }
    singletLast=sc;
  }

  if(p.xCorrSparseArray==NULL){
    xCorrSparseArray=NULL;
  } else {
    xCorrSparseArray = (kSparseMatrix *)calloc((size_t)xCorrSparseArraySize, (size_t)sizeof(kSparseMatrix));
    for(j=0;j<xCorrSparseArraySize;j++) xCorrSparseArray[j]=p.xCorrSparseArray[j];
  }

  kojakBins=p.kojakBins;
  if(p.kojakSparseArray==NULL){
    kojakSparseArray=NULL;
  } else {
    for(j=0;j<kojakBins;j++){
      if(p.kojakSparseArray[j]==NULL){
        kojakSparseArray[j]=NULL;
      } else {
        kojakSparseArray[j] = new char[(int)invBinSize+1];
        for(i=0;i<(unsigned int)invBinSize+1;i++) kojakSparseArray[j][i]=p.kojakSparseArray[j][i];
      }
    }
  }

  singletBins=p.singletBins;
  if(singletBins==0) singletList=NULL;
  else {
    singletList = new list<kSingletScoreCard*>*[singletBins];
    for(j=0;j<singletBins;j++){
      if(p.singletList[j]==NULL) singletList[j]=NULL;
      else {
        singletList[j] = new list<kSingletScoreCard*>;
        list<kSingletScoreCard*>::iterator it = p.singletList[j]->begin();
        while(it!=p.singletList[j]->end()){
          singletList[j]->emplace_back(*it);
          it++;
        }
      }
    }
  }


}
  
KSpectrum::~KSpectrum(){
  delete spec;
  delete precursor;
  delete singlets;
  if(xCorrSparseArray!=NULL) free(xCorrSparseArray);

  while(singletFirst!=NULL){
    kSingletScoreCard* tmp=singletFirst;
    singletFirst=singletFirst->next;
    delete tmp;
  }
  singletLast=NULL;

  int j;
  if(kojakSparseArray!=NULL){
    for(j=0;j<kojakBins;j++){
      if(kojakSparseArray[j]!=NULL) delete [] kojakSparseArray[j];
    }
    delete [] kojakSparseArray;
  }

  if (singletList != NULL){
    for (j = 0; j<singletBins; j++){
      if (singletList[j] != NULL) delete singletList[j];
    }
    delete[] singletList;
  }
}


/*============================
  Operators
============================*/
KSpectrum& KSpectrum::operator=(const KSpectrum& p){
  if(this!=&p){
    unsigned int i;
    int j;
    delete spec;
    spec = new vector<kSpecPoint>;
    for(i=0;i<p.spec->size();i++) spec->push_back(p.spec->at(i));
    for(i=0;i<20;i++) topHit[i]=p.topHit[i];
    delete precursor;
    precursor = new vector<kPrecursor>;
    for(i=0;i<p.precursor->size();i++) precursor->push_back(p.precursor->at(i));
    delete singlets;
    singlets = new vector<KTopPeps>;
    for (i = 0; i<p.singlets->size(); i++) singlets->push_back(p.singlets->at(i));

    binOffset = p.binOffset;
    binSize = p.binSize;
    instrumentPrecursor = p.instrumentPrecursor;
    invBinSize = p.invBinSize;
    charge = p.charge;
    maxIntensity = p.maxIntensity;
    mz = p.charge;
    scanNumber = p.scanNumber;
    rTime = p.rTime;
    xCorrArraySize = p.xCorrArraySize;
    xCorrSparseArraySize = p.xCorrSparseArraySize;
    lowScore = p.lowScore;

    for (i = 0; i<HISTOSZ; i++) histogram[i] = p.histogram[i];
    histogramCount = p.histogramCount;
    histoMaxIndex = p.histoMaxIndex;

    for (i = 0; i<HISTOSZ; i++) histogramSinglet[i] = p.histogramSinglet[i];
    histogramSingletCount = p.histogramSingletCount;

    //diagnostics - probably temporary
    tmpIntercept = p.tmpIntercept;
    tmpSlope = p.tmpSlope;
    tmpIStartCorr = p.tmpIStartCorr;
    tmpINextCorr = p.tmpINextCorr;
    tmpIMaxCorr = p.tmpIMaxCorr;
    tmpRSquare = p.tmpRSquare;
    tmpSingletIntercept = p.tmpSingletIntercept;
    tmpSingletSlope = p.tmpSingletSlope;
    tmpSingletIStartCorr = p.tmpSingletIStartCorr;
    tmpSingletINextCorr = p.tmpSingletINextCorr;
    tmpSingletIMaxCorr = p.tmpSingletIMaxCorr;
    tmpSingletRSquare = p.tmpSingletRSquare;
    tmpSingCount = p.tmpSingCount;
    tmpHistCount = p.tmpHistCount;

    cc = p.cc;
    sc = p.sc;

    singletCount=p.singletCount;
    singletMax=p.singletMax;
    singletFirst=NULL;
    singletLast=NULL;
    kSingletScoreCard* sc=NULL;
    kSingletScoreCard* tmp=p.singletFirst;
    if(tmp!=NULL) {
      singletFirst=new kSingletScoreCard(*tmp);
      sc=singletFirst;
      tmp=tmp->next;
      while(tmp!=NULL){
        sc->next=new kSingletScoreCard(*tmp);
        sc->next->prev=sc;
        sc=sc->next;
        tmp=tmp->next;
      }
      singletLast=sc;
    }

    if(xCorrSparseArray!=NULL) free(xCorrSparseArray);
    if(p.xCorrSparseArray==NULL){
      xCorrSparseArray=NULL;
    } else {
      xCorrSparseArray = (kSparseMatrix *)calloc((size_t)xCorrSparseArraySize, (size_t)sizeof(kSparseMatrix));
      for(j=0;j<xCorrSparseArraySize;j++) xCorrSparseArray[j]=p.xCorrSparseArray[j];
    }
    
    if(kojakSparseArray!=NULL){
      for(j=0;j<kojakBins;j++){
        if(kojakSparseArray[j]!=NULL) delete [] kojakSparseArray[j];
      }
      delete [] kojakSparseArray;
    }
    kojakBins=p.kojakBins;
    if(p.kojakSparseArray==NULL){
      kojakSparseArray=NULL;
    } else {
      for(j=0;j<kojakBins;j++){
        if(p.kojakSparseArray[j]==NULL){
          kojakSparseArray[j]=NULL;
        } else {
          kojakSparseArray[j] = new char[(int)invBinSize+1];
          for(i=0;i<(unsigned int)invBinSize+1;i++) kojakSparseArray[j][i]=p.kojakSparseArray[j][i];
        }
      }
    }

    if (singletList != NULL){
      for (j = 0; j<singletBins; j++){
        if (singletList[j] != NULL) delete singletList[j];
      }
      delete[] singletList;
    }
    singletBins = p.singletBins;
    if (singletBins == 0) singletList = NULL;
    else {
      singletList = new list<kSingletScoreCard*>*[singletBins];
      for (j = 0; j<singletBins; j++){
        if (p.singletList[j] == NULL) singletList[j] = NULL;
        else {
          singletList[j] = new list<kSingletScoreCard*>;
          list<kSingletScoreCard*>::iterator it = p.singletList[j]->begin();
          while (it != p.singletList[j]->end()){
            singletList[j]->emplace_back(*it);
            it++;
          }
        }
      }
    }
  }
  return *this;
}

kSpecPoint& KSpectrum::operator [](const int &i){
  return spec->at(i);
}


/*============================
  Accessors
============================*/
double KSpectrum::getBinOffset(){
  return binOffset;
}

int KSpectrum::getCharge(){
  return charge;
}

bool KSpectrum::getInstrumentPrecursor(){
  return instrumentPrecursor;
}

double KSpectrum::getInvBinSize(){
  return invBinSize;
}

float KSpectrum::getMaxIntensity(){
  return maxIntensity;
}

double KSpectrum::getMZ(){
  return mz;
}

kPrecursor& KSpectrum::getPrecursor(int i){
  return precursor->at(i);
}

kPrecursor* KSpectrum::getPrecursor2(int i){
  return &precursor->at(i);
}

float KSpectrum::getRTime(){
  return rTime;
}

int KSpectrum::getScanNumber(){
  return scanNumber;
}

kScoreCard& KSpectrum::getScoreCard(int i){
  return topHit[i];
}

int KSpectrum::getSingletCount(){
  return singletCount;
}

kSingletScoreCard& KSpectrum::getSingletScoreCard(int i){
  if(i>=singletCount) return *singletLast;
  kSingletScoreCard* sc=singletFirst;
  int j=0;
  while(j<i){
    if(sc->next==NULL) break;
    sc=sc->next;
    j++;
  }
  return *sc;
}

KTopPeps* KSpectrum::getTopPeps(int i){
  return &singlets->at(i);
}

int KSpectrum::size(){
  return (int)spec->size();
}

int KSpectrum::sizePrecursor(){
  return (int)precursor->size();
}


/*============================
  Modifiers
============================*/
void KSpectrum::addPoint(kSpecPoint& s){
  spec->push_back(s);
}

void KSpectrum::addPrecursor(kPrecursor& p, int sz){
  precursor->push_back(p);
  KTopPeps tp;
  tp.singletMax=sz;
  tp.resetSingletList(p.monoMass);
  singlets->push_back(tp);
}

void KSpectrum::clear(){
  spec->clear();
  precursor->clear();
  singlets->clear();
}

void KSpectrum::clearPrecursors(){
  precursor->clear();
  singlets->clear();
}

void KSpectrum::erasePrecursor(int i){
  precursor->erase(precursor->begin()+i);
  singlets->erase(singlets->begin()+i);
}

void KSpectrum::setCharge(int i){
  charge=i;
}

void KSpectrum::setInstrumentPrecursor(bool b){
  instrumentPrecursor=b;
}

void KSpectrum::setMaxIntensity(float f){
  maxIntensity=f;
}

void KSpectrum::setMZ(double d){
  mz=d;
}

void KSpectrum::setRTime(float f){
  rTime=f;
}

void KSpectrum::setScanNumber(int i){
  scanNumber=i;
}

/*============================
  Functions
============================*/
bool KSpectrum::calcEValue(kParams* params, KDecoys& decoys) {
  int i;
  int iLoopCount;
  int iMaxCorr;
  int iStartCorr;
  int iNextCorr;
  double dSlope,dSlopeS;
  double dIntercept,dInterceptS;
  double dRSquare,dRSquareS;
  bool bSkipXL=false;
  bool bSingletFail=false;

  tmpSingCount = histogramSingletCount;
  tmpHistCount = histogramCount;
  if (topHit[0].simpleScore == 0) return true; //no need to do any of this if there are no PSMs...

  /*
  if(histogramSingletCount < DECOY_SIZE) {
    if(!generateSingletDecoys(params, decoys)) {
      bSingletFail=true;
    }
  }

  if(!bSingletFail){
    linearRegression3(dSlopeS, dInterceptS, iMaxCorr, iStartCorr, iNextCorr, dRSquareS);
    tmpSingletIntercept = (float)dInterceptS;  // b
    tmpSingletSlope = (float)dSlopeS;  // m
    tmpSingletIStartCorr = (float)iStartCorr;
    tmpSingletINextCorr = (float)iNextCorr;
    tmpSingletIMaxCorr = (short)iMaxCorr;
    tmpSingletRSquare = dRSquareS;
  } else {
    dInterceptS=0;
    dSlopeS=0;
    dRSquareS=0;
    tmpSingletIntercept = 0;  // b
    tmpSingletSlope = 0;  // m
    tmpSingletIStartCorr = 0;
    tmpSingletINextCorr = 0;
    tmpSingletIMaxCorr = 0;
    tmpSingletRSquare = 0;
  }
  */

  if (histogramCount < DECOY_SIZE) {
    
        //if(!generateXcorrDecoys(params,decoys)) return false;
        //if (!generateXcorrDecoysXL(params, decoys)) return false;
 
    //if(!bSkipXL){ //if singletdecoys were successful, use xl decoy analysis
      //if (!generateXLDecoys(params, decoys)) {
        //if (!generateXcorrDecoys(params, decoys)) return false; 
        //if (!generateXcorrDecoysXL(params, decoys)) return false;
      //}
   // }
    if (!generateXcorrDecoys(params, decoys)) return false;
    //if (!generateXcorrDecoysXL(params, decoys)) return false;
  }

  linearRegression2(dSlope, dIntercept, iMaxCorr, iStartCorr, iNextCorr,dRSquare);
  histoMaxIndex = iMaxCorr;

  //diagnostics - probably temporary
  tmpIntercept = (float)dIntercept;  // b
  tmpSlope = (float)dSlope;  // m
  tmpIStartCorr = (float)iStartCorr;
  tmpINextCorr = (float)iNextCorr;
  tmpIMaxCorr = (short)iMaxCorr;
  tmpRSquare = dRSquare;

  dSlope *= 10.0;
  dSlopeS *=10.0;

  /*
  iLoopCount = 20; //score all e-values among top hits?
  for (i = 0; i<iLoopCount; i++) {
    if (topHit[i].simpleScore == 0) break; //score all e-values among top hits?
    if (dSlope >= 0.0) {
      topHit[i].eVal = 999.0;
      topHit[i].eVal1 = 999.0;
      topHit[i].eVal2 = 999.0;
    } else {
      topHit[i].eVal = pow(10.0, dSlope * topHit[i].simpleScore + dIntercept);
      topHit[i].eVal1 = pow(10.0, dSlope * topHit[i].score1 + dIntercept);
      topHit[i].eVal2 = pow(10.0, dSlope * topHit[i].score2 + dIntercept);
      if (topHit[i].eVal>999.0) topHit[i].eVal=999.0;
      if (topHit[i].eVal1>999.0) topHit[i].eVal1 = 999.0;
      if (topHit[i].eVal2>999.0) topHit[i].eVal2 = 999.0;
    }
  }
  */
  iLoopCount = 20; //score all e-values among top hits?
  double topScore=topHit[0].simpleScore;
  for (i = 0; i<iLoopCount; i++) {
    if (topHit[i].simpleScore == 0) break; //score all e-values among top hits?
    if (dSlope >= 0.0) {
      topHit[i].eVal = 999.0;
    } else {
      topHit[i].eVal = pow(10.0, dSlope * topHit[i].simpleScore + dIntercept);
      if (topHit[i].eVal>999.0) topHit[i].eVal = 999.0;
    }
    //score individual peptides
    if(topHit[i].score2>0){
      if(topHit[i].simpleScore==topScore){ //only do this for top hits right now: it is slow...
        //if (scanNumber == 16660) cout << "E1" << "\t" << topHit[i].score1 << "\t" << topHit[i].mass1 << endl;
        topHit[i].eVal1 = generateSingletDecoys2(params,decoys,topHit[i].score1,topHit[i].mass1,(int)topHit[i].precursor,topHit[i].score2);
        //if(topHit[i].eVal1>999.0) topHit[i].eVal1=999.0;
        //if (scanNumber == 16660) cout << "E2" << "\t" << topHit[i].score2 << "\t" << topHit[i].mass2 << endl;
        topHit[i].eVal2 = generateSingletDecoys2(params, decoys, topHit[i].score2, topHit[i].mass2, (int)topHit[i].precursor,topHit[i].score1);
        //if (topHit[i].eVal2>999.0) topHit[i].eVal2 = 999.0;
      }
      /*
      if(dSlopeS>=0.0){
        topHit[i].eVal1 = 999.0;
        topHit[i].eVal2 = 999.0;
      } else {
        topHit[i].eVal1 = pow(10.0, dSlopeS * topHit[i].score1 + dInterceptS);
        topHit[i].eVal2 = pow(10.0, dSlopeS * topHit[i].score2 + dInterceptS);
        if (topHit[i].eVal1>999.0) topHit[i].eVal1 = 999.0;
        if (topHit[i].eVal2>999.0) topHit[i].eVal2 = 999.0;
      }
      */
    } else {
      if (dSlope >= 0.0) {
        topHit[i].eVal1 = 999.0;
      } else {
        topHit[i].eVal1 = pow(10.0, dSlope * topHit[i].score1 + dIntercept);
        if (topHit[i].eVal1>999.0) topHit[i].eVal1 = 999.0;
      }
      topHit[i].eVal2 = 999.0;
    }
  }
  return true;
}

void KSpectrum::checkScore(kScoreCard& s){
  unsigned int i;
  unsigned int j;

  //edge case for "reversible" cross-links: check if already matches top hit identically
  //note that such duplications still occur below the top score, but shouldn't influence the final result to the user
  int k=0;
  while(k<20 && s.simpleScore==topHit[k].simpleScore){
    if(s.pep1==topHit[k].pep1 && s.pep2==topHit[k].pep2 && s.k1==topHit[k].k1 && s.k2==topHit[k].k2){
      if(s.mods1->size()==topHit[k].mods1->size() && s.mods2->size()==topHit[k].mods2->size()){
        for(i=0;i<s.mods1->size();i++){
          if(s.mods1->at(i).mass!=topHit[k].mods1->at(i).mass || s.mods1->at(i).pos!=topHit[k].mods1->at(i).pos) break;
        }
        for(j=0;j<s.mods2->size();j++){
          if(s.mods2->at(j).mass!=topHit[k].mods2->at(j).mass || s.mods2->at(j).pos!=topHit[k].mods2->at(j).pos) break;
        }
        if(i==s.mods1->size() && j==s.mods2->size()) return;
      }
    }
    k++;
  }

  for(i=0;i<20;i++){
    if(s.simpleScore > topHit[i].simpleScore) {
      for(j=19;j>i;j--) {
        topHit[j]=topHit[j-1];
      }
      topHit[i] = s;
      lowScore=topHit[19].simpleScore;
      return;
    }
  }
}

//This function is now deprecated...
void KSpectrum::checkSingletScore(kSingletScoreCard& s){

  cout << "CHECK SS" << endl;
  printf("CheckSS\n");

  kSingletScoreCard* sc;
  kSingletScoreCard* cur;
  size_t ind;
  
  //If list is empty, add the score card
  if(singletCount==0){
    singletFirst=new kSingletScoreCard(s);
    singletLast=singletFirst;
    singletCount++;

    //add to singlet list
    ind = (int)(s.mass/10);
    if(singletList[ind]==NULL) singletList[ind] = new list<kSingletScoreCard*>;
    singletList[ind]->emplace_back(singletFirst);

    return;
  }  

  //check if we can just add to the end
  if(s.simpleScore<=singletLast->simpleScore){
    //check if we need to store the singlet
    //if(singletCount>=singletMax) return;

    singletLast->next=new kSingletScoreCard(s);
    singletLast->next->prev=singletLast;
    singletLast=singletLast->next;
    singletCount++;

    if(singletCount>singletMax) cout << "Appended max" << endl;

    //add to singlet list
    ind = (int)(s.mass / 10);
    if (singletList[ind] == NULL) singletList[ind] = new list<kSingletScoreCard*>;
    singletList[ind]->emplace_back(singletLast);

    return;
  }
  
  //check if it goes in the front
  if(s.simpleScore>=singletFirst->simpleScore){
    singletFirst->prev=new kSingletScoreCard(s);
    singletFirst->prev->next=singletFirst;
    singletFirst=singletFirst->prev;

    //add to singlet list
    ind = (int)(s.mass / 10);
    if (singletList[ind] == NULL) singletList[ind] = new list<kSingletScoreCard*>;
    singletList[ind]->emplace_back(singletFirst);
    singletCount++;

    if (singletCount>singletMax) {
      cout << "Maxxed out: " << singletCount << "\t" << singletLast->simpleScore << "\t" << singletLast->prev->simpleScore << endl;
    }

    if(singletCount>singletMax){
      int i = singletCount;
      cur = singletLast;
      while (i>singletMax){ //step to singletMax position
        cur = cur->prev;
        i--;
      }
      while (cur->next!=NULL && cur->next->simpleScore==cur->simpleScore){ //step to first instance of score lower than singletMax
        cur=cur->next;
      }

      //delete everything hereafter
      while(cur->next!=NULL){ 
        sc=cur->next;
        cur->next=sc->next;
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
      return;
    }
  }


  //scan to find insertion point
  cout << "WTF" << endl;
  cur = singletFirst->next;
  int i=1;
  while(s.simpleScore < cur->simpleScore){
    i++;
    cur=cur->next;
  }

  sc=new kSingletScoreCard(s);
  sc->prev=cur->prev;
  sc->next=cur;
  cur->prev->next=sc;
  cur->prev=sc;
  if(sc->prev==NULL) singletFirst=sc;

  //add to singlet list
  ind = (int)(s.mass / 10);
  if (singletList[ind] == NULL) singletList[ind] = new list<kSingletScoreCard*>;
  singletList[ind]->emplace_back(sc);
  singletCount++;

  if (singletCount>singletMax) {
    cout << "Middle Maxxed out: " << singletCount << "\t" << singletLast->simpleScore << "\t" << singletLast->prev->simpleScore << endl;
  }

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

//from Comet
// Make synthetic decoy spectra to fill out correlation histogram by going
// through each candidate peptide and rotating spectra in m/z space.
bool KSpectrum::generateSingletDecoys(kParams* params, KDecoys& decoys) {
  int i;
  int n;
  int j;
  int k;
  int maxZ;
  int z;
  int r;
  int key;
  int pos;
  int xlSite;
  int xlLen;
  int badPre=0;
  double dBion;
  double dYion;
  double dXcorr;
  double dFragmentIonMass = 0.0;
  double targetMass;
  double diffMass;

  // DECOY_SIZE is the minimum # of decoys required or else this function isn't
  // called.  So need to generate iLoopMax more xcorr scores for the histogram.
  int iLoopMax = DECOY_SIZE - histogramSingletCount;
  int seed = rand() % 3000;
  int decoyIndex;

  j = 0;
  for (i = 0; i<iLoopMax; i++) { // iterate through required # decoys
    dXcorr = 0.0;
    decoyIndex = (seed + i) % 3000;

    //grab precursor at random
    badPre=0;
    r = rand() % precursor->size();
    while(precursor->at(r).monoMass<params->minPepMass*2+params->xLink->at(0).mass) {
      badPre++;
      if(badPre==50) return false;
      r = rand() % precursor->size();
    }
    maxZ = precursor->at(r).charge;
    if (maxZ>4) maxZ = 4;
    n = (int)(precursor->at(r).monoMass / 2 - params->xLink->at(0).mass-params->minPepMass); //note, might have multiple xlink masses
    if (n == 0) targetMass = precursor->at(r).monoMass / 2;
    else targetMass = precursor->at(r).monoMass / 2 + rand() % n;
    diffMass = precursor->at(r).monoMass-targetMass;

    //find link site - somewhat wasted cycles
    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {
      if (decoys.decoyIons[decoyIndex].pdIonsN[j]>targetMass) break;
    }
    xlSite=rand()%(j-1);
    xlLen=j-1;

    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions
      dBion = decoys.decoyIons[decoyIndex].pdIonsN[j];
      dYion = decoys.decoyIons[decoyIndex].pdIonsC[j];
      if(j>=xlSite) dBion+=diffMass;
      if(j>=xlLen-xlSite) dYion+=diffMass;
      if (dBion>targetMass && dYion>targetMass) break; //stop when fragment ion masses exceed precursor mass

      for (n = 0; n<6; n++) {
        if (!params->ionSeries[n]) continue;
        switch (n) {
        case 0: dFragmentIonMass = dBion - 27.9949141; break;
        case 1: dFragmentIonMass = dBion; break;
        case 2: dFragmentIonMass = dBion + 17.026547; break;
        case 3: dFragmentIonMass = dYion + 25.9792649; break;
        case 4: dFragmentIonMass = dYion; break;
        case 5: dFragmentIonMass = dYion - 16.0187224; break;
        default: break;
        }
        if (dFragmentIonMass>targetMass) continue;

        for (z = 1; z<maxZ; z++) {
          mz = (dFragmentIonMass + (z - 1)*1.007276466) / z;
          mz = params->binSize * (int)(mz*invBinSize + params->binOffset);
          key = (int)mz;
          if (key >= kojakBins) break;
          if (kojakSparseArray[key] == NULL) continue;
          pos = (int)((mz - key)*invBinSize);
          dXcorr += kojakSparseArray[key][pos];
        }
      }
    }

    if (dXcorr <= 0.0) dXcorr = 0.0;
    k = (int)(dXcorr*0.05 + 0.5);  // 0.05=0.005*10; see KAnalysis::kojakScoring
    if (k < 0) k = 0;
    else if (k >= HISTOSZ) k = HISTOSZ - 1;
    histogramSinglet[k]++;
    histogramSingletCount++;
  }

  return true;
}

//from Comet
// Make synthetic decoy spectra to fill out correlation histogram by going
// through each candidate peptide and rotating spectra in m/z space.
double KSpectrum::generateSingletDecoys2(kParams* params, KDecoys& decoys, double xcorr, double mass, int preIndex,double score2) {
  int i;
  int n;
  int j;
  int k;
  int maxZ;
  int z;
  int key;
  int pos;
  int xlSite;
  int xlLen;
  double dBion;
  double dYion;
  double dXcorr;
  double dFragmentIonMass = 0.0;
  double diffMass;
  double preMass;

  int tempHistogram[HISTOSZ];
  for(i=0;i<HISTOSZ;i++) tempHistogram[i]=0;

  int seed = rand() % DECOY_SIZE;
  int decoyIndex;

  //compute modification mass
  maxZ = precursor->at(preIndex).charge;
  if (maxZ>4) maxZ = 4;
  preMass = precursor->at(preIndex).monoMass;
  diffMass=preMass-mass;

  for (i = 0; i<DECOY_SIZE - 1; i++) { // iterate through required # decoys
    dXcorr = 0.0;
    decoyIndex = (seed + i) % DECOY_SIZE;
 
    //find link site - somewhat wasted cycles
    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {
      if (decoys.decoyIons[decoyIndex].pdIonsN[j]>mass) break;
    }
    if(j<1) return 999.0;
    else if(j==1) xlSite=0;
    else xlSite = rand() % (j - 1);
    xlLen = j - 1;

    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions
      dBion = decoys.decoyIons[decoyIndex].pdIonsN[j];
      dYion = decoys.decoyIons[decoyIndex].pdIonsC[j];
      if (j >= xlSite) dBion += diffMass;
      if (j >= xlLen - xlSite) dYion += diffMass;
      if (dBion>preMass && dYion>preMass) break; //stop when fragment ion masses exceed precursor mass

      for (n = 0; n<6; n++) {
        if (!params->ionSeries[n]) continue;
        switch (n) {
        case 0: dFragmentIonMass = dBion - 27.9949141; break;
        case 1: dFragmentIonMass = dBion; break;
        case 2: dFragmentIonMass = dBion + 17.026547; break;
        case 3: dFragmentIonMass = dYion + 25.9792649; break;
        case 4: dFragmentIonMass = dYion; break;
        case 5: dFragmentIonMass = dYion - 16.0187224; break;
        default: break;
        }
        if (dFragmentIonMass>preMass) continue;

        for (z = 1; z<maxZ; z++) {
          mz = (dFragmentIonMass + (z - 1)*1.007276466) / z;
          mz = params->binSize * (int)(mz*invBinSize + params->binOffset);
          key = (int)mz;
          if (key >= kojakBins) break;
          if (kojakSparseArray[key] == NULL) continue;
          pos = (int)((mz - key)*invBinSize);
          dXcorr += kojakSparseArray[key][pos];
        }
      }
    }

    if (dXcorr <= 0.0) dXcorr = 0.0;
    k = (int)(dXcorr*0.05+score2*10 + 0.5);  // 0.05=0.005*10; see KAnalysis::kojakScoring
    if (k < 0) k = 0;
    else if (k >= HISTOSZ) k = HISTOSZ - 1;
    tempHistogram[k]++;
  }

  //add this peptide's xcorr? is this necessary?
  k = (int)((xcorr+score2)*10 + 0.5);
  if (k < 0) k = 0;
  else if (k >= HISTOSZ) k = HISTOSZ - 1;
  tempHistogram[k]++;

  //Do linear regression and compute e-value
  double dSlope,dIntercept,dRSquare;
  int iMaxCorr,iStartCorr,iNextCorr;
  linearRegression4(tempHistogram, dSlope, dIntercept, iMaxCorr, iStartCorr, iNextCorr, dRSquare);
  double eVal = pow(10.0, dSlope * 10 * (xcorr+score2) + dIntercept);
  if(eVal>999.0) eVal=999.0;

  /*
  if(scanNumber==16660){
    cout <<"\n16660" << endl;
    cout << dSlope << "\t" << dIntercept << "\t" << dRSquare << endl;
    for(i=0;i<HISTOSZ;i++){
      cout << i << "\t" << tempHistogram[i] << endl;
    }
    cout << xcorr << "\t" << xcorr+score2 << "\t" << eVal << endl;
  }
  if(scanNumber>17000) exit(1);
  */

  return eVal;
}

bool KSpectrum::generateXLDecoys(kParams* params, KDecoys& decoys) {
  int i;
  int n;
  int j;
  int k;
  int maxZ;
  int z;
  int r,r2;
  int key;
  int pos;
  int xlSite;
  int xlLen;
  int badPre;
  double dSum;
  double d;
  double dBion;
  double dYion;
  double dXcorr;
  double dFragmentIonMass = 0.0;
  double targetMass;
  double diffMass;

  // DECOY_SIZE is the minimum # of decoys required or else this function isn't
  // called.  So need to generate iLoopMax more xcorr scores for the histogram.
  int iLoopMax = DECOY_SIZE - histogramCount;
  int seed = rand() % 3000;
  int decoyIndex;

  j = 0;
  for (i = 0; i<iLoopMax; i++) { // iterate through required # decoys
    dXcorr = 0.0;
    decoyIndex = (seed + i) % 3000;

    //grab precursor at random and xl at random
    badPre=0;
    r = rand() % precursor->size();
    r2 = rand() % params->xLink->size();
    while (precursor->at(r).monoMass<params->minPepMass * 2 + params->xLink->at(r2).mass) {
      badPre++;
      r = rand() % precursor->size();
      r2 = rand() % params->xLink->size();
      if(badPre==50) return false;
    }
    maxZ = precursor->at(r).charge;
    if (maxZ>4) maxZ = 4;
    targetMass = (precursor->at(r).monoMass-params->xLink->at(r2).mass-2*params->minPepMass)*((double)rand()/RAND_MAX)+params->minPepMass;
    diffMass = precursor->at(r).monoMass - targetMass;

    //find link site - somewhat wasted cycles
    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {
      if (decoys.decoyIons[decoyIndex].pdIonsN[j]>targetMass) break;
    }
    xlSite = rand() % (j - 1);
    xlLen = j - 1;

    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions
      dBion = decoys.decoyIons[decoyIndex].pdIonsN[j];
      dYion = decoys.decoyIons[decoyIndex].pdIonsC[j];
      if (j >= xlSite) dBion += diffMass;
      if (j >= xlLen - xlSite) dYion += diffMass;
      if (dBion>targetMass && dYion>targetMass) break; //stop when fragment ion masses exceed precursor mass

      for (n = 0; n<6; n++) {
        if (!params->ionSeries[n]) continue;
        switch (n) {
        case 0: dFragmentIonMass = dBion - 27.9949141; break;
        case 1: dFragmentIonMass = dBion; break;
        case 2: dFragmentIonMass = dBion + 17.026547; break;
        case 3: dFragmentIonMass = dYion + 25.9792649; break;
        case 4: dFragmentIonMass = dYion; break;
        case 5: dFragmentIonMass = dYion - 16.0187224; break;
        default: break;
        }
        if (dFragmentIonMass>targetMass) continue;

        for (z = 1; z<maxZ; z++) {
          mz = (dFragmentIonMass + (z - 1)*1.007276466) / z;
          mz = params->binSize * (int)(mz*invBinSize + params->binOffset);
          key = (int)mz;
          if (key >= kojakBins) break;
          if (kojakSparseArray[key] == NULL) continue;
          pos = (int)((mz - key)*invBinSize);
          dXcorr += kojakSparseArray[key][pos];
        }
      }
    }

    if (dXcorr <= 0.0) dXcorr = 0.0;
    dXcorr*=0.05;

    //find complement peptide randomly from proportion of singlet scores
    d=(double)rand()/RAND_MAX;
    dSum=0;
    for(n=0;n<HISTOSZ;n++){
      dSum+=(double)histogramSinglet[n]/histogramSingletCount;
      if(d<=dSum){
        dXcorr+=n;
        break;
      }
    }
    k = (int)(dXcorr + 0.5);
    if (k < 0) k = 0;
    else if (k >= HISTOSZ) k = HISTOSZ - 1;
    histogram[k]++;
    histogramCount++;
  }

  return true;
}

//from Comet
// Make synthetic decoy spectra to fill out correlation histogram by going
// through each candidate peptide and rotating spectra in m/z space.
bool KSpectrum::generateXcorrDecoys(kParams* params, KDecoys& decoys) {
  int i;
  int n;
  int j;
  int k;
  int maxZ;
  int z;
  int r;
  int key;
  int pos;
  double dBion;
  double dYion;
  double dXcorr;
  double dFragmentIonMass = 0.0;
  int myCount=0;

  tmpSingCount = histogramSingletCount;
  tmpHistCount = histogramCount;

  // DECOY_SIZE is the minimum # of decoys required or else this function isn't
  // called.  So need to generate iLoopMax more xcorr scores for the histogram.
  int iLoopMax = DECOY_SIZE - histogramCount;
  //int seed = rand() % DECOY_SIZE;
  int seed = (scanNumber*histogramCount);
  if(seed<0) seed=-seed;
  seed = seed % DECOY_SIZE; //don't always start at the top, but not random either; remains reproducible across threads
  //if(scanNumber==12896) cout << "SEED: " << seed << endl;
  int decoyIndex;

  size_t maxPre = precursor->size();
  r=0;
  for (i = 0; i<iLoopMax; i++) { // iterate through required # decoys
    dXcorr = 0.0;
    decoyIndex = (seed + i) % DECOY_SIZE;

    //grab precursor at random
    //r = rand() % precursor->size();
    //iterate over precursors
    r++;
    if(r>=(int)precursor->size()) r=0;
    maxZ = precursor->at(r).charge;
    if (maxZ>4) maxZ = 4;

    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions
      dBion = decoys.decoyIons[decoyIndex].pdIonsN[j];
      dYion = decoys.decoyIons[decoyIndex].pdIonsC[j];
      if (dBion>precursor->at(r).monoMass && dYion>precursor->at(r).monoMass) break; //stop when fragment ion masses exceed precursor mass

      for (n = 0; n<6; n++) {
        if (!params->ionSeries[n]) continue;
        switch (n) {
        case 0: dFragmentIonMass = dBion - 27.9949141; break;
        case 1: dFragmentIonMass = dBion; break;
        case 2: dFragmentIonMass = dBion + 17.026547; break;
        case 3: dFragmentIonMass = dYion + 25.9792649; break;
        case 4: dFragmentIonMass = dYion; break;
        case 5: dFragmentIonMass = dYion - 16.0187224; break;
        default: break;
        }
        if(dFragmentIonMass<0) {
          cout << "\ndFragIonMass = " << dFragmentIonMass << endl;
          cout << "j = " << j << endl;
          cout << "decoyIndex = " << decoyIndex << endl;
          cout << "dBion = " << dBion << endl;
        }
        if (dFragmentIonMass>precursor->at(r).monoMass) continue;

        for (z = 1; z<maxZ; z++) {
          mz = (dFragmentIonMass + (z - 1)*1.007276466) / z;
          mz = params->binSize * (int)(mz*invBinSize + params->binOffset);
          key = (int)mz;
          if (key >= kojakBins) break;
          if (kojakSparseArray[key] == NULL) continue;
          pos = (int)((mz - key)*invBinSize);
          dXcorr += kojakSparseArray[key][pos];
        }
      }
    }

    if (dXcorr <= 0.0) dXcorr = 0.0;
    k = (int)(dXcorr*0.05 + 0.5);  // 0.05=0.005*10; see KAnalysis::kojakScoring
    if (k < 0) k = 0;
    else if (k >= HISTOSZ) k = HISTOSZ - 1;
    histogram[k]++;
    histogramCount++;
    myCount++;
    //if (scanNumber == 12896) cout << "myCount: " << myCount-1 << "\t" << k << "\t" << dXcorr*0.05 << endl;
  }
  return true;
}

bool KSpectrum::generateXcorrDecoysXL(kParams* params, KDecoys& decoys) {
  int i;
  int n;
  int j;
  int k;
  int maxZ;
  int z;
  int r;
  int key;
  int pos;
  double dBion;
  double dYion;
  double dXcorr;
  double dFragmentIonMass = 0.0;
  int myCount = 0;

  double p1Mass;
  double p2Mass;
  double totalMass;
  double targetMass;
  double pBindA;
  double pBindB;


  //tmpSingCount = histogramSingletCount;
  //tmpHistCount = histogramCount;

  // DECOY_SIZE is the minimum # of decoys required or else this function isn't
  // called.  So need to generate iLoopMax more xcorr scores for the histogram.
  int iLoopMax = DECOY_SIZE - histogramCount;
  int seed = (scanNumber*histogramCount);
  if (seed<0) seed = -seed;
  seed = seed % DECOY_SIZE; //don't always start at the top, but not random either; remains reproducible across threads

  int seed2 = (histogramCount*histogramCount);
  if (seed2<0) seed2 = -seed2;
  seed2 = seed2 % DECOY_SIZE;

  int decoyIndex;

  size_t maxPre = precursor->size();
  r = 0;
  for (i = 0; i<iLoopMax; i++) { // iterate through required # decoys
    dXcorr = 0.0;
    decoyIndex = (seed + i) % DECOY_SIZE;
    //cout << "Decoy Index 1: " << decoyIndex << endl;

    //iterate over precursors -- does this matter???
    r++;
    if (r >= (int)precursor->size()) r = 0;
    maxZ = precursor->at(r).charge;
    if (maxZ>4) maxZ = 4;

    //determine size of each peptide
    totalMass = precursor->at(r).monoMass - params->xLink->at(0).mass; //note only using first linker mass here. Maybe randomize?
    targetMass = totalMass - 2*params->minPepMass;
    p1Mass = targetMass*((double)rand() / RAND_MAX)+params->minPepMass;
    p2Mass = totalMass-p1Mass;

    //cout << precursor->at(r).monoMass << " = " << params->xLink->at(0).mass << " + " << p1Mass << " + " << p2Mass << endl;

    //score p1
    //determine attachment point (by mass)
    pBindA = (p1Mass - 100)*((double)rand() / RAND_MAX)+50;
    pBindB = p1Mass-pBindA;
    targetMass = precursor->at(r).monoMass - p1Mass;
    //cout << "p1: " << pBindA << "\t" << pBindB << "\t" << targetMass << endl;
    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions
      dBion = decoys.decoyIons[decoyIndex].pdIonsN[j];
      dYion = decoys.decoyIons[decoyIndex].pdIonsC[j];
      if(dBion>pBindA) dBion+=targetMass;
      if (dYion>pBindB) dYion+=targetMass;
      if (dBion>precursor->at(r).monoMass && dYion>precursor->at(r).monoMass) break; //stop when fragment ion masses exceed precursor mass

      for (n = 0; n<6; n++) {
        if (!params->ionSeries[n]) continue;
        switch (n) {
        case 0: dFragmentIonMass = dBion - 27.9949141; break;
        case 1: dFragmentIonMass = dBion; break;
        case 2: dFragmentIonMass = dBion + 17.026547; break;
        case 3: dFragmentIonMass = dYion + 25.9792649; break;
        case 4: dFragmentIonMass = dYion; break;
        case 5: dFragmentIonMass = dYion - 16.0187224; break;
        default: break;
        }
        if (dFragmentIonMass>precursor->at(r).monoMass) continue;

        for (z = 1; z<maxZ; z++) {
          mz = (dFragmentIonMass + (z - 1)*1.007276466) / z;
          mz = params->binSize * (int)(mz*invBinSize + params->binOffset);
          key = (int)mz;
          if (key >= kojakBins) break;
          if (kojakSparseArray[key] == NULL) continue;
          pos = (int)((mz - key)*invBinSize);
          dXcorr += kojakSparseArray[key][pos];
        }
      }
    }
    //cout << "Did p1" << endl;

    //get a different decoy peptide
    decoyIndex = (seed2 + i) % DECOY_SIZE;
    //cout << "Decoy Index 2: " << decoyIndex << endl;

    //score p2
    //determine attachment point (by mass)
    pBindA = (p2Mass - 100)*((double)rand() / RAND_MAX) + 50;
    pBindB = p2Mass - pBindA;
    targetMass = precursor->at(r).monoMass - p2Mass;
    //cout << "p2: " << pBindA << "\t" << pBindB << "\t" << targetMass << endl;
    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions
      dBion = decoys.decoyIons[decoyIndex].pdIonsN[j];
      dYion = decoys.decoyIons[decoyIndex].pdIonsC[j];
      if (dBion>pBindA) dBion += targetMass;
      if (dYion>pBindB) dYion += targetMass;
      if (dBion>precursor->at(r).monoMass && dYion>precursor->at(r).monoMass) break; //stop when fragment ion masses exceed precursor mass

      for (n = 0; n<6; n++) {
        if (!params->ionSeries[n]) continue;
        switch (n) {
        case 0: dFragmentIonMass = dBion - 27.9949141; break;
        case 1: dFragmentIonMass = dBion; break;
        case 2: dFragmentIonMass = dBion + 17.026547; break;
        case 3: dFragmentIonMass = dYion + 25.9792649; break;
        case 4: dFragmentIonMass = dYion; break;
        case 5: dFragmentIonMass = dYion - 16.0187224; break;
        default: break;
        }
        if (dFragmentIonMass>precursor->at(r).monoMass) continue;

        for (z = 1; z<maxZ; z++) {
          mz = (dFragmentIonMass + (z - 1)*1.007276466) / z;
          mz = params->binSize * (int)(mz*invBinSize + params->binOffset);
          key = (int)mz;
          if (key >= kojakBins) break;
          if (kojakSparseArray[key] == NULL) continue;
          pos = (int)((mz - key)*invBinSize);
          dXcorr += kojakSparseArray[key][pos];
        }
      }
    }
    //cout << "Did p2" << endl;

    if (dXcorr <= 0.0) dXcorr = 0.0;
    k = (int)(dXcorr*0.05 + 0.5);  // 0.05=0.005*10; see KAnalysis::kojakScoring
    if (k < 0) k = 0;
    else if (k >= HISTOSZ) k = HISTOSZ - 1;
    histogram[k]++;
    histogramCount++;
    myCount++;
    //if (scanNumber == 12896) cout << "myCount: " << myCount-1 << "\t" << k << "\t" << dXcorr*0.05 << endl;
  }
  return true;
}

//from Comet
void KSpectrum::linearRegression(double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared) {
  double Sx, Sxy;      // Sum of square distances.
  double Mx, My;       // means
  double dx, dy;
  double b, a;
  double SumX, SumY;   // Sum of X and Y values to calculate mean.
  double SST,SSR;
  double rsq;
  double bestRSQ;
  double bestSlope;
  double bestInt;

  double dCummulative[HISTOSZ];  // Cummulative frequency at each xcorr value.

  int i;
  int bestNC;
  int bestStart;
  int iNextCorr;    // 2nd best xcorr index
  int iMaxCorr = 0;   // max xcorr index
  int iStartCorr;
  int iNumPoints;

  // Find maximum correlation score index.
  for (i = HISTOSZ - 2; i >= 0; i--) {
    if (histogram[i] > 0)  break;
  }
  iMaxCorr = i;
  if(iMaxCorr<2) {
    slope=0;
    intercept=0;
    iMaxXcorr=0;
    iStartXcorr=0;
    iNextXcorr=0;
    rSquared=0;
    return;
  }

  //get last datapoint before first zero value; I think this is best indicator of boundaries
  iNextCorr = 0;
  for (i = 0; i<iMaxCorr; i++){
    if (histogram[i]==0) break;
    iNextCorr = i;
  }

  //get the next non-zero value below iMaxCorr
  /* MH: version 2,replaced original.
  iNextCorr=0;
  for(i=iMaxCorr-1;i>=0;i--){
    if(histogram[i]>0){
      iNextCorr=i;
      break;
    }
  }
  */
  /* original from comet
  //MH: this also causes problems if there are early gaps in the scores.
  iNextCorr = 0;
  for (i = 0; i<iMaxCorr; i++)  {
    if (histogram[i] == 0)  {
      // register iNextCorr if there's a histo value of 0 consecutively
      if (histogram[i + 1] == 0 || i + 1 == iMaxCorr) {
        if (i>0) iNextCorr = i - 1;
        break;
      }
    }
  }
  */

  /* MH: This does nothing!
  if (i == iMaxCorr) {
    cout << scanNumber << "\tNOT POSSIBLE: " << i << "\t" << iMaxCorr;
    iNextCorr = iMaxCorr;
    if (iMaxCorr>12) iNextCorr = iMaxCorr - 2;
  }
  */

  // Create cummulative distribution function from iNextCorr down, skipping the outliers.
  //MH: its the skipping the outliers that is causing problems, I think...
  /*
  dCummulative[iNextCorr] = histogram[iNextCorr];
  for (i = iNextCorr - 1; i >= 0; i--) {
    dCummulative[i] = dCummulative[i + 1] + histogram[i];
    //if (histogram[i + 1] == 0) dCummulative[i + 1] = 0.0;
  }
  */
  //More aggressive version summing everything below the max
  dCummulative[iMaxCorr-1] = histogram[iMaxCorr-1];
  for (i = iMaxCorr - 2; i >= 0; i--) {
    dCummulative[i] = dCummulative[i + 1] + histogram[i];
    //if (histogram[i + 1] == 0) dCummulative[i + 1] = 0.0;
  }

  // log10
  for (i = iMaxCorr-1; i >= 0; i--)  {
    histogram[i] = (int)dCummulative[i];  // First store cummulative in histogram. MH:...and stomp all over the original...hard to troubleshoot later
    dCummulative[i] = log10(dCummulative[i]);
  }

  // log10
  /*MH old version
  for (i = iNextCorr; i >= 0; i--)  {
    histogram[i] = (int)dCummulative[i];  // First store cummulative in histogram. MH:...and stomp all over the original...hard to troubleshoot later
    dCummulative[i] = log10(dCummulative[i]);
  }
  */

  iStartCorr = 0;
  //MH: futzing around here
  //if (iNextCorr >= 30) iStartCorr = (int)(iNextCorr - iNextCorr*0.25);
  //else if (iNextCorr >= 15) iStartCorr = (int)(iNextCorr - iNextCorr*0.5);

  bool bRight=false;
  //if((histogram[iStartCorr]-histogram[iStartCorr+1]) < (histogram[iNextCorr-1]-histogram[iNextCorr])) bRight=false; //get rid of outlier with least slope?
  //else bRight=true;
  bestRSQ=0;
  bestNC=0;
  bestSlope=0;
  bestInt=0;
  rsq = Mx = My = a = b = 0.0;
  while (iStartCorr < iNextCorr-1) {
    Sx = Sxy = SumX = SumY = 0.0;
    iNumPoints = 0;

    // Calculate means.
    for (i = iStartCorr; i <= iNextCorr; i++) {
      if (histogram[i] > 0) {
        SumY += dCummulative[i];
        SumX += i;
        iNumPoints++;
      }
    }
    if (iNumPoints > 0) {
      Mx = SumX / iNumPoints;
      My = SumY / iNumPoints;
    } else {
      Mx = My = 0.0;
    }

    // Calculate sum of squares.
    SST=0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      //if (dCummulative[i] > 0) {  //MH: this will allow skipping of data points!!!
        dx = i - Mx;
        dy = dCummulative[i] - My;
        Sx += dx*dx;
        Sxy += dx*dy;
        SST += dy*dy;
      //}
    }
    b = Sxy / Sx;
    //if (Sx > 0) b = Sxy / Sx;   // slope //MH: don't see how Sx could be zero or less given the rules above...
    //else b = 0;
    a = My - b*Mx;  // y-intercept

    //MH: compute R2
    SSR=0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      //if (dCummulative[i] > 0) {
        dy=dCummulative[i]-(b*i+a);
        SSR+=(dy*dy);
      //}
    }
    rsq = 1-SSR/SST;

    if(rsq>bestRSQ){
      bestRSQ=rsq;
      bestNC=iNextCorr;
      bestSlope=b;
      bestInt=a;
      bestStart=iStartCorr;
      if(rsq>0.95) break; //stop when fit is good enough? Don't trim too many outliers?
      if(dCummulative[iStartCorr]-dCummulative[iStartCorr+1] > dCummulative[iNextCorr-1]-dCummulative[iNextCorr]){
      //if(bRight){
        iNextCorr--;
        continue;
      } else {
        //iNextCorr=bestNC;
        iStartCorr++;
        continue;
      }
    } else {
      break;

      if(bRight) {
        if(iNextCorr>10 && dCummulative[iNextCorr]<10){
          iNextCorr--;
          continue;
        }
        bRight=false;
        iNextCorr=bestNC;
        iStartCorr++;
        continue;
      } else {
        iStartCorr--;
        break;
      }
    }
  }

  slope = bestSlope;
  intercept = bestInt;
  iMaxXcorr = iMaxCorr;
  iStartXcorr = bestStart;
  iNextXcorr = bestNC;
  rSquared = bestRSQ;
}

//from Comet
void KSpectrum::linearRegression2(double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared) {
  double Sx, Sxy;      // Sum of square distances.
  double Mx, My;       // means
  double dx, dy;
  double b, a;
  double SumX, SumY;   // Sum of X and Y values to calculate mean.
  double SST, SSR;
  double rsq;
  double bestRSQ;
  double bestSlope;
  double bestInt;

  double dCummulative[HISTOSZ];  // Cummulative frequency at each xcorr value.

  int i;
  int bestNC;
  int bestStart;
  int iNextCorr;    // 2nd best xcorr index
  int iMaxCorr = 0;   // max xcorr index
  int iStartCorr;
  int iNumPoints;

  // Find maximum correlation score index.
  for (i = HISTOSZ - 2; i >= 0; i--) {
    if (histogram[i] > 0)  break;
  }
  iMaxCorr = i;

  //bail now if there is no width to the distribution
  if (iMaxCorr<3) {
    slope = 0;
    intercept = 0;
    iMaxXcorr = 0;
    iStartXcorr = 0;
    iNextXcorr = 0;
    rSquared = 0;
    return;
  }

  //More aggressive version summing everything below the max
  dCummulative[iMaxCorr - 1] = histogram[iMaxCorr - 1];
  for (i = iMaxCorr - 2; i >= 0; i--) {
    dCummulative[i] = dCummulative[i + 1] + histogram[i];
  }

  //get middle-ish datapoint as seed. Using count/10.
  for (i = 0; i<iMaxCorr; i++){
    if (dCummulative[i] < histogramCount/10) break;  
  }
  if(i>=(iMaxCorr-1)) iNextCorr=iMaxCorr-2;
  else iNextCorr = i;

  // log10...and stomp all over the original...hard to troubleshoot later
  for (i = iMaxCorr - 1; i >= 0; i--)  {
    histogram[i] = (int)dCummulative[i];
    dCummulative[i] = log10(dCummulative[i]);
  }

  iStartCorr = iNextCorr-1;
  iNextCorr++;

  bool bRight = false; // which direction to add datapoint from
  bestRSQ = 0;
  bestNC = 0;
  bestSlope = 0;
  bestInt = 0;
  rsq = Mx = My = a = b = 0.0;
  while (true) {
    Sx = Sxy = SumX = SumY = 0.0;
    iNumPoints = 0;

    // Calculate means.
    for (i = iStartCorr; i <= iNextCorr; i++) {
      if (histogram[i] > 0) {
        SumY += dCummulative[i];
        SumX += i;
        iNumPoints++;
      }
    }
    if (iNumPoints > 0) {
      Mx = SumX / iNumPoints;
      My = SumY / iNumPoints;
    } else {
      Mx = My = 0.0;
    }

    // Calculate sum of squares.
    SST = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dx = i - Mx;
      dy = dCummulative[i] - My;
      Sx += dx*dx;
      Sxy += dx*dy;
      SST += dy*dy;
    }
    b = Sxy / Sx;
    a = My - b*Mx;  // y-intercept

    //MH: compute R2
    SSR = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dy = dCummulative[i] - (b*i + a);
      SSR += (dy*dy);
    }
    rsq = 1 - SSR / SST;

    if (rsq>0.95 || rsq>bestRSQ){
      if(rsq>bestRSQ || iNextCorr-iStartCorr+1<8){ //keep better RSQ only if more than 8 datapoints, otherwise keep every RSQ below 8 datapoints
        bestRSQ = rsq;
        bestNC = iNextCorr;
        bestSlope = b;
        bestInt = a;
        bestStart = iStartCorr;
      }
      if(bRight){
        if (iNextCorr<(iMaxCorr - 1)) iNextCorr++;
        else if (iStartCorr>0) iStartCorr--;
        else break;
      } else {
        if(iStartCorr>0) iStartCorr--;
        else if(iNextCorr<(iMaxCorr-1)) iNextCorr++;
        else break;
      }
      bRight=!bRight;
    } else {
      break;
    }
  }

  slope = bestSlope;
  intercept = bestInt;
  iMaxXcorr = iMaxCorr;
  iStartXcorr = bestStart;
  iNextXcorr = bestNC;
  rSquared = bestRSQ;
}

//from Comet. THIS IS TEMPORARY! for testing singlets only...
void KSpectrum::linearRegression3(double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared) {
  double Sx, Sxy;      // Sum of square distances.
  double Mx, My;       // means
  double dx, dy;
  double b, a;
  double SumX, SumY;   // Sum of X and Y values to calculate mean.
  double SST, SSR;
  double rsq;
  double bestRSQ;
  double bestSlope;
  double bestInt;

  double dCummulative[HISTOSZ];  // Cummulative frequency at each xcorr value.

  int i;
  int bestNC;
  int bestStart;
  int iNextCorr;    // 2nd best xcorr index
  int iMaxCorr = 0;   // max xcorr index
  int iStartCorr;
  int iNumPoints;

  // Find maximum correlation score index.
  for (i = HISTOSZ - 2; i >= 0; i--) {
    if (histogramSinglet[i] > 0)  break;
  }
  iMaxCorr = i;

  //bail now if there is no width to the distribution
  if (iMaxCorr<3) {
    slope = 0;
    intercept = 0;
    iMaxXcorr = 0;
    iStartXcorr = 0;
    iNextXcorr = 0;
    rSquared = 0;
    return;
  }

  //More aggressive version summing everything below the max
  dCummulative[iMaxCorr - 1] = histogramSinglet[iMaxCorr - 1];
  for (i = iMaxCorr - 2; i >= 0; i--) {
    dCummulative[i] = dCummulative[i + 1] + histogramSinglet[i];
  }

  //get middle-ish datapoint as seed. Using count/10.
  for (i = 0; i<iMaxCorr; i++){
    if (dCummulative[i] < histogramSingletCount / 10) break;
  }
  if (i >= (iMaxCorr - 1)) iNextCorr = iMaxCorr - 2;
  else iNextCorr = i;

  // log10...and stomp all over the original...hard to troubleshoot later
  for (i = iMaxCorr - 1; i >= 0; i--)  {
    histogramSinglet[i] = (int)dCummulative[i];
    dCummulative[i] = log10(dCummulative[i]);
  }

  iStartCorr = iNextCorr - 1;
  iNextCorr++;

  bool bRight = false; // which direction to add datapoint from
  bestRSQ = 0;
  bestNC = 0;
  bestSlope = 0;
  bestInt = 0;
  rsq = Mx = My = a = b = 0.0;
  while (true) {
    Sx = Sxy = SumX = SumY = 0.0;
    iNumPoints = 0;

    // Calculate means.
    for (i = iStartCorr; i <= iNextCorr; i++) {
      if (histogramSinglet[i] > 0) {
        SumY += dCummulative[i];
        SumX += i;
        iNumPoints++;
      }
    }
    if (iNumPoints > 0) {
      Mx = SumX / iNumPoints;
      My = SumY / iNumPoints;
    } else {
      Mx = My = 0.0;
    }

    // Calculate sum of squares.
    SST = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dx = i - Mx;
      dy = dCummulative[i] - My;
      Sx += dx*dx;
      Sxy += dx*dy;
      SST += dy*dy;
    }
    b = Sxy / Sx;
    a = My - b*Mx;  // y-intercept

    //MH: compute R2
    SSR = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dy = dCummulative[i] - (b*i + a);
      SSR += (dy*dy);
    }
    rsq = 1 - SSR / SST;

    if (rsq>0.95 || rsq>bestRSQ){
      if (rsq>bestRSQ || iNextCorr - iStartCorr + 1<8){ //keep better RSQ only if more than 8 datapoints, otherwise keep every RSQ below 8 datapoints
        bestRSQ = rsq;
        bestNC = iNextCorr;
        bestSlope = b;
        bestInt = a;
        bestStart = iStartCorr;
      }
      if (bRight){
        if (iNextCorr<(iMaxCorr - 1)) iNextCorr++;
        else if (iStartCorr>0) iStartCorr--;
        else break;
      } else {
        if (iStartCorr>0) iStartCorr--;
        else if (iNextCorr<(iMaxCorr - 1)) iNextCorr++;
        else break;
      }
      bRight = !bRight;
    } else {
      break;
    }
  }

  slope = bestSlope;
  intercept = bestInt;
  iMaxXcorr = iMaxCorr;
  iStartXcorr = bestStart;
  iNextXcorr = bestNC;
  rSquared = bestRSQ;
}

//from Comet. THIS IS TEMPORARY! for testing singlets only...
void KSpectrum::linearRegression4(int* histo, double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared) {
  double Sx, Sxy;      // Sum of square distances.
  double Mx, My;       // means
  double dx, dy;
  double b, a;
  double SumX, SumY;   // Sum of X and Y values to calculate mean.
  double SST, SSR;
  double rsq;
  double bestRSQ;
  double bestSlope;
  double bestInt;

  double dCummulative[HISTOSZ];  // Cummulative frequency at each xcorr value.

  int i;
  int bestNC;
  int bestStart;
  int iNextCorr;    // 2nd best xcorr index
  int iMaxCorr = 0;   // max xcorr index
  int iStartCorr;
  int iNumPoints;

  // Find maximum correlation score index.
  for (i = HISTOSZ - 2; i >= 0; i--) {
    if (histo[i] > 0)  break;
  }
  iMaxCorr = i;

  //bail now if there is no width to the distribution
  if (iMaxCorr<3) {
    slope = 0;
    intercept = 0;
    iMaxXcorr = 0;
    iStartXcorr = 0;
    iNextXcorr = 0;
    rSquared = 0;
    return;
  }

  //More aggressive version summing everything below the max
  dCummulative[iMaxCorr - 1] = histo[iMaxCorr - 1];
  for (i = iMaxCorr - 2; i >= 0; i--) {
    dCummulative[i] = dCummulative[i + 1] + histo[i];
  }

  //get middle-ish datapoint as seed. Using count/10.
  for (i = 0; i<iMaxCorr; i++){
    if (dCummulative[i] < DECOY_SIZE / 10) break;
  }
  if (i >= (iMaxCorr - 1)) iNextCorr = iMaxCorr - 2;
  else iNextCorr = i;

  // log10...and stomp all over the original...hard to troubleshoot later
  for (i = iMaxCorr - 1; i >= 0; i--)  {
    histo[i] = (int)dCummulative[i];
    dCummulative[i] = log10(dCummulative[i]);
  }

  iStartCorr = iNextCorr - 1;
  iNextCorr++;

  bool bRight = false; // which direction to add datapoint from
  bestRSQ = 0;
  bestNC = 0;
  bestSlope = 0;
  bestInt = 0;
  rsq = Mx = My = a = b = 0.0;
  while (true) {
    Sx = Sxy = SumX = SumY = 0.0;
    iNumPoints = 0;

    // Calculate means.
    for (i = iStartCorr; i <= iNextCorr; i++) {
      if (histo[i] > 0) {
        SumY += dCummulative[i];
        SumX += i;
        iNumPoints++;
      }
    }
    if (iNumPoints > 0) {
      Mx = SumX / iNumPoints;
      My = SumY / iNumPoints;
    } else {
      Mx = My = 0.0;
    }

    // Calculate sum of squares.
    SST = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dx = i - Mx;
      dy = dCummulative[i] - My;
      Sx += dx*dx;
      Sxy += dx*dy;
      SST += dy*dy;
    }
    b = Sxy / Sx;
    a = My - b*Mx;  // y-intercept

    //MH: compute R2
    SSR = 0;
    for (i = iStartCorr; i <= iNextCorr; i++) {
      dy = dCummulative[i] - (b*i + a);
      SSR += (dy*dy);
    }
    rsq = 1 - SSR / SST;

    if (rsq>0.95 || rsq>bestRSQ){
      if (rsq>bestRSQ || iNextCorr - iStartCorr + 1<8){ //keep better RSQ only if more than 8 datapoints, otherwise keep every RSQ below 8 datapoints
        bestRSQ = rsq;
        bestNC = iNextCorr;
        bestSlope = b;
        bestInt = a;
        bestStart = iStartCorr;
      }
      if (bRight){
        if (iNextCorr<(iMaxCorr - 1)) iNextCorr++;
        else if (iStartCorr>0) iStartCorr--;
        else break;
      } else {
        if (iStartCorr>0) iStartCorr--;
        else if (iNextCorr<(iMaxCorr - 1)) iNextCorr++;
        else break;
      }
      bRight = !bRight;
    } else {
      break;
    }
  }

  slope = bestSlope;
  intercept = bestInt;
  iMaxXcorr = iMaxCorr;
  iStartXcorr = bestStart;
  iNextXcorr = bestNC;
  rSquared = bestRSQ;
}

void KSpectrum::resetSingletList(){
  size_t j;
  double max;
  if (singletList != NULL){
    for (j = 0; j<singletBins; j++){
      if (singletList[j] != NULL) delete singletList[j];
    }
    delete[] singletList;
  }
  max=precursor->at(0).monoMass;
  for(j=1;j<precursor->size();j++){
    if (precursor->at(j).monoMass>max) max = precursor->at(j).monoMass;
  }
  singletBins=(int)(max/10+1);
  singletList = new list<kSingletScoreCard*>*[singletBins];
  for (j = 0; j<singletBins; j++) singletList[j] = NULL;
}

void KSpectrum::sortMZ(){
  qsort(&spec->at(0),spec->size(),sizeof(kSpecPoint),compareMZ);
}

void KSpectrum::xCorrScore(bool b){
  if(b) CometXCorr();
  else  kojakXCorr();
}


/*============================
  Private Functions
============================*/
void KSpectrum::CometXCorr(){
  int i;
  int j;
  int iTmp;
  double dTmp;
  double dSum;
  double *pdTempRawData;
  double *pdTmpFastXcorrData;
  float  *pfFastXcorrData;
  kPreprocessStruct pPre;

  pPre.iHighestIon = 0;
  pPre.dHighestIntensity = 0;

  BinIons(&pPre);

  pdTempRawData = (double *)calloc((size_t)xCorrArraySize, (size_t)sizeof(double));
  if (pdTempRawData == NULL) {
    fprintf(stderr, " Error - calloc(pdTempRawData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  pdTmpFastXcorrData = (double *)calloc((size_t)xCorrArraySize, (size_t)sizeof(double));
  if (pdTmpFastXcorrData == NULL) {
    fprintf(stderr, " Error - calloc(pdTmpFastXcorrData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  pfFastXcorrData = (float *)calloc((size_t)xCorrArraySize, (size_t)sizeof(float));
  if (pfFastXcorrData == NULL) {
    fprintf(stderr, " Error - calloc(pfFastXcorrData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  // Create data for correlation analysis.
  MakeCorrData(pdTempRawData, &pPre, 50.0);

  // Make fast xcorr spectrum.
  dSum=0.0;
  for (i=0; i<75; i++) dSum += pPre.pdCorrelationData[i].intensity;
  for (i=75; i < xCorrArraySize +75; i++) {
    if (i<xCorrArraySize) dSum += pPre.pdCorrelationData[i].intensity;
    if (i>=151) dSum -= pPre.pdCorrelationData[i-151].intensity;
    pdTmpFastXcorrData[i-75] = (dSum - pPre.pdCorrelationData[i-75].intensity)* 0.0066666667;
  }

  xCorrSparseArraySize=1;
  for (i=0; i<xCorrArraySize; i++) {
    dTmp = pPre.pdCorrelationData[i].intensity - pdTmpFastXcorrData[i];
    pfFastXcorrData[i] = (float)dTmp;

    // Add flanking peaks if used
    iTmp = i-1;
    if (iTmp >= 0) pfFastXcorrData[i] += (float) ((pPre.pdCorrelationData[iTmp].intensity - pdTmpFastXcorrData[iTmp])*0.5);

    iTmp = i+1;
    if (iTmp < xCorrArraySize) pfFastXcorrData[i] += (float) ((pPre.pdCorrelationData[iTmp].intensity - pdTmpFastXcorrData[iTmp])*0.5);

    //MH: Count number of sparse entries needed
    if (i>0 && pfFastXcorrData[i] != pfFastXcorrData[i-1]) xCorrSparseArraySize++;

  }
  pfFastXcorrData[0] = 0.0;

  free(pPre.pdCorrelationData);
  free(pdTmpFastXcorrData);

  //MH: Add one more slot for the last bin
  xCorrSparseArraySize++;

  //MH: Fill sparse matrix
  if(xCorrSparseArray!=NULL) free(xCorrSparseArray);
  xCorrSparseArray = (kSparseMatrix *)calloc((size_t)xCorrSparseArraySize, (size_t)sizeof(kSparseMatrix));
  if (xCorrSparseArray == NULL) {
    fprintf(stderr, " Error - calloc(pScoring->pSparseFastXcorrData[%d]).\n\n", xCorrSparseArraySize);
    exit(1);
  }
  xCorrSparseArray[0].bin=0;
  xCorrSparseArray[0].fIntensity=0;
  j=1;

  for (i=1; i<xCorrArraySize; i++){

    if (pfFastXcorrData[i] != pfFastXcorrData[i-1]){
      xCorrSparseArray[j].bin = i;
      xCorrSparseArray[j++].fIntensity = pfFastXcorrData[i];
    }
  }
  xCorrSparseArray[j].bin=i;
  xCorrSparseArray[j].fIntensity=0;

  free(pfFastXcorrData);
  free(pdTempRawData);

}

void KSpectrum::kojakXCorr(){
  int i;
  int j;
  int iTmp;
  double dTmp;
  double dSum;
  double *pdTempRawData;
  double *pdTmpFastXcorrData;
  float  *pfFastXcorrData;
  kPreprocessStruct pPre;

  pPre.iHighestIon = 0;
  pPre.dHighestIntensity = 0;

  BinIons(&pPre);
  //cout << scanNumber << ": " << kojakBins << "\t" << xCorrArraySize << "\t" << invBinSize << "\t" << (int)invBinSize+1 << endl;
  kojakSparseArray=new char*[kojakBins];
  for(i=0;i<kojakBins;i++) kojakSparseArray[i]=NULL;

  pdTempRawData = (double *)calloc((size_t)xCorrArraySize, (size_t)sizeof(double));
  if (pdTempRawData == NULL) {
    fprintf(stderr, " Error - calloc(pdTempRawData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  pdTmpFastXcorrData = (double *)calloc((size_t)xCorrArraySize, (size_t)sizeof(double));
  if (pdTmpFastXcorrData == NULL) {
    fprintf(stderr, " Error - calloc(pdTmpFastXcorrData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  pfFastXcorrData = (float *)calloc((size_t)xCorrArraySize, (size_t)sizeof(float));
  if (pfFastXcorrData == NULL) {
    fprintf(stderr, " Error - calloc(pfFastXcorrData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  // Create data for correlation analysis.
  MakeCorrData(pdTempRawData, &pPre, 50.0);

  // Make fast xcorr spectrum.
  dSum=0.0;
  for (i=0; i<75; i++) dSum += pPre.pdCorrelationData[i].intensity;
  for (i=75; i < xCorrArraySize +75; i++) {
    if (i<xCorrArraySize) dSum += pPre.pdCorrelationData[i].intensity;
    if (i>=151) dSum -= pPre.pdCorrelationData[i-151].intensity;
    pdTmpFastXcorrData[i-75] = (dSum - pPre.pdCorrelationData[i-75].intensity)* 0.0066666667;
  }

  xCorrSparseArraySize=1;
  for (i=0; i<xCorrArraySize; i++) {
    dTmp = pPre.pdCorrelationData[i].intensity - pdTmpFastXcorrData[i];
    pfFastXcorrData[i] = (float)dTmp;

    // Add flanking peaks if used
    iTmp = i-1;
    if (iTmp >= 0) pfFastXcorrData[i] += (float) ((pPre.pdCorrelationData[iTmp].intensity - pdTmpFastXcorrData[iTmp])*0.5);

    iTmp = i+1;
    if (iTmp < xCorrArraySize) pfFastXcorrData[i] += (float) ((pPre.pdCorrelationData[iTmp].intensity - pdTmpFastXcorrData[iTmp])*0.5);

  }
  free(pdTmpFastXcorrData);

  //MH: Fill sparse matrix
  for(i=0;i<xCorrArraySize;i++){
    if(pfFastXcorrData[i]>0.5 || pfFastXcorrData[i]<-0.5){

      //Fill in missing masses as a result of adding flanking peaks
      if(pPre.pdCorrelationData[i].mass==0){
        j=1;
        while(true){
          if( (i+j)<xCorrArraySize){
            if(pPre.pdCorrelationData[i+j].mass>0){
              pPre.pdCorrelationData[i].mass=pPre.pdCorrelationData[i+j].mass-j*binSize;
              break;
            }
          }
          if( (i-j)>-1){
            if(pPre.pdCorrelationData[i-j].mass>0){
              pPre.pdCorrelationData[i].mass=pPre.pdCorrelationData[i-j].mass+j*binSize;
              break;
            }
          }
          j++;
        }
      }

      //convert i to sparse array key
      //dTmp=pPre.pdCorrelationData[i].mass+binSize*binOffset;
      dTmp=binSize*i;
      iTmp=(int)dTmp;
      //cout << i << "\t" << pfFastXcorrData[i] << "\t" << dTmp << "\t" << iTmp << endl;
      if(kojakSparseArray[iTmp]==NULL) {
        kojakSparseArray[iTmp]=new char[(int)invBinSize+1];
        for(j=0;j<(int)invBinSize+1;j++) kojakSparseArray[iTmp][j]=0;
      }
      j=(int)((dTmp-iTmp)*invBinSize/*+0.5*/);
      //cout << (dTmp-iTmp) << "\t" << (dTmp-iTmp)*invBinSize/*+0.5*/ << endl;
      //cout << j << endl;
      //if( j>(int)invBinSize) {
      //  cout << "ERROR!" << endl;
      //  exit(0);
      //}
      if(pfFastXcorrData[i]>127) kojakSparseArray[iTmp][j]=127;
      else if(pfFastXcorrData[i]<-128) kojakSparseArray[iTmp][j]=-128;
      else if(pfFastXcorrData[i]>0) kojakSparseArray[iTmp][j]=(char)(pfFastXcorrData[i]+0.5);
      else kojakSparseArray[iTmp][j]=(char)(pfFastXcorrData[i]-0.5);
      //cout << i << "\t" << iTmp << "\t" << j << "\t" << (int)kojakSparseArray[iTmp][j] << endl;
    }
  }

  /*
  if(scanNumber==11368){
    for(i=0;i<kojakBins;i++){
      if(kojakSparseArray[i]==NULL) {
        cout << i << "\tNULL" << endl;
        continue;
      }
      for(j=0;j<(int)invBinSize+1;j++){
        cout << i << "\t" << j << "\t" << (int)kojakSparseArray[i][j] << endl;
      }
    }
  }
  */

  //exit(1);
  free(pPre.pdCorrelationData);
  free(pfFastXcorrData);
  free(pdTempRawData);

}

void KSpectrum::BinIons(kPreprocessStruct *pPre) {
  int i;
  unsigned int j;
  double dPrecursor;
  double dIon;
  double dIntensity;

  // Just need to pad iArraySize by 75.
  dPrecursor=0;
  for(j=0;j<precursor->size();j++){
    if(precursor->at(j).monoMass>dPrecursor) dPrecursor=precursor->at(j).monoMass;
  }
  xCorrArraySize = (int)((dPrecursor + 100.0) / binSize);
  kojakBins = (int)(spec->at(spec->size()-1).mass+100.0);

  pPre->pdCorrelationData = (kSpecPoint *)calloc(xCorrArraySize, (size_t)sizeof(kSpecPoint));
  if (pPre->pdCorrelationData == NULL) {
    fprintf(stderr, " Error - calloc(pdCorrelationData[%d]).\n\n", xCorrArraySize);
    exit(1);
  }

  i = 0;
  while(true) {
    if (i >= (int)spec->size()) break;

    dIon = spec->at(i).mass;
    dIntensity = spec->at(i).intensity;   
    i++;

    if (dIntensity > 0.0) {
      if (dIon < (dPrecursor + 50.0)) {

        //#define BIN(dMass) (int)(dMass*invBinSize + binOffset)
        int iBinIon = (int)(dIon*invBinSize+binOffset);
        dIntensity = sqrt(dIntensity);
        if (iBinIon > pPre->iHighestIon) pPre->iHighestIon = iBinIon;

        if ((iBinIon < xCorrArraySize) && (dIntensity > pPre->pdCorrelationData[iBinIon].intensity)) {
          if (dIntensity > pPre->pdCorrelationData[iBinIon].intensity) {
            pPre->pdCorrelationData[iBinIon].intensity = (float)dIntensity;
            pPre->pdCorrelationData[iBinIon].mass = dIon;
          }
          if (pPre->pdCorrelationData[iBinIon].intensity > pPre->dHighestIntensity) pPre->dHighestIntensity = pPre->pdCorrelationData[iBinIon].intensity;    
        }
      }
    }
  }

  //Clear spectrum data that we no longer need
  spec->clear();

}

// pdTempRawData now holds raw data, pdCorrelationData is windowed data.
void KSpectrum::MakeCorrData(double *pdTempRawData, kPreprocessStruct *pPre, double scale){
  int  i;
  int  ii;
  int  iBin;
  int  iWindowSize;
  int  iNumWindows=10;
  double dMaxWindowInten;
  double dMaxOverallInten;
  double dTmp1;
  double dTmp2;

  dMaxOverallInten = 0.0;

  // Normalize maximum intensity to 100.
  dTmp1 = 1.0;
  if (pPre->dHighestIntensity > 0.000001) dTmp1 = 100.0 / pPre->dHighestIntensity;

  for (i=0; i < xCorrArraySize; i++) {
    pdTempRawData[i] = pPre->pdCorrelationData[i].intensity*dTmp1;
    pPre->pdCorrelationData[i].intensity=0.0;
    if (dMaxOverallInten < pdTempRawData[i]) dMaxOverallInten = pdTempRawData[i];
  }

  iWindowSize = (int) ceil( (double)(pPre->iHighestIon)/iNumWindows);

  for (i=0; i<iNumWindows; i++){
    dMaxWindowInten = 0.0;
    for (ii=0; ii<iWindowSize; ii++) {   // Find max inten. in window.
      iBin = i*iWindowSize+ii;
      if (iBin < xCorrArraySize) {
        if (pdTempRawData[iBin] > dMaxWindowInten)dMaxWindowInten = pdTempRawData[iBin];
      }
    }

    if (dMaxWindowInten > 0.0) {
      dTmp1 = scale / dMaxWindowInten;
      dTmp2 = 0.05 * dMaxOverallInten;

      for (ii=0; ii<iWindowSize; ii++){    // Normalize to max inten. in window.      
        iBin = i*iWindowSize+ii;
        if (iBin < xCorrArraySize){
          if (pdTempRawData[iBin] > dTmp2) pPre->pdCorrelationData[iBin].intensity = (float)(pdTempRawData[iBin]*dTmp1);
        }
      }
    }
  }

}


/*============================
  Utilities
============================*/
int KSpectrum::compareIntensity(const void *p1, const void *p2){
  const kSpecPoint d1 = *(kSpecPoint *)p1;
  const kSpecPoint d2 = *(kSpecPoint *)p2;
  if(d1.intensity<d2.intensity) return -1;
  else if(d1.intensity>d2.intensity) return 1;
  else return 0;
}

int KSpectrum::compareMZ(const void *p1, const void *p2){
  const kSpecPoint d1 = *(kSpecPoint *)p1;
  const kSpecPoint d2 = *(kSpecPoint *)p2;
  if(d1.mass<d2.mass) return -1;
  else if(d1.mass>d2.mass) return 1;
  else return 0;
}
