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

  nativeID.clear();

  int j;
  for (j = 0; j<HISTOSZ; j++) histogram[j] = 0;
  histogramCount = 0;
  histoMaxIndex = 0;

  //for (j = 0; j<HISTOSZ; j++) histogramSinglet[j] = 0;
  //histogramSingletCount = 0;

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
  //tmpSingCount=0;
  tmpHistCount=0;

  cc=0;
  sc=0;

  peakCounts=0;
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
  nativeID=p.nativeID;

  for (i = 0; i<HISTOSZ; i++) histogram[i] = p.histogram[i];
  histogramCount = p.histogramCount;
  histoMaxIndex = p.histoMaxIndex;

  //for (i = 0; i<HISTOSZ; i++) histogramSinglet[i] = p.histogramSinglet[i];
  //histogramSingletCount = p.histogramSingletCount;

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
  //tmpSingCount = p.tmpSingCount;
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
    nativeID = p.nativeID;

    for (i = 0; i<HISTOSZ; i++) histogram[i] = p.histogram[i];
    histogramCount = p.histogramCount;
    histoMaxIndex = p.histoMaxIndex;

    //for (i = 0; i<HISTOSZ; i++) histogramSinglet[i] = p.histogramSinglet[i];
    //histogramSingletCount = p.histogramSingletCount;

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
    //tmpSingCount = p.tmpSingCount;
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

string KSpectrum::getNativeID(){
  return nativeID;
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

void KSpectrum::setNativeID(string s){
  nativeID=s;
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
bool KSpectrum::calcEValue(kParams* params, KDecoys& decoys, KDatabase& db) {
  int i;
  int iLoopCount;
  int iMaxCorr;
  int iStartCorr;
  int iNextCorr;
  double dSlope;
  double dIntercept;
  double dRSquare;
  bool bSkipXL=false;
  bool bSingletFail=false;

  //tmpSingCount = histogramSingletCount;
  tmpHistCount = histogramCount;
  if (topHit[0].simpleScore == 0) return true; //no need to do any of this if there are no PSMs...

  //precompute which ion series to use
  decoyIonSz=0;
  for (i = 0; i<6; i++){
    if (params->ionSeries[i]) {
      if (i<3) {
        decoyIons[decoyIonSz].b = true;
        if (i == 0) decoyIons[decoyIonSz++].mass = -27.9949141;
        else if (i == 1) decoyIons[decoyIonSz++].mass = 0;
        else decoyIons[decoyIonSz++].mass = 17.026547;
      } else {
        decoyIons[decoyIonSz].b = false;
        if (i == 3) decoyIons[decoyIonSz++].mass = 25.9792649;
        else if (i == 4) decoyIons[decoyIonSz++].mass = 0;
        else decoyIons[decoyIonSz++].mass = -16.0187224;
      }
    }
  }

  if (histogramCount < decoys.decoySize) {
    if (!generateXcorrDecoys(params, decoys)) return false;
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

  //reorder top scoring peptide so that ties always appear in the same order instead
  //of the order in which the search threads finished (which can change from run to run).
  string dStr = params->decoy;
  refreshScore(db,dStr);

  iLoopCount = 20; //score all e-values among top hits?
  double topScore=topHit[0].simpleScore;
  for (i = 0; i<iLoopCount; i++) {
    if (topHit[i].simpleScore == 0) break; //score all e-values among top hits?
    if (dSlope >= 0.0) {
      topHit[i].eVal = 1e12;
    } else {
      topHit[i].eVal = pow(10.0, dSlope * topHit[i].simpleScore + dIntercept);
      if (topHit[i].eVal>1e12) topHit[i].eVal = 1e12;
    }

    //score individual peptides
    if(topHit[i].score2>0){
      if(topHit[i].simpleScore==topScore){ //only do this for top hits right now: it is slow...
        if(i>0){
          //check if we've computed these already - happens with one of the peptides in ties.
          //Note: there are slight differences if the alternate peptide in a tie score has a slightly different mass, resulting in a different decoy
          //distribution should the order of peptides change in the next run.
          if(topHit[i].score1==topHit[i-1].score1) topHit[i].eVal1=topHit[i-1].eVal1;
          else topHit[i].eVal1 = generateSingletDecoys2(params, decoys, topHit[i].score1, topHit[i].mass1, (int)topHit[i].precursor);
          if(topHit[i].score2==topHit[i-1].score2) topHit[i].eVal2=topHit[i-1].eVal2;
          else topHit[i].eVal2 = generateSingletDecoys2(params, decoys, topHit[i].score2, topHit[i].mass2, (int)topHit[i].precursor);
        } else {
          topHit[i].eVal1 = generateSingletDecoys2(params,decoys,topHit[i].score1,topHit[i].mass1,(int)topHit[i].precursor);
          topHit[i].eVal2 = generateSingletDecoys2(params, decoys, topHit[i].score2, topHit[i].mass2, (int)topHit[i].precursor);
        }
      }
    } else {
      if (dSlope >= 0.0) {
        topHit[i].eVal1 = 1e12;
      } else {
        topHit[i].eVal1 = pow(10.0, dSlope * topHit[i].score1 + dIntercept);
        if (topHit[i].eVal1>1e12) topHit[i].eVal1 = 1e12;
      }
      topHit[i].eVal2 = 1e12;
    }
  }

  return true;
}

bool KSpectrum::checkDecoy(KDatabase& db, string& dStr, kScoreCard& hit){
  bool bDecoy = false;
  size_t i;
  kPeptide pep;
  pep = db.getPeptide(hit.pep1);
  for (i = 0; i<pep.map->size(); i++){
    if (db[pep.map->at(i).index].name.find(dStr) != string::npos) bDecoy = true;
  }
  if (!bDecoy && hit.pep2>-1){ //only check second peptide if necessary
    pep = db.getPeptide(hit.pep2);
    for (i = 0; i<pep.map->size(); i++){
      if (db[pep.map->at(i).index].name.find(dStr) != string::npos) bDecoy = true;
    }
  }
  return bDecoy;
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


//from Comet
// Make synthetic decoy spectra to fill out correlation histogram by going
// through each candidate peptide and rotating spectra in m/z space.
double KSpectrum::generateSingletDecoys2(kParams* params, KDecoys& decoys, double xcorr, double mass, int preIndex) {
  int i;
  int n;
  int j;
  int k;
  int maxZ;
  int z;
  int key;
  int pos;
  int xlSite=0;
  int xlLen;
  double dXcorr;
  double dFragmentIonMass = 0.0;
  double diffMass;
  double preMass;

  //int tempHistogram[HISTOSZ];
  for(i=0;i<HISTOSZ;i++) tempHistogram[i]=0;

  int seed = (int)(scanNumber*mass); 
  if (seed<0) seed = -seed;
  seed = seed % DECOY_SIZE; //don't always start at the top, but not random either; remains reproducible across threads
  int decoyIndex;

  //compute modification mass
  maxZ = precursor->at(preIndex).charge;
  if (maxZ>4) maxZ = 4;
  preMass = precursor->at(preIndex).monoMass;
  diffMass=preMass-mass;

  //Does this function need as many DECOY_SIZE as the other? Can this be shortened?
  for (i = 0; i<decoys.decoySize; i++) { // iterate through required # decoys
    dXcorr = 0.0;
    decoyIndex = (seed + i) % DECOY_SIZE;
 
    //find link site - somewhat wasted cycles
    for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {
      if (decoys.decoyIons[decoyIndex].pdIonsN[j]>mass) break;
    }
    if(j<1) return 1e12;

    xlSite++;
    if(xlSite>=(j-1)) xlSite=0;
    xlLen = j;

    for (n = 0; n<decoyIonSz; n++) { //iterate over each ion series
      for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions
      
        if (decoyIons[n].b) {
          dFragmentIonMass = decoys.decoyIons[decoyIndex].pdIonsN[j] + decoyIons[n].mass;
          if (j >= xlSite) dFragmentIonMass+=diffMass;
        } else {
          dFragmentIonMass = decoys.decoyIons[decoyIndex].pdIonsC[j] + decoyIons[n].mass;
          if (j >= xlLen - xlSite) dFragmentIonMass += diffMass;
        }
        if (dFragmentIonMass>preMass) break;

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

    //score the cleavage product ions
    if(params->cleavageProducts->size()>0){
      for (n = 0; n<decoyIonSz; n++) { //iterate over each ion series
        for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions

          if (decoyIons[n].b) {
            if (j < xlSite) continue;
          } else {
            if (j < xlLen - xlSite) continue;
          }

          for (size_t a = 0; a<params->cleavageProducts->size();a++){
            if (decoyIons[n].b) {
              dFragmentIonMass = decoys.decoyIons[decoyIndex].pdIonsN[j] + decoyIons[n].mass + params->cleavageProducts->at(a);
            } else {
              dFragmentIonMass = decoys.decoyIons[decoyIndex].pdIonsC[j] + decoyIons[n].mass + params->cleavageProducts->at(a);
            }
            if (dFragmentIonMass>preMass) break;

            mz = dFragmentIonMass + 1.007276466; //only check 1+ ions
            mz = params->binSize * (int)(mz*invBinSize + params->binOffset);
            key = (int)mz;
            if (key >= kojakBins) break;
            if (kojakSparseArray[key] == NULL) continue;
            pos = (int)((mz - key)*invBinSize);
            dXcorr += kojakSparseArray[key][pos];
          }
        }
      }
    }

    if (dXcorr <= 0.0) dXcorr = 0.0;
    k = (int)(dXcorr*0.05+ 0.5);  // 0.05=0.005*10; see KAnalysis::kojakScoring
    if (k < 0) k = 0;
    else if (k >= HISTOSZ) k = HISTOSZ - 1;
    tempHistogram[k]++;
  }

  //Do linear regression and compute e-value
  double dSlope,dIntercept,dRSquare;
  int iMaxCorr,iStartCorr,iNextCorr;
  linearRegression4(tempHistogram, dSlope, dIntercept, iMaxCorr, iStartCorr, iNextCorr, dRSquare);
  double eVal = pow(10.0, dSlope * 10 * xcorr + dIntercept);
  if(eVal>1e12) eVal=1e12;

  return eVal;
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
  double dXcorr;
  double dFragmentIonMass = 0.0;
  int myCount=0;

  //tmpSingCount = histogramSingletCount;
  tmpHistCount = histogramCount;

  // DECOY_SIZE is the minimum # of decoys required or else this function isn't
  // called.  So need to generate iLoopMax more xcorr scores for the histogram.
  int iLoopMax = decoys.decoySize - histogramCount;
  int seed = (scanNumber*histogramCount);
  if(seed<0) seed=-seed;
  seed = seed % decoys.decoySize; //don't always start at the top, but not random either; remains reproducible across threads
  int decoyIndex;

  size_t maxPre = precursor->size();
  r=0;

  for (i = 0; i<iLoopMax; i++) { // iterate through required # decoys
    dXcorr = 0.0;
    decoyIndex = (seed + i) % decoys.decoySize;

    //iterate over precursors
    r++;
    if(r>=(int)precursor->size()) r=0;
    maxZ = precursor->at(r).charge;
    if (maxZ>4) maxZ = 4;

    for (n = 0; n<decoyIonSz; n++) { //iterate over each ion series
      for (j = 0; j<MAX_DECOY_PEP_LEN; j++) {  // iterate through decoy fragment ions
      
        if (decoyIons[n].b) dFragmentIonMass = decoys.decoyIons[decoyIndex].pdIonsN[j] + decoyIons[n].mass;
        else dFragmentIonMass = decoys.decoyIons[decoyIndex].pdIonsC[j] + decoyIons[n].mass;
        if (dFragmentIonMass>precursor->at(r).monoMass) break;

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
    
  }
  return true;
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

  //double dCummulative[HISTOSZ];  // Cummulative frequency at each xcorr value.

  int i;
  int bestNC;
  int bestStart;
  int iNextCorr;    // 2nd best xcorr index
  int iMaxCorr = 0;   // max xcorr index
  int iStartCorr;
  int iNumPoints;

  //for diagnostics
  for(i=0;i<HISTOSZ;i++) histogramO[i]=histogram[i];

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

  //double dCummulative[HISTOSZ];  // Cummulative frequency at each xcorr value.

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
  if (iMaxCorr<2) {
    slope = 0;
    intercept = 0;
    iMaxXcorr = 0;
    iStartXcorr = 0;
    iNextXcorr = 0;
    rSquared = 0;
    return;
  }

  //More aggressive version summing everything starting at the max score
  dCummulative[iMaxCorr] = histo[iMaxCorr];
  for (i = iMaxCorr - 1; i >= 0; i--) {
    dCummulative[i] = dCummulative[i + 1] + histo[i];
  }
  //edge case where all decoys are the top score (essentially the other peptide score)
  //This should no longer occur with the 2.0.0alpha15 fixes.
  if(dCummulative[0]<1){ 
    slope = 0;
    intercept = 0;
    iMaxXcorr = 0;
    iStartXcorr = 0;
    iNextXcorr = 0;
    rSquared = 0;
    return;
  }

  // log10...and stomp all over the original...hard to troubleshoot later
  for (i = iMaxCorr; i >= 0; i--)  {
    histo[i] = (int)dCummulative[i];
    dCummulative[i] = log10(dCummulative[i]);
  }

  //avoid edge effects from sampling one really high decoy score.
  while(iMaxCorr>2 && dCummulative[iMaxCorr]<0.01) iMaxCorr--;
  iNextCorr=iMaxCorr;
  iStartCorr = iNextCorr/2;
  if(iNextCorr-iStartCorr<2) iStartCorr=iNextCorr-2;

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
      if (rsq>bestRSQ || iNextCorr - iStartCorr + 1<4){ //keep better RSQ only if more than 4 datapoints, otherwise keep every RSQ below 8 datapoints
        bestRSQ = rsq;
        bestNC = iNextCorr;
        bestSlope = b;
        bestInt = a;
        bestStart = iStartCorr;
      }
      if (bRight){
        if (iNextCorr<iMaxCorr) iNextCorr++;
        else if (iStartCorr>0) iStartCorr--;
        else break;
      } else {
        if (iStartCorr>0) iStartCorr--;
        else if (iNextCorr<iMaxCorr) iNextCorr++;
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

//TODO: improve this function to account for modifications when checking peptide sequence...
void KSpectrum::refreshScore(KDatabase& db, string& dStr){
  //skip any lists that are empty
  if(topHit[0].simpleScore==0) return;

  //if not a tie, we're done now.
  if(topHit[0].simpleScore>topHit[1].simpleScore) return;

  int max;
  for(max=2;max<20;max++){
    if(topHit[max].simpleScore<topHit[0].simpleScore) break;
  }

  //sort from 0 to max
  string pep1, pep2;
  int lType1,lType2;
  bool decoy1,decoy2;
  for(int a=0;a<max-1;a++){
    for(int b=a+1;b<max;b++){

      //sort according to decoy state, linker type (single, loop, cross), peptide sequence, link pos, second peptide, second link pos
      decoy1=checkDecoy(db,dStr,topHit[a]);
      decoy2=checkDecoy(db,dStr,topHit[b]);
      if(decoy2<decoy1) goto swap_hit;  //check decoy status first
      else if(decoy1<decoy2) continue;

      //check link type next
      lType1=0;
      if(topHit[a].k1>=0 && topHit[a].k2>=0) lType1=1;
      if(topHit[a].pep2>=0) lType1=2;
      lType2 = 0;
      if (topHit[b].k1 >= 0 && topHit[b].k2 >= 0) lType2 = 1;
      if (topHit[b].pep2 >= 0) lType2 = 2;
      if(lType2<lType1) goto swap_hit;
      else if(lType1<lType2) continue;

      //check peptide sequence (first only)
      db.getPeptideSeq(topHit[a].pep1,pep1);
      db.getPeptideSeq(topHit[b].pep1,pep2);
      int c;
      c=pep2.compare(pep1);
      if (c<0) goto swap_hit;
      else if(c>0) continue;

      //Check modifications
      if(topHit[b].mods1->size()<topHit[a].mods1->size()) goto swap_hit;
      if(topHit[a].mods1->size()<topHit[b].mods1->size()) continue;
      for(size_t d=0;d<topHit[a].mods1->size();d++){
        if(topHit[b].mods1->at(d).pos<topHit[a].mods1->at(d).pos) goto swap_hit;
        if(topHit[b].mods1->at(d).mass<topHit[a].mods1->at(d).mass) goto swap_hit;
      }

      //check link position (first peptide)
      if(topHit[b].k1<topHit[a].k1) goto swap_hit;
      else if(lType1==1 && lType2==1 && topHit[a].k1==topHit[b].k1){
        if (topHit[b].k2<topHit[a].k2) goto swap_hit; //check 2nd pos on loop links
        continue;
      }

      //check second peptide (cross-links only)
      if(lType1==2 && lType2==2){
        db.getPeptideSeq(topHit[a].pep2, pep1);
        db.getPeptideSeq(topHit[b].pep2, pep2);
        int c; 
        c = pep2.compare(pep1);
        if (c<0) goto swap_hit;
        else if(c==0){

          //check modifications
          if (topHit[b].mods2->size()<topHit[a].mods2->size()) goto swap_hit;
          if (topHit[a].mods2->size()<topHit[b].mods2->size()) continue;
          for (size_t d = 0; d<topHit[a].mods2->size(); d++){
            if (topHit[b].mods2->at(d).pos<topHit[a].mods2->at(d).pos) goto swap_hit;
            if (topHit[b].mods2->at(d).mass<topHit[a].mods2->at(d).mass) goto swap_hit;
          }

          //check link position
          if (topHit[b].k2<topHit[a].k2) goto swap_hit;
        }
      }

      continue; //already in order

      //the swap code
      swap_hit:
      kScoreCard tmp = topHit[a];
      topHit[a] = topHit[b];
      topHit[b] = tmp;
    }
  }

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


/*============================
  Private Functions
============================*/
void KSpectrum::kojakXCorr(double* pdTempRawData, double* pdTmpFastXcorrData, float* pfFastXcorrData, kPreprocessStruct*& pPre){
  int i;
  int j;
  int iTmp;
  double dTmp;
  double dSum;

  pPre->iHighestIon=0;
  pPre->dHighestIntensity=0;
  BinIons(pPre);

  memset(pdTempRawData,0,xCorrArraySize*sizeof(double));
  memset(pdTmpFastXcorrData, 0, xCorrArraySize*sizeof(double));
  memset(pfFastXcorrData, 0, xCorrArraySize*sizeof(float));
  kojakSparseArray=new char*[kojakBins];
  for(i=0;i<kojakBins;i++) kojakSparseArray[i]=NULL;


  // Create data for correlation analysis.
  MakeCorrData(pdTempRawData, pPre, 50.0);

  // Make fast xcorr spectrum.
  kSpecPoint *pdCorrelationData = pPre->pdCorrelationData;
  dSum=0.0;
  for (i=0; i<75; i++) dSum += pdCorrelationData[i].intensity;
  for (i=75; i < xCorrArraySize +75; i++) {
    if (i<xCorrArraySize && pdCorrelationData[i].intensity>0) dSum += pdCorrelationData[i].intensity;
    if (i >= 151 && pdCorrelationData[i - 151].intensity>0) dSum -= pdCorrelationData[i - 151].intensity;
    pdTmpFastXcorrData[i-75] = (dSum - pdCorrelationData[i-75].intensity)* 0.0066666667;
  }

  xCorrSparseArraySize=1;

  double dTmp0 = pdCorrelationData[0].intensity - pdTmpFastXcorrData[0];
  double dTmp1 = pdCorrelationData[1].intensity - pdTmpFastXcorrData[1];
  double dTmp2 = pdCorrelationData[2].intensity - pdTmpFastXcorrData[2];
  pfFastXcorrData[0] = (float)(dTmp0+dTmp1*0.5);
  pfFastXcorrData[1] = (float)(dTmp1 + (dTmp0 + dTmp2)*0.5);
  for(i=2;i<xCorrArraySize-1;i++){
    dTmp0=dTmp1;
    dTmp1=dTmp2;
    dTmp2 = pdCorrelationData[i+1].intensity - pdTmpFastXcorrData[i+1];
    pfFastXcorrData[i]=(float)(dTmp1+(dTmp0+dTmp2)*0.5);
  }
  pfFastXcorrData[xCorrArraySize-1] = (float)(dTmp2 + dTmp1*0.5);

  //MH: Fill sparse matrix
  for(i=0;i<xCorrArraySize;i++){
    if(pfFastXcorrData[i]>0.5 || pfFastXcorrData[i]<-0.5){

      dTmp=binSize*i;
      iTmp=(int)dTmp;
      if(kojakSparseArray[iTmp]==NULL) {
        kojakSparseArray[iTmp]=new char[(int)invBinSize+1];
        for(j=0;j<(int)invBinSize+1;j++) kojakSparseArray[iTmp][j]=0;
      }
      j=(int)((dTmp-iTmp)*invBinSize/*+0.5*/);

      if(pfFastXcorrData[i]>127) kojakSparseArray[iTmp][j]=127;
      else if(pfFastXcorrData[i]<-128) kojakSparseArray[iTmp][j]=-128;
      else if(pfFastXcorrData[i]>0) kojakSparseArray[iTmp][j]=(char)(pfFastXcorrData[i]+0.5);
      else kojakSparseArray[iTmp][j]=(char)(pfFastXcorrData[i]-0.5);
    }
  }

}

void KSpectrum::BinIons(kPreprocessStruct *pPre) {
  int i;
  unsigned int j;
  double dPrecursor;
  double dIon;
  double dIntensity;
  kSpecPoint *pdCorrelationData = pPre->pdCorrelationData;

  // Just need to pad iArraySize by 75.
  dPrecursor=0;
  for(j=0;j<precursor->size();j++){
    if(precursor->at(j).monoMass>dPrecursor) dPrecursor=precursor->at(j).monoMass;
  }
  xCorrArraySize = (int)((spec->at(spec->size() - 1).mass + 100.0) / binSize);
  kojakBins = (int)(spec->at(spec->size()-1).mass+100.0);

  memset(pdCorrelationData, 0, xCorrArraySize*sizeof(kSpecPoint));

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

        if ((iBinIon < xCorrArraySize) && (dIntensity > pdCorrelationData[iBinIon].intensity)) {
          if (dIntensity > pdCorrelationData[iBinIon].intensity) {
            pdCorrelationData[iBinIon].intensity = (float)dIntensity;
            pdCorrelationData[iBinIon].mass = dIon;
          }
          if (pdCorrelationData[iBinIon].intensity > pPre->dHighestIntensity) pPre->dHighestIntensity = pdCorrelationData[iBinIon].intensity;    
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
  int  iNumWindows=10;
  int  iWindowSize = (int)ceil((double)(pPre->iHighestIon) / iNumWindows);
  double dMaxWindowInten[10];
  double dMaxOverallInten;
  double dTmp1;
  double dTmp2;

  kSpecPoint *pdCorrelationData = pPre->pdCorrelationData;
  memset(&dMaxWindowInten,0,10*sizeof(double));

  dMaxOverallInten = 0.0;
  // Normalize maximum intensity to 100.
  dTmp1 = 1.0;
  if (pPre->dHighestIntensity > 0.000001) dTmp1 = 100.0 / pPre->dHighestIntensity;

  int x=0;
  int c=0;
  for (i=0; i < xCorrArraySize; i++) {
    dTmp2 = pdCorrelationData[i].intensity*dTmp1;
    pdTempRawData[i] = dTmp2;
    pdCorrelationData[i].intensity = 0;
    if(x<iNumWindows){
      if(dMaxWindowInten[x]<dTmp2) dMaxWindowInten[x]=dTmp2;
      c++;
      if (c == iWindowSize){
        c=0;
        x++;
      }
    }
  }

  dMaxOverallInten=100;

  dTmp2 = 0.05 * dMaxOverallInten;
  for (i=0; i<iNumWindows; i++){
    if (dMaxWindowInten[i] > 0.0) {
      dTmp1 = scale / dMaxWindowInten[i];

      for (ii=0; ii<iWindowSize; ii++){    // Normalize to max inten. in window.      
        iBin = i*iWindowSize+ii;
        if (iBin < xCorrArraySize){
          if (pdTempRawData[iBin] > dTmp2) pdCorrelationData[iBin].intensity = (float)(pdTempRawData[iBin]*dTmp1);
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
