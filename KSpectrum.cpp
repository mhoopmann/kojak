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
}

KSpectrum::KSpectrum(const KSpectrum& p){
  unsigned int i;
  int j;
  spec = new vector<kSpecPoint>;
  for(i=0;i<p.spec->size();i++) spec->push_back(p.spec->at(i));
  for(i=0;i<20;i++) topHit[i]=p.topHit[i];
  precursor = new vector<kPrecursor>;
  for(i=0;i<p.precursor->size();i++) precursor->push_back(p.precursor->at(i));

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

}
  
KSpectrum::~KSpectrum(){
  delete spec;
  delete precursor;
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

void KSpectrum::addPrecursor(kPrecursor& p){
  precursor->push_back(p);
}

void KSpectrum::clear(){
  spec->clear();
  precursor->clear();
}

void KSpectrum::erasePrecursor(int i){
  precursor->erase(precursor->begin()+i);
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
void KSpectrum::checkScore(kScoreCard& s){
  unsigned int i;
  unsigned int j;

  //char str[256];
  //sprintf(str,"%d.txt",scanNumber);
  //FILE* f=fopen(str,"at");
  //fprintf(f,"%.4f\n",s.simpleScore);
  //fclose(f);

  //edge case for similar precursors: check if already matches top hit identically
  if(s.simpleScore==topHit[0].simpleScore && s.pep1==topHit[0].pep1 && s.pep2==topHit[0].pep2 && s.k1==topHit[0].k1 && s.k2==topHit[0].k2){
    if(s.mods1->size()==topHit[0].mods1->size() && s.mods2->size()==topHit[0].mods2->size()){
      for(i=0;i<s.mods1->size();i++){
        if(s.mods1->at(i).mass!=topHit[0].mods1->at(i).mass || s.mods1->at(i).pos!=topHit[0].mods1->at(i).pos) break;
      }
      for(j=0;j<s.mods2->size();j++){
        if(s.mods2->at(j).mass!=topHit[0].mods2->at(j).mass || s.mods2->at(j).pos!=topHit[0].mods2->at(j).pos) break;
      }
      if(i==s.mods1->size() && j==s.mods2->size()) return;
    }
  }

  for(i=0;i<20;i++){
    if(s.simpleScore > topHit[i].simpleScore) {
      for(j=19;j>i;j--) {
        topHit[j]=topHit[j-1];
      }
      topHit[i] = s;
      return;
    }
  }
}

void KSpectrum::checkSingletScore(kSingletScoreCard& s){

  kSingletScoreCard* sc;
  kSingletScoreCard* cur;

  //if(scanNumber==7653) cout << "SCAN NUMBER: " << scanNumber << "\t" << singletCount << endl;
  
  //If list is empty, add the score card
  if(singletCount==0){
    //if(scanNumber==7653) cout << "First time: " << s.simpleScore << endl;
    singletFirst=new kSingletScoreCard(s);
    singletLast=singletFirst;
    singletCount++;
    //if(scanNumber==7653) cout << "New Count: " << singletCount << endl;
    return;
  }  

  //check if we can just add to the end
  if(s.simpleScore<singletLast->simpleScore){
    //if(scanNumber==7653) cout << "Tail end " << s.simpleScore << " < " << singletLast->simpleScore << endl;
    //check if we need to store the singlet
    if(singletCount==singletMax) return;

    singletLast->next=new kSingletScoreCard(s);
    singletLast->next->prev=singletLast;
    singletLast=singletLast->next;
    singletCount++;
    return;
  }
  
  //check if it goes in the front
  if(s.simpleScore>=singletFirst->simpleScore){
    //if(scanNumber==7653) cout << "First " << s.simpleScore << " >= " << singletFirst->simpleScore << endl;

    singletFirst->prev=new kSingletScoreCard(s);
    singletFirst->prev->next=singletFirst;
    singletFirst=singletFirst->prev;
    if(singletCount<singletMax) {
      singletCount++;
      //if(scanNumber==7653) cout << "New count: " << singletCount  << endl;
    } else {
      //if(scanNumber==7653) cout << "Delteting" << endl;
      cur=singletLast;
      singletLast=singletLast->prev;
      singletLast->next=NULL;
      delete cur;
    }
    return;
  }


  //scan to find insertion point
  
  //if(scanNumber==7653) {
  //  cout << "The List: " << endl;
  //  cur = singletFirst;
  //  while(cur!=NULL) {
  //    cout << "\t" << cur->simpleScore << "\t" << cur << endl;
  //    cur=cur->next;
  //  }
  //}
  

  //if(scanNumber==7653) cout << "Iterating from " << singletCount << endl;
  cur = singletFirst->next;
  int i=1;
  while(s.simpleScore < cur->simpleScore){
    i++;
    //if(scanNumber==7653) cout << i << "\t" << s.simpleScore << "\t" << cur->simpleScore << endl;
    cur=cur->next;
    //if(cur==NULL) cout << "We have a problem" << endl;
  }

  //if(scanNumber==7653) cout << "Adding " << s.simpleScore << " before: " << i+1 << "\t" << cur->simpleScore << endl;
  sc=new kSingletScoreCard(s);
  sc->prev=cur->prev;
  sc->next=cur;
  cur->prev->next=sc;
  cur->prev=sc;
  if(sc->prev==NULL) singletFirst=sc;
  if(singletCount<singletMax) {
    singletCount++;
    //if(scanNumber==7653) cout << "New count: " << singletCount  << endl;
  } else {
    //if(scanNumber==7653) cout << "Delteting" << endl;
    cur=singletLast;
    singletLast=singletLast->prev;
    singletLast->next=NULL;
    delete cur;
  }

  //if(scanNumber==7653) cout << "Done!" << endl;

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
