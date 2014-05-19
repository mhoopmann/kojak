#include "KSpectrum.h"

/*============================
  Constructors & Destructors
============================*/
KSpectrum::KSpectrum(int i){
  binSize=0.03;
  invBinSize=1.0/binSize;
  charge = 0; 
  maxIntensity=0;
  mz = 0;
  precursor = new vector<kPrecursor>;
  spec = new vector<specPoint>;
  scanNumber = 0;
  rTime = 0;
  xCorrArraySize=0;
  xCorrSparseArraySize=0;
  xCorrSparseArray=NULL;
  topCount=i;
  topSinglet = new kSingletScoreCard[topCount];
}

KSpectrum::KSpectrum(const KSpectrum& p){
  unsigned int i;
  int j;
  spec = new vector<specPoint>;
  for(i=0;i<p.spec->size();i++) spec->push_back(p.spec->at(i));
  for(i=0;i<20;i++) topHit[i]=p.topHit[i];
  topCount=p.topCount;
  topSinglet = new kSingletScoreCard[topCount];
  for(j=0;j<topCount;j++) topSinglet[j]=p.topSinglet[j];
  precursor = new vector<kPrecursor>;
  for(i=0;i<p.precursor->size();i++) precursor->push_back(p.precursor->at(i));

  binSize = p.binSize;
  invBinSize = p.invBinSize;
  charge= p.charge;
  maxIntensity = p.maxIntensity;
  mz = p.mz;
  scanNumber = p.scanNumber;
  rTime = p.rTime;
  xCorrArraySize = p.xCorrArraySize;
  xCorrSparseArraySize = p.xCorrSparseArraySize;

  if(p.xCorrSparseArray==NULL){
    xCorrSparseArray=NULL;
  } else {
    xCorrSparseArray = (kSparseMatrix *)calloc((size_t)xCorrSparseArraySize, (size_t)sizeof(kSparseMatrix));
    for(j=0;j<xCorrSparseArraySize;j++) xCorrSparseArray[j]=p.xCorrSparseArray[j];
  }
}
  
KSpectrum::~KSpectrum(){
  delete spec;
  delete precursor;
  delete [] topSinglet;
  if(xCorrSparseArray!=NULL) free(xCorrSparseArray);
}


/*============================
  Operators
============================*/
KSpectrum& KSpectrum::operator=(const KSpectrum& p){
  if(this!=&p){
    unsigned int i;
    int j;
    delete spec;
    spec = new vector<specPoint>;
    for(i=0;i<p.spec->size();i++) spec->push_back(p.spec->at(i));
    for(i=0;i<20;i++) topHit[i]=p.topHit[i];
    topCount=p.topCount;
    delete [] topSinglet;
    topSinglet = new kSingletScoreCard[p.topCount];
    for(j=0;j<p.topCount;j++) topSinglet[j]=p.topSinglet[j];
    delete precursor;
    precursor = new vector<kPrecursor>;
    for(i=0;i<p.precursor->size();i++) precursor->push_back(p.precursor->at(i));

    binSize = p.binSize;
    invBinSize = p.invBinSize;
    charge = p.charge;
    maxIntensity = p.maxIntensity;
    mz = p.charge;
    scanNumber = p.scanNumber;
    rTime = p.rTime;
    xCorrArraySize = p.xCorrArraySize;
    xCorrSparseArraySize = p.xCorrSparseArraySize;

    if(xCorrSparseArray!=NULL) free(xCorrSparseArray);
    if(p.xCorrSparseArray==NULL){
      xCorrSparseArray=NULL;
    } else {
      xCorrSparseArray = (kSparseMatrix *)calloc((size_t)xCorrSparseArraySize, (size_t)sizeof(kSparseMatrix));
      for(j=0;j<xCorrSparseArraySize;j++) xCorrSparseArray[j]=p.xCorrSparseArray[j];
    }
  }
  return *this;
}

specPoint& KSpectrum::operator [](const int &i){
  return spec->at(i);
}


/*============================
  Accessors
============================*/
int KSpectrum::getCharge(){
  return charge;
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

float KSpectrum::getRTime(){
  return rTime;
}

int KSpectrum::getScanNumber(){
  return scanNumber;
}

kScoreCard& KSpectrum::getScoreCard(int i){
  return topHit[i];
}

kSingletScoreCard& KSpectrum::getSingletScoreCard(int i){
  return topSinglet[i];
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
void KSpectrum::addPoint(specPoint& s){
  spec->push_back(s);
}

void KSpectrum::addPrecursor(kPrecursor& p){
  precursor->push_back(p);
}

void KSpectrum::clear(){
  spec->clear();
  precursor->clear();
}

void KSpectrum::setCharge(int i){
  charge=i;
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

  for(i=0;i<20;i++){
    if(s.simpleScore > topHit[i].simpleScore) {
      for(j=19;j>i;j--) topHit[j]=topHit[j-1];
      topHit[i] = s;
      return;
    }
  }
}

void KSpectrum::checkSingletScore(kSingletScoreCard& s){
  int i;
  int j;

  for(i=0;i<topCount;i++){
    if(s.simpleScore > topSinglet[i].simpleScore) {
      for(j=topCount-1;j>i;j--) topSinglet[j]=topSinglet[j-1];
      topSinglet[i] = s;
      return;
    }
  }
}

void KSpectrum::sortMZ(){
  qsort(&spec->at(0),spec->size(),sizeof(specPoint),compareMZ);
}

void KSpectrum::xCorrScore(){
  CometXCorr();
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
  MakeCorrData(pdTempRawData, &pPre);

  // Make fast xcorr spectrum.
  dSum=0.0;
  for (i=0; i<75; i++) dSum += pPre.pdCorrelationData[i];
  for (i=75; i < xCorrArraySize +75; i++) {
    if (i<xCorrArraySize) dSum += pPre.pdCorrelationData[i];
    if (i>=151) dSum -= pPre.pdCorrelationData[i-151];
    pdTmpFastXcorrData[i-75] = (dSum - pPre.pdCorrelationData[i-75])* 0.0066666667;
  }

  xCorrSparseArraySize=1;
  for (i=0; i<xCorrArraySize; i++) {
    dTmp = pPre.pdCorrelationData[i] - pdTmpFastXcorrData[i];
    pfFastXcorrData[i] = (float)dTmp;

    // Add flanking peaks if used
    iTmp = i-1;
    if (iTmp >= 0) pfFastXcorrData[i] += (float) ((pPre.pdCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp])*0.5);

    iTmp = i+1;
    if (iTmp < xCorrArraySize) pfFastXcorrData[i] += (float) ((pPre.pdCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp])*0.5);

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

  pPre->pdCorrelationData = (double *)calloc(xCorrArraySize, (size_t)sizeof(double));
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

        //#define BIN(dMass) (int)(dMass*invBinSize + 1.0)
        int iBinIon = (int)(dIon*invBinSize+1.0); //Why +1.0?
        dIntensity = sqrt(dIntensity);
        if (iBinIon > pPre->iHighestIon) pPre->iHighestIon = iBinIon;

        if ((iBinIon < xCorrArraySize) && (dIntensity > pPre->pdCorrelationData[iBinIon])) {
          if (dIntensity > pPre->pdCorrelationData[iBinIon]) pPre->pdCorrelationData[iBinIon] = dIntensity;
          if (pPre->pdCorrelationData[iBinIon] > pPre->dHighestIntensity) pPre->dHighestIntensity = pPre->pdCorrelationData[iBinIon];    
        }
      }
    }
  }

  //Clear spectrum data that we no longer need
  spec->clear();

}

// pdTempRawData now holds raw data, pdCorrelationData is windowed data.
void KSpectrum::MakeCorrData(double *pdTempRawData, kPreprocessStruct *pPre){
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
    pdTempRawData[i] = pPre->pdCorrelationData[i]*dTmp1;
    pPre->pdCorrelationData[i]=0.0;
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
      dTmp1 = 50.0 / dMaxWindowInten;
      dTmp2 = 0.05 * dMaxOverallInten;

      for (ii=0; ii<iWindowSize; ii++){    // Normalize to max inten. in window.      
        iBin = i*iWindowSize+ii;
        if (iBin < xCorrArraySize){
          if (pdTempRawData[iBin] > dTmp2) pPre->pdCorrelationData[iBin] = pdTempRawData[iBin]*dTmp1;
        }
      }
    }
  }

}


/*============================
  Utilities
============================*/
int KSpectrum::compareIntensity(const void *p1, const void *p2){
  const specPoint d1 = *(specPoint *)p1;
  const specPoint d2 = *(specPoint *)p2;
  if(d1.intensity<d2.intensity) return -1;
  else if(d1.intensity>d2.intensity) return 1;
  else return 0;
}

int KSpectrum::compareMZ(const void *p1, const void *p2){
  const specPoint d1 = *(specPoint *)p1;
  const specPoint d2 = *(specPoint *)p2;
  if(d1.mass<d2.mass) return -1;
  else if(d1.mass>d2.mass) return 1;
  else return 0;
}
