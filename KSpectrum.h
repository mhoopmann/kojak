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

#ifndef _KSPECTRUM_H
#define _KSPECTRUM_H

#include <cmath>
#include <list>
#include <vector>
#include "KStructs.h"
#include "KTopPeps.h"
#include "CometDecoys.h"

#define HISTOSZ 152

using namespace std;

class KSpectrum {

public:

  //Constructors & Destructors
  KSpectrum(const int& i, const double& bs, const double& os);
  KSpectrum(const KSpectrum& p);
  ~KSpectrum();

  //Operators
  KSpectrum&  operator=(const KSpectrum& p);
  kSpecPoint&  operator[](const int& i);

  //Data Members
  kSparseMatrix*  xCorrSparseArray;
  int             xCorrSparseArraySize;
  char**          kojakSparseArray;
  int             kojakBins;
  
  int singletBins;
  list<kSingletScoreCard*>**  singletList;
  float lowScore;

  int cc;
  int sc;

  //diagnostics - probably temporary
  float tmpIntercept;
  float tmpSlope;
  float tmpIStartCorr;
  float tmpINextCorr;
  double tmpRSquare;
  int tmpIMaxCorr;
  float tmpSingletIntercept;
  float tmpSingletSlope;
  float tmpSingletIStartCorr;
  float tmpSingletINextCorr;
  double tmpSingletRSquare;
  int tmpSingletIMaxCorr;
  int tmpSingCount;
  int tmpHistCount;


  //Accessors
  double              getBinOffset          ();
  int                 getCharge             ();
  bool                getInstrumentPrecursor();
  double              getInvBinSize         ();
  float               getMaxIntensity       ();
  double              getMZ                 ();
  kPrecursor&         getPrecursor          (int i);
  kPrecursor*         getPrecursor2         (int i);
  float               getRTime              ();
  int                 getScanNumber         ();
  kScoreCard&         getScoreCard          (int i);
  int                 getSingletCount       ();
  kSingletScoreCard&  getSingletScoreCard   (int i);
  KTopPeps*           getTopPeps            (int index);
  int                 size                  ();
  int                 sizePrecursor         ();
  
  kSingletScoreCard*    singletFirst;   //pointer to start of linked list
  kSingletScoreCard*    singletLast;    //pointer to end of linked list
  int                   singletMax;

  int histogram[HISTOSZ];
  int histogramCount;
  int histoMaxIndex;
  int histogramSinglet[HISTOSZ];
  int histogramSingletCount;

  //Modifiers
  void addPoint               (kSpecPoint& s);
  void addPrecursor           (kPrecursor& p, int sz);
  void clear                  ();
  void erasePrecursor         (int i);
  void setCharge              (int i);
  void setInstrumentPrecursor (bool b);
  void setMaxIntensity        (float f);
  void setMZ                  (double d);
  void setPrecursor           (double d, int i);
  void setRTime               (float f);
  void setScanNumber          (int i);

  //Functions
  bool  calcEValue          (kParams* params, KDecoys& decoys);
  void  checkScore          (kScoreCard& s);
  void  checkSingletScore   (kSingletScoreCard& s);
  bool  generateSingletDecoys(kParams* params, KDecoys& decoys);
  bool generateXLDecoys      (kParams* params, KDecoys& decoys);
  bool  generateXcorrDecoys (kParams* params, KDecoys& decoys);
  void  linearRegression    (double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared);
  void  linearRegression2   (double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared);
  void  linearRegression3   (double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared);
  void  resetSingletList    ();
  void  sortMZ              ();
  void  xCorrScore          (bool b);

private:

  //Data members
  double                binOffset;
  double                binSize;
  int                   charge;
  bool                  instrumentPrecursor;
  double                invBinSize;
  float                 maxIntensity;
  double                mz;
  vector<kPrecursor>*   precursor;
  float                 rTime;
  int                   scanNumber;
  int                   singletCount;
  
  
  vector<KTopPeps>*     singlets;
  vector<kSpecPoint>*   spec;
  kScoreCard            topHit[20];
  int                   xCorrArraySize;

  //Functions
  void BinIons      (kPreprocessStruct *pPre);
  void CometXCorr   ();
  void MakeCorrData (double *pdTempRawData, kPreprocessStruct *pPre, double scale);

  void kojakXCorr   ();

  //Utilities
  static int compareIntensity (const void *p1,const void *p2);
  static int compareMZ        (const void *p1,const void *p2);



};

#endif
