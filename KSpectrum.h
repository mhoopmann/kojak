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

#include <algorithm>
#include <cmath>
#include <list>
#include <vector>
#include "KDB.h"
#include "KStructs.h"
#include "KTopPeps.h"
#include "CometDecoys.h"

#define HISTOSZ 152

//a small structure for defining which ions (a,b,c,x,y,z) to compute for decoys
typedef struct sDecoyIons{
  bool b;
  double mass;
} sDecoyIons;

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
  //int tmpSingCount;
  int tmpHistCount;
  int peakCounts;


  //Accessors
  double              getBinOffset          ();
  int                 getCharge             ();
  bool                getInstrumentPrecursor();
  double              getInvBinSize         ();
  float               getMaxIntensity       ();
  double              getMZ                 ();
  std::string         getNativeID           ();
  kPrecursor&         getPrecursor          (int i);
  kPrecursor*         getPrecursor2         (int i);
  float               getRTime              ();
  int                 getScanNumber         ();
  kScoreCard&         getScoreCard          (int i);
  int                 getSingletCount       ();
  KTopPeps*           getTopPeps            (int index);
  int                 size                  ();
  int                 sizePrecursor         ();
  
  int                   singletMax;

  int histogram[HISTOSZ];
  int histogramCount;
  int histoMaxIndex;
  double dCummulative[HISTOSZ];
  int tempHistogram[HISTOSZ];

  //For diagnostics only
  int histogramO[HISTOSZ];

  //Modifiers
  void addPoint               (kSpecPoint& s);
  void addPrecursor           (kPrecursor& p, int sz);
  void clear                  ();
  void erasePrecursor         (int i);
  void setCharge              (int i);
  void setInstrumentPrecursor (bool b);
  void setMaxIntensity        (float f);
  void setMZ                  (double d);
  void setNativeID            (std::string s);
  void setPrecursor           (double d, int i);
  void setRTime               (float f);
  void setScanNumber          (int i);

  //Functions
  bool  calcEValue          (kParams* params, KDecoys& decoys, KDatabase& db);
  void  clearPrecursors     ();
  bool  checkDecoy          (KDatabase& db, std::string& dStr, kScoreCard& hit);  //To be run AFTER analysis completes.
  void  checkScore          (kScoreCard& s, ePSMList listID=listSingle);
  double  generateSingletDecoys2(kParams* params, KDecoys& decoys, double xcorr, double mass, int preIndex);
  bool  generateXcorrDecoys (kParams* params, KDecoys& decoys);
  void  linearRegression2   (double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared);
  void  linearRegression4   (int* histo, double& slope, double& intercept, int&  iMaxXcorr, int& iStartXcorr, int& iNextXcorr, double& rSquared);
  void  refreshScore        (KDatabase& db, std::string& dStr);  //To be run AFTER analysis completes. Looks at top scores, if a tie, make sure decoys are listed second (to help TPP analysis)
  void  sortMZ              ();
  void  sortIntensityRev    () { sort(spec->begin(),spec->end(),compareIntensityRev); }
  void  kojakXCorr          (double* pdTempRawData, double* pdTmpFastXcorrData, float* pfFastXcorrData, kPreprocessStruct*& pPre);

  kScoreCard            topSingle;
  kScoreCard            topLoop;
  kScoreCard            topXL;

private:

  //Data members
  double                binOffset;
  double                binSize;
  int                   charge;
  bool                  instrumentPrecursor;
  double                invBinSize;
  float                 maxIntensity;
  double                mz;
  std::string           nativeID;
  std::vector<kPrecursor>*   precursor;
  float                 rTime;
  int                   scanNumber;
  int                   singletCount;

  int decoyIonSz;
  sDecoyIons decoyIons[6];
  
  std::vector<KTopPeps>*     singlets;
  std::vector<kSpecPoint>*   spec;
  kScoreCard            topHit[20];

  int                   xCorrArraySize;

  //Functions
  void BinIons        (kPreprocessStruct *pPre);
  void MakeCorrData   (double *pdTempRawData, kPreprocessStruct *pPre, double scale);

  //Utilities
  static int compareIntensity (const void *p1,const void *p2);
  static bool compareIntensityRev (const kSpecPoint& p1, const kSpecPoint& p2) { return p1.intensity>p2.intensity;}
  static int compareMZ        (const void *p1,const void *p2);

};

#endif
