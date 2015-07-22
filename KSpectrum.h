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
#include <vector>
#include "KStructs.h"

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
  int                 size                  ();
  int                 sizePrecursor         ();

  //Modifiers
  void addPoint               (kSpecPoint& s);
  void addPrecursor           (kPrecursor& p);
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
  void  checkScore        (kScoreCard& s);
  void  checkSingletScore (kSingletScoreCard& s);
  void  sortMZ            ();
  void  xCorrScore        (bool b);

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
  kSingletScoreCard*    singletFirst;   //pointer to start of linked list
  kSingletScoreCard*    singletLast;    //pointer to end of linked list
  int                   singletMax;
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
