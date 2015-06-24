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

#ifndef _KPRECURSOR_H
#define _KPRECURSOR_H

#include "KStructs.h"
#include "KSpectrum.h"

#include "CHardklor2.h"
#include "CModelLibrary.h"
#include "CHardklorSetting.h"
#include "CHardklorVariant.h"
#include "MSObject.h"
#include "MSReader.h"
#include <cmath>
#include <deque>

#define GAUSSCONST 5.5451774444795623

typedef struct kRTProfile{
  int   scan;
  float rt;
} kRTProfile;

typedef struct kScanBin{
  int   index;
  float intensity;
} kScanBin;

class KPrecursor {
public:

  KPrecursor(kParams* p);
  ~KPrecursor();

  bool  estimatePrecursor (KSpectrum& s);
  int   getSpecRange      (KSpectrum& pls);
  bool  setFile           ();

  vector<int> preCharges;

private:

  //Functions
  void    averageScansCentroid  (vector<Spectrum*>& s, Spectrum& avg, double min, double max);
  void    centroid              (Spectrum& s, Spectrum& out, double resolution, int instrument=0);
  int     findPeak              (Spectrum& s, double mass);
  int     findPeak              (Spectrum& s, double mass, double prec);
  double  polynomialBestFit     (vector<double>& x, vector<double>& y, vector<double>& coeff, int degree=2);

  //Data Members
  char      fileName[256];
  int       lastScan;
  MSObject  mso;
  MSReader  msr;
  MSReader  msw;
  Spectrum  avg;
  Spectrum  cent;
  Spectrum  spec;
  
  CAveragine*      averagine;
	CMercury8*       mercury;
  CModelLibrary*   models;
  CHardklor2*      h;
  CHardklor*       hO;
  CHardklorSetting hs;
  CHardklorSetting hs2;
  CHardklorSetting hs4;
  CHardklorVariant hv;

  deque<Spectrum>*  centBuf;

  kParams*   params;

  //Utilities
  static int    comparePLPScanNum (const void *p1, const void *p2);
  static int    compareScanBinRev (const void *p1, const void *p2);

};

#endif
