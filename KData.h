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

#ifndef _KDATA_H
#define _KDATA_H

#include "KDB.h"
#include "KIons.h"
#include "KParams.h"
#include "KPrecursor.h"
#include "KSpectrum.h"
#include "MSReader.h"
#include "mzIMLTools.h"
#include "pepXMLWriter.h"
#include <iostream>

/*
#ifdef _MSC_VER
#include <direct.h>
#define getcwd _getcwd
#define slashdir '\\'
#else
#define slashdir '/'
#endif
*/

class KData {
public:

  KData();
  KData(kParams* p);
  ~KData();

  KSpectrum& operator[ ](const int& i);
  KSpectrum& at(const int& i);
  KSpectrum* getSpectrum(const int& i);

  void      buildXLTable      ();
  bool      checkLink         (char p1Site, char p2Site, int linkIndex);
  void      diagSinglet       ();
  bool      getBoundaries     (double mass1, double mass2, vector<int>& index, bool* buffer);
  bool      getBoundaries2    (double mass, double prec, vector<int>& index, bool* buffer);
  int       getCounterMotif   (int motifIndex, int counterIndex);
  kLinker&  getLink           (int i);
  double    getMaxMass        ();
  double    getMinMass        ();
  kXLMotif& getMotif          (int motifIndex);
  int       getMotifCount     ();
  int       getXLIndex        (int motifIndex, int xlIndex);
  char**    getXLTable        ();
  bool      mapPrecursors     ();
  void      outputDiagnostics (FILE* f, KSpectrum& s, KDatabase& db);
  bool      outputIntermediate(KDatabase& db);
  bool      outputMzID        (CMzIdentML& m, KDatabase& db, KParams& par, kResults& r);
  bool      outputPepXML      (PXWSpectrumQuery& p, KDatabase& db, kResults& r);
  bool      outputPercolator  (FILE* f, KDatabase& db, kResults& r, int count);
  bool      outputResults     (KDatabase& db, KParams& par);
  void      readLinkers       (char* fn);
  bool      readSpectra       ();
  void      setLinker         (kLinker x);
  void      setVersion        (const char* v);
  int       size              ();
  int       sizeLink          ();
  void      xCorr             (bool b);

private:

  //Data Members
  bool* bScans;
  char               version[32];
  char**             xlTable;
  vector<KSpectrum>  spec;
  vector<kLinker>    link;  //just cross-links, not mono-links
  vector<kMass>      massList;
  kParams*           params;
  KIons              aa;
  kXLMotif           motifs[20]; //lets put a cap on this for now
  int                motifCount;
  kXLTarget          xlTargets[128][5]; //capping analysis at 5 crosslinkers for now

  //Utilities
  void        centroid          (Spectrum& s, Spectrum& out, double resolution, int instrument=0);
  void        collapseSpectrum  (Spectrum& s);
  static int  compareInt        (const void *p1, const void *p2);
  static int  compareMassList   (const void *p1, const void *p2);
  int         getCharge         (Spectrum& s, int index, int next);
  double      polynomialBestFit (vector<double>& x, vector<double>& y, vector<double>& coeff, int degree=2);
  bool        processPath       (const char* in_path, char* out_path);
  string      processPeptide    (kPeptide& pep, vector<kPepMod>* mod, KDatabase& db);

};


#endif
