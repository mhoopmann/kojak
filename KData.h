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
#include "KLog.h"
#include "KParams.h"
#include "KSpectrum.h"
#include "MSReader.h"
#include "mzIMLTools.h"
#include "pepXMLWriter.h"
#include "NeoPepXMLParser.h"
#include "Threading.h"
#include "ThreadPool.h"
#include <iostream>

#include "CHardklor2.h"
#include "CModelLibrary.h"
#include "CHardklorSetting.h"
#include "CHardklorVariant.h"

#define GAUSSCONST 5.5451774444795623

typedef struct kScanBin{
  int   index;
  float intensity;
} kScanBin;

//=============================
// Structures for threading
//=============================
typedef struct kMS2struct{
  MSToolkit::Spectrum* s;
  KSpectrum* pls;
  int state;
  bool thread;
  kMS2struct(MSToolkit::Spectrum* sp, int topCount, double binSize, double binOffset){
    s = sp;
    pls = new KSpectrum(topCount, binSize, binOffset);
    state = 0;
    thread = false;
  }
  ~kMS2struct(){
    delete s;
    pls = NULL;  //this needs to be deleted elsewhere
    //delete pls;
  }
} kMS2struct;

class KData {
public:

  KData();
  KData(kParams* p);
  ~KData();

  KSpectrum& operator[ ](const int& i);
  KSpectrum& at(const int& i);
  KSpectrum* getSpectrum(const int& i);

  void      addProteins       (void* sh, KDatabase& db, int pIndex, bool xl, int linkA, int linkB);
  void      addSearchScore    (CnpxSearchHit& sh, std::string name, double value, std::string fmt);
  void      addSearchScore    (CnpxSearchHit& sh, std::string name, int value);
  void      addXlinkScore     (CnpxLinkedPeptide& lp, std::string name, double value, std::string fmt);
  void      addXlinkScore     (CnpxLinkedPeptide& lp, std::string name, int value);
  void      addXlinkScore     (CnpxXLink& lp, std::string name, double value, std::string fmt);
  void      addXlinkScore     (CnpxXLink& lp, std::string name, int value);
  void      buildXLTable      ();
  bool      checkLink         (char p1Site, char p2Site, int linkIndex);
  void      diagSinglet       ();
  bool      getBoundaries     (double mass1, double mass2, std::vector<int>& index, bool* buffer);
  bool      getBoundaries2    (double mass, double prec, std::vector<int>& index, bool* buffer);
  int       getCounterMotif   (int motifIndex, int counterIndex);
  kLinker&  getLink           (int i);
  double    getMaxMass        ();
  double    getMinMass        ();
  kXLMotif& getMotif          (int motifIndex);
  int       getMotifCount     ();
  int       getXLIndex        (int motifIndex, int xlIndex);
  char**    getXLTable        ();
  CnpxModificationInfo  makeModificationInfo(std::vector<kPepMod>& mods, std::string peptide, bool n15, bool nTerm, bool cTerm);
  void      outputDiagnostics (FILE* f, KSpectrum& s, KDatabase& db);
  bool      outputIntermediate(KDatabase& db);
  bool      outputMzID        (CMzIdentML& m, KDatabase& db, KParams& par, kResults& r);
  bool      outputNeoPepXML   (CnpxSpectrumQuery& p, KDatabase& db, kResults& r);
  bool      outputPepXML      (PXWSpectrumQuery& p, KDatabase& db, kResults& r);
  bool      outputPercolator  (FILE* f, KDatabase& db, kResults& r, int count);
  bool      outputResults     (KDatabase& db, KParams& par);
  void      readLinkers       (char* fn);
  bool      readSpectra       ();
  void      setLinker         (kLinker x);
  void      setLog            (KLog* c);
  void      setVersion        (const char* v);
  int       size              ();
  int       sizeLink          ();


  void report();

private:

  //--Data Members --//
  //Variables for peptide spectrum analysis
  bool*                   bScans;
  char                    version[32];
  char**                  xlTable;
  std::vector<KSpectrum*> spec;
  std::vector<kLinker>    link;  //just cross-links, not mono-links
  std::vector<kMass>      massList;
  static kParams*         params;
  KIons                   aa;
  kXLMotif                motifs[20]; //lets put a cap on this for now
  int                     motifCount;
  kXLTarget               xlTargets[128][5]; //capping analysis at 5 crosslinkers for now
  static KLog*            klog;

  //Common memory to be shared by all threads during spectral processing
  static bool* memoryPool;          
  static double** tempRawData;
  static double** tmpFastXcorrData;
  static float**  fastXcorrData;
  static Mutex    mutexMemoryPool;
  static kPreprocessStruct** preProcess;

  //Optimized file loading structures for spectral processing
  static std::deque<MSToolkit::Spectrum*> dMS1;
  static std::vector<MSToolkit::Spectrum*> vMS1Buffer;
  static Mutex mutexLockMS1;
  static CHardklor2** h;
  static CHardklor**  hO;
  static CHardklorSetting hs;
  static CHardklorSetting hs2;
  static CHardklorSetting hs4;
  static Mutex* mutexHardklor;
  static CAveragine** averagine;
  static CMercury8** mercury;
  static CModelLibrary* models;
  static bool* bHardklor;
  static int maxPrecursorMass;

  //Utilities
  static void centroid(MSToolkit::Spectrum* s, KSpectrum* out, double resolution, int instrument = 0);
  static void collapseSpectrum(KSpectrum& s);
  static int  compareInt        (const void *p1, const void *p2);
  static int  compareMassList   (const void *p1, const void *p2);
  static int  getCharge(KSpectrum& s, int index, int next);
  static double polynomialBestFit (std::vector<double>& x, std::vector<double>& y, std::vector<double>& coeff, int degree=2);
  bool        processPath       (const char* in_path, char* out_path);
  std::string processPeptide    (kPeptide& pep, std::vector<kPepMod>* mod, KDatabase& db);
  void        processProtein    (int pepIndex, int site, char linkSite, std::string& prot, std::string& sites, bool& decoy, KDatabase& db);
  void        writeMzIDDatabase (CMzIdentML& m, KDatabase& db);
  bool        writeMzIDEnzyme   (pxwBasicXMLTag t, CEnzymes& e);
  void        writeMzIDPE       (CMzIdentML& m, CSpectrumIdentificationItem& m_sii, int pepID, KDatabase& db);
  std::string writeMzIDSIP      (CMzIdentML& m, std::string& sRef, KParams& par);

  //spectral processing functions
  static void averageScansCentroid(std::vector<MSToolkit::Spectrum*>& s, MSToolkit::Spectrum& avg, double min, double max);
  static int  findPeak(MSToolkit::Spectrum* s, double mass);
  static int  findPeak(MSToolkit::Spectrum* s, double mass, double prec);
  static void formatMS2(MSToolkit::Spectrum* s, KSpectrum* pls);
  void initHardklor();
  void memoryAllocate();
  void memoryFree();
  static void processMS2(kMS2struct* s);
  static int  processPrecursor(kMS2struct* s, int tIndex);
  void releaseHardklor();

  //static functions
  static bool compareSpecPoint(const kSpecPoint& p1, const kSpecPoint& p2){ return p1.mass<p2.mass; }
  static int compareScanBinRev2(const void *p1, const void *p2);

 
};


#endif
