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

#ifndef _KANALYSIS_H
#define _KANALYSIS_H

#include "KDB.h"
#include "KData.h"
#include "KLog.h"
#include "KIons.h"
#include "Threading.h"
#include "ThreadPool.h"

//=============================
// Structures for threading
//=============================
struct kAnalysisStruct {
  bool*       bKIonsMem;    //Pointer to the memory manager array to mark memory is in use
  Mutex*      mutex;        //Pointer to a mutex for protecting memory
  kPeptide*   pep;
  int         pepIndex;
  kAnalysisStruct(Mutex* m, kPeptide* p, int i){
    mutex=m;
    pep=p;
    pepIndex=i;
  }
  ~kAnalysisStruct(){
    //Mark that memory is not being used, but do not delete it here.
    Threading::LockMutex(*mutex);
    if(bKIonsMem!=NULL) *bKIonsMem=false;
    bKIonsMem=NULL;
    Threading::UnlockMutex(*mutex);
    mutex=NULL;   //release mutex
    pep=NULL;
  }
};

struct kAnalysisNCStruct {
  bool*               bKIonsMem;    //Pointer to the memory manager array to mark memory is in use
  Mutex*              mutex;        //Pointer to a mutex for protecting memory
  std::vector<kPeptideB>*  pep;
  int                 pepIndex;
  kAnalysisNCStruct(Mutex* m, std::vector<kPeptideB>* p, int i){
    mutex=m;
    pep=p;
    pepIndex=i;
  }
  ~kAnalysisNCStruct(){
    //Mark that memory is not being used, but do not delete it here.
    Threading::LockMutex(*mutex);
    if(bKIonsMem!=NULL) *bKIonsMem=false;
    bKIonsMem=NULL;
    Threading::UnlockMutex(*mutex);
    mutex=NULL;   //release mutex
    pep=NULL;
  }
};

struct kAnalysisRelStruct {
  bool*       bKIonsMem;    //Pointer to the memory manager array to mark memory is in use
  Mutex*      mutex;        //Pointer to a mutex for protecting memory
  KSpectrum*  spec;
  kAnalysisRelStruct(Mutex* m, KSpectrum* s){
    mutex=m;
    spec=s;
  }
  ~kAnalysisRelStruct(){
    //Mark that memory is not being used, but do not delete it here.
    Threading::LockMutex(*mutex);
    if(bKIonsMem!=NULL) *bKIonsMem=false;
    bKIonsMem=NULL;
    Threading::UnlockMutex(*mutex);
    mutex=NULL;   //release mutex
    spec=NULL;
  }
};

class KAnalysis{
public:

  //Constructors & Destructors
  KAnalysis  (kParams& p, KDatabase* d, KData* dat);
  ~KAnalysis ();

  //Master Functions
  bool doPeptideAnalysis ();
  bool doEValueAnalysis  ();

  void setLog(KLog* c);
  //bool doPeptideAnalysisNC ();
  //__int64 xCorrCount;

private:

  //Thread-start functions
  static void analyzePeptideProc (kAnalysisStruct* s); 
  static void analyzeEValueProc  (KSpectrum* s);

  //Analysis functions
  static bool analyzePeptide(kPeptide* p, int pepIndex, int iIndex);

  //Private Functions
  bool         allocateMemory          (int threads);
  static bool  analyzeSinglets         (kPeptide& pep, int index, double lowLinkMass, double highLinkMass, int iIndex);
  static bool  analyzeSingletsNC       (kPeptide& pep, int index, int iIndex);
  static void  checkXLMotif            (int motifA, char* motifB, std::vector<int>& v);
  void         deallocateMemory        (int threads);
  static int   findMass                (kSingletScoreCardPlus* s, int sz, double mass);
  static void  scoreSingletSpectra     (int index, int sIndex, double mass, int len, int pep, char k, double minMass, int iIndex);
  static void  scoreSpectra            (std::vector<int>& index, int sIndex, double modMass, int pep1, int pep2, int k1, int k2, int link, int iIndex, char linkSite1, char linkSite2);
  static float kojakScoring            (int specIndex, double modMass, int sIndex, int iIndex, int& match, int& conFrag, int z = 0);
  static void  setBinList              (kMatchSet* m, int iIndex, int charge, double preMass, kPepMod* mods, char modLen);

  //Data Members
  static bool*      bKIonsManager;
  static KDatabase* db;
  static double     highLinkMass;
  static KIons*     ions;
  static double     lowLinkMass;
  static double     maxMass;
  static double     minMass;
  static kParams    params;
  static KData*     spec;
  static char**     xlTable;
  static bool**     scanBuffer;

  static int        numIonSeries;

  static int* pepMassSize;
  static double**   pepMass;
  static bool** pepBin;
  static int* pepBinSize;
  static void makePepLists();
  static bool findCompMass(int motif, double low, double high);
  static int skipCount;
  static int nonSkipCount;

  static bool* soloLoop;
  static bool firstPass;

  static KDecoys decoys;
  static KLog* klog;

  static bool scoreSingletSpectra2(int index, int sIndex, double mass, double xlMass, int counterMotif, int len, int pep, char k, double minMass, int iIndex, char linkSite, int linkIndex);

  static Mutex  mutexKIonsManager; 
  static Mutex* mutexSpecScore; //these signal PSM list reads/additions/deleteions
  static Mutex** mutexSingletScore; //these signal singlet list reads/additions/deletions

  //Utilities
  static int compareD           (const void *p1,const void *p2);
  static int comparePeptideBMass(const void *p1,const void *p2);
  static int compareSSCPlus     (const void *p1,const void *p2);
  
};

#endif
