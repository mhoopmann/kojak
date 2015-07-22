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
  bool        crossLink;
  kAnalysisStruct(Mutex* m, kPeptide* p, int i, bool b){
    mutex=m;
    pep=p;
    pepIndex=i;
    crossLink=b;
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
  vector<kPeptideB>*  pep;
  int                 pepIndex;
  kAnalysisNCStruct(Mutex* m, vector<kPeptideB>* p, int i){
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
  KSpectrum*  spec;
  kAnalysisRelStruct(KSpectrum* s){
    spec=s;
  }
  ~kAnalysisRelStruct(){
    spec=NULL;
  }
};

class KAnalysis{
public:

  //Constructors & Destructors
  KAnalysis  (kParams& p, KDatabase* d, KData* dat);
  ~KAnalysis ();

  //Master Functions
  bool doPeptideAnalysis   (bool crossLink);
  //bool doPeptideAnalysisNC ();
  bool doRelaxedAnalysis   ();
  
  //__int64 xCorrCount;

private:

  //Thread-start functions
  static void analyzePeptideProc(kAnalysisStruct* s); 
  //static void analyzePeptideNCProc(kAnalysisNCStruct* s);  
  static void analyzeRelaxedProc(kAnalysisRelStruct* s);

  //Analysis functions
  static bool analyzePeptide(kPeptide* p, int pepIndex, int iIndex, bool crossLink);
  //static void analyzePeptideNC(vector<kPeptideB>* p, int pIndex, int iIndex);
  static void analyzeRelaxed(KSpectrum* sp);

  //Private Functions
  bool         allocateMemory          (int threads);
  static bool  analyzeSinglets         (kPeptide& pep, int index, double lowLinkMass, double highLinkMass, int iIndex);
  static bool  analyzeSingletsNoLysine (kPeptide& pep, int sIndex, int index, bool linkable, int iIndex);
  void         deallocateMemory        ();
  static int   findMass                (kSingletScoreCardPlus* s, int sz, double mass);
  //static void  scoreNCSpectra          (vector<int>& index, double mass, bool linkable1, bool linkable2, int pep1, int pep2, int iIndex);
  static void  scoreSingletSpectra     (int index, int sIndex, double mass, int len, int pep, char k, bool linkable, double minMass, int iIndex);
  static void  scoreSpectra            (vector<int>& index, int sIndex, double modMass, bool linkable, int pep1, int pep2, int k1, int k2, int link, int iIndex);
  static float xCorrScoring            (KSpectrum& s, double modMass, int sIndex, int iIndex);
  static float kojakScoring            (int specIndex, double modMass, int sIndex, int iIndex);

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

  static int        numIonSeries;

  static Mutex  mutexKIonsManager; 
  static Mutex* mutexSpecScore; 

  //Utilities
  static int compareD           (const void *p1,const void *p2);
  static int comparePeptideBMass(const void *p1,const void *p2);
  static int compareSSCPlus     (const void *p1,const void *p2);
  
};

#endif
