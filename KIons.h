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

#ifndef _KIONS_H
#define _KIONS_H

#include "KStructs.h"
#include "KIonSet.h"
#include <vector>

using namespace std;

typedef struct kModPos{
  int     pos;        //position of aa to add
  int     mod;        //index of mods at this position
  int     modCount;   //total mods prior to this position
  double  mass;       //total mass prior to this position
  double  modMass;    //total mod mass prior to this position
} kModPos;

typedef struct kModType{
  bool    xl;
  double  mass;
} kModType;

//max 10 mods per amino acid
typedef struct kMod{
  int       count;
  kModType  mod[10];
} kMod;

class KIons {
public:
  KIons();
  ~KIons();

  //Functions
  void      addFixedMod       (char mod, double mass);
  void      addMod            (char mod, bool xl, double mass);
  void      buildIons         ();
  void      buildLoopIons     (double linkMass, int link1, int link2);
  void      buildSingletIons  (int link);
  double    getAAMass         (char aa, bool n15=false);
  double    getFixedModMass   (char aa);
  void      modIonsRec        (int start, int link, int index, int depth, bool xl);
  void      modIonsRec2       (int start, int link, int index, int depth, bool xl);
  void      modLoopIonsRec    (int start, int link, int link2, int index, int depth, bool xl);
  void      modLoopIonsRec2   (int start, int link, int link2, int index, int depth, bool xl);
  void      reset             ();

  //Accessors
  KIonSet&  operator[ ]   (const int& i);
  KIonSet*  at            (const int& i);
  int       getIonCount   ();
  double    getModMass    (int index);
  int       getModMassSize();
  double*   getMods       ();
  void      getPeptide    (bool bPepOne, char* seq);
  int       getPeptideLen ();
  void      getPeptideMods(vector<kPepMod>& v);
  int       size          ();

  //Modifiers
  void  setAAMass       (char aa, double mass, bool n15=false);
  void  setMaxModCount  (int i);
  void  setModFlags     (bool monoMods, bool difMods);
  void  setPeptide      (bool bPepOne, char* seq, int len, double mass, bool nTerm, bool cTerm, bool n15);

  //Data Members
  double* modList;
  bool site[128]; //possible sites of linkage based on parameters
  //bool linkN;
  //bool linkC;

private:

  void addModIonSet(int index, char aa, int pos, int modIndex, int loopPos=-1);
  void buildSeries(int setNum);
  void clearSeries();
  
  double  aaMass[128];
  double  aaMassn15[128];
  double  aaFixedModMass[128];
  kMod    aaMod[128];   //inefficient memory usage, but not by much in the grand scheme.
  double  protFixedModMassC;
  double  protFixedModMassN;
  kMod  protModMassC;
  kMod  protModMassN;

  bool  monoModsOnXL;
  bool  diffModsOnXL;

  double  pep1Mass;
  double  pep2Mass;
  int     modIndex;

  int ionCount;
  int maxModCount;
  int pep1Len;
  int pep2Len;

  char* pep1;
  char* pep2;

  bool n15Pep1; //if peptide is of the 15N-labeled variety;
  bool n15Pep2;

  bool nPep1; //peptide has protein n-terminus
  bool nPep2; //peptide has protein n-terminus
  bool cPep1; //peptide has protein c-terminus
  bool cPep2; //peptide has protein c-terminus

  vector<kModPos> modQueue;
  vector<double>  modMassArray;
  vector<KIonSet> sets;

  //Utilities
  static int compareD(const void *p1,const void *p2);

};


#endif