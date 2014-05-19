#ifndef _KIONS_H
#define _KIONS_H

#include "KStructs.h"
#include <vector>

using namespace std;

typedef struct kModPos{
  int     pos;        //position of aa to add
  int     mod;        //index of mods at this position
  int     modCount;   //total mods prior to this position
  double  mass;       //total mass prior to this position
  double  modMass;    //total mod mass prior to this position
} kModPos;

//max 10 mods per amino acid
typedef struct kMod{
  int     count;
  double  mass[10];
} kMod;


class KIons {
public:
  KIons();
  ~KIons();

  //Functions
  void    addFixedMod       (char mod, double mass);
  void    addMod            (char mod, double mass);
  void    buildIons         (double linkMass, int link);
  void    buildLoopIons     (double linkMass, int link1, int link2);
  double  buildModIons      (double linkMass, int link);
  void    buildNCIons       ();
  void    buildSingletIons  (int link);
  void    buildXIons        (double linkMass, int link1, int link2);

  //Accessors
  int       getIonCount   ();
  double    getModMass    (int index);
  int       getModMassSize();
  double*   getMods       ();
  void      getPeptide    (bool bPepOne, char* seq);
  void      getPeptideMods(vector<kPepMod>& v);

  //Modifiers
  void  setMaxModCount  (int i);
  void  setPeptide      (bool bPepOne, char* seq, int len, double mass);

  //Data Members
  double  bIons[6][512];
  double  yIons[6][512];
  double* modList;

private:

  void modMasses(double mass, int index, int count);

  
  double    tyIons[512];
  
  double  aaMass[128];
  kMod    aaMod[128];   //inefficient memory usage, but not by much in the grand scheme.

  double  pep1Mass;
  double  pep2Mass;
  int     modIndex;

  int ionCount;
  int maxModCount;
  int pep1Len;
  int pep2Len;

  char* pep1;
  char* pep2;

  vector<kModPos> modQueue;
  vector<double>  modMassArray;

  //Utilities
  static int compareD(const void *p1,const void *p2);

};

#endif
