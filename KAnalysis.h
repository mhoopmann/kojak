#ifndef _KANALYSIS_H
#define _KANALYSIS_H

#include "KDB.h"
#include "KData.h"
#include "KIons.h"

class KAnalysis{
public:

  //Constructors & Destructors
  KAnalysis  (kParams& p);
  ~KAnalysis ();

  //User Functions
  bool analyzePeptides            (KDatabase& db, KData& d, bool crossLink);
  bool analyzePeptidesNonCovalent (KDatabase& db, KData& d);
  bool analyzeRelaxed             (KDatabase& db, KData& d);
  bool analyzeSinglets            (KDatabase& db, KData& d, kPeptide& pep, int index, double lowLinkMass, double highLinkMass);
  bool analyzeSingletsNoLysine    (KDatabase& db, KData& d, kPeptide& pep, int index, bool linkable);

  //__int64 xCorrCount;


private:

  //Private Functions
  int   findMass            (kSingletScoreCardPlus* s, double mass);
  void  scoreNCSpectra      (vector<int>& index, KData& d, double mass, bool linkable1, bool linkable2, int pep1, int pep2);
  void  scoreSingletSpectra (int index, KData& d, double mass, int len, int pep, int k, bool linkable, double minMass);
  void  scoreSpectra        (vector<int>& index, KData& d, double mass, vector<kPepMod>* v, bool linkable, int pep1, int pep2, int k1, int k2, int link);
  float simpleScoring       (KSpectrum& s);
  float xCorrScoring        (KSpectrum& s);
  float xCorrScoring2       (KSpectrum& s, double modMass);

  //Data Members
  KIons      ions;
  kParams    params;

  //Utilities
  static int compareD           (const void *p1,const void *p2);
  static int comparePeptideBMass(const void *p1,const void *p2);
  static int compareSSCPlus     (const void *p1,const void *p2);
  
};

#endif
