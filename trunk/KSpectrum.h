#ifndef _KSPECTRUM_H
#define _KSPECTRUM_H

#include <cmath>
#include <vector>
#include "KStructs.h"

using namespace std;

typedef struct specPoint{
  double mass;
  float intensity;
} specPoint;

class KSpectrum {

public:

  //Constructors & Destructors
  KSpectrum(int i);
  KSpectrum(const KSpectrum& p);
  ~KSpectrum();

  //Operators
  KSpectrum&  operator=(const KSpectrum& p);
  specPoint&  operator[](const int& i);

  //Data Members
  kSparseMatrix*       xCorrSparseArray;
  
  //Accessors
  int                 getCharge           ();
  double              getInvBinSize       ();
  float               getMaxIntensity     ();
  double              getMZ               ();
  kPrecursor&         getPrecursor        (int i);
  float               getRTime            ();
  int                 getScanNumber       ();
  kScoreCard&         getScoreCard        (int i);
  kSingletScoreCard&  getSingletScoreCard (int i);
  int                 size                ();
  int                 sizePrecursor       ();

  //Modifiers
  void addPoint         (specPoint& s);
  void addPrecursor     (kPrecursor& p);
  void clear            ();
  void setCharge        (int i);
  void setMaxIntensity  (float f);
  void setMZ            (double d);
  void setPrecursor     (double d, int i);
  void setRTime         (float f);
  void setScanNumber    (int i);

  //Functions
  void  checkScore        (kScoreCard& s);
  void  checkSingletScore (kSingletScoreCard& s);
  void  sortMZ            ();
  void  xCorrScore        ();

private:

  //Data members
  double                binSize;
  int                   charge;
  double                invBinSize;
  float                 maxIntensity;
  double                mz;
  vector<kPrecursor>*   precursor;
  float                 rTime;
  int                   scanNumber;
  vector<specPoint>*    spec;
  kScoreCard            topHit[20];
  kSingletScoreCard*    topSinglet;
  int                   topCount;
  int                   xCorrArraySize;
  int                   xCorrSparseArraySize;

  //Functions
  void BinIons      (kPreprocessStruct *pPre);
  void CometXCorr   ();
  void MakeCorrData (double *pdTempRawData, kPreprocessStruct *pPre);

  //Utilities
  static int compareIntensity (const void *p1,const void *p2);
  static int compareMZ        (const void *p1,const void *p2);



};

#endif
