#ifndef _KPRECURSOR_H
#define _KPRECURSOR_H

#include "KStructs.h"
#include "KSpectrum.h"

#include "CHardklor.h"
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


class KPrecursor {
public:

  KPrecursor(kParams* p);
  ~KPrecursor();

  bool  estimatePrecursor (KSpectrum& s);
  bool  getSpecRange      (KSpectrum& s);
  bool  setFile           ();

private:

  //Functions
  void    averageScansBin       (vector<Spectrum*>& s, Spectrum& avg, double min, double max);
  void    averageScansCentroid  (vector<Spectrum*>& s, Spectrum& avg, double min, double max);
  void    centroid              (Spectrum& s, Spectrum& out, double resolution, int instrument=0);
  int     findPeak              (Spectrum& s, double mass);
  int     findPeak              (Spectrum& s, double mass, double prec);
  double  polynomialBestFit     (vector<double>& x, vector<double>& y, vector<double>& coeff, int degree=2);

  //Data Members
  char      fileName[256];
  MSObject  mso;
  MSReader  msr;
  MSReader  msw;
  Spectrum  avg;
  Spectrum  cent;
  Spectrum  spec;
  
  CAveragine*      averagine;
	CMercury8*       mercury;
  CHardklor*       h;
  CHardklorSetting hs;
  CHardklorSetting hs2;
  CHardklorSetting hs4;
  CHardklorVariant hv;

  deque<Spectrum>*  scanBuf;
  deque<Spectrum>*  centBuf;

  kParams*   params;

  //Utilities
  static int    comparePLPScanNum(const void *p1, const void *p2);

};

#endif
