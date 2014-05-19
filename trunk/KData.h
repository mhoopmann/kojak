#ifndef _KDATA_H
#define _KDATA_H

#include "KDB.h"
#include "KPrecursor.h"
#include "KSpectrum.h"
#include "MSReader.h"
#include <iostream>

class KData {
public:

  KData();
  KData(kParams* p);
  ~KData();

  KSpectrum& operator[ ](const int& i);

  bool      getBoundaries     (double mass1, double mass2, vector<int>& index);
  bool      getBoundaries2    (double mass, double prec, vector<int>& index);
  kLinker&  getLink           (int i);
  double    getMaxMass        ();
  double    getMinMass        ();
  bool      mapPrecursors     ();
  bool      outputPercolator  (char* out, char* tag, KDatabase& db, int pLen);
  bool      outputResults     (char* out, KDatabase& db, bool bEnrich);
  void      readLinkers       (char* fn);
  bool      readSpectra2      ();
  void      setLinker         (double mass, int type);
  void      setVersion        (char* v);
  int       size              ();
  int       sizeLink          ();
  void      xCorr             ();

private:

  //Data Members
  char               version[32];
  vector<KSpectrum>  spec;
  vector<kLinker>    link;
  vector<kMass>      massList;
  kParams*           params;

  //Utilities
  void        centroid          (Spectrum& s, Spectrum& out, double resolution, int instrument=0);
  void        collapseSpectrum  (Spectrum& s);
  static int  compareInt        (const void *p1, const void *p2);
  static int  compareMassList   (const void *p1, const void *p2);
  int         getCharge         (Spectrum& s, int index, int next);
  double      polynomialBestFit (vector<double>& x, vector<double>& y, vector<double>& coeff, int degree=2);

};


#endif
