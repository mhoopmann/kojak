#ifndef _KIONSET_H
#define _KIONSET_H

class KIonSet{
public:
  KIonSet(int sz, double m);
  KIonSet(const KIonSet& k);
  ~KIonSet();

  KIonSet& operator=(const KIonSet& k);

  double* bIons[6];
  double* yIons[6];
  double* mods;
  double  mass;
  double  difMass;
  int     len;
};

#endif