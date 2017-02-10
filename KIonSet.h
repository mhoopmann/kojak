#ifndef _KIONSET_H
#define _KIONSET_H

class KIonSet{
public:
  KIonSet(int sz, double m);
  KIonSet(const KIonSet& k);
  ~KIonSet();

  KIonSet& operator=(const KIonSet& k);

  double** aIons;
  double** bIons;
  double** cIons;
  double** xIons;
  double** yIons;
  double** zIons;
  double* mods;
  double  mass;
  double  difMass;
  int     len;
  bool modNTerm;
  bool modCTerm;
};

#endif