#ifndef _KIONSET_H
#define _KIONSET_H

typedef struct kISValue{
  double mz;
  int key;
  int pos;
  kISValue(){
    mz=0;
    key=0;
    pos=0;
  }
} kISValue;

class KIonSet{
public:
  KIonSet();
  KIonSet(const KIonSet& k);
  ~KIonSet();

  KIonSet& operator=(const KIonSet& k);

  void makeIndex(double binSize, double binOffset, bool a, bool b, bool c, bool x, bool y, bool z);
  void setIons(int sz, double m);

  kISValue** aIons;
  kISValue** bIons;
  kISValue** cIons;
  kISValue** xIons;
  kISValue** yIons;
  kISValue** zIons;
  double* mods;
  double  mass;
  double  difMass;
  double nTermMass;
  double cTermMass;
  int     len;
  bool    index;

private:
  void freeMem();
};

#endif