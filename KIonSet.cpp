#include "KIonSet.h"

KIonSet::KIonSet(int sz, double m){
  int i,j;
  len=sz;
  mass=m;
  difMass=0;
  for(j=0;j<6;j++){
    bIons[j]=new double[len];
    yIons[j]=new double[len];
    for(i=0;i<len;i++){
      bIons[j][i]=0;
      yIons[j][i]=0;
    }
  }
  mods=new double[len];
  for(i=0;i<len;i++) mods[i]=0;
}

KIonSet::KIonSet(const KIonSet& k){
  int i,j;
  len=k.len;
  mass=k.mass;
  difMass=k.difMass;
  for(j=0;j<6;j++){
    bIons[j]=new double[len];
    yIons[j]=new double[len];
    for(i=0;i<len;i++){
      bIons[j][i]=k.bIons[j][i];
      yIons[j][i]=k.yIons[j][i];
    }
  }
  mods=new double[len];
  for(i=0;i<len;i++) mods[i]=k.mods[i];
}
  
KIonSet::~KIonSet(){
  for(int j=0;j<6;j++){
    delete [] bIons[j];
    delete [] yIons[j];
  }
  delete [] mods;
}
  
KIonSet& KIonSet::operator=(const KIonSet& k){
  if(this!=&k){
    int i,j;
    len=k.len;
    mass=k.mass;
    difMass=k.difMass;
    for(j=0;j<6;j++){
      delete [] bIons[j];
      delete [] yIons[j];
      bIons[j]=new double[len];
      yIons[j]=new double[len];
      for(i=0;i<len;i++){
        bIons[j][i]=k.bIons[j][i];
        yIons[j][i]=k.yIons[j][i];
      }
    }
    delete [] mods;
    mods=new double[len];
    for(i=0;i<len;i++) mods[i]=k.mods[i];
  }
  return *this;
}