#include "KIonSet.h"

KIonSet::KIonSet(int sz, double m){
  int i,j;
  len=sz;
  mass=m;
  difMass=0;
  aIons=new double*[4];
  bIons=new double*[4];
  cIons=new double*[4];
  xIons=new double*[4];
  yIons=new double*[4];
  zIons=new double*[4];
  for(j=0;j<4;j++){
    aIons[j]=new double[len];
    bIons[j]=new double[len];
    cIons[j]=new double[len];
    xIons[j]=new double[len];
    yIons[j]=new double[len];
    zIons[j]=new double[len];
    for(i=0;i<len;i++){
      aIons[j][i]=0;
      bIons[j][i]=0;
      cIons[j][i]=0;
      xIons[j][i]=0;
      yIons[j][i]=0;
      zIons[j][i]=0;
    }
  }
  mods=new double[len];
  for(i=0;i<len;i++) mods[i]=0;
  modNTerm=false;
  modCTerm=false;
}

KIonSet::KIonSet(const KIonSet& k){
  int i,j;
  len=k.len;
  mass=k.mass;
  difMass=k.difMass;
  aIons=new double*[4];
  bIons=new double*[4];
  cIons=new double*[4];
  xIons=new double*[4];
  yIons=new double*[4];
  zIons=new double*[4];
  for(j=0;j<4;j++){
    aIons[j]=new double[len];
    bIons[j]=new double[len];
    cIons[j]=new double[len];
    xIons[j]=new double[len];
    yIons[j]=new double[len];
    zIons[j]=new double[len];
    for(i=0;i<len;i++){
      aIons[j][i]=k.aIons[j][i];
      bIons[j][i]=k.bIons[j][i];
      cIons[j][i]=k.cIons[j][i];
      xIons[j][i]=k.xIons[j][i];
      yIons[j][i]=k.yIons[j][i];
      zIons[j][i]=k.zIons[j][i];
    }
  }
  mods=new double[len];
  for(i=0;i<len;i++) mods[i]=k.mods[i];
  modNTerm=k.modNTerm;
  modCTerm=k.modCTerm;
}
  
KIonSet::~KIonSet(){
  for(int j=0;j<4;j++){
    delete [] aIons[j];
    delete [] bIons[j];
    delete [] cIons[j];
    delete [] xIons[j];
    delete [] yIons[j];
    delete [] zIons[j];
  }
  delete [] aIons;
  delete [] bIons;
  delete [] cIons;
  delete [] xIons;
  delete [] yIons;
  delete [] zIons;
  delete [] mods;
}
  
KIonSet& KIonSet::operator=(const KIonSet& k){
  if(this!=&k){
    int i,j;
    len=k.len;
    mass=k.mass;
    difMass=k.difMass;
    
    for(j=0;j<4;j++){
      delete [] aIons[j];
      delete [] bIons[j];
      delete [] cIons[j];
      delete [] xIons[j];
      delete [] yIons[j];
      delete [] zIons[j];
      aIons[j] = new double[len];
      bIons[j] = new double[len];
      cIons[j] = new double[len];
      xIons[j] = new double[len];
      yIons[j] = new double[len];
      zIons[j] = new double[len];
      for(i=0;i<len;i++){
        aIons[j][i]=k.aIons[j][i];
        bIons[j][i]=k.bIons[j][i];
        cIons[j][i]=k.cIons[j][i];
        xIons[j][i]=k.xIons[j][i];
        yIons[j][i]=k.yIons[j][i];
        zIons[j][i]=k.zIons[j][i];
      }
    }

    delete [] mods;
    mods=new double[len];
    for(i=0;i<len;i++) mods[i]=k.mods[i];
    modNTerm = k.modNTerm;
    modCTerm = k.modCTerm;
  }
  return *this;
}