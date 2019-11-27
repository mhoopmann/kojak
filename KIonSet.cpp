#include "KIonSet.h"
#include <iostream>
KIonSet::KIonSet(){
  int i,j;
  len=0;  //zero length arrays
  mass=0;
  difMass=0;
  aIons = new kISValue*[4];
  bIons = new kISValue*[4];
  cIons = new kISValue*[4];
  xIons = new kISValue*[4];
  yIons = new kISValue*[4];
  zIons = new kISValue*[4];
  for(j=0;j<4;j++){
    aIons[j] = new kISValue[len];
    bIons[j] = new kISValue[len];
    cIons[j] = new kISValue[len];
    xIons[j] = new kISValue[len];
    yIons[j] = new kISValue[len];
    zIons[j] = new kISValue[len];
    /*
    for(i=0;i<len;i++){
      aIons[j][i]=0;
      bIons[j][i]=0;
      cIons[j][i]=0;
      xIons[j][i]=0;
      yIons[j][i]=0;
      zIons[j][i]=0;
    }
    */
  }
  mods=new double[len];
  for(i=0;i<len;i++) mods[i]=0;
  nTermMass=0;
  cTermMass=0;
  index=false;
}

KIonSet::KIonSet(const KIonSet& k){
  int i,j;
  len=k.len;
  mass=k.mass;
  difMass=k.difMass;
  aIons = new kISValue*[4];
  bIons = new kISValue*[4];
  cIons = new kISValue*[4];
  xIons = new kISValue*[4];
  yIons = new kISValue*[4];
  zIons = new kISValue*[4];
  for(j=0;j<4;j++){
    aIons[j] = new kISValue[len];
    bIons[j] = new kISValue[len];
    cIons[j] = new kISValue[len];
    xIons[j] = new kISValue[len];
    yIons[j] = new kISValue[len];
    zIons[j] = new kISValue[len];
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
  nTermMass=k.nTermMass;
  cTermMass=k.cTermMass;
  index=k.index;
}
  
KIonSet::~KIonSet(){
  freeMem();
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
      aIons[j] = new kISValue[len];
      bIons[j] = new kISValue[len];
      cIons[j] = new kISValue[len];
      xIons[j] = new kISValue[len];
      yIons[j] = new kISValue[len];
      zIons[j] = new kISValue[len];
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
    nTermMass = k.nTermMass;
    cTermMass = k.cTermMass;
    index=k.index;
  }
  return *this;
}

void KIonSet::freeMem(){
  for (int j = 0; j<4; j++){
    delete[] aIons[j];
    delete[] bIons[j];
    delete[] cIons[j];
    delete[] xIons[j];
    delete[] yIons[j];
    delete[] zIons[j];
  }
  delete[] aIons;
  delete[] bIons;
  delete[] cIons;
  delete[] xIons;
  delete[] yIons;
  delete[] zIons;
  delete[] mods;
}

void KIonSet::makeIndex(double binSize, double binOffset, bool a, bool b, bool c, bool x, bool y, bool z){
  int i,j;
  double mz;
  double invBinSize = 1.0/binSize;
  for(i=0;i<len;i++){
    for (j = 1; j<4; j++){
      if(a && aIons[j][i].mz>0){
        mz = binSize * (int)(aIons[j][i].mz * invBinSize + binOffset);
        aIons[j][i].key = (int)mz;
        aIons[j][i].pos = (int)((mz - aIons[j][i].key)*invBinSize);
      }
      if (b && bIons[j][i].mz>0){
        mz = binSize * (int)(bIons[j][i].mz * invBinSize + binOffset);
        bIons[j][i].key = (int)mz;
        bIons[j][i].pos = (int)((mz - bIons[j][i].key)*invBinSize);
      }
      if (c && cIons[j][i].mz>0){
        mz = binSize * (int)(cIons[j][i].mz * invBinSize + binOffset);
        cIons[j][i].key = (int)mz;
        cIons[j][i].pos = (int)((mz - cIons[j][i].key)*invBinSize);
      }

      if(x && xIons[j][i].mz>0){
        mz = binSize * (int)(xIons[j][i].mz * invBinSize + binOffset);
        xIons[j][i].key = (int)mz;
        xIons[j][i].pos = (int)((mz - xIons[j][i].key)*invBinSize);
      }
      if (y && yIons[j][i].mz>0){
        mz = binSize * (int)(yIons[j][i].mz * invBinSize + binOffset);
        yIons[j][i].key = (int)mz;
        yIons[j][i].pos = (int)((mz - yIons[j][i].key)*invBinSize);
      }
      if (z && zIons[j][i].mz>0){
        mz = binSize * (int)(zIons[j][i].mz * invBinSize + binOffset);
        zIons[j][i].key = (int)mz;
        zIons[j][i].pos = (int)((mz - zIons[j][i].key)*invBinSize);
      }
    }
  }
  index=true;
}

void KIonSet::setIons(int sz, double m){
  freeMem();

  int i, j;
  len = sz;
  mass = m;
  difMass = 0;
  aIons = new kISValue*[4];
  bIons = new kISValue*[4];
  cIons = new kISValue*[4];
  xIons = new kISValue*[4];
  yIons = new kISValue*[4];
  zIons = new kISValue*[4];
  for (j = 0; j<4; j++){
    aIons[j] = new kISValue[len];
    bIons[j] = new kISValue[len];
    cIons[j] = new kISValue[len];
    xIons[j] = new kISValue[len];
    yIons[j] = new kISValue[len];
    zIons[j] = new kISValue[len];
  }
  mods = new double[len];
  for (i = 0; i<len; i++) mods[i] = 0;
  nTermMass = 0;
  cTermMass = 0;
}