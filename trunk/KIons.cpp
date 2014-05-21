/*
Copyright 2014, Michael R. Hoopmann, Institute for Systems Biology

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include "KIons.h"

KIons::KIons(){
  int i;
  for(i=0;i<128;i++){
    aaMass[i]=0;
    aaMod[i].count=0;
  }
  aaMass['A']=71.0371103;
  aaMass['C']=103.0091803;
  aaMass['D']=115.0269385;
  aaMass['E']=129.0425877;
  aaMass['F']=147.0684087;
  aaMass['G']=57.0214611;
  aaMass['H']=137.0589059;
  aaMass['I']=113.0840579;
  aaMass['K']=128.0949557;
  aaMass['L']=113.0840579;
  aaMass['M']=131.0404787;
  aaMass['N']=114.0429222;
  aaMass['P']=97.0527595;
  aaMass['Q']=128.0585714;
  aaMass['R']=156.1011021;
  aaMass['S']=87.0320244;
  aaMass['T']=101.0476736;
  aaMass['V']=99.0684087;
  aaMass['W']=186.0793065;
  aaMass['Y']=163.0633228;
  modList=NULL;
  pep1=NULL;
  pep2=NULL;
  maxModCount=0;
  ionCount=0;
}

KIons::~KIons(){
  if(modList!=NULL) delete [] modList;
  pep1=NULL;
  pep2=NULL;
}

void KIons::addFixedMod(char mod, double mass){
  aaMass[mod]+=mass;
}

void KIons::addMod(char mod, double mass){
  aaMod[mod].mass[aaMod[mod].count++]=mass;
}

void KIons::buildIons(double linkMass, int link){
	
  int b;
	int i;
  int j;
  int y;
	double mMass;
  double totalMass;
  double fragMass;

  //set up boundaries
  ionCount=pep1Len-1;
  b=0;
  y=ionCount-1;

  totalMass=pep1Mass+linkMass;

	//b- & y-ions from first peptide
  mMass=0;
	for(i=0;i<ionCount;i++){
    mMass+=aaMass[pep1[i]];
    if(i>=link) fragMass=mMass+linkMass;
    else fragMass = mMass;
    bIons[0][b] = fragMass;
    yIons[0][y] = (totalMass-fragMass);
    b++;
    y--;
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<6;j++){
      bIons[j][i] = ((bIons[0][i]+1.007276466*j)/j);
      yIons[j][i] = ((yIons[0][i]+1.007276466*j)/j);
    }
  }

}

void KIons::buildLoopIons(double linkMass, int link1, int link2){
	
  int b;
	int i;
  int j;
  int y;
  int len;
	double mMass;
  double totalMass;
  double fragMass;

  //set up boundaries
  len = pep1Len-1;
  ionCount=pep1Len-1-(link2-link1);
  b=0;
  y=ionCount-1;

  totalMass=pep1Mass+linkMass;

	//b- & y-ions from first peptide
  mMass=0;
	for(i=0;i<len;i++){
    mMass+=aaMass[pep1[i]];
    if(i>=link1 && i<link2) continue;
    else if(i>=link2) fragMass=mMass+linkMass;
    else fragMass = mMass;
    bIons[0][b] = fragMass;
    yIons[0][y] = (totalMass-fragMass);
    b++;
    y--;
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<6;j++){
      bIons[j][i] = ((bIons[0][i]+1.007276466*j)/j);
      yIons[j][i] = ((yIons[0][i]+1.007276466*j)/j);
    }
  }

}

double KIons::buildModIons(double linkMass, int link){

  int b;
  int i;
  int j;
  int y;
  int modCount;
	double mMass;
  double fragMass;
  double modMass;
  double totalMass;
  bool bMatch;
  kModPos mp;

  ionCount=pep1Len-1;
  totalMass=pep1Mass+linkMass;

  while(modQueue.size()>0){
    
    bMatch=true;

    mp=modQueue[modQueue.size()-1];
    modQueue.pop_back();

    modMass=mp.modMass;
    modCount=mp.modCount;
    mMass=mp.mass;

    if(mp.pos==-1) {
      modList[0]=0;
      if(aaMod['@'].count>mp.mod){
        modList[0]=aaMod['@'].mass[mp.mod];
        modMass+=aaMod['@'].mass[mp.mod++];
        modCount++;
        modQueue.push_back(mp);
      }
      mp.pos=0;
    }

    //set position for each ion series
    b=mp.pos;
    y=ionCount-mp.pos-1;

  	//calc b- and y-ion masses;
    while(b<ionCount){
      modList[b+1]=0;
      mp.mass=mMass;
      mp.modMass=modMass;
      mp.modCount=modCount;
      mMass+=aaMass[pep1[b]];
      if(b==mp.pos && aaMod[pep1[b]].count>mp.mod){
        modList[b+1]=aaMod[pep1[b]].mass[mp.mod];
        modMass+=aaMod[pep1[b]].mass[mp.mod++];
        mp.pos=b;
        modQueue.push_back(mp);
        modCount++;
        if(modCount>maxModCount){
          bMatch=false;
          break;
        }
      } else if(b>mp.pos && aaMod[pep1[b]].count>0){
        
        mp.mod=1;
        mp.pos=b;
        modQueue.push_back(mp);

        modList[b+1]=aaMod[pep1[b]].mass[0];
        modMass+=aaMod[pep1[b]].mass[0];
        modCount++;
        if(modCount>maxModCount){
          bMatch=false;
          break;
        }
      }

      //record the b- and y-ions making note if we include the linker mass
      fragMass=mMass+modMass;
		  if(b>=link) fragMass+=linkMass;
      bIons[0][b] = fragMass;
      tyIons[y] = (totalMass-fragMass);

      b++;
      y--; 
	  }

    if(!bMatch) continue;
    if(modMass==0) return 0;

    for(i=0;i<ionCount;i++){
      yIons[0][i]=tyIons[i]+modMass;
      for(j=1;j<6;j++){
        bIons[j][i] = ((bIons[0][i]+1.007276466*j)/j);
        yIons[j][i] = ((yIons[0][i]+1.007276466*j)/j);
      }
    }
    return totalMass+modMass;

  }

  return 0;

}

void KIons::buildNCIons(){
	
  int b;
	int i;
  int j;
  int y;
  int p1;
  int p2;
	double mMass;
  double totalMass;

  double b1[512];
  double b2[512];
  double y1[512];
  double y2[512];

  //set up boundaries
  ionCount=pep1Len-1;
  b=0;
  y=ionCount-1;

  totalMass=pep1Mass;

	//b- & y-ions from first peptide
  mMass=0;
	for(i=0;i<ionCount;i++){
    mMass+=aaMass[pep1[i]];
    b1[b] = mMass;
    y1[y] = (totalMass-mMass);
    b++;
    y--;
  }

  //If the two peptides are different
  p1=ionCount;
  p2=pep2Len-1;
  b=0;
  y=p2-1;

  totalMass=pep2Mass;

  //b- & y-ions from first peptide
  mMass=0;
  for(i=0;i<p2;i++){
    mMass+=aaMass[pep2[i]];
    b2[b] = mMass;
    y2[y] = (totalMass-mMass);
    b++;
    y--;
  }

  //Merge two series in order
  i=0;
  j=0;
  ionCount=0;
  while(i<p1 && j<p2){
    if(b1[i]<b2[j]) bIons[0][ionCount++]=b1[i++];
    else bIons[0][ionCount++]=b2[j++];
  }
  while(i<p1) bIons[0][ionCount++]=b1[i++];
  while(j<p2) bIons[0][ionCount++]=b2[j++];
  i=0;
  j=0;
  b=0;
  while(i<p1 && j<p2){
    if(y1[i]<y2[j]) yIons[0][b++]=y1[i++];
    else yIons[0][b++]=y2[j++];
  }
  while(i<p1) yIons[0][b++]=y1[i++];
  while(j<p2) yIons[0][b++]=y2[j++];

  //propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<6;j++){
      bIons[j][i] = ((bIons[0][i]+1.007276466*j)/j);
      yIons[j][i] = ((yIons[0][i]+1.007276466*j)/j);
    }
  }

}

void KIons::buildSingletIons(int link){

  int b;
	int i;
  int j;
  int y;
	double mMass;

  //set up boundaries
  ionCount=pep1Len-1;
  b=0;
  y=ionCount-1;

	//b- & y-ions from first peptide
  mMass=0;
	for(i=0;i<ionCount;i++){
    mMass+=aaMass[pep1[i]];
    if(i>=link) { //negative mass indicates that the ion needs a modmass as well
      yIons[0][y--] = pep1Mass-mMass;
      bIons[0][b++] = -mMass;
    } else {
      yIons[0][y--] = -(pep1Mass-mMass);
      bIons[0][b++] = mMass;
    }
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<6;j++){
      if(bIons[0][i]<0) bIons[j][i] = ((bIons[0][i]-1.007276466*j)/j);
      else bIons[j][i] = ((bIons[0][i]+1.007276466*j)/j);
      if(yIons[0][i]<0) yIons[j][i] = ((yIons[0][i]-1.007276466*j)/j);
      else yIons[j][i] = ((yIons[0][i]+1.007276466*j)/j);
    }
  }

}

void KIons::buildXIons(double linkMass, int link1, int link2){
	
  int b;
	int i;
  int j;
  int y;
	double mMass1;
  double mMass2;
  double totalMass;
  double fragMass1;
  double fragMass2;

  //Need to check total ion count and return a false if it exceeds array sizes

	//Calc totalMass
	totalMass = pep1Mass + pep2Mass + linkMass;

  //clear ion series for each charge state
  ionCount = pep1Len+pep2Len-2;
  b=0;
  y=ionCount-1;

  i=0;
  j=0;
  mMass1=aaMass[pep1[i]];
  mMass2=aaMass[pep2[j]];
  if(i>=link1) fragMass1=mMass1+linkMass+pep2Mass;
  else fragMass1=mMass1;
  if(j>=link2) fragMass2=mMass2+linkMass+pep1Mass;
  else fragMass2=mMass2;
  
  while(i<pep1Len || j<pep2Len){
    if(fragMass1<fragMass2){
      bIons[0][b] = fragMass1;
      yIons[0][y] = (totalMass-fragMass1);
      i++;
      b++;
      y--;
      if(i<pep1Len){
        mMass1+=aaMass[pep1[i]];
        if(i>=link1) fragMass1=mMass1+linkMass+pep2Mass;
        else fragMass1=mMass1;
      } else {
        fragMass1+=1000.0;
      }
    } else {
      bIons[0][b] = fragMass2;
      yIons[0][y] = (totalMass-fragMass2);
      j++;
      b++;
      y--;

      if(j<pep2Len){
        mMass2+=aaMass[pep2[j]];
        if(j>=link2) fragMass2=mMass2+linkMass+pep1Mass;
        else fragMass2=mMass2;
      } else {
        fragMass2+=1000.0;
      }
    }
  }

  for(i=0;i<ionCount;i++){
    for(j=1;j<6;j++){
      bIons[j][i] = ((bIons[0][i]+1.007276466*j)/j);
      yIons[j][i] = ((yIons[0][i]+1.007276466*j)/j);
    }
  }

}


double KIons::getModMass(int index){
  return modMassArray[index];
}

int KIons::getModMassSize(){
  return modMassArray.size();
}

double* KIons::getMods(){
  return modList;
}

int KIons::getIonCount (){
  return ionCount;
}

void KIons::getPeptide(bool bPepOne, char *seq){
  if(bPepOne){
    strncpy(seq,pep1,pep1Len);
    seq[pep1Len]='\0';
  } else {
    strncpy(seq,pep2,pep2Len);
    seq[pep1Len]='\0';
  }
}

void KIons::getPeptideMods(vector<kPepMod>& v){
  unsigned int i;
  kPepMod m;
  v.clear();
  for(i=0;i<modQueue.size();i++){
    if(modQueue[i].pos==-1) {
      m.pos=0; //convert n-terminus to first AA
      m.mass=modList[m.pos];
    } else {
      m.pos=(char)modQueue[i].pos;
      m.mass=modList[m.pos+1];
    }
    v.push_back(m);
  }
}

void KIons::modMasses(double mass, int index, int count){

  int i,j;
  for(i=index;i<128;i++){
    if(aaMod[i].count==0) continue;
    if(i==64 && count>0) continue;
    for(j=0;j<aaMod[i].count;j++){
      modMassArray.push_back(mass+aaMod[i].mass[j]);
      //printf("%.6lf\n",mass+aaMod[i].mass[j]);
      if(count+1<maxModCount) modMasses(mass+aaMod[i].mass[j],i,count+1);
    }
  }
}

void KIons::setMaxModCount(int i){
  maxModCount=i;
}

void KIons::setPeptide(bool bPepOne, char* seq, int len, double mass){
  if(bPepOne){
    pep1=seq;
    pep1Len=len;
    pep1Mass=mass;

    //Find better location for this block?
    if(modList!=NULL) delete [] modList;
    modList = new double[len];
    for(int i=0;i<len;i++) modList[i]=0;

    kModPos mp;
    modQueue.clear();
    mp.mass=0;
    mp.mod=0;
    mp.pos=-1;
    mp.modMass=0;
    mp.modCount=0;
    modQueue.push_back(mp);

  } else {
    pep2=seq;
    pep2Len=len;
    pep2Mass=mass;
  }

}

/*============================
  Utilities
============================*/
int KIons::compareD(const void *p1, const void *p2){
  const double d1 = *(double *)p1;
  const double d2 = *(double *)p2;
  if(d1<d2) return -1;
  else if(d1>d2) return 1;
  else return 0;
}

