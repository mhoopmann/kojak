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
  diffModsOnXL=false;
  monoModsOnXL=false;
}

KIons::~KIons(){
  if(modList!=NULL) delete [] modList;
  pep1=NULL;
  pep2=NULL;
}

KIonSet& KIons::operator[ ](const int& i){
  return sets[i];
}

KIonSet* KIons::at(const int& i){
  return &sets[i];
}


void KIons::addFixedMod(char mod, double mass){
  aaMass[mod]+=mass;
}

void KIons::addMod(char mod, bool xl, double mass){
  aaMod[mod].mod[aaMod[mod].count].xl=xl;
  aaMod[mod].mod[aaMod[mod].count++].mass=mass;
}

void KIons::buildIons(){
	
  int b;
	int i;
  int j;
  int y;
  double fragMass;

  //set up boundaries
  ionCount=sets[0].len;
  b=0;
  y=ionCount-1;

	//b- & y-ions from first peptide
  fragMass=0;
	for(i=0;i<ionCount;i++){
    fragMass+=aaMass[pep1[i]];
    sets[0].aIons[0][b] = fragMass - 27.9949141;
    sets[0].bIons[0][b] = fragMass;
    sets[0].cIons[0][b] = fragMass + 17.026547;
    sets[0].xIons[0][y] = (pep1Mass-fragMass) + 25.9792649;
    sets[0].yIons[0][y] = (pep1Mass-fragMass);
    sets[0].zIons[0][y] = (pep1Mass-fragMass) - 16.0187224; //z. not z
    b++;
    y--;
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<6;j++){
      sets[0].aIons[j][i] = ((sets[0].aIons[0][i]+1.007276466*j)/j);
      sets[0].bIons[j][i] = ((sets[0].bIons[0][i]+1.007276466*j)/j);
      sets[0].cIons[j][i] = ((sets[0].cIons[0][i]+1.007276466*j)/j);
      sets[0].xIons[j][i] = ((sets[0].xIons[0][i]+1.007276466*j)/j);
      sets[0].yIons[j][i] = ((sets[0].yIons[0][i]+1.007276466*j)/j);
      sets[0].zIons[j][i] = ((sets[0].zIons[0][i]+1.007276466*j)/j);
    }
  }

}

void KIons::buildLoopIons(double linkMass, int link1, int link2){
	
  int a1,a2;
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
  if(link2>link1) {//Link2 must be larger than link1
    a1=link1;
    a2=link2;
  } else {
    a1=link2;
    a2=link1;
  }
  ionCount=pep1Len-1-(a2-a1); 
  b=0;
  y=ionCount-1;

  //update sets mass:
  sets[0].mass+=linkMass;
  totalMass=pep1Mass+linkMass;

	//b- & y-ions from first peptide
  mMass=0;
	for(i=0;i<len;i++){
    mMass+=aaMass[pep1[i]];
    if(i>=a1 && i<a2) continue;
    else if(i>=a2) fragMass=mMass+linkMass;
    else fragMass = mMass;

    sets[0].aIons[0][b] = fragMass - 27.9949141;
    sets[0].bIons[0][b] = fragMass;
    sets[0].cIons[0][b] = fragMass + 17.026547;
    sets[0].xIons[0][y] = (totalMass-fragMass) + 25.9792649;
    sets[0].yIons[0][y] = (totalMass-fragMass);
    sets[0].zIons[0][y] = (totalMass-fragMass) - 16.0187224; //z. not z

    b++;
    y--;
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<6;j++){
      sets[0].aIons[j][i] = ((sets[0].aIons[0][i]+1.007276466*j)/j);
      sets[0].bIons[j][i] = ((sets[0].bIons[0][i]+1.007276466*j)/j);
      sets[0].cIons[j][i] = ((sets[0].cIons[0][i]+1.007276466*j)/j);
      sets[0].xIons[j][i] = ((sets[0].xIons[0][i]+1.007276466*j)/j);
      sets[0].yIons[j][i] = ((sets[0].yIons[0][i]+1.007276466*j)/j);
      sets[0].zIons[j][i] = ((sets[0].zIons[0][i]+1.007276466*j)/j);
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
  ionCount=sets[0].len; //pep1Len-1;
  b=0;
  y=ionCount-1;

	//b- & y-ions from first peptide
  mMass=0;
	for(i=0;i<ionCount;i++){
    mMass+=aaMass[pep1[i]];
    if(i>=link) { //negative mass indicates that the ion needs a modmass as well
      sets[0].xIons[0][y] = pep1Mass-mMass + 25.9792649;
      sets[0].yIons[0][y] = pep1Mass-mMass;
      sets[0].zIons[0][y] = pep1Mass-mMass - 16.0187224;   //z.
      sets[0].aIons[0][b] = -(mMass - 27.9949141);
      sets[0].bIons[0][b] = -mMass;
      sets[0].cIons[0][b] = -(mMass + 17.026547);
    } else {
      sets[0].xIons[0][y] = -(pep1Mass-mMass + 25.9792649);
      sets[0].yIons[0][y] = -(pep1Mass-mMass);
      sets[0].zIons[0][y] = -(pep1Mass-mMass - 16.0187224);
      sets[0].aIons[0][b] = mMass - 27.9949141;
      sets[0].bIons[0][b] = mMass;
      sets[0].cIons[0][b] = mMass + 17.026547;
    }
    y--;
    b++;
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<6;j++){
      if(sets[0].aIons[0][i]<0) sets[0].aIons[j][i] = ((sets[0].aIons[0][i]-1.007276466*j)/j);
      else sets[0].aIons[j][i] = ((sets[0].aIons[0][i]+1.007276466*j)/j);
      if(sets[0].bIons[0][i]<0) sets[0].bIons[j][i] = ((sets[0].bIons[0][i]-1.007276466*j)/j);
      else sets[0].bIons[j][i] = ((sets[0].bIons[0][i]+1.007276466*j)/j);
      if(sets[0].cIons[0][i]<0) sets[0].cIons[j][i] = ((sets[0].cIons[0][i]-1.007276466*j)/j);
      else sets[0].cIons[j][i] = ((sets[0].cIons[0][i]+1.007276466*j)/j);
      if(sets[0].xIons[0][i]<0) sets[0].xIons[j][i] = ((sets[0].xIons[0][i]-1.007276466*j)/j);
      else sets[0].xIons[j][i] = ((sets[0].xIons[0][i]+1.007276466*j)/j);
      if(sets[0].yIons[0][i]<0) sets[0].yIons[j][i] = ((sets[0].yIons[0][i]-1.007276466*j)/j);
      else sets[0].yIons[j][i] = ((sets[0].yIons[0][i]+1.007276466*j)/j);
      if(sets[0].zIons[0][i]<0) sets[0].zIons[j][i] = ((sets[0].zIons[0][i]-1.007276466*j)/j);
      else sets[0].zIons[j][i] = ((sets[0].zIons[0][i]+1.007276466*j)/j);
    }
  }

}

void KIons::modIonsRec(int start, int link, int index, int depth, bool xl){
  int i,j,k,n;

  for(i=start;i<pep1Len-1;i++){

    //don't modify site where the cross-linker is bound.
    if(i==link) continue;

    //Check if amino acid is on the modification list
    for(j=0;j<aaMod[pep1[i]].count;j++){

      //skip mods not allowed on this peptide
      if(xl && !aaMod[pep1[i]].mod[j].xl && !diffModsOnXL) continue;
      if(xl && aaMod[pep1[i]].mod[j].xl && !monoModsOnXL) continue;

      //Add masses
      KIonSet s=sets[index];

      for(k=i;k<ionCount;k++){
        for(n=1;n<6;n++){
          if(s.aIons[0][k]<0) s.aIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.aIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
          if(s.bIons[0][k]<0) s.bIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.bIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
          if(s.cIons[0][k]<0) s.cIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.cIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
        }
      }
      for(k=ionCount-i;k<ionCount;k++){
        for(n=1;n<6;n++){
          if(s.xIons[0][k]<0) s.xIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.xIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
          if(s.yIons[0][k]<0) s.yIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.yIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
          if(s.zIons[0][k]<0) s.zIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.zIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
        }
      }
      s.mods[i]=aaMod[pep1[i]].mod[j].mass;
      s.mass+=aaMod[pep1[i]].mod[j].mass;
      s.difMass+=aaMod[pep1[i]].mod[j].mass;

      //Add to list
      sets.push_back(s);

      //solve another one
      if(depth+1<maxModCount) modIonsRec(i+1,link,(int)(sets.size())-1,depth+1,xl);
    }
  }

}

void KIons::modLoopIonsRec(int start, int link, int link2, int index, int depth, bool xl){
  int i,j,k,n;
  int a,b; //first and second link positions
  int pos;

  //always order the links from first to second position
  if(link<link2){
    a=link;
    b=link2;
  } else {
    a=link2;
    b=link;
  }

  for(i=start;i<pep1Len-1;i++){

    //don't modify site where the cross-linker is bound.
    if(i==a || i==b) continue;

    if(i<a) pos=i;
    else if(i<b) pos=a;
    else pos=i-b+a;

    //Check if amino acid is on the modification list
    for(j=0;j<aaMod[pep1[i]].count;j++){

      //skip mods not allowed on this peptide
      if(xl && !aaMod[pep1[i]].mod[j].xl) continue;

      //Add masses
      KIonSet s=sets[index];
      for(k=pos;k<ionCount;k++){
        for(n=1;n<6;n++){
          if(s.aIons[0][k]<0) s.aIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.aIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
          if(s.bIons[0][k]<0) s.bIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.bIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
          if(s.cIons[0][k]<0) s.cIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.cIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
        }
      }
      for(k=ionCount-pos;k<ionCount;k++){
        for(n=1;n<6;n++){
          if(s.xIons[0][k]<0) s.xIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.xIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
          if(s.yIons[0][k]<0) s.yIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.yIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
          if(s.zIons[0][k]<0) s.zIons[n][k] -= (aaMod[pep1[i]].mod[j].mass/n);
          else s.zIons[n][k] += (aaMod[pep1[i]].mod[j].mass/n);
        }
      }
      s.mods[i]=aaMod[pep1[i]].mod[j].mass;
      s.mass+=aaMod[pep1[i]].mod[j].mass;
      s.difMass+=aaMod[pep1[i]].mod[j].mass;

      //Add to list
      sets.push_back(s);

      //solve another one
      if(depth+1<maxModCount) modLoopIonsRec(i+1,link,link2,(int)(sets.size())-1,depth+1,xl);
    }
  }

}

void KIons::reset(){
  sets.clear();
  KIonSet k(pep1Len,pep1Mass);
  sets.push_back(k);
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

int KIons::getPeptideLen(){
  return pep1Len;
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

void KIons::setMaxModCount(int i){
  maxModCount=i;
}

void KIons::setModFlags(bool monoMods, bool difMods){
  monoModsOnXL=monoMods;
  diffModsOnXL=difMods;
}

void KIons::setPeptide(bool bPepOne, char* seq, int len, double mass){
  if(bPepOne){
    pep1=seq;
    pep1Len=len;
    pep1Mass=mass;

    sets.clear();
    KIonSet k(pep1Len,pep1Mass);
    sets.push_back(k);

  } else {
    pep2=seq;
    pep2Len=len;
    pep2Mass=mass;
  }

}

int KIons::size(){
  return (int)sets.size();
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

