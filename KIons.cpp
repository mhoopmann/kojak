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

using namespace std;

KIons::KIons(){
  int i;
  for(i=0;i<128;i++){
    aaMass[i]=0;
    aaMassn15[i]=0;
    aaFixedModMass[i]=0;
    aaMod[i].count=0;
    site[i]=false;
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
  aaMass['U']=150.9536303;
  aaMass['W']=186.0793065;
  aaMass['Y']=163.0633228;
  aaMass['c']=0;
  aaMass['n']=0;
  aaMass['$']=0;
  aaMass['%']=0;

  aaMassn15['A'] = 72.0341452;
  aaMassn15['C'] = 104.0062152;
  aaMassn15['D'] = 116.0239734;
  aaMassn15['E'] = 130.0396226;
  aaMassn15['F'] = 148.0654436;
  aaMassn15['G'] = 58.018496;
  aaMassn15['H'] = 140.0500106;
  aaMassn15['I'] = 114.0810928;
  aaMassn15['K'] = 130.0890255;
  aaMassn15['L'] = 114.0810928;
  aaMassn15['M'] = 132.0375136;
  aaMassn15['N'] = 116.036992;
  aaMassn15['P'] = 98.0497944;
  aaMassn15['Q'] = 130.0526412;
  aaMassn15['R'] = 160.0892417;
  aaMassn15['S'] = 88.0290593;
  aaMassn15['T'] = 102.0447085;
  aaMassn15['U'] = 151.9506652;
  aaMassn15['V'] = 100.0654436;
  aaMassn15['W'] = 188.0733763;
  aaMassn15['Y'] = 164.0603577;
  aaMassn15['c'] = 0;
  aaMassn15['n'] = 0;
  aaMassn15['$'] = 0;
  aaMassn15['%'] = 0;

  modList=NULL;
  pep1=NULL;
  pep2=NULL;
  maxModCount=0;
  ionCount=0;
  diffModsOnXL=false;
  monoModsOnXL=false;

  seriesA=false;
  seriesB=false;
  seriesC=false;
  seriesX=false;
  seriesY=false;
  seriesZ=false;
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
  aaMassn15[mod]+=mass;
  aaFixedModMass[mod]=mass;
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
  if(n15Pep1){
    fragMass = aaMassn15['n'];
    if (nPep1) fragMass += aaMassn15['$'];
  } else {
    fragMass=aaMass['n'];
    if (nPep1) fragMass+=aaMass['$'];
  }
	for(i=0;i<ionCount;i++){
    if (n15Pep1) fragMass += aaMassn15[pep1[i]];
    else fragMass+=aaMass[pep1[i]];
    if (i == ionCount - 1) {
      if (n15Pep1) {
        fragMass += aaMassn15['c'];
        if (cPep1) fragMass += aaMassn15['%'];
      } else {
        fragMass += aaMass['c'];
        if (cPep1) fragMass += aaMass['%'];
      }
    }
    if(seriesA) sets[0].aIons[0][b].mz = fragMass - 27.9949141;
    if (seriesB) sets[0].bIons[0][b].mz = fragMass;
    if (seriesC) sets[0].cIons[0][b].mz = fragMass + 17.026547;
    if (seriesX) sets[0].xIons[0][y].mz = (pep1Mass - fragMass) + 25.9792649;
    if (seriesY) sets[0].yIons[0][y].mz = (pep1Mass - fragMass);
    if (seriesZ) sets[0].zIons[0][y].mz = (pep1Mass - fragMass) - 16.0187224; //z. not z
    b++;
    y--;
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<4;j++){
      if(seriesA) sets[0].aIons[j][i].mz = ((sets[0].aIons[0][i].mz + 1.007276466*j) / j);
      if (seriesB) sets[0].bIons[j][i].mz = ((sets[0].bIons[0][i].mz + 1.007276466*j) / j);
      if (seriesC) sets[0].cIons[j][i].mz = ((sets[0].cIons[0][i].mz + 1.007276466*j) / j);
      if (seriesX) sets[0].xIons[j][i].mz = ((sets[0].xIons[0][i].mz + 1.007276466*j) / j);
      if (seriesY) sets[0].yIons[j][i].mz = ((sets[0].yIons[0][i].mz + 1.007276466*j) / j);
      if (seriesZ) sets[0].zIons[j][i].mz = ((sets[0].zIons[0][i].mz + 1.007276466*j) / j);
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
  if(n15Pep1){
    mMass = aaMassn15['n'];
    if (nPep1) mMass += aaMassn15['$'];
  } else {
    mMass = aaMass['n'];
    if (nPep1) mMass += aaMass['$'];
  }
	for(i=0;i<len;i++){
    if (n15Pep1) mMass += aaMassn15[pep1[i]];
    else mMass+=aaMass[pep1[i]];
    if (i == len - 1) {
      if(n15Pep1) {
        mMass += aaMassn15['c'];
        if (cPep1) mMass += aaMassn15['%'];
      } else {
        mMass += aaMass['c'];
        if (cPep1) mMass += aaMass['%'];
      }
    }
    if(i>=a1 && i<a2) continue;
    else if(i>=a2) fragMass=mMass+linkMass;
    else fragMass = mMass;

    if (seriesA) sets[0].aIons[0][b].mz = fragMass - 27.9949141;
    if (seriesB) sets[0].bIons[0][b].mz = fragMass;
    if (seriesC) sets[0].cIons[0][b].mz = fragMass + 17.026547;
    if (seriesX) sets[0].xIons[0][y].mz = (totalMass - fragMass) + 25.9792649;
    if (seriesY) sets[0].yIons[0][y].mz = (totalMass - fragMass);
    if (seriesZ) sets[0].zIons[0][y].mz = (totalMass - fragMass) - 16.0187224; //z. not z

    b++;
    y--;
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<4;j++){
      if (seriesA) sets[0].aIons[j][i].mz = ((sets[0].aIons[0][i].mz + 1.007276466*j) / j);
      if (seriesB) sets[0].bIons[j][i].mz = ((sets[0].bIons[0][i].mz + 1.007276466*j) / j);
      if (seriesC) sets[0].cIons[j][i].mz = ((sets[0].cIons[0][i].mz + 1.007276466*j) / j);
      if (seriesX) sets[0].xIons[j][i].mz = ((sets[0].xIons[0][i].mz + 1.007276466*j) / j);
      if (seriesY) sets[0].yIons[j][i].mz = ((sets[0].yIons[0][i].mz + 1.007276466*j) / j);
      if (seriesZ) sets[0].zIons[j][i].mz = ((sets[0].zIons[0][i].mz + 1.007276466*j) / j);
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
  if(n15Pep1){
    mMass = aaMassn15['n'];
    if (nPep1) mMass += aaMassn15['$'];
  } else {
    mMass=aaMass['n'];
    if (nPep1) mMass += aaMass['$'];
  }
	for(i=0;i<ionCount;i++){
    if (n15Pep1) mMass += aaMassn15[pep1[i]];
    else mMass+=aaMass[pep1[i]];
    if (i == ionCount - 1) {
      if(n15Pep1){
        mMass += aaMassn15['n'];
        if (cPep1) mMass += aaMassn15['%'];
      } else {
        mMass += aaMass['n'];
        if (cPep1) mMass += aaMass['%'];
      }
    }
    if(i>=link) { //negative mass indicates that the ion needs a modmass as well
      if (seriesX) sets[0].xIons[0][y].mz = pep1Mass - mMass + 25.9792649;
      if (seriesY) sets[0].yIons[0][y].mz = pep1Mass - mMass;
      if (seriesZ) sets[0].zIons[0][y].mz = pep1Mass - mMass - 16.0187224;   //z.
      if (seriesA) sets[0].aIons[0][b].mz = -(mMass - 27.9949141);
      if (seriesB) sets[0].bIons[0][b].mz = -mMass;
      if (seriesC) sets[0].cIons[0][b].mz = -(mMass + 17.026547);
    } else {
      if (seriesX) sets[0].xIons[0][y].mz = -(pep1Mass - mMass + 25.9792649);
      if (seriesY) sets[0].yIons[0][y].mz = -(pep1Mass - mMass);
      if (seriesZ) sets[0].zIons[0][y].mz = -(pep1Mass - mMass - 16.0187224);
      if (seriesA) sets[0].aIons[0][b].mz = mMass - 27.9949141;
      if (seriesB) sets[0].bIons[0][b].mz = mMass;
      if (seriesC) sets[0].cIons[0][b].mz = mMass + 17.026547;
    }
    y--;
    b++;
	}

  //Propagate remaining charge states
  for(i=0;i<ionCount;i++){
    for(j=1;j<4;j++){
      if (seriesA){
      if (sets[0].aIons[0][i].mz<0) sets[0].aIons[j][i].mz = ((sets[0].aIons[0][i].mz - 1.007276466*j) / j);
      else sets[0].aIons[j][i].mz = ((sets[0].aIons[0][i].mz + 1.007276466*j) / j);
      }
      if (seriesB){
      if (sets[0].bIons[0][i].mz<0) sets[0].bIons[j][i].mz = ((sets[0].bIons[0][i].mz - 1.007276466*j) / j);
      else sets[0].bIons[j][i].mz = ((sets[0].bIons[0][i].mz + 1.007276466*j) / j);
      }
      if (seriesC){
      if (sets[0].cIons[0][i].mz<0) sets[0].cIons[j][i].mz = ((sets[0].cIons[0][i].mz - 1.007276466*j) / j);
      else sets[0].cIons[j][i].mz = ((sets[0].cIons[0][i].mz + 1.007276466*j) / j);
      }
      if (seriesX){
      if (sets[0].xIons[0][i].mz<0) sets[0].xIons[j][i].mz = ((sets[0].xIons[0][i].mz - 1.007276466*j) / j);
      else sets[0].xIons[j][i].mz = ((sets[0].xIons[0][i].mz + 1.007276466*j) / j);
      }
      if (seriesY){
      if (sets[0].yIons[0][i].mz<0) sets[0].yIons[j][i].mz = ((sets[0].yIons[0][i].mz - 1.007276466*j) / j);
      else sets[0].yIons[j][i].mz = ((sets[0].yIons[0][i].mz + 1.007276466*j) / j);
      }
      if (seriesZ){
      if (sets[0].zIons[0][i].mz<0) sets[0].zIons[j][i].mz = ((sets[0].zIons[0][i].mz - 1.007276466*j) / j);
      else sets[0].zIons[j][i].mz = ((sets[0].zIons[0][i].mz + 1.007276466*j) / j);
      }
    }
  }

}

double KIons::getAAMass(char aa, bool n15){
  if(n15) return aaMassn15[aa];
  return aaMass[aa];
}
double KIons::getFixedModMass(char aa){
  return aaFixedModMass[aa];
}

void KIons::addModIonSet(int index, char aa, int pos, int modIndex, int loopPos){
  int k, n;

  //Copy ionset, then add masses to it
  sets.push_back(sets[index]);

  for (k = pos; k<ionCount; k++){
    for (n = 1; n<4; n++){
      if (seriesA){
      if (sets.back().aIons[0][k].mz<0) sets.back().aIons[n][k].mz -= (aaMod[aa].mod[modIndex].mass / n);
      else sets.back().aIons[n][k].mz += (aaMod[aa].mod[modIndex].mass / n);
      }
      if (seriesB){
      if (sets.back().bIons[0][k].mz<0) sets.back().bIons[n][k].mz -= (aaMod[aa].mod[modIndex].mass / n);
      else sets.back().bIons[n][k].mz += (aaMod[aa].mod[modIndex].mass / n);
      }
      if (seriesC){
      if (sets.back().cIons[0][k].mz<0) sets.back().cIons[n][k].mz -= (aaMod[aa].mod[modIndex].mass / n);
      else sets.back().cIons[n][k].mz += (aaMod[aa].mod[modIndex].mass / n);
      }
    }
  }
  for (k = ionCount - pos; k<ionCount; k++){
    for (n = 1; n<4; n++){
      if (seriesX){
      if (sets.back().xIons[0][k].mz<0) sets.back().xIons[n][k].mz -= (aaMod[aa].mod[modIndex].mass / n);
      else sets.back().xIons[n][k].mz += (aaMod[aa].mod[modIndex].mass / n);
      }
      if (seriesY){
      if (sets.back().yIons[0][k].mz<0) sets.back().yIons[n][k].mz -= (aaMod[aa].mod[modIndex].mass / n);
      else sets.back().yIons[n][k].mz += (aaMod[aa].mod[modIndex].mass / n);
      }
      if (seriesZ){
      if (sets.back().zIons[0][k].mz<0) sets.back().zIons[n][k].mz -= (aaMod[aa].mod[modIndex].mass / n);
      else sets.back().zIons[n][k].mz += (aaMod[aa].mod[modIndex].mass / n);
      }
    }
  }
  if (aa == 'n' || aa == '$') {
    sets.back().nTermMass += aaMod[aa].mod[modIndex].mass;
  }
  if (aa == 'c' || aa == '%') {
    sets.back().cTermMass += aaMod[aa].mod[modIndex].mass;
  }
  if (loopPos>-1) sets.back().mods[loopPos] = aaMod[aa].mod[modIndex].mass;
  else sets.back().mods[pos] = aaMod[aa].mod[modIndex].mass;
  sets.back().mass += aaMod[aa].mod[modIndex].mass;
  sets.back().difMass += aaMod[aa].mod[modIndex].mass;

}

void KIons::makeIonIndex(double binSize, double binOffset){
  cout << "KIons::makeIonIndex - get rid of this" << endl;
  /*
  for(size_t i=0;i<sets.size();i++){
    sets[i].makeIndex(binSize,binOffset,seriesA,seriesB,seriesC,seriesX,seriesY,seriesZ);
  }
  */
}

void KIons::modIonsRec(int start, int link, int index, int depth, bool xl){
  int i,j;

  for(i=start;i<pep1Len;i++){

    //don't modify site where the cross-linker is bound.
    if (i == link) continue;

    //Check if amino acid is on the modification list
    for(j=0;j<aaMod[pep1[i]].count;j++){

      //skip mods not allowed on this peptide
      if(xl && !aaMod[pep1[i]].mod[j].xl && !diffModsOnXL) continue;
      if(xl && aaMod[pep1[i]].mod[j].xl && !monoModsOnXL) continue;

      //skip mods if it is xl mod on a cut site
      if (aaMod[pep1[i]].mod[j].xl && i == pep1Len - 1 && !cPep1) continue;

      //Add masses
      addModIonSet(index,pep1[i],i,j);

      //solve another one
      if(depth+1<maxModCount) modIonsRec(i+1,link,(int)(sets.size())-1,depth+1,xl);
    }

    /*
    //Special case for peptide n-terminus
    if (i == 0 && link!=-1) {
      for (j = 0; j<aaMod['n'].count; j++){
        //skip mods not allowed on this peptide
        if (xl && !aaMod['n'].mod[j].xl && !diffModsOnXL) continue;
        if (xl && aaMod['n'].mod[j].xl && !monoModsOnXL) continue;
        //Add masses
        addModIonSet(index, 'n', i, j);
        //solve another one
        if (depth + 1<maxModCount) modIonsRec(i + 1, link, (int)(sets.size()) - 1, depth + 1, xl);
      }

      //Special case for protein n-terminus
      if (nPep1) {
        for (j = 0; j<aaMod['$'].count; j++){
          //skip mods not allowed on this peptide
          if (xl && !aaMod['$'].mod[j].xl && !diffModsOnXL) continue;
          if (xl && aaMod['$'].mod[j].xl && !monoModsOnXL) continue;
          //Add masses
          addModIonSet(index,'$', i, j);
          //solve another one
          if (depth + 1<maxModCount) modIonsRec(i + 1, link, (int)(sets.size()) - 1, depth + 1, xl);
        }
      }
    }

    //Special case for peptide c-terminus
    if (i == pep1Len - 1 && link!=-2){
      for (j = 0; j<aaMod['c'].count; j++){
        //skip mods not allowed on this peptide
        if (xl && !aaMod['c'].mod[j].xl && !diffModsOnXL) continue;
        if (xl && aaMod['c'].mod[j].xl && !monoModsOnXL) continue;
        //Add masses
        addModIonSet(index, 'c', i, j);
        //solve another one
        if (depth + 1<maxModCount) modIonsRec(i + 1, link, (int)(sets.size()) - 1, depth + 1, xl);
      }

      //Special case for protein c-terminus
      if (cPep1){
        for (j = 0; j<aaMod['%'].count; j++){
          //skip mods not allowed on this peptide
          if (xl && !aaMod['%'].mod[j].xl && !diffModsOnXL) continue;
          if (xl && aaMod['%'].mod[j].xl && !monoModsOnXL) continue;
          //Add masses
          addModIonSet(index, '%', i, j);
          //solve another one
          if (depth + 1<maxModCount) modIonsRec(i + 1, link, (int)(sets.size()) - 1, depth + 1, xl);
        }
      }
    }
    */
  
  }

}

void KIons::modIonsRec2(int start, int link, int index, int depth, bool xl){
  int j;

  //if n-terminus can be linked, proceed as if it is linked
  if (link == 0 && nPep1 && site['n']){
    modIonsRec(start, -1, index, depth, xl);
  } else if (link == pep1Len - 1 && cPep1 && site['c']){ //if c-terminus can be linked, proceed as if it is linked
    modIonsRec(start, -1, index, depth, xl);
  }

  //now proceed with linking on an amino acid
  if(link==0 && !site[pep1[link]]) return; //if this is not a linkable side chain, stop here.
  else if (link == pep1Len - 1 && site[pep1[pep1Len-1]]) return;
  modIonsRec(start, link, index, depth, xl);
  //if(link<0) return; //stop now if peptide cannot be linked, as modifications have been exhausted

  //if at n-terminus and aa is not linkable, stop now; n-terminus cannot be modified because it is holding the linker
  if(link==0 && !site[pep1[0]]) return; 
  else if(link==pep1Len-1 && !site[pep1[pep1Len-1]]) return; //ditto for c-terminus

  //From here, proceed as if the link is on the amino acid with a terminal modification
  //Special case for protein n-terminus
  if (nPep1) {
    for (j = 0; j<aaMod['$'].count; j++){
      //skip mods not allowed on this peptide
      if (xl && !aaMod['$'].mod[j].xl && !diffModsOnXL) continue;
      if (xl && aaMod['$'].mod[j].xl && !monoModsOnXL) continue;
      //Add masses
      addModIonSet(index, '$', 0, j);
      //solve another one
      if (depth + 1<maxModCount) modIonsRec(start, link, (int)(sets.size()) - 1, depth + 1, xl);
    }
  } else {    //Check peptide n-terminus mods
    for (j = 0; j<aaMod['n'].count; j++){
      //skip mods not allowed on this peptide
      if (xl && !aaMod['n'].mod[j].xl && !diffModsOnXL) continue;
      if (xl && aaMod['n'].mod[j].xl && !monoModsOnXL) continue;
      //Add masses
      addModIonSet(index, 'n', 0, j);
      //solve another one
      if (depth + 1<maxModCount) modIonsRec(start, link, (int)(sets.size()) - 1, depth + 1, xl);
    }
  }

  //Special case for protein c-terminus
  if (cPep1){
    for (j = 0; j<aaMod['%'].count; j++){
      //skip mods not allowed on this peptide
      if (xl && !aaMod['%'].mod[j].xl && !diffModsOnXL) continue;
      if (xl && aaMod['%'].mod[j].xl && !monoModsOnXL) continue;
      //Add masses
      addModIonSet(index, '%', pep1Len-1, j);
      //solve another one
      if (depth + 1<maxModCount) modIonsRec(start, link, (int)(sets.size()) - 1, depth + 1, xl);
    }
  } else { //Special case for peptide c-terminus
    for (j = 0; j<aaMod['c'].count; j++){
      //skip mods not allowed on this peptide
      if (xl && !aaMod['c'].mod[j].xl && !diffModsOnXL) continue;
      if (xl && aaMod['c'].mod[j].xl && !monoModsOnXL) continue;
      //Add masses
      addModIonSet(index, 'c', pep1Len-1, j);
      //solve another one
      if (depth + 1<maxModCount) modIonsRec(start, link, (int)(sets.size()) - 1, depth + 1, xl);
    }
  }

}

void KIons::modLoopIonsRec(int start, int link, int link2, int index, int depth, bool xl){
  int i,j;
  int pos;
  int trueLink=link;
  if(trueLink==-1) trueLink=0;
  int trueLink2=link2;
  if(trueLink2==-1) trueLink2=pep1Len-1;

  for(i=start;i<pep1Len;i++){

    //don't modify site where the cross-linker is bound.
    if(i==link || i==link2) continue;

    //mod position applied to left or right side of link.
    //if (i<link) pos = i;
    //else if (i<link2) pos = link;
    //else pos = i - link2 + link;
    if(i<trueLink) pos=i;
    else if(i<trueLink2) pos=trueLink;
    else pos=i-trueLink2+trueLink;
    
    //Check if amino acid is on the modification list
    for(j=0;j<aaMod[pep1[i]].count;j++){

      //skip mods not allowed on this peptide
      if (xl && !aaMod[pep1[i]].mod[j].xl && !diffModsOnXL) continue;
      if (xl && aaMod[pep1[i]].mod[j].xl && !monoModsOnXL) continue;

      //skip mod if it is xl on a cut site
      if (aaMod[pep1[i]].mod[j].xl && i == pep1Len - 1 && !cPep1) continue;

      //Add masses
      addModIonSet(index, pep1[i], pos,j,i);

      //solve another one
      if(depth+1<maxModCount) modLoopIonsRec(i+1,link,link2,(int)(sets.size())-1,depth+1,xl);
    }
    
    /*
    //Special case for peptide n-terminus
    if (i == 0) {
      for (j = 0; j<aaMod['n'].count; j++){
        //skip mods not allowed on this peptide
        if (xl && !aaMod['n'].mod[j].xl && !diffModsOnXL) continue;
        if (xl && aaMod['n'].mod[j].xl && !monoModsOnXL) continue;
        //Add masses
        addModIonSet(index, 'n', pos, j,i);
        //solve another one
        if (depth + 1<maxModCount) modLoopIonsRec(i + 1, link, link2,(int)(sets.size()) - 1, depth + 1, xl);
      }

      //Special case for protein n-terminus
      if (nPep1) {
        for (j = 0; j<aaMod['$'].count; j++){
          //skip mods not allowed on this peptide
          if (xl && !aaMod['$'].mod[j].xl && !diffModsOnXL) continue;
          if (xl && aaMod['$'].mod[j].xl && !monoModsOnXL) continue;
          //Add masses
          addModIonSet(index, '$', pos, j, i);
          //solve another one
          if (depth + 1<maxModCount) modLoopIonsRec(i + 1, link, link2, (int)(sets.size()) - 1, depth + 1, xl);
        }
      }
    }
    
    //Special case for peptide c-terinus
    if (i == pep1Len - 1){
      for (j = 0; j<aaMod['c'].count; j++){
        //skip mods not allowed on this peptide
        if (xl && !aaMod['c'].mod[j].xl && !diffModsOnXL) continue;
        if (xl && aaMod['c'].mod[j].xl && !monoModsOnXL) continue;
        //Add masses
        addModIonSet(index, 'c', pos, j,i);
        //solve another one
        if (depth + 1<maxModCount) modLoopIonsRec(i + 1, link, link2,(int)(sets.size()) - 1, depth + 1, xl);
      }

      //Special case for protein c-terinus
      if (cPep1){
        for (j = 0; j<aaMod['%'].count; j++){
          //skip mods not allowed on this peptide
          if (xl && !aaMod['%'].mod[j].xl && !diffModsOnXL) continue;
          if (xl && aaMod['%'].mod[j].xl && !monoModsOnXL) continue;
          //Add masses
          addModIonSet(index, '%', pos, j, i);
          //solve another one
          if (depth + 1<maxModCount) modLoopIonsRec(i + 1, link, link2, (int)(sets.size()) - 1, depth + 1, xl);
        }
      }
    }
    */

  }

}

void KIons::modLoopIonsRec2(int start, int link, int link2, int index, int depth, bool xl){
  int j;
  bool skipN=false;
  bool skipC=false;

  //small note to self, this doesn't check the possibility, for example, of a c-terminal lysine linked to the c-terminus.
  //ditto for the n-terminus and n-terminal lysine.

  //if n-terminus can be linked, proceed as if it is linked
  if (link == 0 && nPep1 && site['n']){
    if (link2 == pep1Len - 1 && cPep1 && site['c']){ //if c-terminus can be linked, proceed as if it is linked
      modLoopIonsRec(start, -1, -1, index, depth, xl);
    } else {
      modLoopIonsRec(start, -1, link2, index, depth, xl);
    }
  } else if (link2==pep1Len-1 && cPep1 && site['c']){ //if c-terminus can be linked, proceed as if it is linked
    modLoopIonsRec(start, link, -1, index, depth, xl);
  }

  //now proceed with linking on an amino acid
  modLoopIonsRec(start, link, link2, index, depth, xl);

  //if at n-terminus and aa is not linkable, stop now; n-terminus cannot be modified because it is holding the linker
  if (link == 0 && !site[pep1[0]]) skipN=true;
  if (link2 == pep1Len - 1 && !site[pep1[pep1Len - 1]]) skipC=true; //ditto for c-terminus

  //From here, proceed as if the link is on the amino acid with a terminal modification
  if(!skipN){
    if (nPep1) { //Special case for protein n-terminus
      for (j = 0; j<aaMod['$'].count; j++){
        //skip mods not allowed on this peptide
        if (xl && !aaMod['$'].mod[j].xl && !diffModsOnXL) continue;
        if (xl && aaMod['$'].mod[j].xl && !monoModsOnXL) continue;
        //Add masses
        addModIonSet(index, '$', 0, j, 0);
        //solve another one
        if (depth + 1<maxModCount) modLoopIonsRec(start, link, link2, (int)(sets.size()) - 1, depth + 1, xl);
      }
    } else { //Special case for peptide n-terminus
      for (j = 0; j<aaMod['n'].count; j++){
        //skip mods not allowed on this peptide
        if (xl && !aaMod['n'].mod[j].xl && !diffModsOnXL) continue;
        if (xl && aaMod['n'].mod[j].xl && !monoModsOnXL) continue;
        //Add masses
        addModIonSet(index, 'n', 0, j, 0);
        //solve another one
        if (depth + 1<maxModCount) modLoopIonsRec(start, link, link2, (int)(sets.size()) - 1, depth + 1, xl);
      }      
    }
  }
  
  if(!skipC){
    if (cPep1){//Special case for protein c-terinus
      for (j = 0; j<aaMod['%'].count; j++){
        //skip mods not allowed on this peptide
        if (xl && !aaMod['%'].mod[j].xl && !diffModsOnXL) continue;
        if (xl && aaMod['%'].mod[j].xl && !monoModsOnXL) continue;
        //Add masses
        addModIonSet(index, '%', pep1Len - 1 - link2 + link, j, pep1Len-1);
        //solve another one
        if (depth + 1<maxModCount) modLoopIonsRec(start, link, link2, (int)(sets.size()) - 1, depth + 1, xl);
      }
    } else { //Special case for peptide c-terinus
      for (j = 0; j<aaMod['c'].count; j++){
        //skip mods not allowed on this peptide
        if (xl && !aaMod['c'].mod[j].xl && !diffModsOnXL) continue;
        if (xl && aaMod['c'].mod[j].xl && !monoModsOnXL) continue;
        //Add masses
        addModIonSet(index, 'c', pep1Len-1-link2+link, j, pep1Len-1);
        //solve another one
        if (depth + 1<maxModCount) modLoopIonsRec(start, link, link2, (int)(sets.size()) - 1, depth + 1, xl);
      }
    }
  }

  //Not searching case where peptide n-terminus and c-terminus are both modified

}

void KIons::reset(){
  sets.clear();
  sets.push_back(ionBlank);
  sets[0].setIons(pep1Len,pep1Mass);
}

double KIons::getModMass(int index){
  return modMassArray[index];
}

int KIons::getModMassSize(){
  return (int)modMassArray.size();
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

void KIons::setAAMass(char aa, double mass, bool n15){
  if(n15) aaMassn15[aa]=mass;
  else aaMass[aa]=mass;
}

void KIons::setMaxModCount(int i){
  maxModCount=i;
}

void KIons::setModFlags(bool monoMods, bool difMods){
  monoModsOnXL=monoMods;
  diffModsOnXL=difMods;
}

void KIons::setPeptide(bool bPepOne, char* seq, int len, double mass, bool nTerm, bool cTerm, bool n15){
  if(bPepOne){
    pep1=seq;
    pep1Len=len;
    pep1Mass=mass;
    nPep1=nTerm;
    cPep1=cTerm;
    n15Pep1=n15;

    sets.clear();
    sets.push_back(ionBlank);
    sets[0].setIons(pep1Len,pep1Mass);

  } else {
    pep2=seq;
    pep2Len=len;
    pep2Mass=mass;
    nPep2=nTerm;
    cPep2=cTerm;
    n15Pep2 = n15;
  }

}

void KIons::setSeries(bool a, bool b, bool c, bool x, bool y, bool z){
  seriesA=a;
  seriesB=b;
  seriesC=c;
  seriesX=x;
  seriesY=y;
  seriesZ=z;
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

