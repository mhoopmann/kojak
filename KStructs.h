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

#ifndef _KSTRUCTS_H
#define _KSTRUCTS_H

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <iostream>
#include <exception>

using namespace std;

//FASTA database structure
typedef struct kDB{
  string name;      //FASTA header
  string sequence;  //FASTA sequence
} kDB;

//structure holds peptide mappings to database
typedef struct kPepMap{
  int index;  //protein index
  unsigned short start;  //first aa
  unsigned short stop;   //last aa
} kPepMap;

//Peptide reference to an entry in pldbDB
typedef struct kPeptide{
  double mass;            //monoisotopic, zero mass
  vector<kPepMap>* map;   //array of mappings where peptides appear in more than one place
  vector<int>* vA;        //List of linkage sites for set A
  vector<int>* vB;        //List of linkage sites for set B
  kPeptide(){
    mass=0;
    map=new vector<kPepMap>;
    vA=new vector<int>;
    vB=new vector<int>;
  }
  kPeptide(const kPeptide& m){
    mass=m.mass;
    map=new vector<kPepMap>;
    vA=new vector<int>;
    vB=new vector<int>;
    for(unsigned int i=0;i<m.map->size();i++) map->push_back(m.map->at(i));
    for(unsigned int i=0;i<m.vA->size();i++) vA->push_back(m.vA->at(i));
    for(unsigned int i=0;i<m.vB->size();i++) vB->push_back(m.vB->at(i));
  }
  ~kPeptide(){
    delete map;
    delete vA;
    delete vB;
  }
  kPeptide& operator=(const kPeptide& m){
    if(this!=&m){
      mass=m.mass;
      delete map;
      delete vA;
      delete vB;
      map=new vector<kPepMap>;
      vA=new vector<int>;
      vB=new vector<int>;
      for(unsigned int i=0;i<m.map->size();i++) map->push_back(m.map->at(i));
      for(unsigned int i=0;i<m.vA->size();i++) vA->push_back(m.vA->at(i));
      for(unsigned int i=0;i<m.vB->size();i++) vB->push_back(m.vB->at(i));
    }
    return (*this);
  }
} kPeptide;

typedef struct kPeptideB{
  double  mass;
  int     index;
  bool    linkable;
} kPeptideB;

//For sorting peptide lists
typedef struct kPepSort{
  int index;        //peptide array index
  string sequence;  //peptide sequence
} kPepSort;

typedef struct kLinker{
  double mass;
  int mono;     //0=cross-link, 1=mono-link
  int siteA;    //number represents site list. 1 for vA, 2 for vB.
  int siteB;    //number represents site list. 1 for vA, 2 for vB.
} kLinker;

typedef struct kMass {
  bool    xl;
  int     index;
  double  mass;
} kMass;

typedef struct kSparseMatrix{
  int   bin;
  float fIntensity;
} kSparseMatrix;

typedef struct kEnzymeRules{
  bool cutC[128];
  bool cutN[128];
  bool exceptC[128];
  bool exceptN[128];
} kEnzymeRules;

typedef struct kParams {
  //int     diagnostic;
  int     dimers;
  int     instrument;     //0=Orbi, 1=FTICR
  int     maxMods;
  int     maxPeaks;
  int     miscleave;
  int     ms1Centroid;
  int     ms2Centroid;
  int     ms1Resolution;
  int     ms2Resolution;
  int     preferPrecursor;
  int     setA;
  int     setB;
  int     specProcess;
  int     threads;
  int     topCount;
  int     truncate;
  bool    diffModsOnXL;
  bool    monoLinksOnXL;
  bool    xcorr;
  double  binOffset;
  double  binSize;
  double  enrichment;
  double  maxPepMass;
  double  minPepMass;
  double  percVersion;
  double  ppmPrecursor;
  char    dbFile[256];
  char    decoy[256];
  char    enzyme[32];
  char    msFile[256];
  char    outFile[256];
  char    percolator[256];
  vector<int>*      diag;
  vector<kLinker>*  xLink;
  vector<kLinker>*  mLink;
  vector<kMass>*    mods;
  vector<kMass>*    fMods;
  kParams(){
    //diagnostic=0;
    dimers=0;
    instrument=1;
    maxMods=0;
    maxPeaks=0;
    miscleave=2;
    ms1Centroid=0;
    ms2Centroid=0;
    ms1Resolution=60000;
    ms2Resolution=15000;
    preferPrecursor=0;
    setA=0;
    setB=0;
    specProcess=0;
    threads=1;
    topCount=250;
    truncate=0;
    diffModsOnXL=false;
    monoLinksOnXL=false;
    xcorr=false;
    binSize=0.03;
    binOffset=0.0;
    enrichment=0;
    maxPepMass=4000.0;
    minPepMass=500.0;
    percVersion=2.04;
    ppmPrecursor=25.0;
    strcpy(decoy,"random");
    dbFile[0]='\0';
    strcpy(enzyme,"[KR]|");
    msFile[0]='\0';
    outFile[0]='\0';
    percolator[0]='\0';
    diag = new vector<int>;
    xLink = new vector<kLinker>;
    mLink = new vector<kLinker>;
    mods = new vector<kMass>;
    fMods = new vector<kMass>;
  }
  kParams(const kParams& p){
    //diagnostic=p.diagnostic;
    dimers=p.dimers;
    instrument=p.instrument;
    maxMods=p.maxMods;
    maxPeaks=p.maxPeaks;
    miscleave=p.miscleave;
    ms1Centroid=p.ms1Centroid;
    ms2Centroid=p.ms2Centroid;
    ms1Resolution=p.ms1Resolution;
    ms2Resolution=p.ms2Resolution;
    preferPrecursor=p.preferPrecursor;
    setA=p.setA;
    setB=p.setB;
    specProcess=p.specProcess;
    threads=p.threads;
    topCount=p.topCount;
    truncate=p.truncate;
    diffModsOnXL=p.diffModsOnXL;
    monoLinksOnXL=p.monoLinksOnXL;
    xcorr=p.xcorr;
    binOffset=p.binOffset;
    binSize=p.binSize;
    enrichment=p.enrichment;
    maxPepMass=p.maxPepMass;
    minPepMass=p.minPepMass;
    percVersion=p.percVersion;
    ppmPrecursor=p.ppmPrecursor;
    strcpy(decoy,p.decoy);
    strcpy(dbFile,p.dbFile);
    strcpy(enzyme,p.enzyme);
    strcpy(msFile,p.msFile);
    strcpy(outFile,p.outFile);
    strcpy(percolator,p.percolator);
    diag = new vector<int>;
    xLink = new vector<kLinker>;
    mLink = new vector<kLinker>;
    mods = new vector<kMass>;
    fMods = new vector<kMass>;
    unsigned int i;
    for(i=0;i<p.diag->size();i++) diag->push_back(p.diag->at(i));
    for(i=0;i<p.xLink->size();i++) xLink->push_back(p.xLink->at(i));
    for(i=0;i<p.mLink->size();i++) mLink->push_back(p.mLink->at(i));
    for(i=0;i<p.mods->size();i++) mods->push_back(p.mods->at(i));
    for(i=0;i<p.fMods->size();i++) fMods->push_back(p.fMods->at(i));
  }
  ~kParams(){
    delete diag;
    delete xLink;
    delete mLink;
    delete mods;
    delete fMods;
  }
  kParams& operator=(const kParams& p){
    if(this!=&p){
      //diagnostic=p.diagnostic;
      dimers=p.dimers;
      instrument=p.instrument;
      maxMods=p.maxMods;
      maxPeaks=p.maxPeaks;
      miscleave=p.miscleave;
      ms1Centroid=p.ms1Centroid;
      ms2Centroid=p.ms2Centroid;
      ms1Resolution=p.ms1Resolution;
      ms2Resolution=p.ms2Resolution;
      preferPrecursor=p.preferPrecursor;
      setA=p.setA;
      setB=p.setB;
      specProcess=p.specProcess;
      threads=p.threads;
      topCount=p.topCount;
      truncate=p.truncate;
      diffModsOnXL=p.diffModsOnXL;
      monoLinksOnXL=p.monoLinksOnXL;
      xcorr=p.xcorr;
      binOffset=p.binOffset;
      binSize=p.binSize;
      enrichment=p.enrichment;
      maxPepMass=p.maxPepMass;
      minPepMass=p.minPepMass;
      percVersion=p.percVersion;
      ppmPrecursor=p.ppmPrecursor;
      strcpy(decoy,p.decoy);
      strcpy(dbFile,p.dbFile);
      strcpy(enzyme,p.enzyme);
      strcpy(msFile,p.msFile);
      strcpy(outFile,p.outFile);
      strcpy(percolator,p.percolator);
      delete diag;
      delete xLink;
      delete mLink;
      delete mods;
      delete fMods;
      diag = new vector<int>;
      xLink = new vector<kLinker>;
      mLink = new vector<kLinker>;
      mods = new vector<kMass>;
      fMods = new vector<kMass>;
      unsigned int i;
      for(i=0;i<p.diag->size();i++) diag->push_back(p.diag->at(i));
      for(i=0;i<p.xLink->size();i++) xLink->push_back(p.xLink->at(i));
      for(i=0;i<p.mLink->size();i++) mLink->push_back(p.mLink->at(i));
      for(i=0;i<p.mods->size();i++) mods->push_back(p.mods->at(i));
      for(i=0;i<p.fMods->size();i++) fMods->push_back(p.fMods->at(i));
    }
    return (*this);
  }
} kParams;

typedef struct kSpecPoint{
  double mass;
  float intensity;
} kSpecPoint;

typedef struct kPreprocessStruct { //adapted from Comet
   int iHighestIon;
   double dHighestIntensity;
   kSpecPoint *pdCorrelationData;
} kPreprocessStruct;

typedef struct kPepMod{
  char pos;
  double mass;
} kPepMod;

typedef struct kScoreCard{
  bool    linkable1;
  bool    linkable2;
  int     k1;
  int     k2;
  int     link;
  int     pep1;
  int     pep2;
  int     rank;
  float   simpleScore;
  double  mass;
  double  scoreDiff;
  vector<kPepMod>* mods1;
  vector<kPepMod>* mods2;
  kScoreCard(){
    linkable1=false;
    linkable2=false;
    k1=0;
    k2=0;
    link=0;
    pep1=0;
    pep2=0;
    rank=0;
    simpleScore=0;
    mass=0;
    scoreDiff=0;
    mods1=new vector<kPepMod>;
    mods2=new vector<kPepMod>;
    //cout << mods1 << endl;
    //cout << mods2 << endl;
  }
  kScoreCard(const kScoreCard& p){
    unsigned int i;
    linkable1=p.linkable1;
    linkable2=p.linkable2;
    k1=p.k1;
    k2=p.k2;
    link=p.link;
    pep1=p.pep1;
    pep2=p.pep2;
    rank=p.rank;
    simpleScore=p.simpleScore;
    mass=p.mass;
    scoreDiff=p.scoreDiff;
    mods1=new vector<kPepMod>;
    for(i=0;i<p.mods1->size();i++) mods1->push_back(p.mods1->at(i));
    mods2=new vector<kPepMod>;
    for(i=0;i<p.mods2->size();i++) mods2->push_back(p.mods2->at(i));
  }
  ~kScoreCard(){
    delete mods1;
    delete mods2;
  }
  kScoreCard& operator=(const kScoreCard& p){
    if(this!=&p){
      unsigned int i;
      linkable1=p.linkable1;
      linkable2=p.linkable2;
      k1=p.k1;
      k2=p.k2;
      link=p.link;
      pep1=p.pep1;
      pep2=p.pep2;
      rank=p.rank;
      simpleScore=p.simpleScore;
      mass=p.mass;
      scoreDiff=p.scoreDiff;
      delete mods1;
      mods1=new vector<kPepMod>;
      for(i=0;i<p.mods1->size();i++) mods1->push_back(p.mods1->at(i));
      delete mods2;
      mods2=new vector<kPepMod>;
      for(i=0;i<p.mods2->size();i++) mods2->push_back(p.mods2->at(i));
    }
    return *this;
  }
} kScoreCard;

typedef struct kSingletScoreCard{
  char                len;
  bool                linkable;
  char                k1;
  int                 pep1;
  float               simpleScore;
  double              mass;
  char                modLen;
  kPepMod*            mods;
  kSingletScoreCard*  next;
  kSingletScoreCard*  prev;
  kSingletScoreCard(){
    len=0;
    linkable=false;
    k1=0;
    pep1=0;
    simpleScore=0;
    mass=0;
    modLen=0;
    mods=NULL;
    next=NULL;
    prev=NULL;
  }
  kSingletScoreCard(const kSingletScoreCard& k){
    len=k.len;
    linkable=k.linkable;
    k1=k.k1;
    pep1=k.pep1;
    simpleScore=k.simpleScore;
    mass=k.mass;
    modLen=k.modLen;
    mods=NULL;
    if(modLen>0){
      mods=new kPepMod[modLen];
      for(char i=0;i<modLen;i++) mods[i]=k.mods[i];
    }
    next=NULL;
    prev=NULL;
  }
  ~kSingletScoreCard(){
    if(mods!=NULL) delete [] mods;
    next=NULL;
    prev=NULL;
  }
  kSingletScoreCard& operator=(const kSingletScoreCard& k){
    if(this!=&k){
      len=k.len;
      linkable=k.linkable;
      k1=k.k1;
      pep1=k.pep1;
      simpleScore=k.simpleScore;
      mass=k.mass;
      modLen=k.modLen;
      if(mods!=NULL) {
        delete [] mods;
        mods=NULL;
      }
      if(modLen>0){
        mods=new kPepMod[modLen];
        for(char i=0;i<modLen;i++) mods[i]=k.mods[i];
      }
      next=NULL;
      prev=NULL;
    }
    return *this;
  }
} kSingletScoreCard;

typedef struct kSingletScoreCardPlus{
  char    len;
  bool    linkable;
  char    k1;
  double  mass;
  int     pep1;
  int     rank;
  float   simpleScore;
  int     target;
} kSingletScoreCardPlus;

typedef struct kPrecursor{ 
  int     charge;
  double  corr;
  char    label;
  double  monoMass;
  kPrecursor(){
    charge=0;
    corr=0;
    label=0;
    monoMass=0;
  }
} kPrecursor;

typedef struct kResults{
  bool    decoy;
  bool    linkable1;
  bool    linkable2;
  int     charge;
  int     link1;
  int     link2;
  int     pep1;
  int     pep2;
  int     rank;
  int     scanNumber;
  int     type;
  double  obsMass;
  double  ppm;
  double  psmMass;
  double  score;
  double  scoreDelta;
  double  scorePepDif;
  string  modPeptide1;
  string  modPeptide2;
  string  peptide1;
  string  peptide2;
} kResults;

#endif
