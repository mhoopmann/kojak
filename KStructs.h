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

//FASTA database structure
typedef struct kDB{
  bool decoy;
  std::string name;      //FASTA header
  std::string sequence;  //FASTA sequence
} kDB;

typedef struct kFile{
  std::string input;
  std::string base;
  std::string ext;
} kFile;

//structure holds peptide mappings to database
typedef struct kPepMap{
  int index;  //protein index
  unsigned short start;  //first aa
  unsigned short stop;   //last aa
} kPepMap;

//Peptide reference to an entry in pldbDB
typedef struct kPeptide{
  bool cTerm;
  bool nTerm;
  bool n15;
  char xlSites;
  double mass;            //monoisotopic, zero mass
  std::vector<kPepMap>* map;   //array of mappings where peptides appear in more than one place
  kPeptide(){
    cTerm=false;
    nTerm=false;
    n15=false;
    xlSites=0;
    mass=0;
    map = new std::vector<kPepMap>;
  }
  kPeptide(const kPeptide& m){
    cTerm=m.cTerm;
    nTerm=m.nTerm;
    n15=m.n15;
    xlSites=m.xlSites;
    mass=m.mass;
    map = new std::vector<kPepMap>(*m.map);
  }
  ~kPeptide(){
    delete map;
  }
  kPeptide& operator=(const kPeptide& m){
    if(this!=&m){
      cTerm = m.cTerm;
      nTerm = m.nTerm;
      n15=m.n15;
      xlSites = m.xlSites;
      mass=m.mass;
      delete map;
      map = new std::vector<kPepMap>(*m.map);
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
  std::string sequence;  //peptide sequence
  bool n15;
} kPepSort;

typedef struct kLinker{
  std::string label;
  std::string motifA;
  std::string motifB;
  double mass;
  int mono;     //0=cross-link, 1=mono-link
  int motifAIndex;    //for reference to motif list
  int motifBIndex;    //for reference to motif list
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
  int     decoySize;
  int     instrument;     //0=Orbi, 1=FTICR
  int     intermediate;
  int     isotopeError;
  int     maxMods;
  int     maxPeaks;
  int     minPeaks;
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
  bool    buildDecoy;
  bool    diffModsOnXL;
  bool    dimers;
  bool    dimersXL;
  bool    exportMzID;
  bool    exportPepXML;
  bool    exportPercolator;
  bool    ionSeries[6];
  bool    monoLinksOnXL;
  bool    precursorRefinement;
  bool    turbo;
  bool    xcorr;
  double  binOffset;
  double  binSize;
  double  enrichment;
  double  maxPepMass;
  double  minPepMass;
  double  minPepScore;
  double  percVersion;
  double  ppmPrecursor;
  double  removePrecursor;
  char    dbFile[256];
  char    decoy[256];
  char    enzyme[32];
  char    enzymeName[64];
  char    ext[32];
  char    inFile[1024];  //true input file with full path
  char    msFile[256];   //input file parameter from confic
  char    n15Label[256];
  char    outFile[1024];  //true output file with full path
  char    resPath[1024];
  std::string fullPath;
  std::vector<kMass>*    aaMass;
  std::vector<double>*   cleavageProducts;
  std::vector<int>*      diag;
  std::vector<kLinker>*  xLink;
  std::vector<kLinker>*  mLink;
  std::vector<kMass>*    mods;
  std::vector<kMass>*    fMods;
  kParams(){
    decoySize=5000;
    instrument=1;
    intermediate=0;
    isotopeError=1;
    maxMods=0;
    maxPeaks=0;
    minPeaks = 12;
    miscleave=2;
    ms1Centroid=0;
    ms2Centroid=0;
    ms1Resolution=60000;
    ms2Resolution=15000;
    preferPrecursor=1;
    setA=0;
    setB=0;
    specProcess=0;
    threads=1;
    topCount=250;
    truncate=0;
    buildDecoy = false;
    diffModsOnXL=false;
    dimers=false;
    dimersXL=true;
    exportMzID=false;
    exportPepXML=false;
    exportPercolator=false;
    ionSeries[0]=false; //a-ions
    ionSeries[1]=true;  //b-ions
    ionSeries[2]=false; //c-ions
    ionSeries[3]=false; //x-ions
    ionSeries[4]=true;  //y-ions
    ionSeries[5]=false; //z-ions
    monoLinksOnXL=false;
    precursorRefinement=true;
    turbo=true;
    xcorr=false;
    binSize=0.03;
    binOffset=0.0;
    enrichment=0;
    maxPepMass=4000.0;
    minPepMass=500.0;
    minPepScore=0.1;
    percVersion=2.04;
    ppmPrecursor=25.0;
    removePrecursor=0;
    strcpy(decoy,"random");
    dbFile[0]='\0';
    strcpy(enzyme,"[KR]|{P}");
    strcpy(enzymeName, "Trypsin");
    ext[0]='\0';
    inFile[0]='\0';
    msFile[0]='\0';
    n15Label[0]='\0';
    outFile[0]='\0';
    resPath[0]='\0';
    fullPath = "";
    aaMass = new std::vector<kMass>;
    cleavageProducts = new std::vector<double>;
    diag = new std::vector<int>;
    xLink = new std::vector<kLinker>;
    mLink = new std::vector<kLinker>;
    mods = new std::vector<kMass>;
    fMods = new std::vector<kMass>;
  }
  kParams(const kParams& p){
    decoySize=p.decoySize;
    instrument=p.instrument;
    intermediate=p.intermediate;
    isotopeError=p.isotopeError;
    maxMods=p.maxMods;
    maxPeaks=p.maxPeaks;
    minPeaks = p.minPeaks;
    miscleave=p.miscleave;
    ms1Centroid=p.ms1Centroid;
    ms2Centroid=p.ms2Centroid;
    ms1Resolution=p.ms1Resolution;
    ms2Resolution=p.ms2Resolution;
    preferPrecursor=p.preferPrecursor;
    removePrecursor=p.removePrecursor;
    setA=p.setA;
    setB=p.setB;
    specProcess=p.specProcess;
    threads=p.threads;
    topCount=p.topCount;
    truncate=p.truncate;
    buildDecoy = p.buildDecoy;
    diffModsOnXL=p.diffModsOnXL;
    dimers=p.dimers;
    dimersXL=p.dimersXL;
    exportMzID=p.exportMzID;
    exportPepXML=p.exportPepXML;
    exportPercolator=p.exportPercolator;
    monoLinksOnXL=p.monoLinksOnXL;
    precursorRefinement=p.precursorRefinement;
    turbo=p.turbo;
    xcorr=p.xcorr;
    binOffset=p.binOffset;
    binSize=p.binSize;
    enrichment=p.enrichment;
    maxPepMass=p.maxPepMass;
    minPepMass=p.minPepMass;
    minPepScore=p.minPepScore;
    percVersion=p.percVersion;
    ppmPrecursor=p.ppmPrecursor;
    strcpy(decoy,p.decoy);
    strcpy(dbFile,p.dbFile);
    strcpy(enzyme,p.enzyme);
    strcpy(enzymeName,p.enzymeName);
    strcpy(ext,p.ext);
    strcpy(inFile, p.inFile);
    strcpy(msFile,p.msFile);
    strcpy(n15Label, p.n15Label);
    strcpy(outFile,p.outFile);
    strcpy(resPath, p.resPath);
    fullPath = p.fullPath;
    aaMass = new std::vector<kMass>(*p.aaMass);
    cleavageProducts = new std::vector<double>(*p.cleavageProducts);
    diag = new std::vector<int>(*p.diag);
    xLink = new std::vector<kLinker>(*p.xLink);
    mLink = new std::vector<kLinker>(*p.mLink);
    mods = new std::vector<kMass>(*p.mods);
    fMods = new std::vector<kMass>(*p.fMods);
    for(unsigned int i=0;i<6;i++) ionSeries[i]=p.ionSeries[i];
  }
  ~kParams(){
    delete aaMass;
    delete cleavageProducts;
    delete diag;
    delete xLink;
    delete mLink;
    delete mods;
    delete fMods;
  }
  kParams& operator=(const kParams& p){
    if(this!=&p){
      decoySize=p.decoySize;
      instrument=p.instrument;
      intermediate=p.intermediate;
      isotopeError = p.isotopeError;
      maxMods=p.maxMods;
      maxPeaks=p.maxPeaks;
      minPeaks = p.minPeaks;
      miscleave=p.miscleave;
      ms1Centroid=p.ms1Centroid;
      ms2Centroid=p.ms2Centroid;
      ms1Resolution=p.ms1Resolution;
      ms2Resolution=p.ms2Resolution;
      preferPrecursor=p.preferPrecursor;
      removePrecursor=p.removePrecursor;
      setA=p.setA;
      setB=p.setB;
      specProcess=p.specProcess;
      threads=p.threads;
      topCount=p.topCount;
      truncate=p.truncate;
      buildDecoy = p.buildDecoy;
      diffModsOnXL=p.diffModsOnXL;
      dimers=p.dimers;
      dimersXL=p.dimersXL;
      exportMzID = p.exportMzID;
      exportPepXML=p.exportPepXML;
      exportPercolator=p.exportPercolator;
      monoLinksOnXL=p.monoLinksOnXL;
      precursorRefinement = p.precursorRefinement;
      turbo = p.turbo;
      xcorr=p.xcorr;
      binOffset=p.binOffset;
      binSize=p.binSize;
      enrichment=p.enrichment;
      maxPepMass=p.maxPepMass;
      minPepMass=p.minPepMass;
      minPepScore=p.minPepScore;
      percVersion=p.percVersion;
      ppmPrecursor=p.ppmPrecursor;
      strcpy(decoy,p.decoy);
      strcpy(dbFile,p.dbFile);
      strcpy(enzyme,p.enzyme);
      strcpy(enzymeName, p.enzymeName);
      strcpy(ext,p.ext);
      strcpy(inFile, p.inFile);
      strcpy(msFile,p.msFile);
      strcpy(n15Label, p.n15Label);
      strcpy(outFile,p.outFile);
      strcpy(resPath, p.resPath);
      fullPath = p.fullPath;
      delete aaMass;
      delete cleavageProducts;
      delete diag;
      delete xLink;
      delete mLink;
      delete mods;
      delete fMods;
      aaMass = new std::vector<kMass>(*p.aaMass);
      cleavageProducts = new std::vector<double>(*p.cleavageProducts);
      diag = new std::vector<int>(*p.diag);
      xLink = new std::vector<kLinker>(*p.xLink);
      mLink = new std::vector<kLinker>(*p.mLink);
      mods = new std::vector<kMass>(*p.mods);
      fMods = new std::vector<kMass>(*p.fMods);
      for(unsigned int i=0;i<6;i++) ionSeries[i]=p.ionSeries[i];
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
  char    precursor;
  char    site1;
  char    site2;
  int     k1;
  int     k2;
  int     link;
  int     pep1;
  int     pep2;
  int     matches1;
  int     matches2;
  int     conFrag1;
  int     conFrag2;
  float   simpleScore;
  double  eVal;
  double  eVal1;
  double  eVal2;
  double  mass;
  double  mass1;
  double  mass2;
  float  score1;
  float  score2;
  float  cpScore1;
  float  cpScore2;
  std::vector<kPepMod>* mods1;
  std::vector<kPepMod>* mods2;
  kScoreCard(){
    linkable1=false;
    linkable2=false;
    precursor=0;
    site1=0;
    site2=0;
    k1=0;
    k2=0;
    link=0;
    pep1=0;
    pep2=0;
    simpleScore=0;
    eVal=0;
    eVal1=0;
    eVal2=0;
    mass=0;
    mass1=0;
    mass2=0;
    score1=0;
    score2=0;
    cpScore1=0;
    cpScore2=0;
    matches1=0;
    matches2=0;
    conFrag1=0;
    conFrag2=0;
    mods1 = new std::vector<kPepMod>;
    mods2 = new std::vector<kPepMod>;
  }
  kScoreCard(const kScoreCard& p){
    linkable1=p.linkable1;
    linkable2=p.linkable2;
    precursor=p.precursor;
    site1 = p.site1;
    site2 = p.site2;
    k1=p.k1;
    k2=p.k2;
    link=p.link;
    pep1=p.pep1;
    pep2=p.pep2;
    simpleScore=p.simpleScore;
    eVal=p.eVal;
    eVal1=p.eVal1;
    eVal2=p.eVal2;
    mass=p.mass;
    mass1=p.mass1;
    mass2=p.mass2;
    score1=p.score1;
    score2=p.score2;
    cpScore1=p.cpScore1;
    cpScore2=p.cpScore2;
    matches1=p.matches1;
    matches2=p.matches2;
    conFrag1=p.conFrag1;
    conFrag2=p.conFrag2;
    mods1 = new std::vector<kPepMod>(*p.mods1);
    mods2 = new std::vector<kPepMod>(*p.mods2);
  }
  ~kScoreCard(){
    delete mods1;
    delete mods2;
  }
  kScoreCard& operator=(const kScoreCard& p){
    if(this!=&p){
      linkable1=p.linkable1;
      linkable2=p.linkable2;
      precursor=p.precursor;
      site1 = p.site1;
      site2 = p.site2;
      k1=p.k1;
      k2=p.k2;
      link=p.link;
      pep1=p.pep1;
      pep2=p.pep2;
      simpleScore=p.simpleScore;
      eVal=p.eVal;
      eVal1 = p.eVal1;
      eVal2 = p.eVal2;
      mass=p.mass;
      mass1=p.mass1;
      mass2=p.mass2;
      score1=p.score1;
      score2=p.score2;
      cpScore1 = p.cpScore1;
      cpScore2 = p.cpScore2;
      matches1 = p.matches1;
      matches2 = p.matches2;
      conFrag1 = p.conFrag1;
      conFrag2 = p.conFrag2;
      delete mods1;
      mods1 = new std::vector<kPepMod>(*p.mods1);
      delete mods2;
      mods2 = new std::vector<kPepMod>(*p.mods2);
    }
    return *this;
  }
} kScoreCard;

typedef struct kSingletScoreCard{
  char                len;
  bool                linkable;
  char                k1;
  int                 conFrag;  //longest run of consecutive fragment ions
  int                 matches;  //total matched fragment ions
  int                 pep1;
  char                pre;
  float               simpleScore;
  float               cpScore;
  double              mass;
  char                modLen;
  char                site;
  kPepMod*            mods;
  kSingletScoreCard*  next;
  kSingletScoreCard*  prev;
  kSingletScoreCard(){
    len=0;
    linkable=false;
    k1=0;
    pep1=0;
    pre=0;
    simpleScore=0;
    cpScore=0;
    mass=0;
    modLen=0;
    site=0;
    matches=0;
    conFrag=0;
    mods=NULL;
    next=NULL;
    prev=NULL;
  }
  kSingletScoreCard(const kSingletScoreCard& k){
    len=k.len;
    linkable=k.linkable;
    k1=k.k1;
    pep1=k.pep1;
    pre=k.pre;
    simpleScore=k.simpleScore;
    cpScore=k.cpScore;
    mass=k.mass;
    modLen=k.modLen;
    site=k.site;
    matches=k.matches;
    conFrag=k.conFrag;
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
      pre=k.pre;
      simpleScore=k.simpleScore;
      cpScore=k.cpScore;
      mass=k.mass;
      modLen=k.modLen;
      site=k.site;
      matches = k.matches;
      conFrag = k.conFrag;
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
  char    motif[20];
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
  bool    decoy1;
  bool    decoy2;
  bool    linkable1;
  bool    linkable2;
  bool    cTerm1; //peptide contains protein c-terminus
  bool    nTerm1; //peptide contains protein n-terminus
  bool    cTerm2; //peptide contains protein c-terminus
  bool    nTerm2; //peptide contains protein n-terminus
  bool    n15Pep1;
  bool    n15Pep2;
  char    linkSite1;  //actual amino acid that is linked
  char    linkSite2;  //actual amino acid that is linked
  int     charge;
  int     link1;
  int     link2;
  int     linkerID;
  int     matches1;
  int     matches2;
  int     conFrag1;
  int     conFrag2;
  int     pep1;
  int     pep2;
  int     rank;
  int     rankA;
  int     rankB;
  int     scanNumber;
  int     type;  //0=single, 1=loop, 2=xl, 3=dimer
  float   rTime;
  double  eVal;
  double  eVal1;
  double  eVal2;
  double  hk;
  double  massA;
  double  massB;
  double  obsMass;
  double  ppm;
  double  psmMass;
  double  score;
  double  scoreA;
  double  scoreB;
  double  scoreDelta;
  double  scorePepDif;
  double  xlMass;
  std::string  baseName;
  std::string  modPeptide1;
  std::string  modPeptide2;
  std::string  peptide1;
  std::string  peptide2;
  std::string  protein1;
  std::string  protein2;
  std::string  protPos1;
  std::string  protPos2;
  std::string  scanID;
  std::string  xlLabel;
  std::vector<kPepMod> mods1;
  std::vector<kPepMod> mods2;
} kResults;

typedef struct kCKey{
  int key;
  int pos;
} kCKey;

typedef struct kMatchSet{
  int sz;
  int charge;
  kCKey** a;
  kCKey** b;
  kCKey** c;
  kCKey** x;
  kCKey** y;
  kCKey** z;
  kMatchSet(){
    sz=0;
    charge=0;
    a=NULL;
    b=NULL;
    c=NULL;
    x=NULL;
    y=NULL;
    z=NULL;
  }
  ~kMatchSet(){
    deAllocate();
  }
  void allocate(int s, int ch){
    deAllocate();
    sz=s;
    charge=ch;
    //kCKey k;
    //k.key=0;
    //k.pos=0;
    a = new kCKey*[charge];
    b = new kCKey*[charge];
    c = new kCKey*[charge];
    x = new kCKey*[charge];
    y = new kCKey*[charge];
    z = new kCKey*[charge];
    for(int i=0;i<charge;i++){
      a[i] = new kCKey[sz];
      b[i] = new kCKey[sz];
      c[i] = new kCKey[sz];
      x[i] = new kCKey[sz];
      y[i] = new kCKey[sz];
      z[i] = new kCKey[sz];
      //for(int j=0;j<sz;j++){
      //  a[i][j]=k;
      //  b[i][j]=k;
      //  c[i][j]=k;
      //  x[i][j]=k;
      //  y[i][j]=k;
      //  z[i][j]=k;
      //}
    }
  }
  void deAllocate(){
    int i;
    if(a!=NULL){
      for(i=0;i<charge;i++) delete [] a[i];
      delete [] a;
    }
    if(b!=NULL){
      for(i=0;i<charge;i++) delete [] b[i];
      delete [] b;
    }
    if(c!=NULL){
      for(i=0;i<charge;i++) delete [] c[i];
      delete [] c;
    }
    if(x!=NULL){
      for(i=0;i<charge;i++) delete [] x[i];
      delete [] x;
    }
    if(y!=NULL){
      for(i=0;i<charge;i++) delete [] y[i];
      delete [] y;
    }
    if(z!=NULL){
      for(i=0;i<charge;i++) delete [] z[i];
      delete [] z;
    }
    a=NULL;
    b=NULL;
    c=NULL;
    x=NULL;
    y=NULL;
    z=NULL;
  }

} kMatchSet;

typedef struct kXLMotif {
  std::string motif;
  int xlIndex[10];
  int counterMotif[10];
} kXLMotif;

typedef struct kXLTarget{
  size_t linkerID;
  bool target[128];
} kXLTarget;

#endif
