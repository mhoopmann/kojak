#ifndef _KSTRUCTS_H
#define _KSTRUCTS_H

#include <string>
#include <vector>

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
  vector<kPepMap>* map;  //array of mappings where peptides appear in more than one place
  vector<int>* vK;
  kPeptide(){
    mass=0;
    map=new vector<kPepMap>;
    vK=new vector<int>;
  }
  kPeptide(const kPeptide& m){
    mass=m.mass;
    map=new vector<kPepMap>;
    vK=new vector<int>;
    for(unsigned int i=0;i<m.map->size();i++) map->push_back(m.map->at(i));
    for(unsigned int i=0;i<m.vK->size();i++) vK->push_back(m.vK->at(i));
  }
  ~kPeptide(){
    delete map;
    delete vK;
  }
  kPeptide& operator=(const kPeptide& m){
    if(this!=&m){
      mass=m.mass;
      delete map;
      delete vK;
      map=new vector<kPepMap>;
      vK=new vector<int>;
      for(unsigned int i=0;i<m.map->size();i++) map->push_back(m.map->at(i));
      for(unsigned int i=0;i<m.vK->size();i++) vK->push_back(m.vK->at(i));
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
} kLinker;

typedef struct kMass {
  int     index;
  double  mass;
} kMass;

typedef struct kSparseMatrix{
  int   bin;
  float fIntensity;
} kSparseMatrix;

typedef struct kParams {
  int     diagnostic;
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
  int     relaxedAnalysis;
  int     specProcess;
  int     topCount;
  int     truncate;
  int     xcorr;
  double  enrichment;
  double  maxPepMass;
  double  minPepMass;
  double  percVersion;
  double  ppmPrecursor;
  double  ppmFragment;
  char    dbFile[256];
  char    decoy[256];
  char    msFile[256];
  char    outFile[256];
  char    percolator[256];
  vector<double>* xLink;
  vector<double>* mLink;
  vector<kMass>* mods;
  vector<kMass>* fMods;
  kParams(){
    diagnostic=0;
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
    relaxedAnalysis=0;
    specProcess=0;
    topCount=250;
    truncate=0;
    xcorr=0;
    enrichment=0;
    maxPepMass=4000.0;
    minPepMass=500.0;
    percVersion=2.04;
    ppmPrecursor=25.0;
    ppmFragment=25.0;
    strcpy(decoy,"random");
    dbFile[0]='\0';
    msFile[0]='\0';
    outFile[0]='\0';
    percolator[0]='\0';
    xLink = new vector<double>;
    mLink = new vector<double>;
    mods = new vector<kMass>;
    fMods = new vector<kMass>;
  }
  kParams(const kParams& p){
    diagnostic=p.diagnostic;
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
    relaxedAnalysis=p.relaxedAnalysis;
    specProcess=p.specProcess;
    topCount=p.topCount;
    truncate=p.truncate;
    xcorr=p.xcorr;
    enrichment=p.enrichment;
    maxPepMass=p.maxPepMass;
    minPepMass=p.minPepMass;
    percVersion=p.percVersion;
    ppmPrecursor=p.ppmPrecursor;
    ppmFragment=p.ppmFragment;
    strcpy(decoy,p.decoy);
    strcpy(dbFile,p.dbFile);
    strcpy(msFile,p.msFile);
    strcpy(outFile,p.outFile);
    strcpy(percolator,p.percolator);
    xLink = new vector<double>;
    mLink = new vector<double>;
    mods = new vector<kMass>;
    fMods = new vector<kMass>;
    unsigned int i;
    for(i=0;i<p.xLink->size();i++) xLink->push_back(p.xLink->at(i));
    for(i=0;i<p.mLink->size();i++) mLink->push_back(p.mLink->at(i));
    for(i=0;i<p.mods->size();i++) mods->push_back(p.mods->at(i));
    for(i=0;i<p.fMods->size();i++) fMods->push_back(p.fMods->at(i));
  }
  ~kParams(){
    delete xLink;
    delete mLink;
    delete mods;
    delete fMods;
  }
  kParams& operator=(const kParams& p){
    if(this!=&p){
      diagnostic=p.diagnostic;
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
      relaxedAnalysis=p.relaxedAnalysis;
      specProcess=p.specProcess;
      topCount=p.topCount;
      truncate=p.truncate;
      xcorr=p.xcorr;
      enrichment=p.enrichment;
      maxPepMass=p.maxPepMass;
      minPepMass=p.minPepMass;
      percVersion=p.percVersion;
      ppmPrecursor=p.ppmPrecursor;
      ppmFragment=p.ppmFragment;
      strcpy(decoy,p.decoy);
      strcpy(dbFile,p.dbFile);
      strcpy(msFile,p.msFile);
      strcpy(outFile,p.outFile);
      strcpy(percolator,p.percolator);
      delete xLink;
      delete mLink;
      delete mods;
      delete fMods;
      xLink = new vector<double>;
      mLink = new vector<double>;
      mods = new vector<kMass>;
      fMods = new vector<kMass>;
      unsigned int i;
      for(i=0;i<p.xLink->size();i++) xLink->push_back(p.xLink->at(i));
      for(i=0;i<p.mLink->size();i++) mLink->push_back(p.mLink->at(i));
      for(i=0;i<p.mods->size();i++) mods->push_back(p.mods->at(i));
      for(i=0;i<p.fMods->size();i++) fMods->push_back(p.fMods->at(i));
    }
    return (*this);
  }
} kParams;

typedef struct kPreprocessStruct { //adapted from Comet
   int iHighestIon;
   double dHighestIntensity;
   double *pdCorrelationData;
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
  vector<kPepMod>* mods;
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
    mods=new vector<kPepMod>;
  }
  kScoreCard(const kScoreCard& p){
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
    mods=new vector<kPepMod>;
    for(unsigned int i=0;i<p.mods->size();i++) mods->push_back(p.mods->at(i));
  }
  ~kScoreCard(){
    delete mods;
  }
  kScoreCard& operator=(const kScoreCard& p){
    if(this!=&p){
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
      delete mods;
      mods=new vector<kPepMod>;
      for(unsigned int i=0;i<p.mods->size();i++) mods->push_back(p.mods->at(i));
    }
    return *this;
  }
} kScoreCard;

typedef struct kSingletScoreCard{
  char  len;
  bool  linkable;
  int   k1;
  int   pep1;
  float simpleScore;
  kSingletScoreCard(){
    len=0;
    linkable=false;
    k1=0;
    pep1=0;
    simpleScore=0;
  }
} kSingletScoreCard;

typedef struct kSingletScoreCardPlus{
  char    len;
  bool    linkable;
  int     k1;
  double  mass;
  int     pep1;
  int     rank;
  float   simpleScore;
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

#endif
