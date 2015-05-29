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

#include "KAnalysis.h"

bool*       KAnalysis::bKIonsManager;
KDatabase*  KAnalysis::db;
double      KAnalysis::highLinkMass;
KIons*      KAnalysis::ions;
double      KAnalysis::lowLinkMass;
double      KAnalysis::maxMass;
double      KAnalysis::minMass;
Mutex       KAnalysis::mutexKIonsManager;
Mutex*      KAnalysis::mutexSpecScore;
kParams     KAnalysis::params;
KData*      KAnalysis::spec;

/*============================
  Constructors & Destructors
============================*/
KAnalysis::KAnalysis(kParams& p, KDatabase* d, KData* dat){
  unsigned int i;
  int j;
  
  //Assign pointers and structures
  params=p;
  db=d;
  spec=dat;

  //Do memory allocations and initialization
  bKIonsManager=NULL;
  ions=NULL;
  allocateMemory(params.threads);
  for(j=0;j<params.threads;j++){
    for(i=0;i<params.fMods->size();i++) ions[j].addFixedMod((char)params.fMods->at(i).index,params.fMods->at(i).mass);
    for(i=0;i<params.mods->size();i++) ions[j].addMod((char)params.mods->at(i).index,params.mods->at(i).xl,params.mods->at(i).mass);
    ions[j].setMaxModCount(params.maxMods);
  }

  //Initalize variables
  maxMass = spec->getMaxMass()+0.25;
  minMass = spec->getMinMass()-0.25;
  for(j=0;j<spec->sizeLink();j++){
    if(spec->getLink(j).mono==0){
      if(lowLinkMass==0) lowLinkMass=spec->getLink(j).mass;
      if(highLinkMass==0) highLinkMass=spec->getLink(j).mass;
      if(spec->getLink(j).mass<lowLinkMass) lowLinkMass=spec->getLink(j).mass;
      if(spec->getLink(j).mass>highLinkMass) highLinkMass=spec->getLink(j).mass;
    }
  }

  //Create mutexes
  Threading::CreateMutex(&mutexKIonsManager);
  mutexSpecScore = new Mutex[spec->size()];
  for(j=0;j<spec->size();j++){
    Threading::CreateMutex(&mutexSpecScore[j]);
  }

  //xCorrCount=0;
}

KAnalysis::~KAnalysis(){
  int i;

  //Destroy mutexes
  Threading::DestroyMutex(mutexKIonsManager);
  for(i=0;i<spec->size();i++){
    Threading::DestroyMutex(mutexSpecScore[i]);
  }
  delete [] mutexSpecScore;
  
  //Deallocate memory and release pointers
  deallocateMemory();
  db=NULL;
  spec=NULL;

}

//============================
//  Public Functions
//============================
bool KAnalysis::doPeptideAnalysis(bool crossLink){
  unsigned int i;
  int iCount;
  int iPercent;
  int iTmp;
  vector<kPeptide>* p;
  vector<int> index;
  vector<kPepMod> mods;

  kScoreCard sc;

  ThreadPool<kAnalysisStruct*>* threadPool = new ThreadPool<kAnalysisStruct*>(analyzePeptideProc,params.threads,params.threads);

  //Set progress meter
  iPercent=0;
  printf("Progress: %2d%%",iPercent);
  fflush(stdout);

  //Set which list of peptides to search (with and without internal lysine)
  p=db->getPeptideList(crossLink);

  //Iterate entire peptide list
  for(i=0;i<p->size();i++){

    //Peptides are sorted by mass. If greater than max mass, stop checking peptides
    if(!crossLink && p->at(i).mass<minMass) continue;
    if(p->at(i).mass>maxMass) break;

    kAnalysisStruct* a = new kAnalysisStruct(&mutexKIonsManager,&p->at(i),i,crossLink);
    threadPool->Launch(a);

  }

  while(iPercent<99){
    //Update progress meter
    iCount=threadPool->NumParamsQueued();
    iTmp=100-(int)((double)iCount/p->size()*100);
    if(iTmp==100) break;
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }
    Threading::ThreadSleep(10);
  }

  while(threadPool->NumActiveThreads()>0){
    Threading::ThreadSleep(10);
  }

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;

  //clean up memory & release pointers
  delete threadPool;
  threadPool=NULL;
  p=NULL;

  return true;
}

//Non-covalent dimers need to be analyzed separately when doing a full search (not the relaxed mode).
//Relaxed mode has non-covalent dimerization built in.
/*
bool KAnalysis::doPeptideAnalysisNC(){
  unsigned int i;
  int k;
  int iPercent;
  int iTmp;
  int iCount;
  kPeptide pep;
  kPeptideB pepB;
  vector<kPeptideB> p;

  //Combine list of linkable and non-linkable peptides
  pepB.linkable=false;
  for(k=0;k<db->getPeptideListSize(false);k++){
    pepB.mass=db->getPeptide(k,false).mass;
    pepB.index=k;
    p.push_back(pepB);
  }
  pepB.linkable=true;
  for(k=0;k<db->getPeptideListSize(true);k++){
    pepB.mass=db->getPeptide(k,true).mass;
    pepB.index=k;
    p.push_back(pepB);
  }

  //sort list by mass
  qsort(&p[0],p.size(),sizeof(kPeptideB),comparePeptideBMass);

  ThreadPool<kAnalysisNCStruct*>* threadPool = new ThreadPool<kAnalysisNCStruct*>(analyzePeptideNCProc,params.threads,params.threads);

  //Set progress meter
  iPercent=0;
  printf("Progress: %2d%%",iPercent);
  fflush(stdout);

  //Iterate entire peptide list
  for(i=0;i<p.size();i++){

    //Peptides are sorted by mass. If greater than max mass, stop checking peptides
    if(p[i].mass>maxMass) break;

    kAnalysisNCStruct* a = new kAnalysisNCStruct(&mutexKIonsManager,&p,i);
    threadPool->Launch(a);

  }

  while(iPercent<99){
    //Update progress meter
    iCount=threadPool->NumParamsQueued();
    iTmp=100-(int)((double)iCount/p.size()*100);
    if(iTmp==100) break;
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }
    Threading::ThreadSleep(10);
  }

  while(threadPool->NumActiveThreads()>0){
    Threading::ThreadSleep(10);
  }

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;

  return true;
}
*/

bool KAnalysis::doRelaxedAnalysis(){
  int i;
  int iCount;
  int iPercent;
  int iTmp;

  ThreadPool<kAnalysisRelStruct*>* threadPool = new ThreadPool<kAnalysisRelStruct*>(analyzeRelaxedProc,params.threads,params.threads);

  //Set progress meter
  iPercent=0;
  printf("Progress: %2d%%",iPercent);
  fflush(stdout);

  for(i=0;i<spec->size();i++){
    kAnalysisRelStruct* a = new kAnalysisRelStruct(&spec->at(i));
    threadPool->Launch(a);
  }

  while(iPercent<99){
    //Update progress meter
    iCount=threadPool->NumParamsQueued();
    iTmp=100-(int)((double)iCount/spec->size()*100);
    if(iTmp==100) break;
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }
    Threading::ThreadSleep(10);
  }

  while(threadPool->NumActiveThreads()>0){
    Threading::ThreadSleep(10);
  }

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;

  return true;
}


//============================
//  Private Functions
//============================

//============================
//  Thread-Start Functions
//============================

//These functions fire off when a thread starts. They pass the variables to for
//each thread-specific analysis to the appropriate function.
void KAnalysis::analyzePeptideProc(kAnalysisStruct* s){
  int i;
  Threading::LockMutex(mutexKIonsManager);
  for(i=0;i<params.threads;i++){
    if(!bKIonsManager[i]){
      bKIonsManager[i]=true;
      break;
    }
  }
  Threading::UnlockMutex(mutexKIonsManager);
  if(i==params.threads){
    cout << "Error in KAnalysis::analyzePeptidesProc" << endl;
    exit(-1);
  }
  s->bKIonsMem = &bKIonsManager[i];
  analyzePeptide(s->pep,s->pepIndex,i,s->crossLink);
  delete s;
  s=NULL;
}

/*
void KAnalysis::analyzePeptideNCProc(kAnalysisNCStruct* s){
  int i;
  Threading::LockMutex(mutexKIonsManager);
  for(i=0;i<params.threads;i++){
    if(!bKIonsManager[i]){
      bKIonsManager[i]=true;
      break;
    }
  }
  Threading::UnlockMutex(mutexKIonsManager);
  if(i==params.threads){
    cout << "Error in KAnalysis::analyzePeptideNCProc" << endl;
    exit(-1);
  }
  s->bKIonsMem = &bKIonsManager[i];
  analyzePeptideNC(s->pep,s->pepIndex,i);
  delete s;
  s=NULL;
}
*/

void KAnalysis::analyzeRelaxedProc(kAnalysisRelStruct* s){
  analyzeRelaxed(s->spec);
  delete s;
  s=NULL;
}

//============================
//  Analysis Functions
//============================

//Analyzes all single peptides. Also analyzes cross-linked peptides when in full search mode, 
//or stage 1 of relaxed mode analysis
bool KAnalysis::analyzePeptide(kPeptide* p, int pepIndex, int iIndex, bool crossLink){

  int j;
  unsigned int k,k2;
  int n;
  double totalMass;
  bool bt;
  vector<int> index;
  vector<kPepMod> mods;
  
  //char str[256];
  //db->getPeptideSeq(p->map->at(0).index,p->map->at(0).start,p->map->at(0).stop,str);
  //cout << "Checking: " << str << endl;

  //Set the peptide, calc the ions, and score it against the spectra
  ions[iIndex].setPeptide(true,&db->at(p->map->at(0).index).sequence[p->map->at(0).start],p->map->at(0).stop-p->map->at(0).start+1,p->mass);
  
  ions[iIndex].buildIons();
  ions[iIndex].modIonsRec(0,-1,0,0,false);

  for(j=0;j<ions[iIndex].size();j++){

    bt=spec->getBoundaries2(ions[iIndex][j].mass,params.ppmPrecursor,index);
    if(bt) scoreSpectra(index,j,ions[iIndex][j].difMass,crossLink,pepIndex,-1,-1,-1,-1,iIndex);
    
    //For searching non-covalent dimers - Might remove altogether. Adds 100% more computation for less than 0.01% more IDs
    if(params.dimers>0) analyzeSingletsNoLysine(*p,j,pepIndex,crossLink,iIndex);
  }

  if(!crossLink) return true;

  //Crosslinked peptides must also search singlets with reciprocol mass on each lysine
  analyzeSinglets(*p,pepIndex,lowLinkMass,highLinkMass,iIndex);

  //also search loop-links
  //check loop-links by iterating through each cross-linker mass
  for(n=0;n<spec->sizeLink();n++){
    if(spec->getLink(n).mono) continue;
    
    totalMass = p->mass+spec->getLink(n).mass;
      
    //skip links that are too large
    if(totalMass>maxMass) continue;

    //if linker moieties are different, mix sites
    if(spec->getLink(n).siteA!=spec->getLink(n).siteB){
      for(k=0;k<p->vA->size();k++){
        for(k2=0;k2<p->vB->size();k2++){
          ions[iIndex].reset();
          ions[iIndex].buildLoopIons(spec->getLink(n).mass,p->vA->at(k)-p->map->at(0).start,p->vB->at(k2)-p->map->at(0).start);
          ions[iIndex].modLoopIonsRec(0,p->vA->at(k)-p->map->at(0).start,p->vB->at(k2)-p->map->at(0).start,0,0,true);
          for(j=0;j<ions[iIndex].size();j++){
            bt=spec->getBoundaries2(ions[iIndex][j].mass,params.ppmPrecursor,index);
            if(bt) scoreSpectra(index,j,0,crossLink,pepIndex,-1,p->vA->at(k)-p->map->at(0).start,p->vB->at(k2)-p->map->at(0).start,n,iIndex);
          }
        }
      }

    //otherwise link across matching sites within peptide
    } else {
      if(spec->getLink(n).siteA==1){
        if(p->vA->size()>0){
          for(k=0;k<p->vA->size()-1;k++){
            for(k2=k+1;k2<p->vA->size();k2++){  
              ions[iIndex].reset();  
              ions[iIndex].buildLoopIons(spec->getLink(n).mass,p->vA->at(k)-p->map->at(0).start,p->vA->at(k2)-p->map->at(0).start);
              ions[iIndex].modLoopIonsRec(0,p->vA->at(k)-p->map->at(0).start,p->vA->at(k2)-p->map->at(0).start,0,0,true);
              for(j=0;j<ions[iIndex].size();j++){
                bt=spec->getBoundaries2(ions[iIndex][j].mass,params.ppmPrecursor,index);
                if(bt) scoreSpectra(index,j,0,crossLink,pepIndex,-1,p->vA->at(k)-p->map->at(0).start,p->vA->at(k2)-p->map->at(0).start,n,iIndex);
              }
            }
          }
        }
      } else {
        if(p->vB->size()>0){
          for(k=0;k<p->vB->size()-1;k++){
            for(k2=k+1;k2<p->vB->size();k2++){  
              ions[iIndex].reset();
              ions[iIndex].buildLoopIons(spec->getLink(n).mass,p->vB->at(k)-p->map->at(0).start,p->vB->at(k2)-p->map->at(0).start);
              ions[iIndex].modLoopIonsRec(0,p->vB->at(k)-p->map->at(0).start,p->vB->at(k2)-p->map->at(0).start,0,0,true);
              for(j=0;j<ions[iIndex].size();j++){
                bt=spec->getBoundaries2(ions[iIndex][j].mass,params.ppmPrecursor,index);
                if(bt) scoreSpectra(index,j,0,crossLink,pepIndex,-1,p->vB->at(k)-p->map->at(0).start,p->vB->at(k2)-p->map->at(0).start,n,iIndex);
              }
            }
          }
        }
      }
    }

  }

  return true;
}

//Analyzes non-covalent dimers when in full search mode.
/*
void KAnalysis::analyzePeptideNC(vector<kPeptideB>* p, int pIndex, int iIndex){
  unsigned int j;
  kPeptide pep;
  double totalMass;
  vector<int> index;
  bool bt;

  //Set the peptide, calc the ions, and score it against the spectra
  pep=db->getPeptide(p->at(pIndex).index,p->at(pIndex).linkable);
  ions[iIndex].setPeptide(true,&db->at(pep.map->at(0).index).sequence[pep.map->at(0).start],pep.map->at(0).stop-pep.map->at(0).start+1,p->at(pIndex).mass);

  //check dimerization with another peptide (including itself)
  for(j=pIndex;j<p->size();j++){

    totalMass=p->at(pIndex).mass+p->at(j).mass;
    if(totalMass>maxMass) break;

    //Get a set of spectra to analyze. Spectra were sorted by precursor mass
    if(totalMass<minMass) bt=false;
    else bt=spec->getBoundaries2(totalMass,params.ppmPrecursor,index);

    if(bt){
      //set second peptide
      pep=db->getPeptide(p->at(j).index,p->at(j).linkable);
      ions[iIndex].setPeptide(false,&db->at(pep.map->at(0).index).sequence[pep.map->at(0).start],pep.map->at(0).stop-pep.map->at(0).start+1,pep.mass);
              
      //Calc the ions, and score it against the spectra   
      ions[iIndex].buildNCIons();
      scoreNCSpectra(index,totalMass,p->at(pIndex).linkable,p->at(j).linkable,p->at(pIndex).index,p->at(j).index, iIndex);

    }

  }
}
*/

//Stage 2 of relaxed mode analysis. Must be performed after analyzePeptides.
void KAnalysis::analyzeRelaxed(KSpectrum* sp){

  int i,j,k,m,n,x;
  unsigned int q;
  int index;
  int count=sp->getSingletCount();

  double ppm;
  double totalMass;

  kSingletScoreCardPlus* s=new kSingletScoreCardPlus[count];
  kSingletScoreCard sc1;
  kScoreCard sc;
  kPeptide pep;
  
  if(sp->getScanNumber()==params.diagnostic){
    FILE* f=fopen("diagnostic.txt","wt");
    fprintf(f,"Scan: %d\n",sp->getScanNumber());
    char strs[256];
    for(k=0;k<count;k++){
      sc1=sp->getSingletScoreCard(k);
      db->getPeptideSeq( db->getPeptideList(sc1.linkable)->at(sc1.pep1).map->at(0).index,db->getPeptideList(sc1.linkable)->at(sc1.pep1).map->at(0).start,db->getPeptideList(sc1.linkable)->at(sc1.pep1).map->at(0).stop,strs);
      for(q=0;q<strlen(strs);q++){
        fprintf(f,"%c",strs[q]);
        for(x=0;x<sc1.modLen;x++){
          if(sc1.mods[x].pos==q) fprintf(f,"[%.2lf]",sc1.mods[x].mass);
        }
        if(q==sc1.k1)fprintf(f,"[x]");
      }
      fprintf(f,"\t%d\t%d\t%.6lf\t%.4lf\t%.4lf\n",sc1.k1,(int)sc1.modLen,sc1.mass,sc1.simpleScore,sc1.simpleScore*sc1.len);
    }
    fclose(f);
  }

  //Make a sortable list of top hits to compare
  for(j=0;j<count;j++){
    sc1=sp->getSingletScoreCard(j);
    s[j].len=sc1.len;
    s[j].k1=sc1.k1;
    s[j].linkable=sc1.linkable;
    s[j].pep1=sc1.pep1;
    s[j].rank=j;
    s[j].simpleScore=sc1.simpleScore;
    if(sc1.simpleScore>0)s[j].mass=sc1.mass;
    else s[j].mass=0;
    
    //Determine if peptide belongs to target or decoy proteins, or both
    k=0; //target protein counts
    n=0; //decoy protein counts
    pep = db->getPeptide(sc1.pep1,sc1.linkable);
    for(i=0;i<(int)pep.map->size();i++) {
      if(db->at(pep.map->at(i).index).name.find(params.decoy)==string::npos) k++;
      else n++;
    }
    if(k>0 && n>0) s[j].target=2;
    else if(k>0) s[j].target=1;
    else s[j].target=0;
  }
  qsort(s,count,sizeof(kSingletScoreCardPlus),compareSSCPlus);

  int xlink[6];
  int peplink[6];
  char cp1;
  char cp2;
  int len1;
  int len2;
  bool bN1;
  bool bN2;
  bool bSkip;
  kPeptide p;

  //Check true cross-links
  for(k=0;k<spec->sizeLink();k++){
    if(spec->getLink(k).mono>0) continue;
    for(j=0;j<6;j++) xlink[j]=0;
    if(spec->getLink(k).siteA==1) xlink[params.setA]++;
    else xlink[params.setB]++;
    if(spec->getLink(k).siteB==1) xlink[params.setA]++;
    else xlink[params.setB]++;

    for(j=0;j<count;j++){
        
      if(s[j].simpleScore<=0) continue; //skips anything we've already considered
      if(!s[j].linkable) continue;
      if(s[j].k1<0) continue;

      p=db->getPeptide(s[j].pep1,true);
      cp1=db->at(p.map->at(0).index).sequence[p.map->at(0).start+s[j].k1];
      len1=db->at(p.map->at(0).index).sequence.size()-1;
      if( (p.map->at(0).start+s[j].k1)<2) bN1=true;
      else bN1=false;

      for(m=0;m<sp->sizePrecursor();m++){
        index=findMass(s,count,sp->getPrecursor(m).monoMass-s[j].mass-spec->getLink(k).mass);
        n=index;
        while(n<count){
          if(s[n].simpleScore<0 || s[n].k1<0){
            n++;
            continue;
          }

          totalMass=s[j].mass+s[n].mass+spec->getLink(k).mass;
          ppm = (totalMass-sp->getPrecursor(m).monoMass)/sp->getPrecursor(m).monoMass*1e6;
          if(fabs( ppm )<=params.ppmPrecursor) {

            //Check to make sure both peptides pair to the cross-linker
            p=db->getPeptide(s[n].pep1,true);
            cp2=db->at(p.map->at(0).index).sequence[p.map->at(0).start+s[n].k1];
            len2=db->at(p.map->at(0).index).sequence.size()-1;
            if( (p.map->at(0).start+s[n].k1)<2) bN2=true;
            else bN2=false;
            
            for(x=0;x<6;x++) peplink[x]=0;

            if(cp1=='K' || bN1) peplink[1]++;
            if(cp1=='D' || cp1=='E' || s[j].k1==len1) peplink[2]++;
            if(cp1=='C') peplink[3]++;
						if(cp1=='Q') peplink[4]++;
            if(cp1=='K' || cp1=='S' || cp1=='T' || cp1=='Y' || bN1) peplink[5]++;

            if(cp2=='K' || bN2) peplink[1]++;
            if(cp2=='D' || cp2=='E' || s[n].k1==len2) peplink[2]++;
            if(cp2=='C') peplink[3]++;
						if(cp2=='Q') peplink[4]++;
            if(cp2=='K' || cp2=='S' || cp2=='T' || cp2=='Y' || bN2) peplink[5]++;

            bSkip=false;
            for(x=0;x<6;x++){
              if(peplink[x]<xlink[x]) {
                bSkip=true;
                break;
              }
            }
            if(bSkip){
              n++;
              continue;
            }

            //Add the cross-link
            sc.simpleScore=s[j].simpleScore*s[j].len+s[n].simpleScore*s[n].len;
            sc.k1=s[j].k1;
            sc.k2=s[n].k1;
            sc.mass=totalMass;
            sc.linkable1=s[j].linkable;
            sc.linkable2=s[n].linkable;
            sc.pep1=s[j].pep1;
            sc.pep2=s[n].pep1;
            sc.link=k;
            sc.rank=s[j].rank+s[n].rank;
            if((s[j].simpleScore*s[j].len) > (s[n].simpleScore*s[n].len)) sc.scoreDiff=sc.simpleScore - s[j].simpleScore*s[j].len;
            else sc.scoreDiff = sc.simpleScore - s[n].simpleScore*s[n].len;
            sc.mods1->clear();
            sc.mods2->clear();
            sc1=sp->getSingletScoreCard(s[j].rank);
            for(i=0;i<sc1.modLen;i++) sc.mods1->push_back(sc1.mods[i]);
            sc1=sp->getSingletScoreCard(s[n].rank);
            for(i=0;i<sc1.modLen;i++) sc.mods2->push_back(sc1.mods[i]);
            sp->checkScore(sc);
          } else if(ppm>params.ppmPrecursor) {
            break;
          }
          n++;
        }
        n=index-1;
        while(n>-1){
          if(s[n].simpleScore<0 || s[n].k1<0){
            n--;
            continue;
          }

          totalMass=s[j].mass+s[n].mass+spec->getLink(k).mass;
          ppm = (totalMass-sp->getPrecursor(m).monoMass)/sp->getPrecursor(m).monoMass*1e6;
          if(fabs( ppm )<=params.ppmPrecursor) {

            //Check to make sure both peptides pair to the cross-linker
            p=db->getPeptide(s[n].pep1,true);
            cp2=db->at(p.map->at(0).index).sequence[p.map->at(0).start+s[n].k1];
            len2=db->at(p.map->at(0).index).sequence.size()-1;
            if( (p.map->at(0).start+s[n].k1)<2) bN2=true;
            else bN2=false;
            
            for(x=0;x<6;x++) peplink[x]=0;

            if(cp1=='K' || bN1) peplink[1]++;
            if(cp1=='D' || cp1=='E' || s[j].k1==len1) peplink[2]++;
            if(cp1=='C') peplink[3]++;
						if(cp1=='Q') peplink[4]++;
            if(cp1=='K' || cp1=='S' || cp1=='T' || cp1=='Y' || bN1) peplink[5]++;

            if(cp2=='K' || bN2) peplink[1]++;
            if(cp2=='D' || cp2=='E' || s[n].k1==len2) peplink[2]++;
            if(cp2=='C') peplink[3]++;
						if(cp2=='Q') peplink[4]++;
            if(cp2=='K' || cp2=='S' || cp2=='T' || cp2=='Y' || bN2) peplink[5]++;

            bSkip=false;
            for(x=0;x<6;x++){
              if(peplink[x]<xlink[x]) {
                bSkip=true;
                break;
              }
            }
            if(bSkip){
              n--;
              continue;
            }

            sc.simpleScore=s[j].simpleScore*s[j].len+s[n].simpleScore*s[n].len;
            sc.k1=s[j].k1;
            sc.k2=s[n].k1;
            sc.mass=totalMass;
            sc.linkable1=s[j].linkable;
            sc.linkable2=s[n].linkable;
            sc.pep1=s[j].pep1;
            sc.pep2=s[n].pep1;
            sc.link=k;
            sc.rank=s[j].rank+s[n].rank;
            if((s[j].simpleScore*s[j].len) > (s[n].simpleScore*s[n].len)) sc.scoreDiff=sc.simpleScore - s[j].simpleScore*s[j].len;
            else sc.scoreDiff = sc.simpleScore - s[n].simpleScore*s[n].len;
            sc.mods1->clear();
            sc.mods2->clear();
            sc1=sp->getSingletScoreCard(s[j].rank);
            for(i=0;i<sc1.modLen;i++) sc.mods1->push_back(sc1.mods[i]);
            sc1=sp->getSingletScoreCard(s[n].rank);
            for(i=0;i<sc1.modLen;i++) sc.mods2->push_back(sc1.mods[i]);
            sp->checkScore(sc);
          } else if(ppm<-params.ppmPrecursor) {
            break;
          }
          n--;
        }
      }
      s[j].simpleScore=-s[j].simpleScore;
    }
  
    //reset scores
    for(j=0;j<count;j++){
      if(s[j].simpleScore<0) s[j].simpleScore=-s[j].simpleScore;
    }

  }

  //reset scores
  for(j=0;j<count;j++){
    if(s[j].simpleScore<0) s[j].simpleScore=-s[j].simpleScore;
  }

  //Check non-covalent dimers
  //Note that this code is out of date and probably does not present scores or mods correctly.
  if(params.dimers==0) {
    delete [] s;
    return;
  }

  for(j=0;j<count;j++){
    if(s[j].simpleScore<=0) continue;
    if(s[j].k1>-1) continue;

    for(m=0;m<sp->sizePrecursor();m++){
      index=findMass(s,count,sp->getPrecursor(m).monoMass-s[j].mass);
      n=index;
      while(n<count){
        if(s[n].simpleScore<=0 || s[n].k1>-1){
          n++;
          continue;
        }
        totalMass=s[j].mass+s[n].mass;
        ppm = (totalMass-sp->getPrecursor(m).monoMass)/sp->getPrecursor(m).monoMass*1e6;
        if(fabs( ppm )<=params.ppmPrecursor) {
          sc.simpleScore=s[j].simpleScore*s[j].len+s[n].simpleScore*s[n].len;
          sc.k1=-1;
          sc.k2=-1;
          sc.mass=totalMass;
          sc.linkable1=s[j].linkable;
          sc.linkable2=s[n].linkable;
          sc.pep1=s[j].pep1;
          sc.pep2=s[n].pep1;
          sc.link=-2;
          sc.rank=s[j].rank+s[n].rank;

          if((s[j].simpleScore*s[j].len) > (s[n].simpleScore*s[n].len)) sc.scoreDiff=sc.simpleScore - s[j].simpleScore*s[j].len;
          else sc.scoreDiff = sc.simpleScore - s[n].simpleScore*s[n].len;
            
          sc.mods1->clear();
          sc.mods2->clear();
          sc1=sp->getSingletScoreCard(s[j].rank);
          for(i=0;i<sc1.modLen;i++) sc.mods1->push_back(sc1.mods[i]);
          sc1=sp->getSingletScoreCard(s[n].rank);
          for(i=0;i<sc1.modLen;i++) sc.mods2->push_back(sc1.mods[i]);

          sp->checkScore(sc);
        } else if(ppm>params.ppmPrecursor) {
          break;
        }
        n++;
      }
      n=index-1;
      while(n>-1){
        if(s[n].simpleScore<=0 || s[n].k1>-1){
          n--;
          continue;
        }
        totalMass=s[j].mass+s[n].mass;
        ppm = (totalMass-sp->getPrecursor(m).monoMass)/sp->getPrecursor(m).monoMass*1e6;
        if(fabs( ppm )<=params.ppmPrecursor) {
          sc.simpleScore=s[j].simpleScore*s[j].len+s[n].simpleScore*s[n].len;
          sc.k1=-1;
          sc.k2=-1;
          sc.mass=totalMass;
          sc.linkable1=s[j].linkable;
          sc.linkable2=s[n].linkable;
          sc.pep1=s[j].pep1;
          sc.pep2=s[n].pep1;
          sc.link=-2;
          sc.rank=s[j].rank+s[n].rank;
          
          if((s[j].simpleScore*s[j].len) > (s[n].simpleScore*s[n].len)) sc.scoreDiff=sc.simpleScore - s[j].simpleScore*s[j].len;
          else sc.scoreDiff = sc.simpleScore - s[n].simpleScore*s[n].len;
            
          sc.mods1->clear();
          sc.mods2->clear();
          sc1=sp->getSingletScoreCard(s[j].rank);
          for(i=0;i<sc1.modLen;i++) sc.mods1->push_back(sc1.mods[i]);
          sc1=sp->getSingletScoreCard(s[n].rank);
          for(i=0;i<sc1.modLen;i++) sc.mods2->push_back(sc1.mods[i]);

          sp->checkScore(sc);
        } else if(ppm<-params.ppmPrecursor) {
          break;
        }
        n--;
      }
    }
    s[j].simpleScore=-s[j].simpleScore;
  }

  delete [] s;
}

bool KAnalysis::analyzeSinglets(kPeptide& pep, int index, double lowLinkMass, double highLinkMass, int iIndex){
  int i;
  unsigned int j;
  unsigned int k;
  int len;
  double minMass;
  double maxMass;
  vector<int> scanIndex;

  //Set Mass boundaries
  minMass=pep.mass+lowLinkMass+params.minPepMass;
  maxMass=pep.mass+highLinkMass+params.maxPepMass;
  minMass-=(minMass/1000000*params.ppmPrecursor);
  maxMass+=(maxMass/1000000*params.ppmPrecursor);

  //Find mod mass as difference between precursor and peptide
  len=(pep.map->at(0).stop-pep.map->at(0).start)+1;
  ions[iIndex].setPeptide(true,&db->at(pep.map->at(0).index).sequence[pep.map->at(0).start],len,pep.mass);
  
  //Iterate every link site as a site for the monolink
  for(k=0;k<pep.vA->size();k++){
  
    //build fragment ions and score against all potential spectra
    ions[iIndex].reset();
    ions[iIndex].buildSingletIons(pep.vA->at(k)-pep.map->at(0).start);
    ions[iIndex].modIonsRec(0,pep.vA->at(k)-pep.map->at(0).start,0,0,true);

    //iterate through all ion sets
    for(i=0;i<ions[iIndex].size();i++){

      //Iterate all spectra from (peptide mass + low linker + minimum mass) to (peptide mass + high linker + maximum mass)
      if(!spec->getBoundaries(minMass+ions[iIndex][i].difMass,maxMass+ions[iIndex][i].difMass,scanIndex)) continue;
      for(j=0;j<scanIndex.size();j++){
        scoreSingletSpectra(scanIndex[j],i,ions[iIndex][i].mass,len,index,(char)(pep.vA->at(k)-pep.map->at(0).start),true,minMass,iIndex);
      }
    }

  }

  //Iterate every alternate link site as a site for the monolink
  for(k=0;k<pep.vB->size();k++){

    //build fragment ions and score against all potential spectra
    ions[iIndex].reset();
    ions[iIndex].buildSingletIons(pep.vB->at(k)-pep.map->at(0).start);
    ions[iIndex].modIonsRec(0,pep.vB->at(k)-pep.map->at(0).start,0,0,true);

    //iterate through all ion sets
    for(i=0;i<ions[iIndex].size();i++){

      //Iterate all spectra from (peptide mass + low linker + minimum mass) to (peptide mass + high linker + maximum mass)
      if(!spec->getBoundaries(minMass+ions[iIndex][i].difMass,maxMass+ions[iIndex][i].difMass,scanIndex)) continue;
      for(j=0;j<scanIndex.size();j++){
        scoreSingletSpectra(scanIndex[j],i,ions[iIndex][i].mass,len,index,pep.vB->at(k)-pep.map->at(0).start,true,minMass,iIndex);
      }
    }

  }

  return true;
}

bool KAnalysis::analyzeSingletsNoLysine(kPeptide& pep, int sIndex, int index, bool linkable, int iIndex){
  unsigned int j;
  double minMass;
  double maxMass;
  vector<int> scanIndex;

  //Set Mass boundaries
  minMass=pep.mass+ions[iIndex][sIndex].difMass+params.minPepMass;
  maxMass=pep.mass+ions[iIndex][sIndex].difMass+params.maxPepMass;
  minMass-=(minMass/1000000*params.ppmPrecursor);
  maxMass+=(maxMass/1000000*params.ppmPrecursor);

  //Iterate all spectra from (peptide mass + minimum mass) to (peptide mass + maximum mass)
  cout << "stray getBoundaries" << endl;
  if(!spec->getBoundaries(minMass,maxMass,scanIndex)) return false;
  for(j=0;j<scanIndex.size();j++){
    scoreSingletSpectra(scanIndex[j],sIndex,ions[iIndex][sIndex].mass,pep.map->at(0).stop-pep.map->at(0).start+1,index,-1,linkable,minMass,iIndex);
  }
  return true;

}


/*============================
  Private Functions
============================*/
bool KAnalysis::allocateMemory(int threads){
  bKIonsManager = new bool[threads];
  ions = new KIons[threads];
  for(int i=0;i<threads;i++) {
    bKIonsManager[i]=false;
    ions[i].setModFlags(params.monoLinksOnXL,params.diffModsOnXL);
  }
  return true;
}

void KAnalysis::deallocateMemory(){
  delete [] bKIonsManager;
  delete [] ions;
}

int KAnalysis::findMass(kSingletScoreCardPlus* s, int sz, double mass){
  int lower=0;
  int mid=sz/2;
  int upper=sz;

  //binary search to closest mass
  while(s[mid].mass!=mass){
		if(lower>=upper) break;
    if(mass<s[mid].mass){
			upper=mid-1;
			mid=(lower+upper)/2;
		} else {
			lower=mid+1;
			mid=(lower+upper)/2;
		}
		if(mid==sz) {
			mid--;
			break;
		}
	}
  return mid;
}

/*
void KAnalysis::scoreNCSpectra(vector<int>& index, double mass, bool linkable1, bool linkable2, int pep1, int pep2, int iIndex){
  unsigned int a;
  kScoreCard sc;
 
  //score spectra
  for(a=0;a<index.size();a++){
    if(params.xcorr==1) sc.simpleScore=xCorrScoring(spec->at(index[a]),iIndex);
    else sc.simpleScore=simpleScoring(spec->at(index[a]),iIndex);
    sc.k1=-1;
    sc.k2=-1;
    sc.mass=mass;
    sc.linkable1=linkable1;
    sc.linkable2=linkable2;
    sc.pep1=pep1;
    sc.pep2=pep2;
    sc.link=-2;
    Threading::LockMutex(mutexSpecScore);
    spec->at(index[a]).checkScore(sc);
    Threading::UnlockMutex(mutexSpecScore);
  }
}
*/

void KAnalysis::scoreSingletSpectra(int index, int sIndex, double mass, int len, int pep, char k, bool linkable, double minMass, int iIndex){
  kSingletScoreCard sc;
  KIonSet* iset;
  kPepMod mod;
  double score1=0;
  double score2=0;
  int i;
  vector<kPepMod> v;

  KSpectrum* s=spec->getSpectrum(index);
  kPrecursor* p;
  int sz=s->sizePrecursor();
  for(i=0;i<sz;i++){
    p=s->getPrecursor2(i);
    if(p->monoMass>minMass){
      if(params.xcorr) score1=xCorrScoring(*s,p->monoMass-mass,sIndex,iIndex);
      else score1=kojakScoring(index,p->monoMass-mass,sIndex,iIndex);
      if(score1>score2) score2=score1;
    }
  }

  //Check the highest score
  sc.len=len;
  sc.simpleScore=(float)score2/len;
  sc.k1=k;
  sc.linkable=linkable;
  sc.pep1=pep;
  sc.mass=mass;
  if(sc.simpleScore>0) {
    v.clear();
    if(sc.mods!=NULL) {
      sc.modLen=0;
      delete [] sc.mods;
      sc.mods=NULL;
    }
    iset=ions[iIndex].at(sIndex);
    if(iset->difMass!=0){
      for(i=0;i<ions[iIndex].getIonCount();i++) {
        if(iset->mods[i]!=0){
          mod.pos=(char)i;
          mod.mass=iset->mods[i];
          v.push_back(mod);
        }
      }
      sc.modLen=(char)v.size();
      sc.mods=new kPepMod[sc.modLen];
      for(i=0;i<(int)sc.modLen;i++) sc.mods[i]=v[i];
    }

    Threading::LockMutex(mutexSpecScore[index]);
    spec->at(index).checkSingletScore(sc);
    Threading::UnlockMutex(mutexSpecScore[index]);

  }

}

void KAnalysis::scoreSpectra(vector<int>& index, int sIndex, double modMass, bool linkable, int pep1, int pep2, int k1, int k2, int link, int iIndex){
  unsigned int a;
  int i;
  kScoreCard sc;
  kPepMod mod;

  //score spectra
  for(a=0;a<index.size();a++){
    sc.mods1->clear();
    sc.mods2->clear();
    if(params.xcorr) sc.simpleScore=xCorrScoring(spec->at(index[a]),modMass,sIndex,iIndex);
    else sc.simpleScore=kojakScoring(index[a],modMass,sIndex,iIndex);
    sc.k1=k1;
    sc.k2=k2;
    sc.mass=ions[iIndex][sIndex].mass;
    sc.linkable1=sc.linkable2=linkable;
    sc.pep1=pep1;
    sc.pep2=pep2;
    sc.link=link;
    if(ions[iIndex][sIndex].difMass!=0){
      for(i=0;i<ions[iIndex].getPeptideLen();i++) {
        if(ions[iIndex][sIndex].mods[i]!=0){
          mod.pos=(char)i;
          mod.mass=ions[iIndex][sIndex].mods[i];
          sc.mods1->push_back(mod);
        }
      }
    }
    Threading::LockMutex(mutexSpecScore[index[a]]);
    spec->at(index[a]).checkScore(sc);
    Threading::UnlockMutex(mutexSpecScore[index[a]]);
  }

}

//An alternative score uses the XCorr metric from the Comet algorithm
//This version allows for fast scoring when the cross-linked mass is added.
float KAnalysis::xCorrScoring(KSpectrum& s, double modMass, int sIndex, int iIndex) { 

  //xCorrCount++;

  double dXcorr;
  double invBinSize=1.0/params.binSize;
  double binOffset=params.binOffset;
  double dif;

  int ionCount=ions[iIndex].getIonCount();
  int k;
  int maxCharge;
  int xx;

  //Grabbing a pointer directly to the ion set is faster than going 
  //through ions[iIndex][sIndex] for all fragment ions.
  KIonSet* ki=ions[iIndex].at(sIndex);

  int i;
  int j;

  unsigned int SpecSize=s.size();

  dXcorr=0.0;
  
  //Get the number of fragment ion series to analyze
  //The number is PrecursorCharge-1
  maxCharge=s.getCharge();
  if(maxCharge>6) maxCharge=6;

  //Iterate all series
  for(k=1;k<maxCharge;k++){

    dif=modMass/k;

    //Ratchet through pfFastXcorrData
    xx=0;
    for(i=0;i<ionCount;i++){
      if(ki->bIons[k][i]<0) j = (int)((dif-ki->bIons[k][i])*invBinSize+binOffset);
      else j = (int)(ki->bIons[k][i]*invBinSize+binOffset);
      while( j >=  s.xCorrSparseArray[xx].bin) xx++;
      dXcorr += s.xCorrSparseArray[xx-1].fIntensity;
    }

    xx=0;
    for(i=0;i<ionCount;i++){
      if(ki->yIons[k][i]<0) j = (int)((dif-ki->yIons[k][i])*invBinSize+binOffset);
      else j = (int)(ki->yIons[k][i]*invBinSize+binOffset);
      while( j >=  s.xCorrSparseArray[xx].bin) xx++;
      dXcorr += s.xCorrSparseArray[xx-1].fIntensity;
    }
  }

  //Scale score appropriately
  if(dXcorr <= 0.0) dXcorr=0.0;
  else dXcorr *= 0.005;

  return float(dXcorr);
}

//An alternative score uses the XCorr metric from the Comet algorithm
//This version allows for fast scoring when the cross-linked mass is added.
float KAnalysis::kojakScoring(int specIndex, double modMass, int sIndex, int iIndex) { 

  KSpectrum* s=spec->getSpectrum(specIndex);
  KIonSet* ki=ions[iIndex].at(sIndex);

  double dXcorr=0.0;
  double invBinSize=s->getInvBinSize();
  double binOffset=params.binOffset;
  double dif;
  double mz;

  int ionCount=ions[iIndex].getIonCount();
  int k;
  int maxCharge=s->getCharge();  

  int i;
  int key;
  int pos;

  unsigned int SpecSize=s->size();

  //The number of fragment ion series to analyze is PrecursorCharge-1
  //However, don't analyze past the 5+ series
  if(maxCharge>6) maxCharge=6;

  //Iterate all series
  for(k=1;k<maxCharge;k++){

    dif=modMass/k;

    //Iterate through pfFastXcorrData
    for(i=0;i<ionCount;i++){

      //get key
      if(ki->bIons[k][i]<0) mz = params.binSize * (int)((dif-ki->bIons[k][i])*invBinSize+binOffset);
      else mz = params.binSize * (int)(ki->bIons[k][i]*invBinSize+binOffset);
      key = (int)mz;
      if(key>=s->kojakBins) break;
      if(s->kojakSparseArray[key]==NULL) continue;
      pos = (int)((mz-key)*invBinSize);
      dXcorr += s->kojakSparseArray[key][pos];
    }

    for(i=0;i<ionCount;i++){
      if(ki->yIons[k][i]<0) mz = params.binSize * (int)((dif-ki->yIons[k][i])*invBinSize+binOffset);
      else mz = params.binSize * (int)(ki->yIons[k][i]*invBinSize+binOffset);
      key = (int)mz;
      if(key>=s->kojakBins) break;
      if(s->kojakSparseArray[key]==NULL) continue;
      pos = (int)((mz-key)*invBinSize);
      dXcorr += s->kojakSparseArray[key][pos];
    }
  }

  //Scale score appropriately
  if(dXcorr <= 0.0) dXcorr=0.0;
  else dXcorr *= 0.005;

  return float(dXcorr);
}


/*============================
  Utilities
============================*/
int KAnalysis::compareD(const void *p1, const void *p2){
  const double d1 = *(double *)p1;
  const double d2 = *(double *)p2;
  if(d1<d2) return -1;
  else if(d1>d2) return 1;
  else return 0;
}

int KAnalysis::comparePeptideBMass(const void *p1, const void *p2){
  const kPeptideB d1 = *(kPeptideB *)p1;
  const kPeptideB d2 = *(kPeptideB *)p2;
  if(d1.mass<d2.mass) return -1;
  else if(d1.mass>d2.mass) return 1;
  else return 0;
}

int KAnalysis::compareSSCPlus(const void *p1, const void *p2){
  const kSingletScoreCardPlus d1 = *(kSingletScoreCardPlus *)p1;
  const kSingletScoreCardPlus d2 = *(kSingletScoreCardPlus *)p2;
  if(d1.mass<d2.mass) return -1;      
  else if(d1.mass>d2.mass) return 1;
  else return 0;
}
