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
char**      KAnalysis::xlTable;
bool**      KAnalysis::scanBuffer;

int         KAnalysis::numIonSeries;

int*        KAnalysis::pepMassSize;
double**    KAnalysis::pepMass;
bool**      KAnalysis::pepBin;
int*        KAnalysis::pepBinSize;

int KAnalysis::skipCount;
int KAnalysis::nonSkipCount;

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
  xlTable = spec->getXLTable();

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

  numIonSeries=0;
  for(i=0;i<6;i++){
    if(params.ionSeries[i]) numIonSeries++;
  }

  //Create mutexes
  Threading::CreateMutex(&mutexKIonsManager);
  mutexSpecScore = new Mutex[spec->size()];
  for(j=0;j<spec->size();j++){
    Threading::CreateMutex(&mutexSpecScore[j]);
  }

  if (params.turbo) {
    cout << "  Turbo mode in use." << endl;
    makePepLists();
    skipCount=0;
    nonSkipCount=0;
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

  if (params.turbo){
    for (i = 0; i < spec->getMotifCount(); i++){
      delete[] pepMass[i];
      delete[] pepBin[i];
    }
    delete[] pepMass;
    delete[] pepMassSize;
    delete[] pepBin;
    delete[] pepBinSize;
  }

  //Deallocate memory and release pointers
  deallocateMemory(params.threads);
  db=NULL;
  spec=NULL;
  xlTable=NULL;
  
}

//============================
//  Public Functions
//============================
bool KAnalysis::doPeptideAnalysis(){
  unsigned int i;
  int iPercent;
  int iTmp;
  vector<kPeptide>* p;
  vector<int> index;
  vector<kPepMod> mods;

  kScoreCard sc;

  ThreadPool<kAnalysisStruct*>* threadPool = new ThreadPool<kAnalysisStruct*>(analyzePeptideProc,params.threads,params.threads,1);

  //Set progress meter
  iPercent=0;
  printf("%2d%%",iPercent);
  fflush(stdout);

  //Set which list of peptides to search (with and without internal lysine)
  p=db->getPeptideList();

  //Iterate entire peptide list
  for(i=0;i<p->size();i++){

    threadPool->WaitForQueuedParams();

    //Peptides are sorted by mass. If greater than max mass, stop checking peptides
    if(p->at(i).xlSites==0 && p->at(i).mass<minMass) continue;
    if(p->at(i).mass>maxMass) break;

    kAnalysisStruct* a = new kAnalysisStruct(&mutexKIonsManager,&p->at(i),i);
    threadPool->Launch(a);

    //Update progress meter
    iTmp=(int)((double)i/p->size()*100);
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }
  }

  threadPool->WaitForQueuedParams();
  threadPool->WaitForThreads();

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;

  //clean up memory & release pointers
  delete threadPool;
  threadPool=NULL;
  p=NULL;

  return true;
}

bool KAnalysis::doRelaxedAnalysis(){
  int i;
  int iPercent;
  int iTmp;

  ThreadPool<kAnalysisRelStruct*>* threadPool = new ThreadPool<kAnalysisRelStruct*>(analyzeRelaxedProc,params.threads,params.threads,1);

  //Set progress meter
  iPercent=0;
  printf("%2d%%",iPercent);
  fflush(stdout);

  for(i=0;i<spec->size();i++){

    threadPool->WaitForQueuedParams();

    kAnalysisRelStruct* a = new kAnalysisRelStruct(&mutexKIonsManager,&spec->at(i));
    threadPool->Launch(a);

    //Update progress meter
    iTmp=(int)((double)i/spec->size()*100);
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }

  }

  threadPool->WaitForQueuedParams();
  threadPool->WaitForThreads();

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;

  //clean up memory & release pointers
  delete threadPool;
  threadPool = NULL;

  /*
  skipCount=0;
  nonSkipCount=0;
  for (i = 0; i < spec->size(); i++){
    skipCount += spec->at(i).sc;
    nonSkipCount += spec->at(i).cc;
  }
  cout << "SkipCount = " << skipCount << endl;
  cout << "CorrCount = " << nonSkipCount << endl;
  */

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
  analyzePeptide(s->pep,s->pepIndex,i);
  delete s;
  s=NULL;
}

void KAnalysis::analyzeRelaxedProc(kAnalysisRelStruct* s){
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
    cout << "Error in KAnalysis::analyzeRelaxedProc" << endl;
    exit(-1);
  }
  s->bKIonsMem = &bKIonsManager[i];
  analyzeRelaxed(s->spec,i);
  delete s;
  s=NULL;
}

//============================
//  Analysis Functions
//============================

//Analyzes all single peptides. Also analyzes cross-linked peptides when in full search mode, 
//or stage 1 of relaxed mode analysis
bool KAnalysis::analyzePeptide(kPeptide* p, int pepIndex, int iIndex){

  int j;
  size_t k,k2,k3;
  bool bt;
  vector<int> index;
  vector<kPepMod> mods;

  //char str[256];
  //db->getPeptideSeq(p->map->at(0).index,p->map->at(0).start,p->map->at(0).stop,str);
  //cout << str << endl;
  //Set the peptide, calc the ions, and score it against the spectra
  ions[iIndex].setPeptide(true,&db->at(p->map->at(0).index).sequence[p->map->at(0).start],p->map->at(0).stop-p->map->at(0).start+1,p->mass,p->nTerm,p->cTerm,p->n15);
  
  ions[iIndex].buildIons();
  ions[iIndex].modIonsRec2(0,-1,0,0,false);

  for(j=0;j<ions[iIndex].size();j++){

    bt=spec->getBoundaries2(ions[iIndex][j].mass,params.ppmPrecursor,index,scanBuffer[iIndex]);
    if(bt) scoreSpectra(index,j,ions[iIndex][j].difMass,pepIndex,-1,-1,-1,-1,iIndex);
    
  }

  if(p->xlSites==0) return true;

  //Crosslinked peptides must also search singlets with reciprocol mass on each lysine
  analyzeSinglets(*p,pepIndex,lowLinkMass,highLinkMass,iIndex);

  if(p->xlSites==1) return true;

  //also search loop-links
  //check loop-links by iterating through each cross-linker mass
  string pepSeq;
  vector<int> xlIndex;
  int x;
  db->getPeptideSeq(*p,pepSeq);

  //iterate over every amino acid (except last - it has nothing to link)
  for (k = 0; k < pepSeq.size()-1; k++){

    //check if aa can be linked
    if (xlTable[pepSeq[k]][0]>-1) {

      //check all possible motifs at the site
      for (x=0;x<20;x++){
        if (xlTable[pepSeq[k]][x]==-1) break;

        //site can be linked, so check remaining sites
        for (k2 = k + 1; k2 < pepSeq.size(); k2++){
          if (k2 == pepSeq.size() - 1){ //handle c-terminus differently
            if (p->cTerm) checkXLMotif(xlTable[pepSeq[k]][x], xlTable['c'],xlIndex);
            else continue;
          } else if (xlTable[pepSeq[k2]][0] == -1) {
            continue;
          } else {
            checkXLMotif((int)xlTable[pepSeq[k]][x], xlTable[pepSeq[k2]],xlIndex);
          }

          if (xlIndex.size()>0){
            for(k3=0;k3<xlIndex.size();k3++){
              if(xlIndex[k3]<0) continue;
              ions[iIndex].reset();
              ions[iIndex].buildLoopIons(spec->getLink(xlIndex[k3]).mass, (int)k, (int)k2);
              ions[iIndex].modLoopIonsRec2(0, (int)k, (int)k2, 0, 0, true);
              for (j = 0; j<ions[iIndex].size(); j++){
                bt = spec->getBoundaries2(ions[iIndex][j].mass, params.ppmPrecursor, index, scanBuffer[iIndex]);
                if (bt) scoreSpectra(index, j, 0, pepIndex, -1, (int)k, (int)k2, xlIndex[k3], iIndex);
              }
            } //k3
          }
        } //k2

      }//x
    }

    //also check the n-terminus
    if (k == 0 && p->nTerm && xlTable['n'][0]>-1) {

      //check all possible motifs at the site
      for (x = 0; x<20; x++){
        if (xlTable['n'][x] == -1) break;

        //site can be linked, so check remaining sites
        for (k2 = k + 1; k2 < pepSeq.size(); k2++){
          if (k2 == pepSeq.size() - 1){ //handle c-terminus differently
            if (p->cTerm) checkXLMotif(xlTable['n'][x], xlTable['c'],xlIndex);
            else continue;
          } else if (xlTable[pepSeq[k2]][0] == -1) {
            continue;
          } else {
            checkXLMotif((int)xlTable['n'][x], xlTable[pepSeq[k2]],xlIndex);
          }

          if (xlIndex.size()>0){
            for (k3 = 0; k3<xlIndex.size(); k3++){
              if(xlIndex[k3]<0) continue;
              ions[iIndex].reset();
              ions[iIndex].buildLoopIons(spec->getLink(xlIndex[k3]).mass, (int)k, (int)k2);
              ions[iIndex].modLoopIonsRec2(0, (int)k, (int)k2, 0, 0, true);
              for (j = 0; j<ions[iIndex].size(); j++){
                bt = spec->getBoundaries2(ions[iIndex][j].mass, params.ppmPrecursor, index, scanBuffer[iIndex]);
                if (bt) scoreSpectra(index, j, 0, pepIndex, -1, (int)k, (int)k2, xlIndex[k3], iIndex);
              }
            } //k3
          }
        } //k2

      }//x
    }

  }//k
  return true;
}

//Stage 2 of relaxed mode analysis. Must be performed after analyzePeptides.
void KAnalysis::analyzeRelaxed(KSpectrum* sp, int iIndex){
  
  int i,j,k,m,n,x;
  int motifCount;
  unsigned int q,d;
  int index;
  int count=sp->getSingletCount();
  string pepSeq;
  char aa;

  double ppm;
  double totalMass;

  kSingletScoreCardPlus* s=new kSingletScoreCardPlus[count];
  kSingletScoreCard sc1;
  kScoreCard sc;
  kPeptide pep;
  
  //<------ Diagnostic output
  for(d=0;d<params.diag->size();d++){
    if(sp->getScanNumber()==params.diag->at(d)){
      char diagStr[256];
      sprintf(diagStr,"diagnostic_%d.txt",params.diag->at(d));
      FILE* f=fopen(diagStr,"wt");
      fprintf(f,"Scan: %d\n",sp->getScanNumber());
      char strs[256];
      for(k=0;k<count;k++){
        sc1=sp->getSingletScoreCard(k);
        db->getPeptideSeq( db->getPeptideList()->at(sc1.pep1).map->at(0).index,db->getPeptideList()->at(sc1.pep1).map->at(0).start,db->getPeptideList()->at(sc1.pep1).map->at(0).stop,strs);
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
      break;
    }
  }
  //<------ End diagnostic

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
    
    pep = db->getPeptide(sc1.pep1);
    db->getPeptideSeq(pep,pepSeq);
    aa=db->at(pep.map->at(0).index).sequence[pep.map->at(0).start + sc1.k1];

    //map the motifs
    motifCount=0;
    s[j].motif[motifCount]=-1;
    for (i = 0; i<20;i++){
      if(xlTable[aa][i]>-1) {
        s[j].motif[motifCount++] = xlTable[aa][i];
        s[j].motif[motifCount]=-1;
      } else {
        break;
      }
    }
    //special case to map n-term or c-term
    for (i = 0; i<pep.map->size();i++){
      //add n-term
      if ((pep.map->at(i).start + sc1.k1) < 2) {
        for (k = 0; k<20; k++){
          if (xlTable['n'][k]==-1) break;
          for (m=0;m<20;m++){
            if (s[j].motif[m]==-1) break;
            if (s[j].motif[m] == xlTable['n'][k]) break;
          }
          if (m == 20) {
            cout << "Max motifs reached in KAnalysis. Please report error." << endl;
            exit(1);
          } else if (m==motifCount){
            s[j].motif[motifCount++] = xlTable['n'][k];
            s[j].motif[motifCount]=-1;
          }
        }
      }
      //add c-term
      n = (int)db->at(pep.map->at(i).index).sequence.size() - 1;
      if (n==sc1.k1) {
        for (k = 0; k<20; k++){
          if (xlTable['c'][k] == -1)break;
          for (m = 0; m<20; m++){
            if (s[j].motif[m] == -1) break;
            if (s[j].motif[m] == xlTable['c'][k]) break;
          }
          if (m == 20) {
            cout << "Max motifs reached in KAnalysis. Please report error." << endl;
            exit(1);
          } else if (m == motifCount){
            s[j].motif[motifCount++] = xlTable['c'][k];
            s[j].motif[motifCount] = -1;
          }
        }
      }
    }

    //Determine if peptide belongs to target or decoy proteins, or both
    k = 0; //target protein counts
    n = 0; //decoy protein counts
    for(i=0;i<(int)pep.map->size();i++) {
      if(db->at(pep.map->at(i).index).name.find(params.decoy)==string::npos) k++;
      else n++;
    }
    if(k>0 && n>0) s[j].target=2;
    else if(k>0) s[j].target=1;
    else s[j].target=0;
  }
  qsort(s,count,sizeof(kSingletScoreCardPlus),compareSSCPlus);

  int cm;
  int counterMotif;
  int xlIndex;
  char cp1;
  char cp2;
  int len1;
  int len2;
  kPeptide p;
  vector<int> matches;
  kMatchSet msTemplate;
  kMatchSet msPartner;
  double dShared;

  //Iterate through sorted list of peptides
  for(j=0;j<count;j++){
        
    if(s[j].simpleScore<=0) continue; //skips anything we've already considered
    if(s[j].k1<0) continue;

    matches.clear();

    //Some basic peptide stats
    p=db->getPeptide(s[j].pep1);
    cp1=db->at(p.map->at(0).index).sequence[p.map->at(0).start+s[j].k1];
    len1=(int)db->at(p.map->at(0).index).sequence.size()-1;

    //iterate through all motifs associated with this peptide
    for (i=0;i<20;i++){
      if (s[j].motif[i]==-1) break;

      //iterate through all countermotifs
      for (cm = 0; cm < 10; cm++){
        counterMotif = spec->getCounterMotif((int)s[j].motif[i],cm);
        if (counterMotif<0) break;

        //Get XL index
        xlIndex = spec->getXLIndex((int)s[j].motif[i],cm);
        if (xlIndex < 0) {
          cout << "Error in xlIndex of analyzedRelaxed() function. Please report." << endl;
          exit(-1);
        }

        //iterate over all precursor masses
        for(m=0;m<sp->sizePrecursor();m++){

          //Grab bin coordinates of all fragment ions
          p=db->getPeptide(s[j].pep1);
          ions[iIndex].setPeptide(true,&db->at(p.map->at(0).index).sequence[p.map->at(0).start],p.map->at(0).stop-p.map->at(0).start+1,p.mass,p.nTerm,p.cTerm,p.n15);
          ions[iIndex].buildSingletIons((int)s[j].k1);
          setBinList(&msTemplate, iIndex, sp->getPrecursor(m).charge, sp->getPrecursor(m).monoMass-s[j].mass, sp->getSingletScoreCard(s[j].rank).mods, sp->getSingletScoreCard(s[j].rank).modLen);

          //get index of paired peptide subset by complement mass
          index=findMass(s,count,sp->getPrecursor(m).monoMass-s[j].mass-spec->getLink(xlIndex).mass);
          n=index;

          //iterate over all possible peptide pairs
          while(n<count){
            if(!params.dimersXL && n==j){ //skip self if not allowed to dimerize
              n++;
              continue;
            }
            if(s[n].simpleScore<0 || s[n].k1<0){ //skip peptides with no score contribution. Not sure if this is necessary anymore
              n++;
              continue;
            }

            //Make sure we didn't already pair these peptides with a different precursor during this pass
            for(x=0;x<(int)matches.size();x++){
              if(matches[x]==n) break;
            }
            if(x<(int)matches.size()) {
              n++;
              continue;
            }

            //compute theoretical mass
            totalMass=s[j].mass+s[n].mass+spec->getLink(xlIndex).mass;
            ppm = (totalMass-sp->getPrecursor(m).monoMass)/sp->getPrecursor(m).monoMass*1e6;

            //if wrong mass, move on to the next peptide or exit loop if next mass is wrong
            if (ppm<-params.ppmPrecursor) {
              n++;
              continue;
            }
            if (ppm>params.ppmPrecursor) break;
          
            //Make sure peptide has valid XL motif
            for (x = 0; x < 20; x++){
              if (s[n].motif[x]==-1) break;
              if (s[n].motif[x]==counterMotif) break;
            }
            if (x == 20) {
              n++;
              continue;
            } else if (s[n].motif[x] == -1){
              n++;
              continue;
            }

            p=db->getPeptide(s[n].pep1);
            cp2=db->at(p.map->at(0).index).sequence[p.map->at(0).start+s[n].k1];
            len2=(int)db->at(p.map->at(0).index).sequence.size()-1;
            
            //Grab bin coordinates of all fragment ions
            ions[iIndex].setPeptide(true, &db->at(p.map->at(0).index).sequence[p.map->at(0).start], p.map->at(0).stop - p.map->at(0).start + 1, p.mass, p.nTerm, p.cTerm,p.n15);
            ions[iIndex].buildSingletIons((int)s[n].k1);
            setBinList(&msPartner, iIndex, sp->getPrecursor(m).charge, sp->getPrecursor(m).monoMass-s[n].mass, sp->getSingletScoreCard(s[n].rank).mods, sp->getSingletScoreCard(s[n].rank).modLen);
            dShared = sharedScore(sp,&msTemplate,&msPartner,sp->getPrecursor(m).charge);

            //Add the cross-link
            sc.simpleScore=(float)(s[j].simpleScore*s[j].len+s[n].simpleScore*s[n].len-dShared);
            sc.k1=s[j].k1;
            sc.k2=s[n].k1;
            sc.mass=totalMass;
            sc.linkable1=s[j].linkable;
            sc.linkable2=s[n].linkable;
            sc.pep1=s[j].pep1;
            sc.pep2=s[n].pep1;
            sc.link=xlIndex;
            sc.rank1=s[j].rank;
            sc.rank2=s[n].rank;
            sc.score1=s[j].simpleScore*s[j].len;
            sc.score2=s[n].simpleScore*s[n].len;
            sc.mass1=s[j].mass;
            sc.mass2=s[n].mass;
            sc.mods1->clear();
            sc.mods2->clear();
            sc1=sp->getSingletScoreCard(s[j].rank);
            for(x=0;x<sc1.modLen;x++) sc.mods1->push_back(sc1.mods[x]);
            sc1=sp->getSingletScoreCard(s[n].rank);
            for(x=0;x<sc1.modLen;x++) sc.mods2->push_back(sc1.mods[x]);
            sp->checkScore(sc);
            matches.push_back(n);

            n++;
          }
      
          //repeat process walking through peptides in other direction
          n=index-1;
          while(n>-1){
            if(!params.dimersXL && n==j){
              n--;
              continue;
            }
            if(s[n].simpleScore<0 || s[n].k1<0){
              n--;
              continue;
            }

            //Make sure we didn't already pair these peptides with a different precursor during this pass
            for(x=0;x<matches.size();x++){
              if(matches[x]==n) break;
            }
            if(x<matches.size()) {
              n--;
              continue;
            }

            //compute theoretical mass
            totalMass = s[j].mass + s[n].mass + spec->getLink(xlIndex).mass;
            ppm = (totalMass - sp->getPrecursor(m).monoMass) / sp->getPrecursor(m).monoMass*1e6;

            //if wrong mass, move on to the next peptide or exit loop if next mass is wrong
            if (ppm>params.ppmPrecursor) {
              n--;
              continue;
            }
            if (ppm<-params.ppmPrecursor) break;

            //Make sure peptide has valid XL motif
            for (x = 0; x < 20; x++){
              if (s[n].motif[x] == -1) break;
              if (s[n].motif[x] == counterMotif) break;
            }
            if (x == 20) {
              n--;
              continue;
            } else if (s[n].motif[x] == -1){
              n--;
              continue;
            }

            p = db->getPeptide(s[n].pep1);
            cp2 = db->at(p.map->at(0).index).sequence[p.map->at(0).start + s[n].k1];
            len2 = (int)db->at(p.map->at(0).index).sequence.size() - 1;

            //Grab bin coordinates of all fragment ions
            ions[iIndex].setPeptide(true, &db->at(p.map->at(0).index).sequence[p.map->at(0).start], p.map->at(0).stop - p.map->at(0).start + 1, p.mass, p.nTerm, p.cTerm,p.n15);
            ions[iIndex].buildSingletIons((int)s[n].k1);
            setBinList(&msPartner, iIndex, sp->getPrecursor(m).charge, sp->getPrecursor(m).monoMass-s[n].mass, sp->getSingletScoreCard(s[n].rank).mods, sp->getSingletScoreCard(s[n].rank).modLen);
            dShared = sharedScore(sp,&msTemplate,&msPartner,sp->getPrecursor(m).charge);
 
            sc.simpleScore=(float)(s[j].simpleScore*s[j].len+s[n].simpleScore*s[n].len-dShared);
            sc.k1=s[j].k1;
            sc.k2=s[n].k1;
            sc.mass=totalMass;
            sc.linkable1=s[j].linkable;
            sc.linkable2=s[n].linkable;
            sc.pep1=s[j].pep1;
            sc.pep2=s[n].pep1;
            sc.link=xlIndex;
            sc.rank1=s[j].rank;
            sc.rank2=s[n].rank;
            sc.score1=s[j].simpleScore*s[j].len;
            sc.score2=s[n].simpleScore*s[n].len;
            sc.mass1=s[j].mass;
            sc.mass2=s[n].mass;
            sc.mods1->clear();
            sc.mods2->clear();
            sc1=sp->getSingletScoreCard(s[j].rank);
            for(x=0;x<sc1.modLen;x++) sc.mods1->push_back(sc1.mods[x]);
            sc1=sp->getSingletScoreCard(s[n].rank);
            for(x=0;x<sc1.modLen;x++) sc.mods2->push_back(sc1.mods[x]);
            sp->checkScore(sc);
            matches.push_back(n);

            n--;
          }
        } //for(m)
      } //for(cm)
    } //for(i)

    //Mark peptide as searched
    s[j].simpleScore=-s[j].simpleScore;
  }
  
  //reset scores
  for(j=0;j<count;j++){
    if(s[j].simpleScore<0) s[j].simpleScore=-s[j].simpleScore;
  }

 
  delete [] s;
  
}

bool KAnalysis::analyzeSinglets(kPeptide& pep, int index, double lowLinkMass, double highLinkMass, int iIndex){

  int i;
  size_t j;
  int k;
  int len;
  char mot[10];
  int m,n,x;
  double minMass;
  double maxMass;
  vector<int> scanIndex;
  string pepSeq;
  bool bSearch;

  //get the peptide sequence
  db->getPeptideSeq(pep,pepSeq);

  //Set Mass boundaries
  minMass=pep.mass+lowLinkMass+params.minPepMass;
  maxMass=pep.mass+highLinkMass+params.maxPepMass;
  minMass-=(minMass/1000000*params.ppmPrecursor);
  maxMass+=(maxMass/1000000*params.ppmPrecursor);

  //Find mod mass as difference between precursor and peptide
  
  len=(pep.map->at(0).stop-pep.map->at(0).start)+1;
  ions[iIndex].setPeptide(true, &db->at(pep.map->at(0).index).sequence[pep.map->at(0).start], len, pep.mass, pep.nTerm, pep.cTerm,pep.n15);
  
  //Iterate every link site
  for(k=0;k<len;k++){
    if (params.turbo){
      m=0;
      if (k == len - 1 && pep.cTerm){ //check if we are at the c-terminus on a c-terminal peptide
        i=0;
        while (xlTable['c'][i]>-1){
          for (n=0;n<m;n++){
            if (xlTable['c'][i]==mot[n]) break;
          }
          if (n==m) mot[m++] = xlTable['c'][i];
          i++;
        }
      } else if (k == len - 1) { //if we are at the c-term of any other peptide, stop the analysis because it cannot be linked.
        continue;
      } else {
        i=0;
        while (xlTable[pepSeq[k]][i]>-1){ //check if amino acid can be linked
          for (n = 0; n<m; n++){
            if (xlTable[pepSeq[k]][i] == mot[n]) break;
          }
          if (n == m) mot[m++] = xlTable[pepSeq[k]][i];
          i++;
        }
        if (m==0 && k == 0 && pep.nTerm){ //if it cannot be linked, but is the n-terminus of a protein, check for a linker
          i = 0;
          while (xlTable['n'][i]>-1){
            for (n = 0; n<m; n++){
              if (xlTable['n'][i] == mot[n]) break;
            }
            if (n == m) mot[m++] = xlTable['n'][i];
            i++;
          }
        }
      }
      if (m==0) continue;
    } else { //shorthand version of above approach when not in turbo mode.
      if (k == len - 1 && pep.cTerm){
        if (xlTable['c'][0] == -1)  continue;
      } else if (k==len-1) {
        continue;
      } else if (xlTable[pepSeq[k]][0]==-1) {
        if (k == 0 && pep.nTerm){
          if(xlTable['n'][0] == -1) continue;
        } else {
          continue;
        }
      }
    }
  
    //build fragment ions and score against all potential spectra
    ions[iIndex].reset();
    ions[iIndex].buildSingletIons(k);
    ions[iIndex].modIonsRec2(0,k,0,0,true);

    //iterate through all ion sets
    for(i=0;i<ions[iIndex].size();i++){

      //Iterate all spectra from (peptide mass + low linker + minimum mass) to (peptide mass + high linker + maximum mass)
      if (!spec->getBoundaries(minMass + ions[iIndex][i].difMass, maxMass + ions[iIndex][i].difMass, scanIndex, scanBuffer[iIndex])) continue;

      for(j=0;j<scanIndex.size();j++){
        if (params.turbo){
          for (n = 0; n < m; n++){
            x=0;
            while (spec->getCounterMotif(mot[n], x)>-1){
              bSearch = scoreSingletSpectra2(scanIndex[j], i, ions[iIndex][i].mass, spec->getLink(spec->getXLIndex((int)mot[n], x)).mass, spec->getCounterMotif((int)mot[n], x), len, index, (char)k, minMass, iIndex);
              if (bSearch) break;
              x++;
            }
            if (bSearch) break;
          }
        } else {
          scoreSingletSpectra(scanIndex[j],i,ions[iIndex][i].mass,len,index,(char)k,minMass,iIndex);
        }
      }
    }

  }

  return true;
}

/*============================
  Private Functions
============================*/
bool KAnalysis::allocateMemory(int threads){
  size_t j,k;
  bKIonsManager = new bool[threads];
  ions = new KIons[threads];
  scanBuffer = new bool*[threads];
  for(int i=0;i<threads;i++) {
    bKIonsManager[i]=false;
    ions[i].setModFlags(params.monoLinksOnXL,params.diffModsOnXL);
    scanBuffer[i] = new bool[spec->size()];
    for(j=0;j<params.xLink->size();j++){
      for(k=0;k<params.xLink->at(j).motifA.size();k++){
        ions[i].site[params.xLink->at(j).motifA[k]]=true;
      }
      for (k = 0; k<params.xLink->at(j).motifB.size(); k++){
        ions[i].site[params.xLink->at(j).motifB[k]] = true;
      }
    }
  }
  return true;
}

void KAnalysis::checkXLMotif(int motifA, char* motifB, vector<int>& v){
  int i, j;
  int cm;
  v.clear();
  for (i = 0; i<10; i++){
    cm = spec->getCounterMotif(motifA, i);
    if (cm<0) return;
    for (j = 0; j<20; j++){
      if (motifB[j]<0) break;
      if (motifB[j] == cm) {
        v.push_back(spec->getXLIndex(motifA, i));
      }
    }
  }
  return;
}

void KAnalysis::deallocateMemory(int threads){
  delete [] bKIonsManager;
  delete [] ions;
  for (int i = 0; i < threads; i++){
    delete[] scanBuffer[i];
  }
  delete[] scanBuffer;
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

void KAnalysis::scoreSingletSpectra(int index, int sIndex, double mass, int len, int pep, char k, double minMass, int iIndex){
  kSingletScoreCard sc;
  KIonSet* iset;
  kPepMod mod;
  double score1=0;
  double score2=0;
  int i;
  //int sk = 0;
  vector<kPepMod> v;

  KSpectrum* s=spec->getSpectrum(index);
  kPrecursor* p;
  int sz=s->sizePrecursor();
  for(i=0;i<sz;i++){
    p=s->getPrecursor2(i);
    if(p->monoMass>minMass){
      //sk++;
      if(params.xcorr) score1=xCorrScoring(*s,p->monoMass-mass,sIndex,iIndex);
      else score1=kojakScoring(index,p->monoMass-mass,sIndex,iIndex);
      if(score1>score2) score2=score1;
    }
  }

  //Check the highest score
  sc.len=len;
  sc.simpleScore=(float)score2/len;
  sc.k1=k;
  sc.linkable=false;
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
          if (i == 0 && iset->modNTerm) mod.term = true;
          else if (i == ions[iIndex].getIonCount() - 1 && iset->modCTerm) mod.term = true;
          else mod.term = false;
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

  //Threading::LockMutex(mutexSpecScore[index]);
  //spec->at(index).cc += sk;
  //Threading::UnlockMutex(mutexSpecScore[index]);

}

bool KAnalysis::scoreSingletSpectra2(int index, int sIndex, double mass, double xlMass, int counterMotif, int len, int pep, char k, double minMass, int iIndex){
  kSingletScoreCard sc;
  KIonSet* iset;
  kPepMod mod;
  double score1 = 0;
  double score2 = 0;
  int bestPre = 0;
  int i;
  //int sk=0;
  //int nsk=0;
  vector<kPepMod> v;

  double low,high,m;
  int lowI,highI;
  bool bScored=false;
  int max = (int)((params.maxPepMass + 1000) / 0.015);

  //cout << "scoreSingletSpectra2" << endl;

  KSpectrum* s = spec->getSpectrum(index);
  kPrecursor* p;
  int sz = s->sizePrecursor();
  //cout << "started loop" << endl;
  for (i = 0; i<sz; i++){
    p = s->getPrecursor2(i);
    if (p->monoMass>minMass){

      low = p->monoMass;
      high = low;
      m = low / 1000000 * params.ppmPrecursor;
      low -= m;
      high += m;
      m = mass + xlMass;
      low -= m;
      high -= m;
      lowI = (int)(low/0.015);
      highI = (int)(high/0.015)+1;

      if (lowI<0) lowI=0;
      if (highI>max) highI=max;

      while (lowI < highI){
        if (pepBin[counterMotif][lowI]) break;
        lowI++;
      }
      if (lowI >= highI){
        //sk++;
        continue;
      }
      
      //nsk++;
      if (params.xcorr) score1 = xCorrScoring(*s, p->monoMass - mass, sIndex, iIndex);
      else score1 = kojakScoring(index, p->monoMass - mass, sIndex, iIndex);
      if (score1>score2) {
        score2 = score1;
        bestPre = i;
      }
      bScored=true;
    }
  }

  //Check the highest score
  sc.len = len;
  sc.simpleScore = (float)score2 / len;
  sc.k1 = k;
  sc.linkable = false;
  sc.pep1 = pep;
  sc.pre = (char)bestPre;
  sc.mass = mass;
  if (sc.simpleScore>0) {
    v.clear();
    if (sc.mods != NULL) {
      sc.modLen = 0;
      delete[] sc.mods;
      sc.mods = NULL;
    }
    iset = ions[iIndex].at(sIndex);
    if (iset->difMass != 0){
      for (i = 0; i<ions[iIndex].getIonCount(); i++) {
        if (iset->mods[i] != 0){
          if (i==0 && iset->modNTerm) mod.term=true;
          else if (i == ions[iIndex].getIonCount()-1 && iset->modCTerm) mod.term=true;
          else mod.term=false;
          mod.pos = (char)i;
          mod.mass = iset->mods[i];
          v.push_back(mod);
        }
      }
      sc.modLen = (char)v.size();
      sc.mods = new kPepMod[sc.modLen];
      for (i = 0; i<(int)sc.modLen; i++) sc.mods[i] = v[i];
    }

    Threading::LockMutex(mutexSpecScore[index]);
    spec->at(index).checkSingletScore(sc);
    Threading::UnlockMutex(mutexSpecScore[index]);

  }

  //Threading::LockMutex(mutexSpecScore[index]);
  //spec->at(index).sc += sk;
  //spec->at(index).cc += nsk;
  //Threading::UnlockMutex(mutexSpecScore[index]);
  return bScored;
}

void KAnalysis::scoreSpectra(vector<int>& index, int sIndex, double modMass, int pep1, int pep2, int k1, int k2, int link, int iIndex){
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
    sc.linkable1=sc.linkable2=false;
    sc.pep1=pep1;
    sc.pep2=pep2;
    sc.link=link;
    if(ions[iIndex][sIndex].difMass!=0){
      for(i=0;i<ions[iIndex].getPeptideLen();i++) {
        if(ions[iIndex][sIndex].mods[i]!=0){
          if (i == 0 && ions[iIndex][sIndex].modNTerm) mod.term = true;
          else if (i == ions[iIndex].getIonCount() - 1 && ions[iIndex][sIndex].modCTerm) mod.term = true;
          else mod.term = false;
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
  int maxCharge=s->getCharge();  

  int i,j,k;
  int key;
  int pos;

  //Assign ion series
  double***  ionSeries;
  ionSeries=new double**[numIonSeries];
  k=0;
  if(params.ionSeries[0]) ionSeries[k++]=ki->aIons;
  if(params.ionSeries[1]) ionSeries[k++]=ki->bIons;
  if(params.ionSeries[2]) ionSeries[k++]=ki->cIons;
  if(params.ionSeries[3]) ionSeries[k++]=ki->xIons;
  if(params.ionSeries[4]) ionSeries[k++]=ki->yIons;
  if(params.ionSeries[5]) ionSeries[k++]=ki->zIons;

  //The number of fragment ion series to analyze is PrecursorCharge-1
  //However, don't analyze past the 5+ series
  if(maxCharge>6) maxCharge=6;

  //Iterate all series
  for(k=1;k<maxCharge;k++){

    dif=modMass/k;

    //Iterate through pfFastXcorrData
    for(j=0;j<numIonSeries;j++){
      for(i=0;i<ionCount;i++){

        //get key
        if(ionSeries[j][k][i]<0) mz = params.binSize * (int)((dif-ionSeries[j][k][i])*invBinSize+binOffset);
        else mz = params.binSize * (int)(ionSeries[j][k][i]*invBinSize+binOffset);
        key = (int)mz;
        if(key>=s->kojakBins) break;
        if(s->kojakSparseArray[key]==NULL) continue;
        pos = (int)((mz-key)*invBinSize);
        dXcorr += s->kojakSparseArray[key][pos];
      }
    }

  }

  //Scale score appropriately
  if(dXcorr <= 0.0) dXcorr=0.0;
  else dXcorr *= 0.005;

  //Clean up memory
  k = 0;
  if (params.ionSeries[0]) ionSeries[k++] = NULL;
  if (params.ionSeries[1]) ionSeries[k++] = NULL;
  if (params.ionSeries[2]) ionSeries[k++] = NULL;
  if (params.ionSeries[3]) ionSeries[k++] = NULL;
  if (params.ionSeries[4]) ionSeries[k++] = NULL;
  if (params.ionSeries[5]) ionSeries[k++] = NULL;

  delete[] ionSeries;

  return float(dXcorr);
}

void KAnalysis::setBinList(kMatchSet* m, int iIndex, int charge, double preMass, kPepMod* mods, char modLen){

  KIonSet* ki=ions[iIndex].at(0);
  double** ionSeries;
  double invBinSize=1.0/params.binSize;
  double mz;
  double dif;
  int key;
  int ionCount=ions[iIndex].getIonCount();
  int i,j;

  //set up all modification mass modifiers
  double* mod;
  double* modRev;
  mod = new double[ionCount];
  modRev = new double[ionCount];
  for(i=0;i<ionCount;i++){
    mod[i]=0;
    modRev[i]=0;
  }
  for(i=0;i<(int)modLen;i++){
    for(j=mods[i].pos;j<ionCount;j++){
      mod[j]+=mods[i].mass;
    }
    for(j=ionCount-mods[i].pos;j<ionCount;j++){
      modRev[j]+=mods[i].mass;
    }
  }
  
  if(charge>6) charge=6;
  m->allocate(ions[iIndex].getIonCount(),charge);

  //populate structure
  for(i=1;i<charge;i++){

    dif=preMass/i;

    //a-ions
    if(params.ionSeries[0]) {
      ionSeries=ki->aIons;
      for(j=0;j<ionCount;j++) {
        if(ionSeries[i][j]<0) mz = params.binSize * (int)((dif-(ionSeries[i][j]-mod[j]/i))*invBinSize+params.binOffset);
        else mz = params.binSize * (int)((ionSeries[i][j]+mod[j]/i)*invBinSize+params.binOffset);
        key = (int)mz;
        m->a[i][j].pos = (int)((mz-key)*invBinSize);
        m->a[i][j].key = key;
      }
    }

    //b-ions
    if(params.ionSeries[1]) {
      ionSeries=ki->bIons;
      for(j=0;j<ionCount;j++) {
        if(ionSeries[i][j]<0) mz = params.binSize * (int)((dif-(ionSeries[i][j]-mod[j]/i))*invBinSize+params.binOffset);
        else mz = params.binSize * (int)((ionSeries[i][j]+mod[j]/i)*invBinSize+params.binOffset);
        key = (int)mz;
        m->b[i][j].pos = (int)((mz-key)*invBinSize);
        m->b[i][j].key = key;
      }
    }

    //c-ions
    if(params.ionSeries[2]) {
      ionSeries=ki->cIons;
      for(j=0;j<ionCount;j++) {
        if(ionSeries[i][j]<0) mz = params.binSize * (int)((dif-(ionSeries[i][j]-mod[j]/i))*invBinSize+params.binOffset);
        else mz = params.binSize * (int)((ionSeries[i][j]+mod[j]/i)*invBinSize+params.binOffset);
        key = (int)mz;
        m->c[i][j].pos = (int)((mz-key)*invBinSize);
        m->c[i][j].key = key;
      }
    }

    //x-ions
    if(params.ionSeries[3]) {
      ionSeries=ki->xIons;
      for(j=0;j<ionCount;j++) {
        if(ionSeries[i][j]<0) mz = params.binSize * (int)((dif-(ionSeries[i][j]-modRev[j]/i))*invBinSize+params.binOffset);
        else mz = params.binSize * (int)((ionSeries[i][j]+modRev[j]/i)*invBinSize+params.binOffset);
        key = (int)mz;
        m->x[i][j].pos = (int)((mz-key)*invBinSize);
        m->x[i][j].key = key;
      }
    }

    //y-ions
    if(params.ionSeries[4]) {
      ionSeries=ki->yIons;
      for(j=0;j<ionCount;j++) {
        if(ionSeries[i][j]<0) mz = params.binSize * (int)((dif-(ionSeries[i][j]-modRev[j]/i))*invBinSize+params.binOffset);
        else mz = params.binSize * (int)((ionSeries[i][j]+modRev[j]/i)*invBinSize+params.binOffset);
        key = (int)mz;
        m->y[i][j].pos = (int)((mz-key)*invBinSize);
        m->y[i][j].key = key;
      }
    }

    //z-ions
    if(params.ionSeries[5]) {
      ionSeries=ki->zIons;
      for(j=0;j<ionCount;j++) {
        if(ionSeries[i][j]<0) mz = params.binSize * (int)((dif-(ionSeries[i][j]-modRev[j]/i))*invBinSize+params.binOffset);
        else mz = params.binSize * (int)((ionSeries[i][j]+modRev[j]/i)*invBinSize+params.binOffset);
        key = (int)mz;
        m->z[i][j].pos = (int)((mz-key)*invBinSize);
        m->z[i][j].key = key;
      }
    }

  }

  delete [] mod;
  delete [] modRev;

}

double KAnalysis::sharedScore(KSpectrum* s, kMatchSet* m1, kMatchSet* m2, int charge){

  //cout << "In sharedScore: " << s->getScanNumber() << endl;
  int i,j,k;
  int maxCharge=charge;
  double dScore=0;

  if(maxCharge>6) maxCharge=6;

  /*
  if(s->getScanNumber()==62111){
    cout << "Pep 1: " << m1->sz << endl;
    for(i=0;i<m1->sz;i++) {
      cout << m1->b[1][i].key << "\t" << m1->b[1][i].pos << endl;
      cout << m1->b[2][i].key << "\t" << m1->b[2][i].pos << endl;
    }
    for(i=0;i<m1->sz;i++) {
      cout << m1->y[1][i].key << "\t" << m1->y[1][i].pos << endl;
      cout << m1->y[2][i].key << "\t" << m1->y[2][i].pos << endl;
    }
    cout << "Pep 2: " << m2->sz << endl;
    for(i=0;i<m2->sz;i++) {
      cout << m2->b[1][i].key << "\t" << m2->b[1][i].pos << endl;
      cout << m2->b[2][i].key << "\t" << m2->b[2][i].pos << endl;
    }
    for(i=0;i<m2->sz;i++) {
      cout << m2->y[1][i].key << "\t" << m2->y[1][i].pos << endl;
      cout << m2->y[2][i].key << "\t" << m2->y[2][i].pos << endl;
    }
  }
  */

  //for(i=0;i<s->sizePrecursor();i++){
  //  cout << "precursor: " << i << "\t" << s->getPrecursor(i).monoMass << endl;
  //}
  
  for(i=1;i<maxCharge;i++){
    //cout << "SS, charge: " << i << endl;
    //a-ions
    if(params.ionSeries[0]) {
      j=0;
      k=0;
      while(j<m1->sz && k<m2->sz){
        if(m1->a[i][j].key==m2->a[i][k].key){
          if(m1->a[i][j].pos==m2->a[i][k].pos){
            if(m2->a[i][k].key>=s->kojakBins) break;
            if(s->kojakSparseArray[m2->a[i][k].key]!=NULL){
              dScore+=s->kojakSparseArray[m2->a[i][k].key][m2->a[i][k].pos];
            }
            j++;
            k++;
          } else {
            if(m1->a[i][j].pos<m2->a[i][k].pos) j++;
            else k++;
          }
        } else {
          if(m1->a[i][j].key<m2->a[i][k].key) j++;
          else k++;
        }
      }
    }

    //b-ions
    if(params.ionSeries[1]) {
      //cout << "SS, B" << endl;
      j=0;
      k=0;
      while(j<m1->sz && k<m2->sz){
        //cout << "J: " << j << " " << m1->b[i][j].key << "," << m1->b[i][j].pos;
        //cout << "\tK: " << k << " " << m2->b[i][k].key << "," << m2->b[i][k].pos << endl;
        if(m1->b[i][j].key==m2->b[i][k].key){
          if(m1->b[i][j].pos==m2->b[i][k].pos){
            //cout << "Match! max = " << s->kojakBins << endl;
            if(m2->b[i][k].key>=s->kojakBins) break;
            if(s->kojakSparseArray[m2->b[i][k].key]!=NULL){
              dScore+=s->kojakSparseArray[m2->b[i][k].key][m2->b[i][k].pos];
            }
            j++;
            k++;
          } else {
            if(m1->b[i][j].pos<m2->b[i][k].pos) j++;
            else k++;
          }
        } else {
          if(m1->b[i][j].key<m2->b[i][k].key) j++;
          else k++;
        }
      }
    }

    //c-ions
    if(params.ionSeries[2]) {
      j=0;
      k=0;
      while(j<m1->sz && k<m2->sz){
        if(m1->c[i][j].key==m2->c[i][k].key){
          if(m1->c[i][j].pos==m2->c[i][k].pos){
            if(m2->c[i][k].key>=s->kojakBins) break;
            if(s->kojakSparseArray[m2->c[i][k].key]!=NULL){
              dScore+=s->kojakSparseArray[m2->c[i][k].key][m2->c[i][k].pos];
            }
            j++;
            k++;
          } else {
            if(m1->c[i][j].pos<m2->c[i][k].pos) j++;
            else k++;
          }
        } else {
          if(m1->c[i][j].key<m2->c[i][k].key) j++;
          else k++;
        }
      }
    }

    //x-ions
    if(params.ionSeries[3]) {
      j=0;
      k=0;
      while(j<m1->sz && k<m2->sz){
        if(m1->x[i][j].key==m2->x[i][k].key){
          if(m1->x[i][j].pos==m2->x[i][k].pos){
            if(m2->x[i][k].key>=s->kojakBins) break;
            if(s->kojakSparseArray[m2->x[i][k].key]!=NULL){
              dScore+=s->kojakSparseArray[m2->x[i][k].key][m2->x[i][k].pos];
            }
            j++;
            k++;
          } else {
            if(m1->x[i][j].pos<m2->x[i][k].pos) j++;
            else k++;
          }
        } else {
          if(m1->x[i][j].key<m2->x[i][k].key) j++;
          else k++;
        }
      }
    }

    //y-ions
    if(params.ionSeries[4]) {
      //cout << "SS, Y" << endl;
      j=0;
      k=0;
      while(j<m1->sz && k<m2->sz){
        //cout << "J: " << j << " " << m1->y[i][j].key << "," << m1->y[i][j].pos;
        //cout << "\tK: " << k << " " << m2->y[i][k].key << "," << m2->y[i][k].pos << endl;
        if(m1->y[i][j].key==m2->y[i][k].key){
          if(m1->y[i][j].pos==m2->y[i][k].pos){
            //cout << "Match! max = " << s->kojakBins << endl;
            if(m2->y[i][k].key>=s->kojakBins) break;
            if(s->kojakSparseArray[m2->y[i][k].key]!=NULL){
              dScore+=s->kojakSparseArray[m2->y[i][k].key][m2->y[i][k].pos];
            }
            j++;
            k++;
          } else {
            if(m1->y[i][j].pos<m2->y[i][k].pos) j++;
            else k++;
          }
        } else {
          if(m1->y[i][j].key<m2->y[i][k].key) j++;
          else k++;
        }
      }
    }

    //z-ions
    if(params.ionSeries[5]) {
      j=0;
      k=0;
      while(j<m1->sz && k<m2->sz){
        if(m1->z[i][j].key==m2->z[i][k].key){
          if(m1->z[i][j].pos==m2->z[i][k].pos){
            if(m2->z[i][k].key>=s->kojakBins) break;
            if(s->kojakSparseArray[m2->z[i][k].key]!=NULL){
              dScore+=s->kojakSparseArray[m2->z[i][k].key][m2->z[i][k].pos];
            }
            j++;
            k++;
          } else {
            if(m1->z[i][j].pos<m2->z[i][k].pos) j++;
            else k++;
          }
        } else {
          if(m1->z[i][j].key<m2->z[i][k].key) j++;
          else k++;
        }
      }
    }

  }

  if(dScore <= 0.0) dScore=0.0;
  else dScore *= 0.005;

  //cout << "Out sharedScore" << endl;
  return dScore;

}

void KAnalysis::makePepLists(){
  int j;
  int m;
  bool bMatch;
  size_t i,k;
  vector<double> v;
  vector<kPeptide>* p=db->getPeptideList();
  string pep;
  
  pepMass = new double*[spec->getMotifCount()];
  pepMassSize = new int[spec->getMotifCount()];
  pepBin = new bool*[spec->getMotifCount()];
  pepBinSize = new int[spec->getMotifCount()];

  //iterate through each motif
  for (j = 0; j<spec->getMotifCount(); j++){
    v.clear();

    //iterate through peptide list
    for (i = 0; i < p->size(); i++){
      //cout << i << " of " << p->size() << endl;

      //skip peptides that cannot be linked
      if (p->at(i).xlSites==0) continue;

      //determine which motifs are allowed for this peptide
      db->getPeptideSeq(p->at(i),pep);
      bMatch=false;
      for (k = 0; k < pep.size(); k++){
        if (xlTable[pep[k]][0]>-1) {
          for (m = 0; m < 20; m++){
            if (xlTable[pep[k]][m]==-1) break;
            if (xlTable[pep[k]][m] == (char)j) {
              bMatch=true;
              break;
            }
          }
        }
        if (bMatch) break; //if we match our desired motif, stop searching
      }

      //skip peptide if it doesn't have the desired motif anywhere
      if (!bMatch) continue;

      //add peptide mass to our list
      v.push_back(p->at(i).mass);

      //add all possible modification masses too
      if (!params.diffModsOnXL && !params.monoLinksOnXL) continue;

      //These should be modified for quick mass calculations without ion series expansion
      ions[0].setPeptide(true, &db->at(p->at(i).map->at(0).index).sequence[p->at(i).map->at(0).start], p->at(i).map->at(0).stop - p->at(i).map->at(0).start + 1, p->at(i).mass, p->at(i).nTerm, p->at(i).cTerm,p->at(i).n15);
      ions[0].buildIons();
      ions[0].modIonsRec(0, -1, 0, 0, false);
      for (m = 1; m<ions[0].size(); m++) v.push_back(ions[0][m].mass);
      
    } //i

    //copy our list to memory
    pepMassSize[j] = (int)v.size();
    pepMass[j] = new double[v.size()];
    for (i = 0; i<v.size(); i++) pepMass[j][i] = v[i];

    //sort list
    qsort(pepMass[j], pepMassSize[j], sizeof(double),compareD);

    pepBinSize[j] = (int)((params.maxPepMass+1000)/0.015);
    pepBin[j] = new bool[pepBinSize[j]];
    for (i = 0; i<pepBinSize[j]; i++) pepBin[j][i]=false;
    for (i = 0; i<pepMassSize[j]; i++) pepBin[j][(int)(pepMass[j][i]/0.015)] = true;

  }

  p = NULL;

}

bool KAnalysis::findCompMass(int motif, double low, double high){
  int sz = pepMassSize[motif];
  int lower = 0;
  int mid = sz / 2;
  int upper = sz;

  //binary search to boundaries
  while (lower < upper){
    //cout << lower << " " << mid << " " << upper << "\t" << pepMass[motif][mid] << "\t" << low << "-" << high << endl;
    if (pepMass[motif][mid]<low){
      lower = mid + 1;
      mid = (lower + upper) / 2;
    } else if (pepMass[motif][mid]>high) {
      upper = mid - 1;
      mid = (lower + upper) / 2;
    } else {
      return true;
    }
  }
  //exit(1);
  return false;
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
