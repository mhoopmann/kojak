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

using namespace std;

bool*       KAnalysis::bKIonsManager;
KDatabase*  KAnalysis::db;
double      KAnalysis::highLinkMass;
KIons*      KAnalysis::ions;
double      KAnalysis::lowLinkMass;
double      KAnalysis::maxMass;
double      KAnalysis::minMass;
Mutex       KAnalysis::mutexKIonsManager;
Mutex*      KAnalysis::mutexSpecScore;
Mutex**     KAnalysis::mutexSingletScore;
kParams     KAnalysis::params;
KData*      KAnalysis::spec;
char**      KAnalysis::xlTable;
bool**      KAnalysis::scanBuffer;

int         KAnalysis::numIonSeries;

int*        KAnalysis::pepMassSize;
double**    KAnalysis::pepMass;
bool**      KAnalysis::pepBin;
int*        KAnalysis::pepBinSize;

bool*       KAnalysis::soloLoop;
bool        KAnalysis::firstPass;

int KAnalysis::skipCount;
int KAnalysis::nonSkipCount;

KDecoys KAnalysis::decoys;
KLog*   KAnalysis::klog;

/*============================
  Constructors & Destructors
============================*/
KAnalysis::KAnalysis(kParams& p, KDatabase* d, KData* dat){
  unsigned int i;
  int j,k;
  
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
    for(i=0;i<params.aaMass->size();i++) ions[j].setAAMass((char)params.aaMass->at(i).index, params.aaMass->at(i).mass, params.aaMass->at(i).xl);
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
  mutexSingletScore = new Mutex*[spec->size()];
  mutexSpecScore = new Mutex[spec->size()];
  for(j=0;j<spec->size();j++){
    Threading::CreateMutex(&mutexSpecScore[j]);
    mutexSingletScore[j] = new Mutex[spec->at(j).sizePrecursor()];
    for(k=0;k<spec->at(j).sizePrecursor();k++){
      Threading::CreateMutex(&mutexSingletScore[j][k]);
    }
  }

  decoys.decoySize=params.decoySize;

  makePepLists();
  skipCount=0;
  nonSkipCount=0;
  klog=NULL;
}

KAnalysis::~KAnalysis(){
  int i,j;

  //Destroy mutexes
  Threading::DestroyMutex(mutexKIonsManager);
  for(i=0;i<spec->size();i++){
    Threading::DestroyMutex(mutexSpecScore[i]);
    for(j=0;j<spec->at(i).sizePrecursor();j++){
      Threading::DestroyMutex(mutexSingletScore[i][j]);
    }
    delete [] mutexSingletScore[i];
  }
  delete [] mutexSingletScore;
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
  klog=NULL;
}

//============================
//  Public Functions
//============================
bool KAnalysis::doPeptideAnalysis(){
  size_t i;
  int iPercent;
  int iTmp;
  vector<kPeptide>* p;
  vector<int> index;
  vector<kPepMod> mods;

  kScoreCard sc;

  firstPass=true;

  ThreadPool<kAnalysisStruct*>* threadPool = new ThreadPool<kAnalysisStruct*>(analyzePeptideProc,params.threads,params.threads,1);

  //Set progress meter
  iPercent=0;
  printf("%2d%%",iPercent);
  fflush(stdout);

  //Set which list of peptides to search (with and without internal lysine)
  p=db->getPeptideList();

  //track non-links and loops
  soloLoop = new bool[p->size()];
  for(i=0;i<p->size();i++) soloLoop[i]=false;

  //get boundaries for first pass
  double lowerBound=(spec->getMinMass()-lowLinkMass)/2-0.25;
  double upperBound=spec->getMaxMass()-highLinkMass-params.minPepMass+0.25;

  //Iterate the peptide for the first pass
  for(i=0;i<p->size();i++){

    if(p->at(i).mass>upperBound) continue;
    if(p->at(i).mass<lowerBound) break;

    threadPool->WaitForQueuedParams();

    kAnalysisStruct* a = new kAnalysisStruct(&mutexKIonsManager,&p->at(i),(int)i);
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


  //Perform the second pass
  firstPass=false;
  if(klog!=NULL) klog->addMessage("Scoring peptides (second pass).",true);
  cout << "  Second pass ... ";

  //get boundary for second pass
  upperBound = (spec->getMaxMass() - lowLinkMass)/2;

  //Set progress meter
  iPercent = 0;
  printf("%2d%%", iPercent);
  fflush(stdout);

  i=p->size() - 1;
  while(true){
    if (p->at(i).mass>upperBound) break;

    threadPool->WaitForQueuedParams();

    kAnalysisStruct* a = new kAnalysisStruct(&mutexKIonsManager, &p->at(i), (int)i);
    threadPool->Launch(a);

    //Update progress meter
    iTmp = (int)((1.0-(double)i / p->size()) * 100);
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }
   
    if(i==0) break;
    i--;
  }
  threadPool->WaitForQueuedParams();
  threadPool->WaitForThreads();

  //Finalize progress meter
  if(iPercent<100) printf("\b\b\b100%%");
  cout << endl;

  //clean up memory & release pointers
  delete [] soloLoop;
  delete threadPool;
  threadPool=NULL;
  p=NULL;
  return true;
}

bool KAnalysis::doEValueAnalysis(){
  int i;
  int iPercent;
  int iTmp;

  ThreadPool<KSpectrum*>* threadPool = new ThreadPool<KSpectrum*>(analyzeEValueProc, params.threads, params.threads, 1);

  //Set progress meter
  iPercent = 0;
  printf("%2d%%", iPercent);
  fflush(stdout);

  //Iterate the peptide for the first pass
  for (i = 0; i<spec->size(); i++){

    threadPool->WaitForQueuedParams();

    KSpectrum* a = &spec->at(i);
    threadPool->Launch(a);

    //Update progress meter
    iTmp = (int)((double)i / spec->size() * 100);
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
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

void KAnalysis::analyzeEValueProc(KSpectrum* s){
  s->calcEValue(&params, decoys,*db);
  s = NULL;
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
  //if(strcmp(str,"VPSKK")==0) cout << str << "\t" << p->mass << endl;
  //Set the peptide, calc the ions, and score it against the spectra
  ions[iIndex].setPeptide(true,&db->at(p->map->at(0).index).sequence[p->map->at(0).start],p->map->at(0).stop-p->map->at(0).start+1,p->mass,p->nTerm,p->cTerm,p->n15);
  
  if(!soloLoop[pepIndex]){ //if we've searched this peptide as solo in the first pass, skip doing so again
    ions[iIndex].buildIons();
    ions[iIndex].modIonsRec2(0,-1,0,0,false);

    for(j=0;j<ions[iIndex].size();j++){
      bt=spec->getBoundaries2(ions[iIndex][j].mass,params.ppmPrecursor,index,scanBuffer[iIndex]);
      if(bt) scoreSpectra(index,j,ions[iIndex][j].difMass,pepIndex,-1,-1,-1,-1,iIndex,-1,-1);
    }

    //search non-covalent dimerization if requested by user
    /* Deprecate this
    if(params.dimers) {
      analyzeSingletsNC(*p, pepIndex, iIndex); //this may need updating since 1.6.0
    }
    */

    if(p->xlSites==0) {
      soloLoop[pepIndex]=true;
      return true;
    }
  } else {
    if (p->xlSites == 0)  return true;
  }

  //Crosslinked peptides must also search singlets with reciprocol mass on each lysine
  analyzeSinglets(*p,pepIndex,lowLinkMass,highLinkMass,iIndex);

  if(p->xlSites==1) {
    soloLoop[pepIndex] = true;
    return true;
  }

  //also search loop-links
  //check loop-links by iterating through each cross-linker mass

  //if we've already searched the loop link, exit now
  if(soloLoop[pepIndex]) return true;

  string pepSeq;
  vector<int> xlIndex;
  int x;
  char site1, site2;
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
            if (p->cTerm) {
              checkXLMotif(xlTable[pepSeq[k]][x], xlTable['c'],xlIndex);
              site1 = pepSeq[k];
              site2='c';
            } else continue;
          } else if (xlTable[pepSeq[k2]][0] == -1) {
            continue;
          } else {
            checkXLMotif((int)xlTable[pepSeq[k]][x], xlTable[pepSeq[k2]],xlIndex);
            site1=pepSeq[k];
            site2=pepSeq[k2];
          }

          if (xlIndex.size()>0){
            for(k3=0;k3<xlIndex.size();k3++){
              if(xlIndex[k3]<0) continue;
              ions[iIndex].reset();
              ions[iIndex].buildLoopIons(spec->getLink(xlIndex[k3]).mass, (int)k, (int)k2);
              ions[iIndex].modLoopIonsRec2(0, (int)k, (int)k2, 0, 0, true);
              for (j = 0; j<ions[iIndex].size(); j++){
                bt = spec->getBoundaries2(ions[iIndex][j].mass, params.ppmPrecursor, index, scanBuffer[iIndex]);
                if (bt) scoreSpectra(index, j, 0, pepIndex, -1, (int)k, (int)k2, xlIndex[k3], iIndex,site1,site2);
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
            if (p->cTerm) {
              checkXLMotif(xlTable['n'][x], xlTable['c'],xlIndex);
              site1='n';
              site2='c';
            } else continue;
          } else if (xlTable[pepSeq[k2]][0] == -1) {
            continue;
          } else {
            checkXLMotif((int)xlTable['n'][x], xlTable[pepSeq[k2]],xlIndex);
            site1 = 'n';
            site2 = pepSeq[k2];
          }

          if (xlIndex.size()>0){
            for (k3 = 0; k3<xlIndex.size(); k3++){
              if(xlIndex[k3]<0) continue;
              ions[iIndex].reset();
              ions[iIndex].buildLoopIons(spec->getLink(xlIndex[k3]).mass, (int)k, (int)k2);
              ions[iIndex].modLoopIonsRec2(0, (int)k, (int)k2, 0, 0, true);
              for (j = 0; j<ions[iIndex].size(); j++){
                bt = spec->getBoundaries2(ions[iIndex][j].mass, params.ppmPrecursor, index, scanBuffer[iIndex]);
                if (bt) scoreSpectra(index, j, 0, pepIndex, -1, (int)k, (int)k2, xlIndex[k3], iIndex,site1,site2);
              }
            } //k3
          }
        } //k2

      }//x
    }

  }//k

  soloLoop[pepIndex] = true;
  return true;
}

bool KAnalysis::analyzeSinglets(kPeptide& pep, int index, double lowLinkMass, double highLinkMass, int iIndex){
  int i;
  size_t j;
  int k;
  int len;
  char mot[10];
  char site[10];
  int m,n,c;
  double minMass;
  double maxMass;
  vector<int> scanIndex;
  string pepSeq;
  bool bSearch;

  int counterMotif;
  int xlIndex;
  double xlMass;

  //get the peptide sequence
  db->getPeptideSeq(pep,pepSeq);

  //Set Mass boundaries
  if(firstPass){
    minMass = pep.mass + lowLinkMass + params.minPepMass;
    maxMass = pep.mass*2 + highLinkMass; //any higher, and this would be the smaller peptide of the cross-link
  } else {
    minMass = pep.mass*2 + lowLinkMass; //any lower, and this would be the larger peptide of the cross-link
    maxMass = pep.mass + highLinkMass + params.maxPepMass;
  }
  minMass-=(minMass/1000000*params.ppmPrecursor);
  maxMass+=(maxMass/1000000*params.ppmPrecursor);

  //Find mod mass as difference between precursor and peptide
  
  len=(pep.map->at(0).stop-pep.map->at(0).start)+1;
  ions[iIndex].setPeptide(true, &db->at(pep.map->at(0).index).sequence[pep.map->at(0).start], len, pep.mass, pep.nTerm, pep.cTerm, pep.n15);
  
  //Iterate every link site
  for(k=0;k<len;k++){
    m=0; //number of motifs (linker-to-site combinations) found in the peptide
    if (k == len - 1 && pep.cTerm){ //check if we are at the c-terminus on a c-terminal peptide
      i=0;
      while (xlTable['c'][i]>-1){ //check if c-terminus can be linked
        for (n=0;n<m;n++){ //if we've seen motif(s) already, check if we are seeing them again
          if (xlTable['c'][i]==mot[n]) break;
        }
        if (n==m) { //we have a new motif
          site[m]='c'; //mark it as c-terminal
          mot[m++] = xlTable['c'][i]; //mark the corresponding site
        }
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
        if (n == m) {
          site[m]=pepSeq[k];
          mot[m++] = xlTable[pepSeq[k]][i];
        }
        i++;
      }
      if (/*m==0 &&*/ k == 0 && pep.nTerm){ // /*if it cannot be linked,*/ but is the n-terminus of a protein, check for a linker
        i = 0;
        while (xlTable['n'][i]>-1){
          for (n = 0; n<m; n++){
            if (xlTable['n'][i] == mot[n]) break;
          }
          if (n == m) {
            site[m] = 'n';
            mot[m++] = xlTable['n'][i];
          }
          i++;
        }
      }
    }
    if (m==0) continue; //no matching motifs, so check next site on peptide
  
    //build fragment ions and score against all potential spectra
    ions[iIndex].reset();
    ions[iIndex].buildSingletIons(k);
    ions[iIndex].modIonsRec2(0,k,0,0,true);
    //ions[iIndex].makeIonIndex(params.binSize, params.binOffset);

    //iterate through all ion sets
    for(i=0;i<ions[iIndex].size();i++){
      for(j=0;j<ions[iIndex][i].len;j++){
        if (ions[iIndex][i].mods[j]==0) continue;
      }

      //Iterate all spectra from (peptide mass + low linker + minimum mass) to (peptide mass + high linker + maximum mass)
      if (!spec->getBoundaries(minMass + ions[iIndex][i].difMass, maxMass + ions[iIndex][i].difMass, scanIndex, scanBuffer[iIndex])) continue;

      //This set of iterations is slow because of the amount of iterating.
      for (n = 0; n < m; n++){ //iterate over sites
        for(c=0;c<10;c++){ //check all crosslinkers that bind to this site
          counterMotif = spec->getCounterMotif(mot[n], c);
          if (counterMotif>-1){ //only check peptide if it has a counterpart at this link site.
            xlIndex = spec->getXLIndex((int)mot[n], c);
            xlMass = spec->getLink(xlIndex).mass;
            for (j = 0; j<scanIndex.size(); j++){ //iterate over all potential spectra
              bSearch = scoreSingletSpectra2(scanIndex[j], i, ions[iIndex][i].mass, xlMass, counterMotif, len, index, (char)k, minMass, iIndex, site[n], xlIndex);
            }
          } else {
            break; //no countermotif means no more crosslinkers
          }
        }
      }
    }
  }
  return true;
}

/* Deprecating
bool KAnalysis::analyzeSingletsNC(kPeptide& pep, int index, int iIndex){
  int len;
  int i;
  size_t j;
  double minMass;
  double maxMass;
  vector<int> scanIndex;

  //Set Mass boundaries
  minMass = pep.mass + params.minPepMass;
  maxMass = pep.mass + params.maxPepMass;
  minMass -= (minMass / 1000000 * params.ppmPrecursor);
  maxMass += (maxMass / 1000000 * params.ppmPrecursor);

  //Find mod mass as difference between precursor and peptide
  len = (pep.map->at(0).stop - pep.map->at(0).start) + 1;
  ions[iIndex].setPeptide(true, &db->at(pep.map->at(0).index).sequence[pep.map->at(0).start], len, pep.mass, pep.nTerm, pep.cTerm, pep.n15);

  //build fragment ions and score against all potential spectra
  ions[iIndex].reset();
  ions[iIndex].buildIons();
  ions[iIndex].modIonsRec2(0, -1, 0, 0, false);
  //ions[iIndex].makeIonIndex(params.binSize, params.binOffset);

  //iterate through all ion sets
  for (i = 0; i<ions[iIndex].size(); i++){
    //Iterate all spectra from (peptide mass + minimum mass) to (peptide mass + maximum mass)
    if (!spec->getBoundaries(minMass + ions[iIndex][i].difMass, maxMass + ions[iIndex][i].difMass, scanIndex, scanBuffer[iIndex])) return false;

    for (j = 0; j<scanIndex.size(); j++){
      scoreSingletSpectra(scanIndex[j], i, ions[iIndex][i].mass, len, index, -1, minMass, iIndex);
    }
  }
  return true;

}
*/

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
    ions[i].setSeries(params.ionSeries[0],params.ionSeries[1],params.ionSeries[2],params.ionSeries[3],params.ionSeries[4], params.ionSeries[5]);
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

//Breakdown of the many parameters:
// index   = spectrum index in data spectra object
// sIndex  = ion set index
// mass    = peptide mass (including all modifications, but not the linker)
// xlMass  = cross-linker mass
// counterMotif = index of complementary linker site (the one on the "unknown" peptide).
// len     = peptide length (amino acids)
// pep     = peptide index in database array
// k       = link site position on peptide (relative to peptide)
// minMass = lowest allowed mass boundary on spectrum (for spectra with multiple possible precursors)
// iIndex  = thread index
// linkSite  = amino acid being linked (or 'n' or 'c')
// linkIndex = motif index of linked amino acid or terminus
bool KAnalysis::scoreSingletSpectra2(int index, int sIndex, double mass, double xlMass, int counterMotif, int len, int pep, char k, double minMass, int iIndex, char linkSite, int linkIndex){ 
  kSingletScoreCard sc;
  kSingletScoreCard* tsc;
  KIonSet* iset;
  kPepMod mod;
  float score = 0;
  int matches;
  int conFrag;
  int i,j,y;
  vector<kPepMod> v;

  double low,high,m;
  int lowI,highI;
  bool bScored=false;
  bool ret;
  int max = (int)((params.maxPepMass + 1000) / 0.015);
  int alpha;

  KSpectrum* s = spec->getSpectrum(index);
  kPrecursor* p;
  KTopPeps* tp;
  kScoreCard protSC;
  size_t x;
  int sz = s->sizePrecursor();

  string seq1,seq2;

  for (i = 0; i<sz; i++){
    p = s->getPrecursor2(i);

    if(!firstPass && (p->monoMass-spec->getLink(linkIndex).mass+0.2)/2>mass){

      //Threading::LockMutex(mutexSingletScore[index][i]);
      tp = s->getTopPeps(i);
      int ind=(int)((p->monoMass-xlMass-mass)/10);
      if(tp->singletList[ind]==NULL) {
        //Threading::UnlockMutex(mutexSingletScore[index][i]);
        continue;
      }

      //TODO: this could be made faster by checking the peptide only once for every partner of the same mass!!!
      list<kSingletScoreCard*>::iterator it = tp->singletList[ind]->begin();
      while(it!=tp->singletList[ind]->end()){
        tsc=*it;
        if(fabs((p->monoMass-tsc->mass-spec->getLink(linkIndex).mass-mass)/p->monoMass*1e6)>params.ppmPrecursor){
          it++;
          continue;
        }
        if(!spec->checkLink(linkSite,tsc->site,linkIndex)){
          it++;
          continue;
        }
        score = kojakScoring(index, p->monoMass - mass, sIndex, iIndex, matches, conFrag, p->charge);
        protSC.simpleScore = tsc->simpleScore + score;
        y = (int)(protSC.simpleScore * 10.0 + 0.5);
        if (y >= HISTOSZ) y = HISTOSZ - 1;
        Threading::LockMutex(mutexSpecScore[index]);  //no matter how low the score, put this test in our histogram.
        s->histogram[y]++;
        s->histogramCount++;
        if (score<params.minPepScore || protSC.simpleScore <= s->lowScore) { //peptide needs a minimum score, and combined score should exceed bottom of best hits
          Threading::UnlockMutex(mutexSpecScore[index]);
          it++;
          continue;
        }
        Threading::UnlockMutex(mutexSpecScore[index]);
       

        protSC.mods1->clear();
        protSC.mods2->clear();

        //alphabetize cross-linked peptides before storing them; this prevents confusion with duplications downstream
        //instead of alphabetical, order them by mass (larger first), this has implications with downstream scoring...
        //resort to alphabetical in case of equal mass
        bool bSecond=false;
        if(mass==tsc->mass){
          if(seq2.size()==0) ret=db->getPeptideSeq(pep,seq2);
          ret=db->getPeptideSeq(tsc->pep1, seq1);
          alpha=seq1.compare(seq2);
          if (alpha>0 || (alpha == 0 && k<tsc->k1)) bSecond=true;
        }
        if (mass>tsc->mass || bSecond){ //pep2 listed first
          protSC.k1 = k;
          protSC.k2 = tsc->k1;
          protSC.site1 = linkSite;
          protSC.site2 = tsc->site;
          protSC.pep1 = pep;
          protSC.pep2 = tsc->pep1;
          protSC.score1 = score;
          protSC.score2 = tsc->simpleScore;
          protSC.mass1 = mass;
          protSC.mass2 = tsc->mass;
          protSC.matches1 = matches;
          protSC.matches2 = tsc->matches;
          protSC.conFrag1 = conFrag;
          protSC.conFrag2 = tsc->conFrag;
          for (x = 0; x<tsc->modLen; x++) protSC.mods2->push_back(tsc->mods[x]);
          ret=true; //indicates where to put mods to pep2
        } else { //pep1 listed first
          protSC.k1 = tsc->k1;
          protSC.k2 = k;
          protSC.site1 = tsc->site;
          protSC.site2 = linkSite;
          protSC.pep1 = tsc->pep1;
          protSC.pep2 = pep;
          protSC.score1 = tsc->simpleScore;
          protSC.score2 = score;
          protSC.mass1 = tsc->mass;
          protSC.mass2 = mass;
          protSC.matches1 = tsc->matches;
          protSC.matches2 = matches;
          protSC.conFrag1 = tsc->conFrag;
          protSC.conFrag2 = conFrag;
          for (x = 0; x<tsc->modLen; x++) protSC.mods1->push_back(tsc->mods[x]);
          ret=false; //indicates where to put mods to pep2
        }
        
        protSC.mass = tsc->mass + spec->getLink(linkIndex).mass + mass;
        protSC.linkable1 = tsc->linkable;
        protSC.linkable2 = false;
        protSC.link = linkIndex;
        protSC.precursor = (char)i;

        //special case for identical peptides
        //if(protSC.pep2==protSC.pep1 && protSC.k1==protSC.k2){ 
        //  protSC.score1/=2;
        //  protSC.score2/=2;
        //  protSC.simpleScore/=2;
        //  it++;
        //  continue; //WARNING...JUST SKIPPING THESE...
        //}

        //index the mods for pep2
        iset = ions[iIndex].at(sIndex);
        if (iset->difMass != 0){
          for (j = 0; j<ions[iIndex].getIonCount(); j++) {
            if (iset->mods[j] != 0){
              if(j==0){
                if(iset->nTermMass!=0){
                  mod.pos=-1;
                  mod.mass=iset->nTermMass;
                  if (ret) protSC.mods1->push_back(mod);
                  else protSC.mods2->push_back(mod);
                }
                if(fabs(iset->mods[j]-iset->nTermMass)<0.0001) continue;
              }
              if (j == ions[iIndex].getIonCount() - 1){
                if (iset->cTermMass != 0){
                  mod.pos = -2;
                  mod.mass = iset->cTermMass;
                  if (ret) protSC.mods1->push_back(mod);
                  else protSC.mods2->push_back(mod);
                }
                if (fabs(iset->mods[j] - iset->cTermMass)<0.0001) continue;
              }
              //if (j == 0 && iset->modNTerm) mod.term = true;
              //else if (j == ions[iIndex].getIonCount() - 1 && iset->modCTerm) mod.term = true;
              //else mod.term = false;
              mod.pos = (char)j;
              mod.mass = iset->mods[j];
              if(ret) protSC.mods1->push_back(mod);
              else protSC.mods2->push_back(mod);
            }
          }
        }
        Threading::LockMutex(mutexSpecScore[index]);
        s->checkScore(protSC);
        Threading::UnlockMutex(mutexSpecScore[index]);
        it++;
      }
    //Threading::UnlockMutex(mutexSingletScore[index][i]);
    }

    if (firstPass && mass>(p->monoMass - spec->getLink(linkIndex).mass) / 2-0.2){
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
      if (lowI >= highI) continue;
      
      score = kojakScoring(index, p->monoMass - mass, sIndex, iIndex, matches, conFrag, p->charge);
      bScored = true;
      //y = (int)(score * 10.0 + 0.5);
      //if (y >= HISTOSZ) y = HISTOSZ - 1;
      //Threading::LockMutex(mutexSpecScore[index]);
      //s->histogramSinglet[y]++;
      //s->histogramSingletCount++;
      //Threading::UnlockMutex(mutexSpecScore[index]);
      if(score<params.minPepScore || score<=0) continue;
      //if(conFrag<2) continue; //FOR TESTING ONLY

      //boost score...FOR TESTING ONLY
      //if(conFrag>2){
      //  score*=(1.0+(double)conFrag/10);
      //}

      Threading::LockMutex(mutexSingletScore[index][i]);
      tp = s->getTopPeps(i);
      if(tp->singletCount>=tp->singletMax && score<tp->singletLast->simpleScore) {
        Threading::UnlockMutex(mutexSingletScore[index][i]);
        continue; //don't bother with the singlet overhead if it won't make the list
      }
      Threading::UnlockMutex(mutexSingletScore[index][i]);

      sc.len = len;
      sc.simpleScore = score;
      sc.k1 = k;
      sc.linkable = false;
      sc.pep1 = pep;
      sc.pre = i;
      sc.mass = mass;
      sc.site = linkSite;
      sc.matches=matches;
      sc.conFrag=conFrag;
      v.clear();
      if (sc.mods != NULL) {
        sc.modLen = 0;
        delete[] sc.mods;
        sc.mods = NULL;
      }
      iset = ions[iIndex].at(sIndex);
      if (iset->difMass != 0){
        for (j = 0; j<ions[iIndex].getIonCount(); j++) {
          if (iset->mods[j] != 0){
            if (j == 0){
              if (iset->nTermMass != 0){
                mod.pos = -1;
                mod.mass = iset->nTermMass;
                v.push_back(mod);
              }
              if (fabs(iset->mods[j] - iset->nTermMass)<0.0001) continue;
            }
            if (j == ions[iIndex].getIonCount() - 1){
              if (iset->cTermMass != 0){
                mod.pos = -2;
                mod.mass = iset->cTermMass;
                v.push_back(mod);
              }
              if (fabs(iset->mods[j] - iset->cTermMass)<0.0001) continue;
            }
            //if (j == 0 && iset->modNTerm) mod.term = true;
            //else if (j == ions[iIndex].getIonCount() - 1 && iset->modCTerm) mod.term = true;
            //else mod.term = false;
            mod.pos = (char)j;
            mod.mass = iset->mods[j];
            v.push_back(mod);
          }
        }
        sc.modLen = (char)v.size();
        sc.mods = new kPepMod[sc.modLen];
        for (j = 0; j<(int)sc.modLen; j++) sc.mods[j] = v[j];
      }

      Threading::LockMutex(mutexSingletScore[index][i]);
      tp = s->getTopPeps(i);
      tp->checkSingletScore(sc);
      Threading::UnlockMutex(mutexSingletScore[index][i]);

      //bScored=true;
    }
  }

  return bScored;
}

void KAnalysis::scoreSpectra(vector<int>& index, int sIndex, double modMass, int pep1, int pep2, int k1, int k2, int link, int iIndex, char linkSite1, char linkSite2){
  unsigned int a;
  int i,z,ps,y;
  kScoreCard sc;
  kPepMod mod;
  int matches;
  int conFrag;
  double mass = ions[iIndex][sIndex].mass;

  //score spectra
  for(a=0;a<index.size();a++){
    
    //find the specific precursor mass in this spectrum to identify the charge state
    z=0;
    for (ps = 0; ps<spec->at(index[a]).sizePrecursor(); ps++){
      kPrecursor* p = spec->at(index[a]).getPrecursor2(ps);
      double ppm = (p->monoMass - mass) / mass*1e6;
      if (ppm<params.ppmPrecursor && ppm>-params.ppmPrecursor){
        z = p->charge;
        break;
      }
    }
    
    sc.simpleScore=kojakScoring(index[a],modMass,sIndex,iIndex, matches, conFrag, z);
    y = (int)(sc.simpleScore * 10.0 + 0.5);
    Threading::LockMutex(mutexSpecScore[index[a]]);
    spec->at(index[a]).histogram[y]++;
    spec->at(index[a]).histogramCount++;
    Threading::UnlockMutex(mutexSpecScore[index[a]]);
    if(sc.simpleScore<0.1)  continue;

    //maybe do all this only if the score is going to make the list? see singlets above
    sc.mods1->clear();
    sc.mods2->clear();
    sc.score1=sc.simpleScore;
    sc.matches1=matches;
    sc.conFrag1=conFrag;
    sc.k1=k1;
    sc.k2=k2;
    sc.site1=linkSite1; //need to be amino acids
    sc.site2=linkSite2; //need to be amino acids
    sc.mass=ions[iIndex][sIndex].mass;
    sc.linkable1=sc.linkable2=false;
    sc.pep1=pep1;
    sc.pep2=pep2;
    sc.link=link;
    sc.precursor=(char)ps;
    if(ions[iIndex][sIndex].difMass!=0){
      for(i=0;i<ions[iIndex].getPeptideLen();i++) {
        if(ions[iIndex][sIndex].mods[i]!=0){
          if (i == 0){
            if (ions[iIndex][sIndex].nTermMass != 0){
              mod.pos = -1;
              mod.mass = ions[iIndex][sIndex].nTermMass;
              sc.mods1->push_back(mod);
            }
            if (fabs(ions[iIndex][sIndex].mods[i] - ions[iIndex][sIndex].nTermMass)<0.0001) continue;
          }
          if (i == ions[iIndex].getIonCount() - 1){
            if (ions[iIndex][sIndex].cTermMass != 0){
              mod.pos = -2;
              mod.mass = ions[iIndex][sIndex].cTermMass;
              sc.mods1->push_back(mod);
            }
            if (fabs(ions[iIndex][sIndex].mods[i] - ions[iIndex][sIndex].cTermMass)<0.0001) continue;
          }
          //if (i == 0 && ions[iIndex][sIndex].modNTerm) mod.term = true;
          //else if (i == ions[iIndex].getIonCount() - 1 && ions[iIndex][sIndex].modCTerm) mod.term = true;
          //else mod.term = false;
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
float KAnalysis::kojakScoring(int specIndex, double modMass, int sIndex, int iIndex, int& match, int& conFrag, int z) { 

  KSpectrum* s=spec->getSpectrum(specIndex);
  KIonSet* ki=ions[iIndex].at(sIndex);
  if (!ki->index) {
    ki->makeIndex(params.binSize, params.binOffset, params.ionSeries[0], params.ionSeries[1], params.ionSeries[2], params.ionSeries[3], params.ionSeries[4], params.ionSeries[5]);
  }

  double dXcorr=0.0;
  double invBinSize=s->getInvBinSize();
  double binOffset=params.binOffset;
  double dif;
  double mz;

  int ionCount=ions[iIndex].getIonCount();
  int maxCharge=z;
  if(maxCharge<1) maxCharge=s->getCharge();  

  int i,j,k;
  int key;
  int pos;
  int con;
  match=0;
  conFrag=0;

  //Assign ion series
  kISValue***  ionSeries;
  ionSeries=new kISValue**[numIonSeries];
  k=0;
  if(params.ionSeries[0]) ionSeries[k++]=ki->aIons;
  if(params.ionSeries[1]) ionSeries[k++]=ki->bIons;
  if(params.ionSeries[2]) ionSeries[k++]=ki->cIons;
  if(params.ionSeries[3]) ionSeries[k++]=ki->xIons;
  if(params.ionSeries[4]) ionSeries[k++]=ki->yIons;
  if(params.ionSeries[5]) ionSeries[k++]=ki->zIons;

  //The number of fragment ion series to analyze is PrecursorCharge-1
  //However, don't analyze past the 3+ series
  if(maxCharge>4) maxCharge=4;

  //Iterate all series
  for(k=1;k<maxCharge;k++){

    dif=modMass/k;

    //Iterate through pfFastXcorrData
    for(j=0;j<numIonSeries;j++){
      con=0;

      for(i=0;i<ionCount;i++){

        //get key -- see if this can be precomputed for the half that doesn't contain the linked peptide
        if(ionSeries[j][k][i].mz<0) {
          mz = params.binSize * (int)((dif-ionSeries[j][k][i].mz)*invBinSize+binOffset);
          key = (int)mz;
          if(key>=s->kojakBins) {
            if (con>conFrag) conFrag = con;
            con = 0;
            break;
          }
          if(s->kojakSparseArray[key]==NULL) {
            if (con>conFrag) conFrag = con;
            con = 0;
            continue;
          }
          pos = (int)((mz-key)*invBinSize);
          dXcorr += s->kojakSparseArray[key][pos];
          if (s->kojakSparseArray[key][pos]>5) {
            match++;
            con++;
          } else {
            if(con>conFrag) conFrag=con;
            con=0;
          }
        } else {
          key=ionSeries[j][k][i].key;
          if (key >= s->kojakBins) {
            if (con>conFrag) conFrag = con;
            con = 0;
            break;
          }
          if (s->kojakSparseArray[key] == NULL) {
            if (con>conFrag) conFrag = con;
            con = 0;
            continue;
          }
          pos = ionSeries[j][k][i].pos;
          dXcorr += s->kojakSparseArray[key][pos];
          if (s->kojakSparseArray[key][pos]>5) {
            match++;
            con++;
          } else {
            if (con>conFrag) conFrag = con;
            con = 0;
          }
        }
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

/* kruft?
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
*/

void KAnalysis::setLog(KLog* c){
  klog=c;
}

//This function determines if a site on a peptide is linkable to another peptide in the database.
//Sites may have multiple partners, such as in dual-linker searches with K-K and K-D/E
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

      //skip peptides that cannot be linked
      if (p->at(i).xlSites==0) continue;

      //determine which motifs are allowed for this peptide
      //check termini first
      bMatch=false;
      if(p->at(i).nTerm){
        for (m = 0; m < 20; m++){
          if (xlTable['n'][m] == -1) break;
          if (xlTable['n'][m] == (char)j) {
            bMatch = true;
            break;
          }
        }
      }
      if(!bMatch && p->at(i).cTerm){
        for (m = 0; m < 20; m++){
          if (xlTable['c'][m] == -1) break;
          if (xlTable['c'][m] == (char)j) {
            bMatch = true;
            break;
          }
        }
      }
      if(!bMatch){
        db->getPeptideSeq(p->at(i),pep);
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
      ions[0].modIonsRec(0, -1, 0, 0, false); //does this need to be modIonsRec2???
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
