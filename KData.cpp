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

#include "KData.h"
#include "Profiler.h"

using namespace std;
using namespace MSToolkit;

extern Profiler prof;

Mutex KData::mutexMemoryPool;
bool* KData::memoryPool;
double** KData::tempRawData;
double** KData::tmpFastXcorrData;
float**  KData::fastXcorrData;
kPreprocessStruct** KData::preProcess;
kParams* KData::params;

/*============================
  Constructors
============================*/
KData::KData(){
  int i,j,k;
  bScans=NULL;
  params=NULL;
  klog=NULL;
  xlTable = new char*[128];
  for (i = 0; i<128; i++) xlTable[i] = new char[20];
  for(i=0;i<128;i++){
    for(j=0;j<5;j++){
      for(k=0;k<128;k++){
        xlTargets[i][j].target[k]=false;
      }
    }
  }
}

KData::KData(kParams* p){
  bScans=NULL;
  klog=NULL;
  params=p;
  size_t i;
  int j,k;
  for(i=0;i<p->fMods->size();i++) aa.addFixedMod((char)p->fMods->at(i).index,p->fMods->at(i).mass);
  xlTable = new char*[128];
  for (j = 0; j<128; j++) xlTable[j] = new char[20];
  for (i = 0; i<128; i++){
    for (j = 0; j<5; j++){
      for (k = 0; k<128; k++){
        xlTargets[i][j].target[k] = false;
      }
    }
  }
}

KData::~KData(){
  params=NULL;
  klog=NULL;
  if(bScans!=NULL) delete[] bScans;
  for (int i = 0; i<128; i++) delete[] xlTable[i];
  delete[] xlTable;
}


/*============================
  Operators
============================*/
KSpectrum& KData::operator [](const int& i){
  return spec[i];
}

//============================
//  Thread-Start Functions
//============================

//These functions fire off when a thread starts. They pass the variables to for
//each thread-specific analysis to the appropriate function.
void KData::xCorrProc(kSpectrumStruct* s){
  int i;
  Threading::LockMutex(mutexMemoryPool);
  for (i = 0; i<params->threads; i++){
    if (!memoryPool[i]){
      memoryPool[i] = true;
      break;
    }
  }
  Threading::UnlockMutex(mutexMemoryPool);
  if (i == params->threads){
    cout << "Error in KData::xCorrProc" << endl;
    exit(-1);
  }
  s->mem=&memoryPool[i];
  s->spec->kojakXCorr(tempRawData[i],tmpFastXcorrData[i],fastXcorrData[i],preProcess[i]);
  delete s;
  //s->xCorrScore(false); //this cruft should be eliminated...xCorrScore() should always use kojak version
  s = NULL;
}


/*============================
  Functions
============================*/
void KData::addProteins(void* sh, KDatabase& db, int pIndex, bool xl, int linkA, int linkB) {
  char c;
  char n;

  kPeptide pep = db.getPeptide(pIndex);
  if(!xl) (*(CnpxSearchHit*)sh).num_tot_proteins = (int)pep.map->size();
  else (*(CnpxLinkedPeptide*)sh).num_tot_proteins = (int)pep.map->size();

  for (size_t j = 0; j<pep.map->size(); j++){
    if (pep.n15 && db[pep.map->at(j).index].name.find(params->n15Label) == string::npos) {
      if (!xl) (*(CnpxSearchHit*)sh).num_tot_proteins--;
      else (*(CnpxLinkedPeptide*)sh).num_tot_proteins--;
      continue;
    }
    if (!pep.n15 && strlen(params->n15Label)>0 && db[pep.map->at(j).index].name.find(params->n15Label) != string::npos) {
      if (!xl) (*(CnpxSearchHit*)sh).num_tot_proteins--;
      else (*(CnpxLinkedPeptide*)sh).num_tot_proteins--;
      continue;
    }
    string protein = "";
    for (size_t i = 0; i<db[pep.map->at(j).index].name.size(); i++){
      if (params->truncate>0 && i == params->truncate) break;
      protein += db[pep.map->at(j).index].name[i];
    }
    
    if (pep.map->at(j).start<1) n = '-';
    else n = db[pep.map->at(j).index].sequence[pep.map->at(j).start - 1];
    if (pep.map->at(j).stop + 1 == db[pep.map->at(j).index].sequence.size()) c = '-';
    else c = db[pep.map->at(j).index].sequence[pep.map->at(j).stop + 1];

    if(j==0){
      if (!xl) {
        (*(CnpxSearchHit*)sh).protein=protein;
        (*(CnpxSearchHit*)sh).peptide_prev_aa=n;
        (*(CnpxSearchHit*)sh).peptide_next_aa=c;
        (*(CnpxSearchHit*)sh).peptide_start_pos = (int)pep.map->at(j).start + 1;
        if (linkA>0) (*(CnpxSearchHit*)sh).protein_link_pos_a = pep.map->at(j).start + linkA;
        if (linkB>0) (*(CnpxSearchHit*)sh).protein_link_pos_b = pep.map->at(j).start + linkB;

      } else {
        (*(CnpxLinkedPeptide*)sh).protein=protein;
        (*(CnpxLinkedPeptide*)sh).peptide_next_aa=c;
        (*(CnpxLinkedPeptide*)sh).peptide_prev_aa=n;
        (*(CnpxLinkedPeptide*)sh).peptide_start_pos = (int)pep.map->at(j).start + 1;
        (*(CnpxLinkedPeptide*)sh).protein_link_pos_a = pep.map->at(j).start + linkA;

      }
    } else {
      CnpxAlternativeProtein ap;
      ap.protein = protein;
      ap.peptide_prev_aa = n;
      ap.peptide_next_aa = c;
      ap.peptide_start_pos = (int)pep.map->at(j).start + 1;
      if (linkA>0) ap.protein_link_pos_a = pep.map->at(j).start + linkA;
      if (linkB>0) ap.protein_link_pos_b = pep.map->at(j).start + linkB;
      if (!xl) (*(CnpxSearchHit*)sh).alternative_protein.push_back(ap);
      else (*(CnpxLinkedPeptide*)sh).alternative_protein.push_back(ap);
    }
  }

}

void KData::addSearchScore(CnpxSearchHit& sh, string name, double value, string fmt){
  char score[32];
  sprintf(score, fmt.c_str(), value);
  CnpxSearchScore ss;
  ss.name = name;
  ss.value = score;
  sh.search_score.push_back(ss);
}

void KData::addSearchScore(CnpxSearchHit& sh, string name, int value){
  char score[32];
  sprintf(score, "%d", value);
  CnpxSearchScore ss;
  ss.name = name;
  ss.value = score;
  sh.search_score.push_back(ss);
}

void KData::addXlinkScore(CnpxLinkedPeptide& lp, string name, double value, string fmt){
  char score[32];
  sprintf(score, fmt.c_str(), value);
  CnpxXLinkScore xls;
  xls.name=name;
  xls.value=score;
  lp.xlink_score.push_back(xls);
}

void KData::addXlinkScore(CnpxLinkedPeptide& lp, string name, int value){
  char score[32];
  sprintf(score, "%d", value);
  CnpxXLinkScore xls;
  xls.name = name;
  xls.value = score;
  lp.xlink_score.push_back(xls);
}

void KData::addXlinkScore(CnpxXLink& xl, string name, double value, string fmt){
  char score[32];
  sprintf(score, fmt.c_str(), value);
  CnpxXLinkScore xls;
  xls.name = name;
  xls.value = score;
  xl.xlink_score.push_back(xls);
}

void KData::addXlinkScore(CnpxXLink& xl, string name, int value){
  char score[32];
  sprintf(score, "%d", value);
  CnpxXLinkScore xls;
  xls.name = name;
  xls.value = score;
  xl.xlink_score.push_back(xls);
}

void KData::buildXLTable(){
  int i, j;
  int xlA, xlB;
  size_t k;

  //initialize/clear values
  motifCount=0;
  for (i = 0; i < 20; i++){
    motifs[i].motif.clear();
    for (int j = 0; j<10; j++){
      motifs[i].xlIndex[j] = -1;
      motifs[i].counterMotif[j] = -1;
    }
  }
  for (i = 0; i < 128; i++){
    for (j = 0; j<20; j++) xlTable[i][j] = -1;
  }

  //Iterate over all cross-linkers
  for (k = 0; k<link.size(); k++){

    //Add motifA, if needed
    for (i = 0; i < motifCount; i++){
      if (motifs[i].motif.compare(link[k].motifA) == 0) break;
    }
    if (i == motifCount) {
      motifs[i].motif = link[k].motifA;
      motifCount++;
    }
    link[k].motifAIndex = i;
    //set index to cross-linker
    for (j = 0; j < 10; j++){
      if (motifs[link[k].motifAIndex].xlIndex[j] == (int)k) break;
      if (motifs[link[k].motifAIndex].xlIndex[j] == -1) {
        motifs[link[k].motifAIndex].xlIndex[j] = (int)k;    
        break;
      }
    }
    xlA = j;

    //Add motifB, if needed
    for (i = 0; i < motifCount; i++){
      if (motifs[i].motif.compare(link[k].motifB) == 0) break;
    }
    if (i == motifCount) { //Add new motif
      motifs[i].motif = link[k].motifB;
      motifCount++;
    }
    link[k].motifBIndex = i;
    //set index
    for (j = 0; j < 10; j++){
      if (motifs[link[k].motifBIndex].xlIndex[j] == (int)k) break;
      if (motifs[link[k].motifBIndex].xlIndex[j] == -1) {
        motifs[link[k].motifBIndex].xlIndex[j] = (int)k;    
        break;
      }
    }
    xlB = j;

    //set counter motifs
    motifs[link[k].motifAIndex].counterMotif[xlA] = link[k].motifBIndex;
    motifs[link[k].motifBIndex].counterMotif[xlB] = link[k].motifAIndex;
  }

  //Build table of linked amino acids
  for (i = 0; i < motifCount; i++){

    //This lookup table is for fast referencing partner link sites when given a site
    for(j=0;j<5;j++){
      if(motifs[i].xlIndex[j]<0) break;
      for (k = 0; k < motifs[motifs[i].counterMotif[j]].motif.size(); k++){
        for(size_t n=0;n<motifs[i].motif.size();n++){
          xlTargets[motifs[i].motif[n]][motifs[i].xlIndex[j]].target[motifs[motifs[i].counterMotif[j]].motif[k]] = true;
        }
      }
    }


    for (k = 0; k < motifs[i].motif.size(); k++){
      for (j = 0; j < 20; j++){
        if (xlTable[motifs[i].motif[k]][j] == -1){
          xlTable[motifs[i].motif[k]][j] = i;
          break;
        }
      }
    }
  }

  //Output for diagnostic purposes
  //for (i = 0; i < motifCount; i++){
  //  cout << "Motif " << i << ": " << &motifs[i].motif[0] << endl;
  //  for(k=0;k<10;k++){
  //    if(motifs[i].counterMotif[k]<0) break;
  //    cout << "  XLink index: " << (int)motifs[i].xlIndex[k] << "\tCounter Motif: " << (int)motifs[i].counterMotif[k] << endl;
  //  }
  //}
  //bool bShow;
  //char cStr[256];
  //string str;
  //for(i=0;i<128;i++){
  //  for(j=0;j<link.size();j++){
  //    bShow=false;
  //    sprintf(cStr,"%c: %s",(char)i,link[j].label);
  //    str=cStr;
  //    for(k=0;k<128;k++){
  //      if(xlTargets[i][j].target[k]) {
  //        str+=" ";
  //        str+=(char)k;
  //        bShow=true;
  //      }
  //    }
  //    if(bShow) cout << str << endl;
  //  }
  //}
  
}

bool KData::checkLink(char p1Site, char p2Site, int linkIndex){
  return xlTargets[p1Site][linkIndex].target[p2Site];
}

bool KData::doXCorr(kParams& params){
  klog->addMessage("Using Kojak modified XCorr scores.", true);
  cout << "  Using Kojak modified XCorr scores." << endl;
  klog->addMessage("Transforming spectra.", true);
  cout << "  Transforming spectra ... ";

  int i;
  int iPercent;
  int iTmp;
  size_t szBuffer=2000;
  size_t index=0;

  memoryAllocate();

  ThreadPool<kSpectrumStruct*>* threadPool = new ThreadPool<kSpectrumStruct*>(xCorrProc, params.threads, params.threads, 1);

  //Set progress meter
  iPercent = 0;
  printf("%2d%%", iPercent);
  fflush(stdout);

  //Iterate the spectra for the first pass
  for (i = 0; i<spec.size(); i++){

    threadPool->WaitForQueuedParams();

    kSpectrumStruct* a=new kSpectrumStruct(&mutexMemoryPool,&spec[i]);
    threadPool->Launch(a);

    //Update progress meter
    iTmp = (int)((double)i / spec.size() * 100);
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
  memoryFree();

  return true;
}

KSpectrum* KData::getSpectrum(const int& i){
  return &spec[i];
}

char** KData::getXLTable(){
  return xlTable;
}

KSpectrum& KData::at(const int& i){
  return spec[i];
}

void KData::diagSinglet(){
  int oddCount = 0;
  double maxScore;
  int bigCount = 0;
  int twoCount = 0;
  double bigScore = 0;
  double bigMass;
  double unknownMass;
  for (size_t b = 0; b<spec.size(); b++){
    //if(spec[b].getScoreCard(0).simpleScore==0) continue;
    maxScore = 0;
    unknownMass = 0;
    for (int q = 0; q<spec[b].sizePrecursor(); q++){
      if (spec[b].getTopPeps(q)->singletCount == 0) continue;
      if (spec[b].getTopPeps(q)->singletFirst->simpleScore>spec[b].getScoreCard(0).simpleScore){
        if (spec[b].getTopPeps(q)->singletFirst->simpleScore>maxScore) {
          maxScore = spec[b].getTopPeps(q)->singletFirst->simpleScore;
          unknownMass = spec[b].getPrecursor(q).monoMass - spec[b].getTopPeps(q)->singletFirst->mass;
        }
      }
    }
    if (maxScore>0) oddCount++;
    if (maxScore>3) bigCount++;
    if (maxScore>3 && spec[b].getScoreCard(0).simpleScore>0 && maxScore>spec[b].getScoreCard(0).simpleScore * 2) twoCount++;
    if (maxScore>bigScore) {
      bigScore = maxScore;
      bigMass = unknownMass;
    }
  }
  cout << " Diagnostics:" << endl;
  cout << "  Odd counts " << oddCount << " of " << spec.size() << endl;
  cout << "  Big counts " << bigCount << " of " << oddCount << endl;
  cout << "  Two fold " << twoCount << " of " << bigCount << endl;
  cout << "  Biggest Score: " << bigScore << "  Mass: " << bigMass << endl;
  cout << endl;
}

bool KData::getBoundaries(double mass1, double mass2, vector<int>& index, bool* buffer){
  int sz=(int)massList.size();

  if(mass1>massList[sz-1].mass) return false;

  int lower=0;
  int mid=sz/2;
  int upper=sz;
  int i;
  int low;
  int high;

  vector<int> v;

  //binary search to closest mass
  while(massList[mid].mass!=mass1){
		if(lower>=upper) break;
    if(mass1<massList[mid].mass){
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

  //Adjust if needed
  if(massList[mid].mass<mass1) mid++;
  if(mid>0 && massList[mid-1].mass>mass1) mid--;
  if(massList[mid].mass>mass2) return false;
  if(mid==sz) return false;
  low=mid;

  //binary search to next mass
  lower=0;
  mid=sz/2;
  upper=sz;
  while(massList[mid].mass!=mass2){
		if(lower>=upper) break;
    if(mass2<massList[mid].mass){
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

  //Adjust if needed
  if(massList[mid].mass>mass2) mid--;
  if(mid<sz-1 && massList[mid+1].mass<mass2) mid++;
  if(mid<0) return false;
  high=mid;

  sz=(int)spec.size();
  memset(buffer,false,sz);
  for (i = low; i <= high; i++) buffer[massList[i].index]=true;
  index.clear();
  for (i = 0; i < sz; i++) {
    if(buffer[i]) index.push_back(i);
  }
  return true;

  /* old method
  for(i=low;i<=high;i++) v.push_back(massList[i].index);

  //Sort indexes and copy to final array, removing duplicates.
  //This may be a potentially slow step and should be profiled
  qsort(&v[0],v.size(),sizeof(int),compareInt);
  index.clear();
  index.push_back(v[0]);
  mid=(int)v.size();
  for(i=1;i<mid;i++){
    if(v[i]!=v[i-1]) index.push_back(v[i]);
  }
	return true;
  */

}

//Get the list of spectrum array indexes to search based on desired mass
bool KData::getBoundaries2(double mass, double prec, vector<int>& index, bool* buffer){
  int sz=(int)massList.size();
  int lower=0;
  int mid=sz/2;
  int upper=sz;
	int i;

  vector<int> v;

  double minMass = mass - (mass/1000000*prec);
  double maxMass = mass + (mass/1000000*prec);

  //binary search to closest mass
  while(massList[mid].mass<minMass || massList[mid].mass>maxMass){
		if(lower>=upper) break;
    if(mass<massList[mid].mass){
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

	//Check that mass is correct
	if(massList[mid].mass<minMass || massList[mid].mass>maxMass) return false;

  v.push_back(massList[mid].index);

	//check left 
  i=mid;
	while(i>0){
		i--;
    if(massList[i].mass<minMass) break;
    v.push_back(massList[i].index);
	}

	//check right
	i=mid;
	while(i<(sz-1)){
		i++;
    if(massList[i].mass>maxMass) break;
    v.push_back(massList[i].index);
	}

  sz = (int)spec.size();
  memset(buffer, false, sz);
  for (i = 0; i < (int)v.size(); i++) buffer[v[i]] = true;
  index.clear();
  for (i = 0; i < sz; i++) {
    if (buffer[i]) index.push_back(i);
  }
  return true;

  //Sort indexes and copy to final array, removing duplicates.
  //This may be a potentially slow step and should be profiled
  qsort(&v[0],v.size(),sizeof(int),compareInt);
  index.clear();
  index.push_back(v[0]);
  mid=(int)v.size();
  for(i=1;i<mid;i++){
    if(v[i]!=v[i-1]) index.push_back(v[i]);
  }

	return true;

}

int KData::getCounterMotif(int motifIndex, int counterIndex){
  if (motifIndex>=motifCount) return -1;
  return motifs[motifIndex].counterMotif[counterIndex];
}

kLinker& KData::getLink(int i){
  return link[i];
}

double KData::getMaxMass(){
  if(massList.size()==0) return 0;
  else return massList[massList.size()-1].mass;
}

double KData::getMinMass(){
  if(massList.size()==0) return 0;
  else return massList[0].mass;
}

int KData::getMotifCount(){
  return motifCount;
}

int KData::getXLIndex(int motifIndex, int xlIndex){
  if (motifIndex >= motifCount) return -1;
  return motifs[motifIndex].xlIndex[xlIndex];
}

CnpxModificationInfo KData::makeModificationInfo(vector<kPepMod>& mods, string peptide, bool n15, bool nTerm, bool cTerm){
  CnpxModificationInfo mi;
  for (size_t i = 0; i<mods.size(); i++){
    if (mods[i].pos == -1) mi.mod_nterm_mass = mods[i].mass;
    else if (mods[i].pos == -2) mi.mod_cterm_mass = mods[i].mass;
    else {
      CnpxModAminoAcidMass maam;
      maam.position=(int)mods[i].pos + 1;
      maam.mass=mods[i].mass + aa.getAAMass(peptide[mods[i].pos], n15);
      maam.variable=mods[i].mass;
      mi.mod_aminoacid_mass.push_back(maam);
    }
  }
  if (nTerm && aa.getFixedModMass('$') != 0)mi.mod_nterm_mass += aa.getFixedModMass('$');
  if (cTerm && aa.getFixedModMass('%') != 0)mi.mod_cterm_mass += aa.getFixedModMass('%');
  mi.mod_nterm_mass += aa.getFixedModMass('n');
  mi.mod_cterm_mass += aa.getFixedModMass('c');
  if (mi.mod_nterm_mass != 0) mi.mod_nterm_mass += 1.00782503;
  if (mi.mod_cterm_mass != 0) mi.mod_cterm_mass += 17.00273963;
  for (size_t i = 0; i<peptide.size(); i++){
    if (aa.getFixedModMass(peptide[i])>0) {
      CnpxModAminoAcidMass maam;
      maam.position=(int)i + 1;
      maam.mass=aa.getAAMass(peptide[i], n15);
      maam.staticMass=aa.getFixedModMass(peptide[i]);
      mi.mod_aminoacid_mass.push_back(maam);
    }
  }
  return mi;
}

//This function tries to assign best possible 18O2 and 18O4 precursor ion mass values
//for all MS2 spectra
bool KData::mapPrecursors(){
  
  int iPercent=0;
  int iTmp;
  
  size_t i;
  int j,k,n;

  KPrecursor pre(params);
  kMass      m;

  int peakCounts=0;
  int specCounts=0;
  int ret;

  int prePre=0;
  int foundPre=0;
  int noPre=0;

  //Open the data file in the precursor mapping object
  //if(!pre.setFile(&p)) return false;

  //Print progress
  if(klog!=NULL) {
    klog->addMessage("Mapping precursors to MS/MS spectra",true);
    pre.setLog(klog);
  }
  printf("  Mapping precursors ... %2d%%",iPercent);
  fflush(stdout);

  //Iterate all MS/MS spectra
  for(i=0;i<spec.size();i++){

    //Update progress
    iTmp=(int)(i*100.0/spec.size());
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }

    bool bAddHardklor=false;
    bool bAddEstimate=false;

    if(params->preferPrecursor==1){
      if(spec[i].sizePrecursor()==0){
        if(spec[i].getCharge()>0) bAddEstimate=true;
        else bAddHardklor=true;
      }
    } else if(params->preferPrecursor==0){
      spec[i].clearPrecursors();
      bAddHardklor=true;
      if (spec[i].getCharge()) bAddEstimate = true;
    } else {
      bAddHardklor=true;
      if(spec[i].getCharge()) bAddEstimate=true;
    }

    if(bAddHardklor){
      //only do Hardklor analysis if data contain precursor scans
      if (params->precursorRefinement) ret = pre.getSpecRange(spec[i]);
    }

    if(bAddEstimate){
      kPrecursor pr;
      pr.monoMass = spec[i].getMZ()*spec[i].getCharge() - 1.007276466*spec[i].getCharge();
      pr.charge = spec[i].getCharge();
      pr.corr = -5;
      spec[i].setCharge(pr.charge);
      spec[i].addPrecursor(pr, params->topCount);
      for (int px = 1; px <= params->isotopeError; px++){
        if (px == 4) break;
        pr.monoMass -= 1.00335483;
        pr.corr -= 0.1;
        spec[i].addPrecursor(pr, params->topCount);
      }
    }

    //Now clean up any duplicate precursors. They should already be in order of priority. Use 5ppm as tolerance
    for(k=0;k<spec[i].sizePrecursor()-1;k++){
      for(n=k+1;n<spec[i].sizePrecursor();n++){
        double m1=spec[i].getPrecursor(k).monoMass;
        double m2=spec[i].getPrecursor(n).monoMass;
        double m=(m1-m2)/m1*1e6;
        if(fabs(m)<5){
          spec[i].erasePrecursor(n);
          n--;
        }
      }
    }

    if(spec[i].sizePrecursor()>0){
      foundPre++;
      specCounts++;
      peakCounts+=spec[i].size();

      //build singletList
      spec[i].resetSingletList();

    }

  }
 

  //Finalize the progress
  printf("\b\b\b100%%");
  cout << endl;

  cout << "  " << specCounts << " spectra with " << peakCounts << " peaks will be analyzed." << endl;
  if (klog != NULL) {
    char tempStr[256];
    sprintf(tempStr,"%d spectra with %d peaks will be analyzed.",specCounts,peakCounts);
    klog->addMessage(string(tempStr),true);
  }

  //Build mass list - this orders all precursor masses, with an index pointing to the actual
  //array position for the spectrum. This is because all spectra will have more than 1
  //precursor mass
  massList.clear();
  for(i=0;i<spec.size();i++){
    m.index=(int)i;
    for(j=0;j<spec[i].sizePrecursor();j++){
      m.mass=spec[i].getPrecursor(j).monoMass;
      massList.push_back(m);
    }
  }

  //sort mass list from low to high
  qsort(&massList[0],massList.size(),sizeof(kMass),compareMassList);

  if(bScans!=NULL) delete[] bScans;
  bScans = new bool[spec.size()];

  return true;
}

void KData::memoryAllocate(){
  //find largest possible array for a spectrum
  int threads=params->threads;
  double xlMass=0;
  for(size_t a=0;a<params->xLink->size();a++){
    if(params->xLink->at(a).mass>xlMass) xlMass=params->xLink->at(a).mass;
  }
  int xCorrArraySize = (int)((params->maxPepMass*2+xlMass + 100.0) / params->binSize);

  //Mark all arrays as available
  memoryPool = new bool[threads];
  for(int a=0;a<threads;a++) memoryPool[a]=false;

  //Allocate arrays
  tempRawData = new double*[threads]();
  for(int a=0;a<threads;a++) tempRawData[a]=new double[xCorrArraySize]();
  
  tmpFastXcorrData = new double*[threads]();
  for(int a=0;a<threads;a++) tmpFastXcorrData[a]=new double[xCorrArraySize]();

  fastXcorrData = new float*[threads]();
  for (int a = 0; a<threads; a++) fastXcorrData[a] = new float[xCorrArraySize]();

  preProcess = new kPreprocessStruct*[threads]();
  for (int a = 0; a<threads; a++) {
    preProcess[a] = new kPreprocessStruct();
    preProcess[a]->pdCorrelationData = new kSpecPoint[xCorrArraySize]();
  }

  //Create mutex
  Threading::CreateMutex(&mutexMemoryPool);
}

void KData::memoryFree(){
  delete [] memoryPool;
  for(int a=0;a<params->threads;a++){
    delete [] tempRawData[a];
    delete [] tmpFastXcorrData[a];
    delete [] fastXcorrData[a];
    delete [] preProcess[a]->pdCorrelationData;
    delete preProcess[a];
  }
  delete [] tempRawData;
  delete [] tmpFastXcorrData;
  delete [] fastXcorrData;
  delete [] preProcess;

  //Destroy mutexes
  Threading::DestroyMutex(mutexMemoryPool);
}

void KData::outputDiagnostics(FILE* f, KSpectrum& s, KDatabase& db){
  size_t i,x;
  int j,k;
  int code;
  char strs[256];
  char st[32];
  string pep1,pep2,tmp;
  kPeptide pep;
  kPrecursor* p;
  KTopPeps* tp;
  kSingletScoreCard* sc;
  kScoreCard psm;
  
  fprintf(f, " <scan id=\"%d\">\n", s.getScanNumber());
  fprintf(f, "  <precursor_list size=\"%d\">\n", s.sizePrecursor());
  
  for (j = 0; j<s.sizePrecursor(); j++) {
    p=s.getPrecursor2(j);
    if(p->corr<-4) code=2;
    else if (p->corr<0)code=3;
    else if(p->corr==0)code=2;
    else code=1;
    fprintf(f, "   <precursor mass=\"%.4lf\" charge=\"%d\" type=\"%d\" hk_corr=\"%.4lf\">\n", p->monoMass, p->charge, code, p->corr);

    tp = s.getTopPeps(j);
    sc = tp->singletFirst;
    k=1;
    while (sc != NULL){
      fprintf(f,"    <peptide rank=\"%d\" sequence=\"",k++);
      db.getPeptideSeq(db.getPeptideList()->at(sc->pep1).map->at(0).index, db.getPeptideList()->at(sc->pep1).map->at(0).start, db.getPeptideList()->at(sc->pep1).map->at(0).stop, strs);
      for (i = 0; i<strlen(strs); i++){
        fprintf(f, "%c", strs[i]);
        for (x = 0; x<sc->modLen; x++){
          if (sc->mods[x].pos == char(i)) fprintf(f, "[%.2lf]", sc->mods[x].mass);
        }
        if (char(i) == sc->k1) fprintf(f, "[x]");
      }
      fprintf(f, "\" link_site=\"%d\" score=\"%.4lf\" matches=\"%d\" longest_run=\"%d\" mass=\"%.4lf\"/>\n", (int)sc->k1+1, sc->simpleScore, sc->matches, sc->conFrag, sc->mass);
      sc = sc->next;
    }
    fprintf(f,"   </precursor>\n");
  }
  fprintf(f,"  </precursor_list>\n");

  k=0;
  for(j=0;j<20;j++){
    if(s.getScoreCard(j).simpleScore>0) k++;
    else break;
  }
  fprintf(f,"  <results_list size=\"%d\">\n",k);
  for (j = 0; j<k; j++){
    fprintf(f,"   <result rank=\"%d\" ",j+1);
    psm = s.getScoreCard(j);
    pep = db.getPeptide(psm.pep1);
    db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, strs);
    pep1.clear();
    if (pep.nTerm && aa.getFixedModMass('$') != 0) {
      sprintf(st, "[%.2lf]", aa.getFixedModMass('$'));
      pep1+=st;
    }
    for (i = 0; i<strlen(strs); i++){
      pep1+=strs[i];
      for (x = 0; x<psm.mods1->size(); x++){
        if (psm.mods1->at(x).pos == (char)i) {
          sprintf(st, "[%.2lf]", psm.mods1->at(x).mass);
          pep1+=st;
        }
      }
      if ((int)i == psm.k1 || (psm.pep2<0 && (int)i == psm.k2)) pep1+="[x]";
    }
    if (pep.cTerm && aa.getFixedModMass('%') != 0) {
      sprintf(st, "[%.2lf]", aa.getFixedModMass('%'));
      pep1+=st;
    }

    pep2.clear();
    if (psm.pep2>-1){
      pep = db.getPeptide(psm.pep2);
      db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, strs);
      if (pep.nTerm && aa.getFixedModMass('$') != 0) {
        sprintf(st, "[%.2lf]", aa.getFixedModMass('$'));
        pep2+=st;
      }
      for (i = 0; i<strlen(strs); i++){
        pep2+=strs[i];
        for (x = 0; x<psm.mods2->size(); x++){
          if (psm.mods2->at(x).pos == (char)i) {
            sprintf(st, "[%.2lf]", psm.mods2->at(x).mass);
            pep2+=st;
          }
        }
        if ((int)i == psm.k2) pep2+="[x]";
      }
      if (pep.cTerm && aa.getFixedModMass('%') != 0) {
        sprintf(st, "[%.2lf]", aa.getFixedModMass('%'));
        pep2+=st;
      }
    } else {
      //fprintf(f2, "\t(%d)", spec[i].getScoreCard(j).k2);
    }
    
    if(pep2.size()>0 && pep1.compare(pep2)>0){
      tmp=pep1;
      pep1=pep2;
      pep2=tmp;
    }

    if (psm.link == -2) pep1+="+";
    else if(pep2.size()>0) pep1+="--";
    pep1+=pep2;
    fprintf(f, "sequence=\"%s\" score=\"%.4lf\" evalue=\"%.3e\" mass=\"%.4lf\"",&pep1[0], psm.simpleScore, psm.eVal, psm.mass);
    if (psm.link<0) fprintf(f, " crosslinker_mass=\"0\"/>\n"); 
    else fprintf(f, " crosslinker_mass=\"%.4lf\"/>\n", link[psm.link].mass);
  }
  fprintf(f, "  </results_list>\n");

  /*
  fprintf(f, "  <histogramSinglet count=\"%d\" intercept=\"%.4f\" slope=\"%.4f\" rsq=\"%.4lf\" start=\"%.4f\" next=\"%.4f\" max=\"%d\">\n", s.histogramSingletCount, s.tmpSingletIntercept, s.tmpSingletSlope, s.tmpSingletRSquare, s.tmpSingletIStartCorr, s.tmpSingletINextCorr, s.tmpSingletIMaxCorr);
  for (j = 0; j<HISTOSZ; j++){
    fprintf(f, "   <bin id=\"%d\" value=\"%d\"/>\n", j, s.histogramSinglet[j]);
  }
  fprintf(f, "  </histogramSinglet>\n");
  */

  fprintf(f, "  <histogram count=\"%d\" intercept=\"%.4f\" slope=\"%.4f\" rsq=\"%.4lf\" start=\"%.4f\" next=\"%.4f\" max=\"%d\">\n", s.histogramCount,s.tmpIntercept,s.tmpSlope,s.tmpRSquare,s.tmpIStartCorr,s.tmpINextCorr,s.tmpIMaxCorr);
  for (j = 0; j<HISTOSZ; j++){
    fprintf(f, "   <bin id=\"%d\" value=\"%d\" score=\"%.1lf\" count=\"%d\"/>\n", j, s.histogram[j],(double)j/10,s.histogramO[j]);
  }
  fprintf(f, "  </histogram>\n");
  

  fprintf(f, " </scan>\n");
}

//Function deprecated. Should be excised.
bool KData::outputIntermediate(KDatabase& db){

  size_t i,n;
  int j, k, x,z;
  char fName[1056];

  FILE* fOut = NULL;

  KTopPeps* tp;

  kPeptide pep;
  kSingletScoreCard* sc;
  char strs[256];
  char strTmp[32];
  string pepSeq;
  string protSeq;

  //Open all the required output files.
  sprintf(fName, "%s.intermediate.xml", params->outFile);
  fOut = fopen(fName, "wt");
  if (fOut == NULL) return false;

  for (i = 0; i<spec.size(); i++) {

    fprintf(fOut, "<spectrum scan=\"%d\" retention_time_sec=\"%.4f\" selected_mz=\"%.8lf\">\n", spec[i].getScanNumber(), spec[i].getRTime(), spec[i].getMZ());
    fprintf(fOut, " <precursorList>\n");
    for (j = 0; j<spec[i].sizePrecursor(); j++){
      fprintf(fOut, "  <precursor mono_mass=\"%.8lf\" charge=\"%d\" corr=\"%.4lf\">\n", spec[i].getPrecursor(j).monoMass, spec[i].getPrecursor(j).charge, spec[i].getPrecursor(j).corr);
      fprintf(fOut, "   <peptideList>\n");
      tp=spec[i].getTopPeps(j);
      sc=tp->singletFirst;
      for (z = 0; z<params->intermediate; z++){
        if (sc==NULL) break;
        db.getPeptideSeq(db.getPeptideList()->at(sc->pep1).map->at(0).index, db.getPeptideList()->at(sc->pep1).map->at(0).start, db.getPeptideList()->at(sc->pep1).map->at(0).stop, strs);
        pepSeq.clear();
        for (k = 0; k<strlen(strs); k++){
          pepSeq += strs[k];
          for (x = 0; x<sc->modLen; x++){
            if (sc->mods[x].pos == k) {
              sprintf(strTmp, "[%.2lf]", sc->mods[x].mass);
              pepSeq+=strTmp;
            }
          }
          if (k == sc->k1) pepSeq+="[x]";
        }

        pep = db.getPeptide(sc->pep1);
        protSeq.clear();
        for (n = 0; n<pep.map->size(); n++){
          protSeq += db[pep.map->at(n).index].name;
          if (n<pep.map->size()-1) protSeq+='-';
        }
        fprintf(fOut, "    <peptide sequence=\"%s\" mass=\"%.8lf\" protein=\"%s\" num_tot_proteins=\"%d\" link_site=\"%d\" score=\"%.4lf\" complement_mass=\"%.8lf\">\n", &pepSeq[0], sc->mass, &protSeq[0], (int)pep.map->size(), sc->k1 + 1, sc->simpleScore,spec[i].getPrecursor((int)sc->pre).monoMass-sc->mass);
        if (sc->modLen>0){
          fprintf(fOut, "     <modificationList>\n");
          for (k = 0; k<sc->modLen; k++){
            fprintf(fOut, "      <modification position=\"%d\" mass=\"%.8lf\"/>\n",sc->mods[k].pos+1,sc->mods[k].mass);
          }
          fprintf(fOut, "     </modificationList>\n");
        }
        fprintf(fOut, "    </peptide>\n");
        sc=sc->next;
      }
      fprintf(fOut, "  </peptideList>\n");
      fprintf(fOut,"  </precursor>\n");
    }
    
    fprintf(fOut, " </precursorList>\n");
    fprintf(fOut, "</spectrum>\n");
  }
  fclose(fOut);
  return true;
}

bool KData::outputMzID(CMzIdentML& m, KDatabase& db, KParams& par, kResults& r){
  size_t i;

  char str[256];
  string cStr;

  //Add/obtain the reference id for the psm data file.
  CSpectraData* m_sd = m.dataCollection.inputs.addSpectraData(params->inFile);
  CSearchDatabase* m_db = m.dataCollection.inputs.addSearchDatabase(params->dbFile);

  //Only one spectrum identification list per kojak search, so here it is.
  CSpectrumIdentificationList* m_sil = &m.dataCollection.analysisData.spectrumIdentificationList[0];
  
  //Only one spectrum identification protocol per kojak search
  CSpectrumIdentificationProtocol* m_sip=&m.analysisProtocolCollection.spectrumIdentificationProtocol[0];
 
  //Create spectrum
  CSpectrumIdentificationResult* m_sir;
  if(false){ //if this is the same as the last spectrum
    m_sir=&m_sil->spectrumIdentificationResult.back();
  } else {
    m_sir=new CSpectrumIdentificationResult();
    sprintf(str, "scan=%d", r.scanNumber);
    m_sir->spectrumID=r.scanID;
    m_sir->name=str;
    m_sir->spectraDataRef=m_sd->id;
    cStr = str;
    sprintf(str, "%s_%d", m_sd->id.c_str(), m_sil->spectrumIdentificationResult.size());
    m_sir->id = str;
    m_sil->spectrumIdentificationResult.push_back(*m_sir);
    m_sir = &m_sil->spectrumIdentificationResult.back();
  }

  //Declare the required classes
  CPeptide m_p,m_p2;

  //Get peptide #1 sequence and modifications
  m_p.peptideSequence.text=r.peptide1;
  for (i = 0; i < r.mods1.size(); i++){
    CModification m_m;
    sCvParam cv;
    m_m.location = (int)r.mods1[i].pos + 1;
    m_m.monoisotopicMassDelta = r.mods1[i].mass;
    if (m_m.location>0) {
      m_m.residues = r.peptide1[r.mods1[i].pos];
      cv = m_sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, m_m.residues);
    } else if(r.mods1[i].pos==-1){
      m_m.residues=".";
      cv = m_sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, m_m.residues,true);
    } else {
      m_m.residues=".";
      cv = m_sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, m_m.residues,false,true);
    }
    m_m.cvParam.push_back(cv);
    m_p.modification.push_back(m_m);
  }
  //add static mods, too
  for (i = 0; i<r.peptide1.size(); i++){
    if (aa.getFixedModMass(r.peptide1[i])>0) {
      CModification m_m;
      sCvParam cv;
      m_m.location = (int)i + 1;
      m_m.monoisotopicMassDelta = aa.getFixedModMass(r.peptide1[i]);
      if (m_m.location>0) m_m.residues = r.peptide1[i];
      else m_m.residues.clear();
      cv = m_sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, m_m.residues, m_m.location == 0);
      m_m.cvParam.push_back(cv);
      m_p.modification.push_back(m_m);
    }
  }

  CSpectrumIdentificationItem m_sii;
  sprintf(str, "%s_%d", m_sir->id.c_str(), (int)m_sir->spectrumIdentificationItem.size());
  m_sii.id = str;
  m_sii.calculatedMassToCharge = (r.psmMass + r.charge*1.007276466) / r.charge;
  m_sii.chargeState = r.charge;
  m_sii.experimentalMassToCharge = (r.obsMass + r.charge*1.007276466) / r.charge;
  m_sii.rank = 1;

  //add scores
  m_sii.addPSMValue("Kojak", "kojak_score", r.score);
  m_sii.addPSMValue("Kojak", "delta_score", r.scoreDelta);
  m_sii.addPSMValue("Kojak", "ppm_error", r.ppm);
  m_sii.addPSMValue("Kojak", "e-value", r.eVal);
  m_sii.addPSMValue("Kojak", "ion_match", r.matches1 + r.matches2);
  
  if(r.type< 2){

    if(r.type==1){ //loop link
      CModification m_m;
      sCvParam cv;
      m_m.location = r.link1;
      m_m.monoisotopicMassDelta = r.xlMass;
      m_m.residues = r.peptide1[r.link1 - 1];
      cv.cvRef = "PSI-MS";
      cv.accession = "MS:1002509";
      cv.name = "cross-link donor";
      m_m.cvParam.push_back(cv);
      m_p.modification.push_back(m_m);

      m_m.clear();
      m_m.location = r.link2;
      m_m.monoisotopicMassDelta = 0;
      m_m.residues = r.peptide1[r.link2 - 1];
      cv.cvRef = "PSI-MS";
      cv.accession = "MS:1002510";
      cv.name = "cross-link acceptor";
      m_m.cvParam.push_back(cv);
      m_p.modification.push_back(m_m);

      m_sii.addPSMValue("Kojak", "link", r.link1, "pep");
      m_sii.addPSMValue("Kojak", "link", r.link2, "pep");
    }

    m_sii.addPSMValue("Kojak", "consecutive_ion_match", r.conFrag1);
    m_sii.peptideRef = m.sequenceCollection.addPeptide(m_p);

    //Add all proteins mapped by this peptide
    writeMzIDPE(m,m_sii,r.pep1,db);
    m_sir->spectrumIdentificationItem.push_back(m_sii);

  } else if(r.type==2){
    CModification m_m;
    sCvParam cv;
    
    //Get peptide #2 sequence and modifications
    m_p2.peptideSequence.text = r.peptide2;
    for (i = 0; i < r.mods2.size(); i++){
      m_m.clear();
      m_m.location = (int)r.mods2[i].pos + 1;
      m_m.monoisotopicMassDelta = r.mods2[i].mass;
      if (m_m.location>0) {
        m_m.residues = r.peptide2[r.mods2[i].pos];
        cv = m_sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, m_m.residues);
      } else if (r.mods2[i].pos == -1){
        m_m.residues = ".";
        cv = m_sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, m_m.residues, true);
      } else {
        m_m.residues = ".";
        cv = m_sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, m_m.residues, false, true);
      }
      m_m.cvParam.push_back(cv);
      m_p2.modification.push_back(m_m);
    }
    //add static mods, too
    for (i = 0; i<r.peptide2.size(); i++){
      if (aa.getFixedModMass(r.peptide2[i])>0) {
        m_m.clear();
        m_m.location = (int)i + 1;
        m_m.monoisotopicMassDelta = aa.getFixedModMass(r.peptide2[i]);
        if (m_m.location>0) m_m.residues = r.peptide2[i];
        else m_m.residues.clear();
        cv = m_sip->modificationParams.back().getModificationCvParam(m_m.monoisotopicMassDelta, m_m.residues, m_m.location == 0);
        m_m.cvParam.push_back(cv);
        m_p2.modification.push_back(m_m);
      }
    }

    //Add cross-linker modifications
    m_m.clear();
    m_m.location=r.link1;
    m_m.monoisotopicMassDelta=r.xlMass;
    m_m.residues=r.peptide1[r.link1-1];
    m_m.cvParam.push_back(m_sip->modificationParams.back().getModificationCvParam(r.xlMass,m_m.residues));
    if(m_m.location==1 && m_m.cvParam[0].accession.compare("MS:1001460")==0){
      m_m.cvParam.clear();
      m_m.location=0;
      m_m.residues.clear();
      m_m.cvParam.push_back(m_sip->modificationParams.back().getModificationCvParam(r.xlMass, ".",true));
    } else if (m_m.cvParam[0].accession.compare("MS:1001460") == 0){
      m_m.cvParam.clear();
      m_m.location+=1;
      m_m.residues.clear();
      m_m.cvParam.push_back(m_sip->modificationParams.back().getModificationCvParam(r.xlMass, ".", false, true));
    }
    m_p.modification.push_back(m_m);

    m_m.clear();
    m_m.location = r.link2;
    m_m.monoisotopicMassDelta = 0;
    m_m.residues = r.peptide2[r.link2 - 1];
    m_m.cvParam.push_back(m_sip->modificationParams.back().getModificationCvParam(0, m_m.residues));
    if (m_m.location == 1 && m_m.cvParam[0].accession.compare("MS:1001460") == 0){
      m_m.cvParam.clear();
      m_m.location = 0;
      m_m.residues.clear();
    } else if (m_m.cvParam[0].accession.compare("MS:1001460") == 0){
      m_m.cvParam.clear();
      m_m.location += 1;
      m_m.residues.clear();
    } else m_m.cvParam.clear();
    m_p2.modification.push_back(m_m);

    //Create identifier for this combination of linked peptides
    string ID=r.peptide1;
    string pRef1,pRef2,xlValue;
    sprintf(str, "%d", r.link1);
    ID+=str;
    for (i = 0; i < r.mods1.size(); i++) {
      sprintf(str, "[%d,%.2lf]", r.mods1[i].pos, r.mods1[i].mass);
      ID += str;
    }
    ID += r.peptide2;
    sprintf(str, "%d", r.link2);
    ID += str;
    for (i = 0; i < r.mods2.size(); i++) {
      sprintf(str, "[%d,%.2lf]", r.mods2[i].pos, r.mods2[i].mass);
      ID += str;
    }
    m.sequenceCollection.addXLPeptides(ID,m_p,m_p2,pRef1,pRef2,xlValue);

    m_sii.peptideRef = pRef1;
    m_sii.addCvParam("MS:1002511", "PSI-MS", "cross-link spectrum identification item", "", "", "", xlValue);
    m_sii.addPSMValue("Kojak", "score", r.scoreA, "pep");
    m_sii.addPSMValue("Kojak", "rank", r.rankA, "pep");
    m_sii.addPSMValue("Kojak", "link", r.link1, "pep");
    m_sii.addPSMValue("Kojak", "e-value", r.eVal1, "pep");
    m_sii.addPSMValue("Kojak", "ion_match", r.matches1, "pep");
    m_sii.addPSMValue("Kojak", "consecutive_ion_match", r.conFrag1, "pep");
    writeMzIDPE(m, m_sii, r.pep1, db);
    m_sir->spectrumIdentificationItem.push_back(m_sii);

    //Add PSM
    CSpectrumIdentificationItem m_sii2;
    sprintf(str, "%s_%d", m_sir->id.c_str(), (int)m_sir->spectrumIdentificationItem.size());
    m_sii2.id = str;
    m_sii2.calculatedMassToCharge = (r.psmMass + r.charge*1.007276466) / r.charge;
    m_sii2.chargeState = r.charge;
    m_sii2.experimentalMassToCharge = (r.obsMass + r.charge*1.007276466) / r.charge;
    m_sii2.rank = 1;

    //add scores
    m_sii2.addPSMValue("Kojak", "kojak_score", r.score);
    m_sii2.addPSMValue("Kojak", "delta_score", r.scoreDelta);
    m_sii2.addPSMValue("Kojak", "ppm_error", r.ppm);
    m_sii2.addPSMValue("Kojak", "e-value", r.eVal);
    m_sii2.addPSMValue("Kojak", "ion_match", r.matches1 + r.matches2);
    //m_sii2.addPSMValue("Kojak", "consecutive_ion_match", (r.conFrag1 + r.conFrag2) / 2);
    m_sii2.addPSMValue("Kojak", "score", r.scoreB, "pep");
    m_sii2.addPSMValue("Kojak", "rank", r.rankB, "pep");
    m_sii2.addPSMValue("Kojak", "link", r.link2, "pep");
    m_sii2.addPSMValue("Kojak", "e-value", r.eVal2, "pep");
    m_sii2.addPSMValue("Kojak", "ion_match", r.matches2, "pep");
    m_sii2.addPSMValue("Kojak", "consecutive_ion_match", r.conFrag2, "pep");

    m_sii2.peptideRef = pRef2;
    m_sii2.addCvParam("MS:1002511", "PSI-MS", "cross-link spectrum identification item", "", "", "", xlValue);
    writeMzIDPE(m, m_sii2, r.pep2, db);
    m_sir->spectrumIdentificationItem.push_back(m_sii2);


  } else {

    //error: dimers not supported in mzID (or perhaps even in Kojak)
  }


  //if (r.type == 2){
  //  mod.location = r.link1;
  //  mod.monoisotopicMassDelta = r.xlMass;
  //  mod.residues = r.peptide1[r.link1 - 1];
  //  mod.cvParam->at(0).cvRef = "XLtmp";
  //  mod.cvParam->push_back(sip->modificationParams.getModificationCvParam(mod.monoisotopicMassDelta, mod.residues));
  //  mods.push_back(mod);
  //  mod.clear();
  //}
  //for (i = 0; i < r.mods1.size(); i++){
  //  mod.location = (int)r.mods1[i].pos + 1;
  //  mod.monoisotopicMassDelta = r.mods1[i].mass;
  //  if (mod.location>0) mod.residues = r.peptide1[r.mods1[i].pos];
  //  else mod.residues.clear();
  //  mod.cvParam->at(0) = sip->modificationParams.getModificationCvParam(mod.monoisotopicMassDelta, mod.residues, mod.location == 0);
  //  mods.push_back(mod);
  //  mod.clear();
  //}
  //add static mods, too
  //for (i = 0; i<r.peptide1.size(); i++){
  //  if (aa.getFixedModMass(r.peptide1[i])>0) {
  //    mod.location = (int)i + 1;
  //    mod.monoisotopicMassDelta = aa.getFixedModMass(r.peptide1[i]);
  //    if (mod.location>0) mod.residues = r.peptide1[i];
  //    else mod.residues.clear();
  //    mod.cvParam->at(0) = sip->modificationParams.getModificationCvParam(mod.monoisotopicMassDelta, mod.residues, mod.location == 0);
  //    mods.push_back(mod);
  //    mod.clear();
  //  }
  //}
  //if (r.type == 2){
  //  vector<CModification> mods2;
  //  mod.location = r.link2;
  //  mod.residues = r.peptide2[r.link2 - 1];
  //  mod.monoisotopicMassDelta = r.xlMass;
  //  mod.cvParam->at(0).cvRef = "XLtmp";
  //  mod.cvParam->push_back(sip->modificationParams.getModificationCvParam(mod.monoisotopicMassDelta, mod.residues));
  //  mods2.push_back(mod);
  //  mod.clear();
  //  for (i = 0; i < r.mods2.size(); i++){
  //    mod.location = (int)r.mods2[i].pos + 1;
  //    mod.monoisotopicMassDelta = r.mods2[i].mass;
  //    if (mod.location>0) mod.residues = r.peptide2[r.mods2[i].pos];
  //    else mod.residues.clear();
  //    mod.cvParam->at(0) = sip->modificationParams.getModificationCvParam(mod.monoisotopicMassDelta, mod.residues, mod.location == 0);
  //    mods2.push_back(mod);
  //    mod.clear();
  //  }
  //  //add static mods, too
  //  for (i = 0; i<r.peptide2.size(); i++){
  //    if (aa.getFixedModMass(r.peptide2[i])>0) {
  //      mod.location = (int)i + 1;
  //      mod.monoisotopicMassDelta = aa.getFixedModMass(r.peptide2[i]);
  //      if (mod.location>0) mod.residues = r.peptide2[i];
  //      else mod.residues.clear();
  //      mod.cvParam->at(0) = sip->modificationParams.getModificationCvParam(mod.monoisotopicMassDelta, mod.residues, mod.location == 0);
  //      mods2.push_back(mod);
  //      mod.clear();
  //    }
  //  }
  //  if (!m.addXLPeptides(r.peptide1, mods, peptide_ref, r.peptide2, mods2, peptide_ref2, value)){
  //    cout << "ERROR adding XLPeptides." << endl;
  //    exit(891);
  //  }
  //} else {
  //  peptide_ref = m.addPeptide(r.peptide1, mods);
  //  peptide_ref2="null";
  //}

  ////Add all proteins mapped by this peptide
  //pepRef.clear();
  //pep = db.getPeptide(r.pep1);
  //for (i = 0; i<pep.map->size(); i++){
  //  if (pep.n15 && db[pep.map->at(i).index].name.find(params->n15Label) == string::npos) {
  //    continue;
  //  }
  //  if (!pep.n15 && strlen(params->n15Label)>0 && db[pep.map->at(i).index].name.find(params->n15Label) != string::npos) {
  //    continue;
  //  }
  //  protein = db[pep.map->at(i).index].name;
  //  proteinDesc = "";
  //  dbSequence_ref = m.addDBSequence(protein, searchDatabase_ref, proteinDesc);

  //  if (pep.map->at(i).start<1) pre = '-';
  //  else pre = db[pep.map->at(i).index].sequence[pep.map->at(i).start - 1];
  //  if ((size_t)pep.map->at(i).stop + 1 == db[pep.map->at(i).index].sequence.size()) post = '-';
  //  else post = db[pep.map->at(i).index].sequence[(size_t)pep.map->at(i).stop + 1];
  //  isDecoy = db[pep.map->at(i).index].decoy;

  //  pepRef.push_back(m.addPeptideEvidence(dbSequence_ref, peptide_ref,(int)pep.map->at(i).start+1,(int)pep.map->at(i).stop+1,pre,post,isDecoy));
  //}

  ////Add PSM
  //sii = sir->addSpectrumIdentificationItem(r.charge, (r.obsMass + 1.007276466*r.charge) / r.charge, 1, pepRef, true, peptide_ref);
  //sii->calculatedMassToCharge = (r.psmMass + 1.007276466*r.charge) / r.charge;
  //if (r.type == 2) sii->addCvParam("MS:1002511", "PSI-MS", "cross-link spectrum identification item", "", "", "", value);
  //if(sir->spectrumIdentificationItem->size()==1) sir->addCvParam("MS:1000894", (double)r.rTime);

  ////if we have 2nd peptide
  //if (peptide_ref2.compare("null") != 0){
  //  //Add all proteins mapped by this peptide
  //  pepRef.clear();
  //  pep = db.getPeptide(r.pep2);
  //  for (i = 0; i<pep.map->size(); i++){
  //    if (pep.n15 && db[pep.map->at(i).index].name.find(params->n15Label) == string::npos) continue;
  //    if (!pep.n15 && strlen(params->n15Label)>0 && db[pep.map->at(i).index].name.find(params->n15Label) != string::npos) continue;
  //    protein = db[pep.map->at(i).index].name;
  //    proteinDesc = "";
  //    dbSequence_ref = m.addDBSequence(protein, searchDatabase_ref, proteinDesc);

  //    if (pep.map->at(i).start<1) pre = '-';
  //    else pre = db[pep.map->at(i).index].sequence[pep.map->at(i).start - 1];
  //    if ((size_t)pep.map->at(i).stop + 1 == db[pep.map->at(i).index].sequence.size()) post = '-';
  //    else post = db[pep.map->at(i).index].sequence[(size_t)pep.map->at(i).stop + 1];
  //    isDecoy = db[pep.map->at(i).index].decoy;

  //    pepRef.push_back(m.addPeptideEvidence(dbSequence_ref, peptide_ref2, pep.map->at(i).start + 1, pep.map->at(i).stop + 1, pre, post, isDecoy));
  //  }

  //  //Add PSM
  //  sii = sir->addSpectrumIdentificationItem(r.charge, (r.obsMass + 1.007276466*r.charge) / r.charge, 1, pepRef, true, peptide_ref2);
  //  sii->calculatedMassToCharge = (r.psmMass + 1.007276466*r.charge) / r.charge;
  //  sii->addCvParam("MS:1002511", "PSI-MS", "cross-link spectrum identification item", "", "", "", value);

  //  //add scores
  //  sii->addPSMValue("Kojak", "kojak_score", r.score);
  //  sii->addPSMValue("Kojak", "delta_score", r.scoreDelta);
  //  sii->addPSMValue("Kojak", "ppm_error", r.ppm);
  //  sii->addPSMValue("Kojak", "e-value", r.eVal);
  //  sii->addPSMValue("Kojak", "ion_match", r.matches1 + r.matches2);
  //  sii->addPSMValue("Kojak", "consecutive_ion_match", (r.conFrag1 + r.conFrag2) / 2);
  //  sii->addPSMValue("Kojak", "score", r.scoreB, "pep");
  //  sii->addPSMValue("Kojak", "rank", r.rankB, "pep");
  //  sii->addPSMValue("Kojak", "link", r.link2, "pep");
  //  sii->addPSMValue("Kojak", "e-value", r.eVal2, "pep");
  //  sii->addPSMValue("Kojak", "ion_match", r.matches2, "pep");
  //  sii->addPSMValue("Kojak", "consecutive_ion_match", r.conFrag2, "pep");
  //}

  //m_sir->spectrumIdentificationItem.push_back(m_sii);
  return true;
}

bool KData::outputNeoPepXML(CnpxSpectrumQuery& p, KDatabase& db, kResults& r){
  CnpxSearchHit sh;

  sh.hit_rank = 1;
  
  sh.calc_neutral_pep_mass = r.psmMass;
  sh.massdiff = r.psmMass - r.obsMass;
  
  if(r.type<2) { //single or loop
    sh.peptide = r.peptide1;
    if(r.type==0) sh.xlink_type="na";
    else sh.xlink_type="loop";

    CnpxModificationInfo mi = makeModificationInfo(r.mods1, r.peptide1, r.n15Pep1, r.nTerm1, r.cTerm1);
    if (!mi.modified_peptide.empty() || mi.mod_cterm_mass != 0 || mi.mod_nterm_mass != 0 || !mi.mod_aminoacid_mass.empty()) sh.modification_info.push_back(mi);

    addProteins(&sh, db, r.pep1, false, r.link1, r.link2);

    if(r.type==1) { //loop link
      CnpxXLink xl;
      xl.identifier = r.xlLabel;
      xl.mass = r.xlMass;
      addXlinkScore(xl, "link", r.link1);
      addXlinkScore(xl, "link", r.link2);
      sh.xlink.push_back(xl);
    }

  } else { //XL
    sh.peptide='-';
    sh.peptide_prev_aa='-';
    sh.peptide_next_aa='-';
    sh.protein='-';
    sh.num_tot_proteins=1;
    sh.xlink_type="xl";

    CnpxXLink xl;
    xl.identifier=r.xlLabel;
    xl.mass=r.xlMass;

    CnpxLinkedPeptide lp;
    lp.designation="alpha";
    lp.peptide=r.peptide1;
    lp.calc_neutral_pep_mass=r.massA;
    lp.complement_mass=r.obsMass-r.massA;

    CnpxModificationInfo mi=makeModificationInfo(r.mods1,r.peptide1,r.n15Pep1,r.nTerm1,r.cTerm1);
    if(!mi.modified_peptide.empty() || mi.mod_cterm_mass!=0 || mi.mod_nterm_mass!=0 || !mi.mod_aminoacid_mass.empty()) lp.modification_info.push_back(mi);

    addProteins(&lp,db,r.pep1,true,r.link1,0);
    
    addXlinkScore(lp,"score",r.scoreA,"%.4lf");
    addXlinkScore(lp,"rank",r.rankA);
    addXlinkScore(lp,"link",r.link1);
    addXlinkScore(lp,"e-value",r.eVal1,"%.3e");
    addXlinkScore(lp,"ion_match",r.matches1);
    addXlinkScore(lp,"consecutive_ion_match",r.conFrag1);
    
    xl.linked_peptide.push_back(lp);
    
    //second peptide
    CnpxLinkedPeptide lp2;
    lp2.designation="beta";
    lp2.peptide=r.peptide2;
    lp2.calc_neutral_pep_mass=r.massB;
    lp2.complement_mass=r.obsMass-r.massB;

    mi = makeModificationInfo(r.mods2, r.peptide2, r.n15Pep2, r.nTerm2, r.cTerm2);
    if (!mi.modified_peptide.empty() || mi.mod_cterm_mass != 0 || mi.mod_nterm_mass != 0 || !mi.mod_aminoacid_mass.empty()) lp2.modification_info.push_back(mi);

    addProteins(&lp2, db, r.pep2, true, r.link2, 0);
    
    addXlinkScore(lp2, "score", r.scoreB, "%.4lf");
    addXlinkScore(lp2, "rank", r.rankB);
    addXlinkScore(lp2, "link", r.link2);
    addXlinkScore(lp2, "e-value", r.eVal2, "%.3e");
    addXlinkScore(lp2, "ion_match", r.matches2);
    addXlinkScore(lp2, "consecutive_ion_match", r.conFrag2);

    xl.linked_peptide.push_back(lp2);

    sh.xlink.push_back(xl);

  } //XL

  addSearchScore(sh, "kojak_score", r.score, "%.4lf");
  addSearchScore(sh, "delta_score", r.scoreDelta, "%.4lf");
  addSearchScore(sh, "ppm_error", r.ppm, "%.4lf");
  addSearchScore(sh, "e-value", r.eVal, "%.3e");
  addSearchScore(sh, "ion_match", r.matches1+r.matches2);
  addSearchScore(sh, "consecutive_ion_match", r.conFrag1+r.conFrag2);


  if(p.search_result.empty()){
    CnpxSearchResult sr;
    sr.search_hit.push_back(sh);
    p.search_result.push_back(sr);
    p.spectrum=r.baseName+"."+to_string(r.scanNumber)+"."+to_string(r.scanNumber)+"."+to_string(r.charge);
    p.start_scan=r.scanNumber;
    p.end_scan=r.scanNumber;
    p.precursor_neutral_mass=r.obsMass;
    p.assumed_charge=r.charge;
    p.index=1;
    p.retention_time_sec=r.rTime*60;
  } else {
    p.search_result[0].search_hit.push_back(sh);
  }

  return true;
}

bool KData::outputPepXML(PXWSpectrumQuery& sq, KDatabase& db, kResults& r){

  unsigned int i;
  unsigned int j;

  char c;
  char n;
  char score[32];

  string peptide;
  string protein;
  string sequence;
  string tStr;

  kPeptide pep;
  kScoreCard sc;
  kScoreCard sc2;

  PXWSearchHit sh;
  PXWSearchHit shB;

  int siteA;
  int siteB;

  sh.hit_rank=1;
  sh.peptide=r.peptide1;
  if(r.type==0){
    sh.calc_neutral_pep_mass=r.psmMass;
    sh.massdiff=r.psmMass-r.obsMass;
  } else if(r.type==1){
    sh.calc_neutral_pep_mass=r.psmMass;
    sh.massdiff=r.psmMass-r.obsMass;
  } else {
    sh.calc_neutral_xl_mass=r.psmMass;
    sh.xl_massdiff=r.psmMass-r.obsMass;
    sh.calc_neutral_pep_mass=r.massA;
    sh.massdiff=r.obsMass-r.massA;
  }
  
  for(i=0;i<r.mods1.size();i++){
    //if(sq.start_scan==31877) cout << r.mods1[i].mass << "\t" << (int)r.mods1[i].pos << "\t" << (int)r.mods1[i].term << endl;
    if (r.mods1[i].pos == -1) sh.modInfo.mod_nterm_mass = r.mods1[i].mass;
    else if (r.mods1[i].pos == -2) sh.modInfo.mod_cterm_mass = r.mods1[i].mass;
    else sh.modInfo.addMod((int)r.mods1[i].pos+1,r.mods1[i].mass+aa.getAAMass(r.peptide1[r.mods1[i].pos],r.n15Pep1),r.mods1[i].mass,true);
  }
  if (r.nTerm1 && aa.getFixedModMass('$')!=0)sh.modInfo.mod_nterm_mass += aa.getFixedModMass('$');
  if (r.cTerm1 && aa.getFixedModMass('%')!=0)sh.modInfo.mod_cterm_mass += aa.getFixedModMass('%');
  sh.modInfo.mod_nterm_mass += aa.getFixedModMass('n');
  sh.modInfo.mod_cterm_mass += aa.getFixedModMass('c');
  if(sh.modInfo.mod_nterm_mass!=0) sh.modInfo.mod_nterm_mass+=1.00782503;
  if(sh.modInfo.mod_cterm_mass!=0) sh.modInfo.mod_cterm_mass+=17.00273963;
  for(i=0;i<r.peptide1.size();i++){
    if(aa.getFixedModMass(r.peptide1[i])>0) {
      sh.modInfo.addMod(i+1,aa.getAAMass(r.peptide1[i],r.n15Pep1),aa.getFixedModMass(r.peptide1[i]),false);
    }
  }

  if(r.type>1){
    shB.hit_rank=1;
    shB.peptide=r.peptide2;
    shB.calc_neutral_pep_mass=r.massB;
    shB.massdiff=r.obsMass-r.massB;

    for(i=0;i<r.mods2.size();i++){
      if (r.mods2[i].pos == -1) shB.modInfo.mod_nterm_mass = r.mods2[i].mass;
      else if (r.mods2[i].pos ==-2) shB.modInfo.mod_cterm_mass = r.mods2[i].mass;
      else shB.modInfo.addMod((int)r.mods2[i].pos+1,r.mods2[i].mass+aa.getAAMass(r.peptide2[r.mods2[i].pos],r.n15Pep2),r.mods2[i].mass,true);
    }
    if (r.nTerm2 && aa.getFixedModMass('$')!=0)shB.modInfo.mod_nterm_mass += aa.getFixedModMass('$');
    if (r.cTerm2 && aa.getFixedModMass('%')!=0)shB.modInfo.mod_cterm_mass += aa.getFixedModMass('%');
    shB.modInfo.mod_nterm_mass += aa.getFixedModMass('n');
    shB.modInfo.mod_cterm_mass += aa.getFixedModMass('c');
    if (shB.modInfo.mod_nterm_mass != 0) shB.modInfo.mod_nterm_mass += 1.00782503;
    if (shB.modInfo.mod_cterm_mass != 0) shB.modInfo.mod_cterm_mass += 17.00273963;
    for(i=0;i<r.peptide2.size();i++){
      if(aa.getFixedModMass(r.peptide2[i])>0) {
        shB.modInfo.addMod(i+1,aa.getAAMass(r.peptide2[i],r.n15Pep2),aa.getFixedModMass(r.peptide2[i]),false);
      }
    }
  }

  if(r.type>0) {
    if(r.type==1) {
      sh.xlink_type="loop";
      sprintf(score,"%d",r.link1);
      sh.addXLScore("link",score);
      sprintf(score,"%d",r.link2);
      sh.addXLScore("link",score);
    } else {
      sh.xlink_type="xl";
      sprintf(score,"%.4lf",r.scoreA);
      sh.addXLScore("score",score);
      sprintf(score,"%d",r.rankA);
      sh.addXLScore("rank",score);
      sprintf(score,"%d",r.link1);
      sh.addXLScore("link",score);
      sprintf(score, "%.3e", r.eVal1);
      sh.addXLScore("e-value", score);
      sprintf(score, "%d", r.matches1);
      sh.addXLScore("ion_match", score);
      sprintf(score, "%d", r.conFrag1);
      sh.addXLScore("consecutive_ion_match", score);

      shB.xlink_type="xl";
      sprintf(score,"%.4lf",r.scoreB);
      shB.addXLScore("score",score);
      sprintf(score,"%d",r.rankB);
      shB.addXLScore("rank",score);
      sprintf(score,"%d",r.link2);
      shB.addXLScore("link",score);
      sprintf(score, "%.3e", r.eVal2);
      shB.addXLScore("e-value", score);
      sprintf(score, "%d", r.matches2);
      shB.addXLScore("ion_match", score);
      sprintf(score, "%d", r.conFrag2);
      shB.addXLScore("consecutive_ion_match", score);
    }
    
  }

  sprintf(score,"%.4lf",r.score);
  sh.addScore("kojak_score",score);
  sprintf(score,"%.4lf",r.scoreDelta);
  sh.addScore("delta_score",score);
  sprintf(score,"%.4lf",r.ppm);
  sh.addScore("ppm_error",score);
  sprintf(score, "%.3e", r.eVal);
  sh.addScore("e-value", score);
  sprintf(score, "%d", r.matches1+r.matches2);
  sh.addScore("ion_match", score);
  if(r.type<2) sprintf(score, "%d", r.conFrag1);
  else sprintf(score, "%d", (r.conFrag1 + r.conFrag2)/2);
  sh.addScore("consecutive_ion_match", score);

  //Get proteins
  pep = db.getPeptide(r.pep1);
  sh.num_tot_proteins=(int)pep.map->size();
  for(j=0;j<pep.map->size();j++){
    if (pep.n15 && db[pep.map->at(j).index].name.find(params->n15Label)==string::npos) {
      sh.num_tot_proteins--;
      continue;
    }
    if (!pep.n15 && strlen(params->n15Label)>0 && db[pep.map->at(j).index].name.find(params->n15Label) != string::npos) {
      sh.num_tot_proteins--;
      continue;
    }
    protein="";
    for(i=0;i<db[pep.map->at(j).index].name.size();i++){
      if(params->truncate>0 && i==params->truncate) break;
      protein+=db[pep.map->at(j).index].name[i];
    }
    if(pep.map->at(j).start<1) n='-';
    else n=db[pep.map->at(j).index].sequence[pep.map->at(j).start-1];
    if(pep.map->at(j).stop+1==db[pep.map->at(j).index].sequence.size()) c='-';
    else c=db[pep.map->at(j).index].sequence[pep.map->at(j).stop+1];
    siteA = pep.map->at(j).start+r.link1;
    if(r.type==1){
      siteB = pep.map->at(j).start+r.link2;
      sh.addProtein(protein, c, n, (int)pep.map->at(j).start + 1,siteA, siteB);
    } else {
      sh.addProtein(protein, c, n, (int)pep.map->at(j).start + 1, siteA);
    }
  }

  if(r.type>1){
    pep = db.getPeptide(r.pep2);
    shB.num_tot_proteins=(int)pep.map->size();
    for(j=0;j<pep.map->size();j++){
      if (pep.n15 && db[pep.map->at(j).index].name.find(params->n15Label) == string::npos) {
        sh.num_tot_proteins--;
        continue;
      }
      if (!pep.n15 && strlen(params->n15Label)>0 && db[pep.map->at(j).index].name.find(params->n15Label) != string::npos) {
        sh.num_tot_proteins--;
        continue;
      }
      protein="";
      for(i=0;i<db[pep.map->at(j).index].name.size();i++){
        if(params->truncate>0 && i==params->truncate) break;
        protein+=db[pep.map->at(j).index].name[i];
      }
      if(pep.map->at(j).start<1) n='-';
      else n=db[pep.map->at(j).index].sequence[pep.map->at(j).start-1];
      if(pep.map->at(j).stop+1==db[pep.map->at(j).index].sequence.size()) c='-';
      else c=db[pep.map->at(j).index].sequence[pep.map->at(j).stop+1];
      siteA = pep.map->at(j).start + r.link2;
      shB.addProtein(protein, c, n, (int)pep.map->at(j).start+1,siteA);
    }
  }

  if(r.type==0) {
    sq.addSearchHit(&sh,NULL,NULL,NULL);
  } else if(r.type==1) {
    sq.addSearchHit(&sh,NULL,&r.xlLabel,&r.xlMass);
  } else {
    sq.addSearchHit(&sh,&shB,&r.xlLabel,&r.xlMass);
  }

  return true;
}

bool KData::outputPercolator(FILE* f, KDatabase& db, kResults& r, int count){

  unsigned int i;
  unsigned int j;

  string peptide;
  string protein;
  string sequence;
  string tStr;
  string p1,p2;

  kPeptide pep;
  kScoreCard sc;
  kScoreCard sc2;

  //Export Results:
  if(r.decoy) fprintf(f,"D-");
  else fprintf(f,"T-");
  fprintf(f,"%s-%d-%.2f",r.baseName.c_str(),r.scanNumber,r.rTime);
  if(count>1) fprintf(f,"-%d",count);
  if(r.decoy) fprintf(f,"\t-1");
  else fprintf(f,"\t1");
  if(params->percVersion>2.04) fprintf(f,"\t%d",r.scanNumber);
  fprintf(f,"\t%.4lf",r.score);
  fprintf(f,"\t%.4lf",r.scoreDelta);
  fprintf(f,"\t%.6lf",-log10(r.eVal));
  if(r.type==2 || r.type==3) {
    if (r.scoreA>r.scoreB) fprintf(f,"\t%.6lf\t%.6lf\t%d\t%d\t%d\t%d\t%d\t%d",-log10(r.eVal1),-log10(r.eVal2),r.matches1+r.matches2,(r.conFrag1+r.conFrag2)/2,r.matches1,r.conFrag1,r.matches2,r.conFrag2);
    else fprintf(f, "\t%.6lf\t%.6lf\t%d\t%d\t%d\t%d\t%d\t%d", -log10(r.eVal2), -log10(r.eVal1), r.matches1 + r.matches2, (r.conFrag1 + r.conFrag2) / 2, r.matches2, r.conFrag2, r.matches1, r.conFrag1);
    fprintf(f,"\t%d\t%.4lf",r.rank,r.scorePepDif);
  } else {
    fprintf(f,"\t%d\t%d",r.matches1,r.conFrag1);
  }
  //if(r.type==1) fprintf(f,"\t1\t0");
  //else if(r.type==2) fprintf(f,"\t0\t1");
  //else if(r.type==3) fprintf(f,"\t0\t0");
  //else fprintf(f,"\t0\t0");
  for(int z=1;z<8;z++){
    if(r.charge==(int)z) fprintf(f,"\t1");
    else fprintf(f,"\t0");
  }
  if(r.charge>7) fprintf(f,"\t1");
  else fprintf(f,"\t0");
  fprintf(f,"\t%.4lf",r.psmMass);
  fprintf(f,"\t%.4lf",r.ppm);
  p1=r.modPeptide1;
  p2=r.modPeptide2;
  if(r.n15Pep1)p1+="-15N";
  if(r.n15Pep2)p2+="-15N";
  if (r.type == 2 || r.type == 3) {
    if (r.peptide1.size()>r.peptide2.size()) fprintf(f,"\t%d\t%d", (int)r.peptide2.size(), (int)r.peptide1.size());
    else fprintf(f,"\t%d\t%d", (int)r.peptide1.size(), (int)r.peptide2.size());
    if(r.type==3) fprintf(f,"\t%d\t-.%s+%s.-",(int)(r.peptide1.size()+r.peptide2.size()),&p1[0],&p2[0]);
    else fprintf(f,"\t%d\t-.%s(%d)--%s(%d).-",(int)(r.peptide1.size()+r.peptide2.size()),&p1[0],r.link1,&p2[0],r.link2);
  } else {
    fprintf(f,"\t%d\t-.%s",(int)r.peptide1.size(),&p1[0]);
    if(r.type==1) fprintf(f,"(%d,%d)-LOOP",r.link1,r.link2);
    fprintf(f,".-");
  }
  

  //export proteins
  pep = db.getPeptide(r.pep1);
  for(j=0;j<pep.map->size();j++){
    protein="";
    for(i=0;i<db[pep.map->at(j).index].name.size();i++){
      if(params->truncate>0 && i==params->truncate) break;
      if(db[pep.map->at(j).index].name[i]==' ') protein+='_';
      else protein+=db[pep.map->at(j).index].name[i];
    }
    fprintf(f,"\t%s",&protein[0]);
  }
  if(r.pep2>=0){
    pep = db.getPeptide(r.pep2);
    for(j=0;j<pep.map->size();j++){
      protein="";
      for(i=0;i<db[pep.map->at(j).index].name.size();i++){
        if(params->truncate>0 && i==params->truncate) break;
        if(db[pep.map->at(j).index].name[i]==' ') protein+='_';
        else protein+=db[pep.map->at(j).index].name[i];
      }
      fprintf(f,"\t%s",&protein[0]);
    }
  }

  fprintf(f,"\n");

  return true;
}

//Function deprecated. Should be excised.
bool KData::outputResults(KDatabase& db, KParams& par){

  size_t i;
  int j,k,n,d;
  char fName[1056];
  char outPath[1056];
  char peptide[256];
  char tmp[16];
  char specID[256];

  kPeptide pep;
  kPeptide pep2;
  kPrecursor precursor;
  kScoreCard tmpSC;
  kScoreCard tmpSC2;

  kEnzymeRules enzyme;
  kResults res;

  NeoPepXMLParser p;
  string pepXML_fileName;
  //PepXMLWriter p;
  //pxwAminoAcidModification aam;
  //pxwTerminalModification tm;
  //pxwMSMSRunSummary rs;
  //pxwSampleEnzyme enz;
  //PXWSearchSummary ss;
  //PXWSpectrumQuery sq;

  CMzIdentML mzID;
  string analysisSoftware_ref;
  string sip_ref;

  bool bBadFiles;
  bool bInter;
  bool bTarget1;
  bool bDupe;
  bool bDiag;

  int scoreIndex;
  int iDupe;

  double topScore;

  string tmpPep1;
  string tmpPep2;
  string outFile;
  string dStr;

  FILE* fOut    = NULL;
  FILE* fIntra  = NULL;
  FILE* fInter  = NULL;
  FILE* fLoop   = NULL;
  FILE* fSingle = NULL;
  FILE* fDimer  = NULL;
  FILE* fDiag   = NULL;

  //Export FASTA database if Kojak generated the decoys.
  if (params->buildDecoy) {
    outFile = params->dbFile;
    i = outFile.find_last_of("/\\");
    if (i != string::npos) outFile = outFile.substr(i + 1);
    outFile = params->fullPath + slashdir + outFile;
    outFile+=".kojak.fasta";
    db.exportDB(outFile);
    strcpy(params->dbFile, outFile.c_str());
    params->buildDecoy=false; //only build the decoy library once if in batch mode
  }

  //Open all the required output files.
  bBadFiles=false;
  sprintf(fName,"%s.kojak.txt",params->outFile);
  fOut=fopen(fName,"wt");
  if(fOut==NULL) bBadFiles=true;
  if(params->exportMzID){
    analysisSoftware_ref = mzID.addAnalysisSoftware("Kojak", version);
    writeMzIDDatabase(mzID,db);
    sip_ref=writeMzIDSIP(mzID,analysisSoftware_ref,par);
    CSpectraData* m_sd = mzID.dataCollection.inputs.addSpectraData(params->inFile);
    CSearchDatabase* m_db = mzID.dataCollection.inputs.addSearchDatabase(params->dbFile);
    CSpectrumIdentificationProtocol* m_sip = mzID.getSpectrumIdentificationProtocol(sip_ref);
    CSpectrumIdentificationList* m_sil = NULL;
    CSpectrumIdentification* si = mzID.addSpectrumIdentification(m_sd->id, m_db->id, m_sip->id, m_sil);
  }
  if(params->exportPercolator) {
    sprintf(fName,"%s.perc.intra.txt",params->outFile);
    fIntra=fopen(fName,"wt");
    if(fIntra==NULL) bBadFiles=true;
    sprintf(fName,"%s.perc.inter.txt",params->outFile);
    fInter=fopen(fName,"wt");
    if(fInter==NULL) bBadFiles=true;
    sprintf(fName,"%s.perc.loop.txt",params->outFile);
    fLoop=fopen(fName,"wt");
    if(fLoop==NULL) bBadFiles=true;
    sprintf(fName,"%s.perc.single.txt",params->outFile);
    fSingle=fopen(fName,"wt");
    if(fSingle==NULL) bBadFiles=true;
    if (params->dimers){
      sprintf(fName, "%s.perc.dimer.txt", params->outFile);
      fDimer = fopen(fName, "wt");
      if (fDimer == NULL) bBadFiles = true;
    }
  }
  if(params->exportPepXML) {
    CnpxMSMSPipelineAnalysis pa;

    char timebuf[80];
    time_t timeNow;
    time(&timeNow);
    strftime(timebuf, 80, "%Y-%m-%dT%H:%M:%S", localtime(&timeNow));
    pa.date.parseDateTime(timebuf);
    
    CnpxMSMSRunSummary rs;
    rs.base_name=params->inFile; 
    rs.base_name = rs.base_name.substr(0, rs.base_name.find_last_of('.'));
    outFile=params->outFile;
    if (rs.base_name[0] == '/'){ //unix
      outFile = outFile.substr(outFile.find_last_of("/") + 1, outFile.size());
    } else { //assuming windows
      outFile = outFile.substr(outFile.find_last_of("\\") + 1, outFile.size());
    }
    rs.raw_data=params->ext;
    rs.raw_data_type="raw";

    //Add the enzyme
    CnpxSampleEnzyme se;
    se.name = params->enzymeName;
    CnpxSpecificity ses;
    enzyme = db.getEnzymeRules();
    for (i = 65; i < 90; i++){
      if (enzyme.cutC[i]) ses.cut += (char)i;
      if (enzyme.exceptN[i]) ses.no_cut += (char)i;
    }
    if (ses.cut.size()>0){
      ses.sense = "C";
    } else {
      ses.sense = "N";
      for (i = 65; i < 90; i++){
        if (enzyme.cutN[i]) ses.cut += (char)i;
        if (enzyme.exceptC[i]) ses.no_cut += (char)i;
      }
    }
    se.specificity.push_back(ses);
    rs.sample_enzyme.push_back(se);
    
    CnpxSearchSummary ss;
    ss.search_engine="Kojak";
    ss.base_name=rs.base_name;
    ss.search_engine_version=version;
    ss.precursor_mass_type="monoisotopic";
    ss.fragment_mass_type="monoisotopic";
    ss.search_id=1;
    processPath(params->dbFile,outPath);

    CnpxSearchDatabase sd;
    sd.local_path=outPath;
    sd.type="AA";
    ss.search_database.push_back(sd);

    CnpxEnzymaticSearchConstraint esc;
    esc.enzyme=params->enzymeName;
    esc.max_num_internal_cleavages=params->miscleave;
    esc.min_number_termini=2;
    ss.enzymatic_search_constraint.push_back(esc);

    //Add modifications
    for (i = 0; i<params->fMods->size(); i++){ 
      if (params->fMods->at(i).index == '$' || params->fMods->at(i).index=='%') { //special case protein termini
        CnpxTerminalModification tm;
        if (params->fMods->at(i).index == '$') tm.terminus="n";
        else tm.terminus="c";
        tm.protein_terminus="Y";
        tm.massdiff = params->fMods->at(i).mass;
        tm.mass = params->fMods->at(i).mass;
        tm.variable = "N";
        ss.terminal_modification.push_back(tm);
      } else if (params->fMods->at(i).index == 'n' || params->fMods->at(i).index == 'c') { //peptide termini
        CnpxTerminalModification tm;
        if (params->fMods->at(i).index == 'n') tm.terminus = "n";
        else tm.terminus = "c";
        tm.protein_terminus = "N";
        tm.massdiff = params->fMods->at(i).mass;
        tm.mass = params->fMods->at(i).mass;
        tm.variable = "N";
        ss.terminal_modification.push_back(tm);
      } else {
        CnpxAminoAcidModification aam;
        aam.aminoacid = (char)params->fMods->at(i).index;
        aam.massdiff = params->fMods->at(i).mass;
        aam.mass = db.getAAMass(params->fMods->at(i).index);
        aam.variable = "N";
        ss.aminoacid_modification.push_back(aam);
      }
    }
    for(i=0;i<params->mods->size();i++){
      if (params->mods->at(i).index == '$' || params->mods->at(i).index == '%') { //special case protein termini
        CnpxTerminalModification tm;
        if (params->mods->at(i).index == '$') tm.terminus="n";
        else tm.terminus="c";
        tm.protein_terminus = "Y";
        tm.massdiff = params->mods->at(i).mass;
        tm.mass = db.getAAMass(params->mods->at(i).index) + params->mods->at(i).mass;
        tm.variable = "Y";
        ss.terminal_modification.push_back(tm);
      } else if (params->mods->at(i).index == 'n' || params->mods->at(i).index == 'c') { //peptide termini
        CnpxTerminalModification tm;
        if (params->mods->at(i).index == 'n') tm.terminus = "n";
        else tm.terminus = "c";
        tm.protein_terminus = "N";
        tm.massdiff = params->mods->at(i).mass;
        tm.mass = db.getAAMass(params->mods->at(i).index) + params->mods->at(i).mass;
        tm.variable = "Y";
        ss.terminal_modification.push_back(tm);
      } else {
        CnpxAminoAcidModification aam;
        aam.aminoacid=(char)params->mods->at(i).index;
        aam.massdiff=params->mods->at(i).mass;
        aam.mass = db.getAAMass(params->mods->at(i).index) + params->mods->at(i).mass;
        aam.variable = "Y";
        ss.aminoacid_modification.push_back(aam);
      }
    }
    for(i=0;i<par.xmlParams.size();i++){
      CnpxParameter px;
      px.name=par.xmlParams[i].name;
      px.value=par.xmlParams[i].value;
      ss.parameter.push_back(px);
    }

    rs.search_summary.push_back(ss);
    pa.msms_run_summary.push_back(rs);
    
    sprintf(fName, "%s.pep.xml", params->outFile);
    pepXML_fileName=fName;
    FILE* ft=fopen(fName,"wt");
    if(ft==NULL) bBadFiles=true;
    else fclose(ft);
    pa.summary_xml = pepXML_fileName;
    
    p.msms_pipeline_analysis.push_back(pa);
  }
  
  if (params->diag->size()>0){ //create diagnostic file if needed
    sprintf(fName, "%s.diag.xml", params->outFile);
    fDiag = fopen(fName, "wt");
    if (fDiag == NULL) bBadFiles = true;
  }

  //check that all output files are valid
  if(bBadFiles){
    if(fOut!=NULL)    fclose(fOut);
    if(fIntra!=NULL)  fclose(fIntra);
    if(fInter!=NULL)  fclose(fInter);
    if(fLoop!=NULL)   fclose(fLoop);
    if(fSingle!=NULL) fclose(fSingle);
    if(fDimer!=NULL)  fclose(fDimer);
    if(fDiag!=NULL)   fclose(fDiag);
    klog->addError("Error exporting results. Please make sure drive is writable.");
    return false;
  }

  //Put the headers on all the files
  fprintf(fOut,"Kojak version %s\n",version);
  fprintf(fOut,"Scan Number\tRet Time\tObs Mass\tCharge\tPSM Mass\tPPM Error\tScore\tdScore\tE-value\tPeptide #1 Score\tPeptide #1 E-value\tPeptide #1\tLinked AA #1\tProtein #1\tProtein #1 Site\tPeptide #2 Score\tPeptide #2 E-value\tPeptide #2\tLinked AA #2\tProtein #2\tProtein #2 Site\tLinker Mass\n");
  if(params->exportPercolator){
    if(params->percVersion>2.04) {
      fprintf(fIntra,"SpecId\tLabel\tscannr\tScore\tdScore\tnegLog10eVal\tnegLog10eValA\tnegLog10eValB\tIonMatch\tConIonMatch\tIonMatchA\tConIonMatchA\tIonMatchB\tConIonMatchB\t");
      fprintf(fInter,"SpecId\tLabel\tscannr\tScore\tdScore\tnegLog10eVal\tnegLog10eValA\tnegLog10eValB\tIonMatch\tConIonMatch\tIonMatchA\tConIonMatchA\tIonMatchB\tConIonMatchB\t");
      fprintf(fLoop,"SpecId\tLabel\tscannr\tScore\tdScore\tnegLog10eVal\tIonMatch\tConIonMatch\t");
      fprintf(fSingle,"SpecId\tLabel\tscannr\tScore\tdScore\tnegLog10eVal\tIonMatch\tConIonMatch\t");
      if (params->dimers) fprintf(fDimer, "SpecId\tLabel\tscannr\tScore\tdScore\tnegLog10eVal\t");
    } else {
      fprintf(fIntra,"SpecId\tLabel\tScore\tdScore\tnegLog10eVal\tnegLog10eValA\tnegLog10eValB\tIonMatch\tConIonMatch\tIonMatchA\tConIonMatchA\tIonMatchB\tConIonMatchB\t");
      fprintf(fInter,"SpecId\tLabel\tScore\tdScore\tnegLog10eVal\tnegLog10eValA\tnegLog10eValB\tIonMatch\tConIonMatch\tIonMatchA\tConIonMatchA\tIonMatchB\tConIonMatchB\t");
      fprintf(fLoop,"SpecId\tLabel\tScore\tdScore\tnegLog10eVal\tIonMatch\tConIonMatch\t");
      fprintf(fSingle,"SpecId\tLabel\tScore\tdScore\tnegLog10eVal\tIonMatch\tConIonMatch\t");
      if (params->dimers) fprintf(fDimer, "SpecId\tLabel\tScore\tdScore\tnegLog10eVal\t");
    }
    fprintf(fIntra,"NormRank\tPPScoreDiff\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\tCharge6\tCharge7\tCharge8Plus\tMass\tPPM\tLenShort\tLenLong\tLenSum\tPeptide\tProteins\n");
    fprintf(fInter,"NormRank\tPPScoreDiff\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\tCharge6\tCharge7\tCharge8Plus\tMass\tPPM\tLenShort\tLenLong\tLenSum\tPeptide\tProteins\n");
    fprintf(fLoop,"Charge1\tCharge2\tCharge3\tCharge4\tCharge5\tCharge6\tCharge7\tCharge8Plus\tMass\tPPM\tLen\tPeptide\tProteins\n");
    fprintf(fSingle,"Charge1\tCharge2\tCharge3\tCharge4\tCharge5\tCharge6\tCharge7\tCharge8Plus\tMass\tPPM\tLen\tPeptide\tProteins\n");
    if (params->dimers) fprintf(fDimer, "NormRank\tPPScoreDiff\\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\tCharge6\tCharge7\tCharge8Plus\tMass\tPPM\tLenShort\tLenLong\tLenSum\tPeptide\tProteins\n");
  }
  if(fDiag!=NULL){
    fprintf(fDiag,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fDiag,"<kojak_analysis date=\"now\">\n");
  }

  res.baseName=params->outFile;
  if (res.baseName[0] == '/'){ //unix
    res.baseName = res.baseName.substr(res.baseName.find_last_of("/") + 1, res.baseName.size());
  } else { //assuming windows
    res.baseName = res.baseName.substr(res.baseName.find_last_of("\\") + 1, res.baseName.size());
  }

  //Output top score for each spectrum
  //Must iterate through all possible precursors for that spectrum
  dStr=params->decoy;
  for(i=0;i<spec.size();i++) {

    //update top hits so that a target result is always first among ties between targets and decoys
    spec[i].refreshScore(db,dStr);

    //Check if we need to output diagnostic information
    bDiag=false;
    if(params->diag->size()>0){
      if(params->diag->at(0)==-1) {
        bDiag=true;
      } else {
        for (d = 0; d<params->diag->size(); d++){
          if (spec[i].getScanNumber() == params->diag->at(d)){
            bDiag=true;
            break;
          }
        }
      }
    }
    if(bDiag) outputDiagnostics(fDiag,spec[i],db);

    scoreIndex=0;
    tmpSC=spec[i].getScoreCard(scoreIndex);
    res.scanNumber=spec[i].getScanNumber();
    res.scanID=spec[i].getNativeID();
    res.rTime=spec[i].getRTime();
    
    CnpxSpectrumQuery sq;
    sq.spectrum = res.baseName + "." + to_string(res.scanNumber) + "." + to_string(res.scanNumber) + "." + to_string(spec[i].getPrecursor(0).charge);
    sq.start_scan = res.scanNumber;
    sq.end_scan = res.scanNumber;
    sq.precursor_neutral_mass = spec[i].getPrecursor(0).monoMass;
    sq.assumed_charge = spec[i].getPrecursor(0).charge;
    sq.index = 1;
    sq.retention_time_sec = spec[i].getRTime() * 60;

    //if there are no matches to the spectrum, return null result and continue
    if(tmpSC.simpleScore==0){
      fprintf(fOut,"%d\t%.4f\t0\t0\t0\t0\t0\t0\t999\t0\t999\t-\t-\t-\t-\t0\t999\t-\t-\t-\t-\t0\n",res.scanNumber,res.rTime);

      if(params->exportPepXML) {
        sq.index = (int)p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.size() + 1;
        p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.push_back(sq);
      }
      continue;
    }

    

    //Export top scoring peptide, plus any ties that occur after it.
    topScore=tmpSC.simpleScore;
    int count=0;
    while(tmpSC.simpleScore==topScore){

      count++;

      //Get precursor ion for the PSM
      precursor=spec[i].getPrecursor((int)tmpSC.precursor);
      res.obsMass = precursor.monoMass;
      res.charge  = precursor.charge;
      res.ppm = (tmpSC.mass - precursor.monoMass) / precursor.monoMass*1e6;
      res.psmMass = tmpSC.mass;
      res.hk = precursor.corr;

      if(params->exportPepXML){
        sq.assumed_charge=res.charge;
        sq.precursor_neutral_mass=res.obsMass;
        sprintf(specID,"%s.%d.%d.%d",outFile.c_str(),res.scanNumber,res.scanNumber,res.charge);
        sq.spectrum=specID;
      }

      //grab the next highest score that matches to the same precursor ion for the delta score
      //do not count ties - look for the first difference
      //if no other match has the same precursor, just take the lowest score in the list
      n=scoreIndex+1;
      while(n<19){
        tmpSC2=spec[i].getScoreCard(n++);
        if(tmpSC2.simpleScore==0) break;
        if(tmpSC2.simpleScore==topScore) continue;
        if(tmpSC2.precursor!=tmpSC.precursor) continue;

        //if peptides and link sites are the same, go to the next one
        //this no longer applies to the top result. duplicates may occur among lower results
        //if(tmpSC.link>-1 && tmpSC2.link>-1 && tmpSC2.pep1==tmpSC.pep1 && tmpSC2.pep2==tmpSC.pep2 && tmpSC2.k1==tmpSC.k1 && tmpSC2.k2==tmpSC.k2){
        //  cout << "Oddity 1: " << spec[i].getScanNumber() << endl;
        //  continue;
        //}
        break;
      }
      res.score       = tmpSC.simpleScore;
      res.scoreDelta  = tmpSC.simpleScore-tmpSC2.simpleScore;
      res.eVal        = tmpSC.eVal;
      res.eVal1       = tmpSC.eVal1;
      res.eVal2       = tmpSC.eVal2;
      res.matches1    = tmpSC.matches1;
      res.matches2    = tmpSC.matches2;
      res.conFrag1    = tmpSC.conFrag1;
      res.conFrag2    = tmpSC.conFrag2;
      if(tmpSC.score1<tmpSC.score2) res.scorePepDif = tmpSC.score1;
      else res.scorePepDif = tmpSC.score2;
      
      KTopPeps* tp = spec[i].getTopPeps((int)tmpSC.precursor);
      kSingletScoreCard* grr;
      grr = tp->singletFirst;
      int rank=1;
      res.rankA=0;
      res.rankB=0;
      while(grr!=NULL){
        if(res.rankA==0 && tmpSC.pep1==grr->pep1 && tmpSC.score1==grr->simpleScore) res.rankA=rank;
        if(res.rankB==0 && tmpSC.score2>grr->simpleScore) res.rankB=rank;
        rank++;
        grr = grr->next;
      }
      if(res.rankB==0) res.rankB=params->topCount;
      res.rank        = res.rankA+res.rankB;
      res.scoreA      = tmpSC.score1;
      res.scoreB      = tmpSC.score2;
      res.massA       = tmpSC.mass1;
      res.massB       = tmpSC.mass2;

      //Get the peptide sequence(s)
      pep = db.getPeptide(tmpSC.pep1);
      db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
      res.peptide1 = peptide;
      res.mods1.clear();
      res.cTerm1 = pep.cTerm;
      res.nTerm1 = pep.nTerm;
      res.linkSite1 = tmpSC.site1;
      if(tmpSC.site2>-1) res.linkSite2=tmpSC.site2; //loop-link
      res.n15Pep1 = pep.n15;
      for(j=0;j<tmpSC.mods1->size();j++) res.mods1.push_back(tmpSC.mods1->at(j));
      res.peptide2 = "";
      if(tmpSC.pep2>=0){
        pep2 = db.getPeptide(tmpSC.pep2);
        db.getPeptideSeq( pep2.map->at(0).index,pep2.map->at(0).start,pep2.map->at(0).stop,peptide);
        res.peptide2 = peptide;
        res.mods2.clear();
        res.cTerm2 = pep2.cTerm;
        res.nTerm2 = pep2.nTerm;
        res.linkSite2 = tmpSC.site2;
        res.n15Pep2 = pep2.n15;
        for(j=0;j<tmpSC.mods2->size();j++) res.mods2.push_back(tmpSC.mods2->at(j));
      }

      //Process the peptide
      res.modPeptide1 = processPeptide(pep,tmpSC.mods1,db);      
      res.modPeptide2 = "";
      if(res.peptide2.size()>0){
        res.modPeptide2=processPeptide(pep2,tmpSC.mods2,db);
      }

      //Get the link positions - relative to the peptide
      res.link1 = tmpSC.k1;
      res.link2 = tmpSC.k2;
      if(res.link1>=0) res.link1++;
      if(res.link2>=0) res.link2++;

      //set link type
      res.type=0;
      if(tmpSC.k1>=0 && tmpSC.k2>=0) res.type=1;
      if(tmpSC.pep1>=0 && tmpSC.pep2>=0) res.type=2;
      if(res.type==2 && tmpSC.k1==-1 && tmpSC.k2==-1) res.type=3;

      if(res.type>0 && res.type!=3) {
        res.xlMass=link[tmpSC.link].mass;
        res.xlLabel=link[tmpSC.link].label;
      }

      //Get the peptide indexes
      res.pep1 = tmpSC.pep1;
      res.pep2 = tmpSC.pep2;
      res.linkable1 = tmpSC.linkable1;
      res.linkable2 = tmpSC.linkable2;

      //Edge case where single peptide is shared between linked and non-linked peptide lists
      //This occurs when the peptide appears multiple times in a database: internally and on
      //the c-terminus for amine reactive cross-linkers, for example.
      bDupe=false;
      if(res.type==0){
        n=scoreIndex+1;
        iDupe=1;
        while(n<19){
          iDupe++;
          tmpSC2=spec[i].getScoreCard(n++);
          if(tmpSC2.simpleScore==0) break;
          if(tmpSC2.simpleScore!=topScore) break;

          //if peptides are the same, but different lists (linked vs. non), use second peptide as location
          if(tmpSC2.linkable1!=tmpSC.linkable1) {
            pep = db.getPeptide(res.pep1);
            db.getPeptideSeq(pep,tmpPep1);
            pep2 = db.getPeptide(tmpSC2.pep1);
            db.getPeptideSeq(pep2,tmpPep2);
            if(tmpPep1.compare(tmpPep2)==0){
              res.pep2=tmpSC2.pep1;
              res.linkable2=tmpSC2.linkable1;
              res.linkSite2=tmpSC2.site1;
              bDupe=true;
              break;
            }
          }
        }
      }

      //Process the protein
      processProtein(res.pep1, res.link1-1, res.linkSite1, res.protein1, res.protPos1, res.decoy1, db);
      if (res.modPeptide2.size()>1) {
        processProtein(res.pep2, res.link2-1, res.linkSite2, res.protein2, res.protPos2, res.decoy2, db);
        if(res.decoy1 || res.decoy2) res.decoy=true;
        else res.decoy=false;
      } else if(res.linkSite2>-1){ //loop link special case.
        processProtein(res.pep1, res.link2 - 1, res.linkSite2, res.protein2, res.protPos2, res.decoy2, db);
        if(!res.decoy1 || !res.decoy2) res.decoy=false;
        else res.decoy=true;
      } else {
        res.decoy=res.decoy1;
      }

      tmpPep1=res.peptide1;
      sprintf(tmp,"(%d)",res.link1);
      tmpPep1+=tmp;
      tmpPep2 = res.peptide2;
      sprintf(tmp, "(%d)", res.link2);
      tmpPep2 += tmp;

      //Export Results:
      fprintf(fOut,"%d",res.scanNumber);
      //fprintf(fOut, "\t%.4lf",res.hk); //this was for diagnostics of hardklor correlation results (or lack of)
      fprintf(fOut,"\t%.4f",res.rTime);
      fprintf(fOut,"\t%.4lf",res.obsMass);
      fprintf(fOut,"\t%d",res.charge);
      fprintf(fOut,"\t%.4lf",res.psmMass);
      fprintf(fOut,"\t%.4lf",res.ppm);
      fprintf(fOut,"\t%.4lf",res.score);
      fprintf(fOut,"\t%.4lf",res.scoreDelta);
      fprintf(fOut,"\t%.3e",res.eVal);
      //fprintf(fOut,"\t%.4lf",res.scorePepDif);
      if (res.scoreA == 0)fprintf(fOut, "\t%.4lf", res.score);
      else fprintf(fOut,"\t%.4lf",res.scoreA);
      fprintf(fOut,"\t%.3e",res.eVal1);
      fprintf(fOut,"\t%s",&res.modPeptide1[0]);
      if(res.n15Pep1) fprintf(fOut,"-15N");
      fprintf(fOut,"\t%d",res.link1);

      //export protein
      fprintf(fOut, "\t%s", res.protein1.c_str());
      /* not sure about this anymore - probably breaking something by removing it
      if(bDupe){
        pep = db.getPeptide(res.pep2);
        for(j=0;j<pep.map->size();j++){
          fprintf(fOut,"%s;",&db[pep.map->at(j).index].name[0]);
          //if(res.link1>=0) fprintf(fOut,"(%d);",pep.map->at(j).start+res.link1); //only non-linked peptides
        }
      }
      */
      if(res.link1>-1) fprintf(fOut,"\t%s",res.protPos1.c_str());
      else fprintf(fOut,"\t-");

      if(res.modPeptide2.size()>1) {
        fprintf(fOut, "\t%.4lf", res.scoreB);
        fprintf(fOut, "\t%.3e", res.eVal2);
        fprintf(fOut,"\t%s",&res.modPeptide2[0]);
        if (res.n15Pep2) fprintf(fOut, "-15N");
        fprintf(fOut,"\t%d",res.link2);
        fprintf(fOut,"\t%s",res.protein2.c_str());
        fprintf(fOut, "\t%s", res.protPos2.c_str());
        if(tmpSC.link>-1)fprintf(fOut,"\t%.4lf",link[tmpSC.link].mass);
        else fprintf(fOut,"\t0");
      } else if(res.link2>-1){
        fprintf(fOut,"\t0\t999\t-\t%d\t-\t%s",res.link2,res.protPos2.c_str());
        fprintf(fOut,"\t%.4lf",link[tmpSC.link].mass);
      } else {
        fprintf(fOut,"\t0\t999\t-\t-1\t-\t-\t0");
      }
      
      fprintf(fOut,"\n");

      if(res.type==2){
        bInter=true;
        pep = db.getPeptide(res.pep1);
        pep2 = db.getPeptide(res.pep2);
        for(j=0;j<pep.map->size();j++){
          for(k=0;k<pep2.map->size();k++){
            if(pep.map->at(j).index==pep2.map->at(k).index){
              bInter=false;
              break;
            }
          }
          if(!bInter) break;
        }
      }

      if(params->exportMzID){
        outputMzID(mzID,db,par,res);
      }
      
      if(params->exportPercolator) {
        switch(res.type){
          case 1:   outputPercolator(fLoop,db,res,count);   break;
          case 2:
            if(bInter)  outputPercolator(fInter,db,res,count);
            else        outputPercolator(fIntra,db,res,count);
            break;
          case 3:   outputPercolator(fDimer, db, res, count);  break;
          default:  outputPercolator(fSingle,db,res,count); break;
        }
      }

      if(params->exportPepXML){
        outputNeoPepXML(sq,db,res);
      }

      //Get the next entry - it must also be exported if it has the same score
      if(bDupe) scoreIndex+=iDupe;
      else scoreIndex++;
      if(scoreIndex>=20) break;
      tmpSC=spec[i].getScoreCard(scoreIndex);
    }

    if(params->exportPepXML) {
      sq.index=(int)p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.size()+1;
      p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.push_back(sq);
      //p.writeSpectrumQuery(sq);
    }

  }

  fclose(fOut);
  if(params->exportPercolator) {
    fclose(fIntra);
    fclose(fInter);
    fclose(fLoop);
    fclose(fSingle);
    if (params->dimers) fclose(fDimer);
  }
  if(params->exportPepXML){
    p.write(pepXML_fileName.c_str(),true);
    //p.closePepXML();
  }
  if (params->exportMzID){
    sprintf(fName, "%s.mzid", params->outFile);
    mzID.writeFile(fName);
  }
  if (fDiag != NULL){
    fprintf(fDiag, "</kojak_analysis>\n");
    fclose(fDiag);
  }

  return true;

}

void KData::readLinkers(char* fn){
  FILE* f;
  kLinker k;
  int ret;

  cout << "Reading linker file...";

  f=fopen(fn,"rt");
  while(!feof(f)){
    ret=fscanf(f,"%lf\t%d\n",&k.mass,&k.mono);
    link.push_back(k);
  }
  fclose(f);
  cout << "Done!" << endl;

  for(unsigned int i=0;i<link.size();i++) cout << "Linker: " << link[i].mass << " is " << link[i].mono << endl;
}

//Reads in raw/mzXML/mzML files. Other formats supported in MSToolkit as well.
bool KData::readSpectra(){

  MSReader   msr;
  Spectrum   s;
  Spectrum   c;
  KSpectrum  pls(params->topCount,params->binSize,params->binOffset);
  kSpecPoint sp;
  float      max;
  kPrecursor pre;

  int totalScans=0;
  int totalPeaks=0;
  int collapsedPeaks=0;
  int iPercent=0;
  int iTmp;

  int i;
  int j;

  char nStr[256];
  string sStr;

  bool doCentroid;

  spec.clear();
  msr.setFilter(MS2);

  //Set progress meter
  printf("%2d%%", iPercent);
  fflush(stdout);

  prof.Init();

  if(!msr.readFile(params->msFile,s)) return false;
  while(s.getScanNumber()>0){

    totalScans++;
    if(s.size()<1) {
      int64 pID=prof.StartTimer("ReadFile");
      msr.readFile(NULL,s);
      prof.StopTimer(pID);
      continue;
    }

    int64 pID=prof.StartTimer("ProcessSpectrum");
    pls.clear();
    pls.setRTime(s.getRTime());
    pls.setScanNumber(s.getScanNumber());
    s.getNativeID(nStr,256);
    sStr=nStr;
    pls.setNativeID(sStr);
    max=0;

    doCentroid=false;
    switch(s.getCentroidStatus()){
    case 0:
      if(params->ms2Centroid) {
        char tmpStr[256];
        sprintf(tmpStr,"Kojak parameter indicates MS/MS data are centroid, but spectrum %d labeled as profile.",s.getScanNumber());
        klog->addError(string(tmpStr));
      } else doCentroid=true;
      break;
    case 1:
      if (!params->ms2Centroid) {
        klog->addWarning(0, "Spectrum is labeled as centroid, but Kojak parameter indicates data are profile. Ignoring Kojak parameter.");
      }
      break;
    default:
      if(!params->ms2Centroid) doCentroid=true;
      break;
    }

    //If not centroided, do so now.
    if(doCentroid){
      centroid(s, c, params->ms2Resolution, params->instrument);
      totalPeaks += c.size();
    } else c=s;

    //remove precursor if requested
    if(params->removePrecursor>0){
      double pMin=s.getMZ()-params->removePrecursor;
      double pMax=s.getMZ()+params->removePrecursor;
      for(i=0;i<c.size();i++){
        if(c[i].mz>pMin && c[i].mz<pMax) c[i].intensity=0;
      }
    }

    //Collapse the isotope peaks
    if (params->specProcess == 1 && c.size()>1) {
      collapseSpectrum(c);
      collapsedPeaks += c.size();
    }

    //If user limits number of peaks to analyze, sort by intensity and take top N
    if (params->maxPeaks>0){
      if (c.size()>1) c.sortIntensityRev();
      if (c.size()<params->maxPeaks) j = c.size();
      else j = params->maxPeaks;
    } else {
      j = c.size();
    }
    for (i = 0; i<j; i++){
      sp.mass = c[i].mz;
      sp.intensity = c[i].intensity;
      pls.addPoint(sp);
      if (sp.intensity>max) max = sp.intensity;
    }
    pls.setMaxIntensity(max);

    //Sort again by MZ, if needed
    if (pls.size()>1 && params->maxPeaks>0) pls.sortMZ();

    //Get any additional information user requested
    pls.setCharge(s.getCharge());
    pls.setMZ(s.getMZ());
    if(params->preferPrecursor>0){
      if(s.getMonoMZ()>0 && s.getCharge()>0){
        pre.monoMass=s.getMonoMZ()*s.getCharge()-s.getCharge()*1.007276466;
        pre.charge=s.getCharge();
        pre.corr=0;
        pls.addPrecursor(pre,params->topCount);
        for(int px=1;px<=params->isotopeError;px++){
          if(px==4) break;
          pre.monoMass -= 1.00335483;
          pre.corr -= 0.1;
          pls.addPrecursor(pre, params->topCount);
         }
        pls.setInstrumentPrecursor(true);
      }
    }
    prof.StopTimer(pID);

    //Add spectrum (if it has enough data points) to data object and read next file
    pID=prof.StartTimer("Add2Array");
    if(pls.size()>params->minPeaks) spec.push_back(pls);
    prof.StopTimer(pID);

    /*
    for(unsigned int d=0;d<params->diag->size();d++){
      if(pls.getScanNumber()==params->diag->at(d)){
        char diagStr[256];
        sprintf(diagStr,"diagnostic_spectrum_%d.txt",params->diag->at(d));
        FILE* f=fopen(diagStr,"wt");
        fprintf(f,"Scan: %d\t%d\n",pls.getScanNumber(),pls.size());
        for(int k=0;k<pls.size();k++) fprintf(f,"%.6lf\t%.0f\n",pls[k].mass,pls[k].intensity);
        fclose(f);
        break;
      }
    }
    */

    //Update progress meter
    iTmp = msr.getPercent();
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }
    pID = prof.StartTimer("ReadFile");
    msr.readFile(NULL,s);
    prof.StopTimer(pID);
  }

  //Finalize progress meter
  if(iPercent<100) printf("\b\b\b100%%");
  cout << endl;

  prof.Release();

  cout << "  " << spec.size() << " total spectra have enough data points (" << params->minPeaks << " peaks) for searching." << endl;
  //cout << totalScans << " total scans were loaded." <<  endl;
  //cout << totalPeaks << " total peaks in original data." << endl;
  //cout << collapsedPeaks << " peaks after collapsing." << endl;
  //cout << finalPeaks << " peaks after top N." << endl;
	return true;
}

void KData::setLinker(kLinker x){
  if(x.mono==0) link.push_back(x);
}

void KData::setLog(KLog* c){
  klog=c;
}

void KData::setVersion(const char* v){
  strcpy(version,v);
}

int KData::size(){
  return (int)spec.size();
}

int KData::sizeLink(){
  return (int)link.size();
}

void KData::xCorr(bool b){
  if(b) {
    klog->addMessage("Using XCorr scores.",true);
    cout << "  Using XCorr scores." << endl;
  } else  {
    klog->addMessage("Using Kojak modified XCorr scores.",true);
    cout << "  Using Kojak modified XCorr scores." << endl;
  }

  klog->addMessage("Transforming spectra.",true);
  cout << "  Transforming spectra ... ";
  int iTmp;
  int iPercent = 0;
  printf("%2d%%", iPercent);
  fflush(stdout);
  for(size_t i=0;i<spec.size();i++) {
    spec[i].xCorrScore(b);

    //Update progress meter
    iTmp = (int)((double)i / spec.size() * 100);
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }

  }

  //Finalize progress meter
  if(iPercent<100) printf("\b\b\b100%%");
  cout << endl;
}

/*============================
  Private Utilities
============================*/
//First derivative method, returns base peak intensity of the set
void KData::centroid(Spectrum& s, Spectrum& out, double resolution, int instrument){
  int i,j;
  float maxIntensity;
  int bestPeak;
  bool bLastPos;

	int nextBest;
	double FWHM;
  double maxMZ = s[s.size()-1].mz+1.0;
	Peak_T centroid;

	vector<double> x;
	vector<double> y;
	vector<double> c;
	int left, right;
	bool bPoly;
	float lastIntensity;

	out.clear();

  bLastPos=false;
	for(i=0;i<s.size()-1;i++){

    if(s[i].intensity<s[i+1].intensity) {
      bLastPos=true;
      continue;
    } else {
      if(bLastPos){
				bLastPos=false;

				//find max and add peak
				maxIntensity=0;
				for(j=i;j<i+1;j++){
				  if (s[j].intensity>maxIntensity){
				    maxIntensity=s[j].intensity;
				    bestPeak = j;
				  }
				}

				//walk left and right to find bounds above half max
				left=right=bestPeak;
				lastIntensity=maxIntensity;
				for(left=bestPeak-1;left>0;left--){
					if(s[left].intensity<(maxIntensity/3) || s[left].intensity>lastIntensity){
						left++;
						break;
					}
					lastIntensity=s[left].intensity;
				}
				lastIntensity=maxIntensity;
				for(right=bestPeak+1;right<s.size()-1;right++){
					if(s[right].intensity<(maxIntensity/3) || s[right].intensity>lastIntensity){
						right--;
						break;
					}
					lastIntensity=s[right].intensity;
				}

				//if we have at least 5 data points, try polynomial fit
				double r2;
				bPoly=false;
				if((right-left+1)>4){
					x.clear();
					y.clear();
					for(j=left;j<=right;j++){
						x.push_back(s[j].mz);
						y.push_back(log(s[j].intensity));
					}
					r2=polynomialBestFit(x,y,c);
					if(r2>0.95){
						bPoly=true;
						centroid.mz=-c[1]/(2*c[2])+c[3];
						centroid.intensity=(float)exp(c[0]-c[2]*(c[1]/(2*c[2]))*(c[1]/(2*c[2])));
					} else {

					}
				}

				if(!bPoly){
					//Best estimate of Gaussian centroid
					//Get 2nd highest point of peak
					if(bestPeak==s.size()) nextBest=bestPeak-1;
					else if(s[bestPeak-1].intensity > s[bestPeak+1].intensity) nextBest=bestPeak-1;
					else nextBest=bestPeak+1;

					//Get FWHM
					switch(instrument){
						case 0: FWHM = s[bestPeak].mz*sqrt(s[bestPeak].mz)/(20*resolution); break;  //Orbitrap
						case 1: FWHM = s[bestPeak].mz*s[bestPeak].mz/(400*resolution); break;				//FTICR
						default: break;
					}

					//Calc centroid MZ (in three lines for easy reading)
					centroid.mz = pow(FWHM,2)*log(s[bestPeak].intensity/s[nextBest].intensity);
					centroid.mz /= GAUSSCONST*(s[bestPeak].mz-s[nextBest].mz);
					centroid.mz += (s[bestPeak].mz+s[nextBest].mz)/2;

					//Calc centroid intensity
					centroid.intensity=(float)(s[bestPeak].intensity/exp(-pow((s[bestPeak].mz-centroid.mz)/FWHM,2)*GAUSSCONST));
				}

				//some peaks are funny shaped and have bad gaussian fit.
				//if error is more than 10%, keep existing intensity
				if( fabs((s[bestPeak].intensity - centroid.intensity) / centroid.intensity * 100) > 10 ||
            //not a good check for infinity
            centroid.intensity>9999999999999.9 ||
            centroid.intensity < 0 ) {
					centroid.intensity=s[bestPeak].intensity;
				}

				//Hack until I put in mass ranges
				if(centroid.mz<0 || centroid.mz>maxMZ) {
					//do nothing if invalid mz
				} else {
					out.add(centroid);
				}
			
      }

    }
  }

}

//Function tries to remove isotopes of signals by stacking the intensities on the monoisotopic peak
//Also creates an equal n+1 peak in case wrong monoisotopic peak was identified.
void KData::collapseSpectrum(Spectrum& s){
  int i,j,k,n;
  int charge,z;
  int maxIndex;
  float max;
  float cutoff;
  vector<int> dist;

  Spectrum s2;

  while(true){
    max=0.1f;
    for(i=0;i<s.size();i++){
      if(s[i].intensity>max){
        max=s[i].intensity;
        maxIndex=i;
      }
    }

    //finish and exit function
    if(max<1) break;

    dist.clear();
    dist.push_back(maxIndex);

    //check right
    j=maxIndex+1;
    while(j<s.size() && (s[j].mz-s[maxIndex].mz)<1.1){
      if(s[j].intensity<1) {
        j++;
        continue;
      }
      charge=getCharge(s,maxIndex,j);

      if(charge==0){
        j++;
        continue;
      }

      //try stepping along at same charge state here out
      //note that if this doesn't work, it doesn't go back and look for a different charge state
      dist.push_back(j);
      k=j;
      n=j+1;
      while(n<s.size() && (s[n].mz-s[k].mz)<1.1){
        if(s[n].intensity<1) {
          n++;
          continue;
        }
        z=getCharge(s,k,n);
        if(z>0 && z<charge) {
          break;
        } else if(z==charge && (s[n].mz-s[k].mz)>(0.99/charge) && (s[n].mz-s[k].mz)<(1.0041/charge)) {
          dist.push_back(n);
          k=n;
          n++;
        } else {
          n++;
        }
      }
      break;
    }

    //if nothing found to the right, quit here?
    if(dist.size()==1){
      s2.add(s[dist[0]]);
      s[dist[0]].intensity=0;
      continue;
    }

    //step to the left
    j=maxIndex-1;
    while(j>=0 && (s[maxIndex].mz-s[j].mz)<1.1){
      if(s[j].intensity<1) {
        j--;
        continue;
      }
      z=getCharge(s,j,maxIndex);
      if(z!=charge){
        j--;
        continue;
      }

      //try stepping along at same charge state here out
      dist.push_back(j);
      k=j;
      n=j-1;
      while(n>=0 && (s[k].mz-s[n].mz)<1.1){
        if(s[n].intensity<1) {
          n--;
          continue;
        }
        z=getCharge(s,n,k);
        //printf("\tleft\t%.6lf\t%.6lf\t%d\n",s[n].mz,s[k].mz-s[n].mz,z);
        if(z>0 && z<charge) {
          break;
        } else if(z==charge && s[k].mz-s[n].mz > 0.99/charge && s[k].mz-s[n].mz < 1.0041/charge) {
          dist.push_back(n);
          k=n;
          n--;
        } else {
          n--;
        }
      }
      break;
    }


    //Only accept size of 2 if charge is 1 or 2
    if(dist.size()==2){
      if(charge<3){
        max=s[dist[0]].intensity+s[dist[1]].intensity;
        s2.add(s[dist[0]].mz,max);
        s[dist[1]].intensity=0;
       // s2.add(s[dist[1]].mz,max);
      } else {
        s2.add(s[dist[0]]);
       // s2.add(s[dist[1]]);
      }
      s[dist[0]].intensity=0;
      //s[dist[1]].intensity=0;
    } else {
      cutoff=max/20;
      max=0;
      j=dist[0];
      k=dist[1];
      for(i=0;i<(int)dist.size();i++) {
        if(dist[i]<j && s[dist[i]].intensity>cutoff){
          k=j;
          j=dist[i];
        }
        if(s[dist[i]].intensity>cutoff){
          max+=s[dist[i]].intensity;
          s[dist[i]].intensity=0;
        }
      }
      s2.add(s[j].mz,max);
      //s2.add(s[k].mz,max);
    }

  }

  s2.sortMZ();
  s.clearPeaks();
  for(i=0;i<s2.size();i++) {
    if(i<s2.size()-1 && s2[i].mz==s2[i+1].mz){
      if(s2[i].intensity>s2[i+1].intensity) s.add(s2[i]);
      else s.add(s2[i+1]);
      i++;
    } else {
      s.add(s2[i]);
    }
  }
  
}

int KData::compareInt(const void *p1, const void *p2){
  int d1 = *(int *)p1;
  int d2 = *(int *)p2;
  if(d1<d2) {
		return -1;
	} else if(d1>d2) {
  	return 1;
  } else {
	  return 0;
  }
}

int KData::compareMassList(const void *p1, const void *p2){
  kMass d1 = *(kMass *)p1;
  kMass d2 = *(kMass *)p2;
  if(d1.mass<d2.mass) {
		return -1;
	} else if(d1.mass>d2.mass) {
  	return 1;
  } else {
	  return 0;
  }
}

int KData::getCharge(Spectrum& s, int index, int next){
  double mass;

  mass=s[next].mz-s[index].mz;
  if(mass>0.99 && mass<1.007) return 1;
  else if(mass>0.495 && mass<0.5035) return 2;
  else if(mass>0.33 && mass<0.335667) return 3;
  else if(mass>0.2475 && mass<0.25175) return 4;
  else if(mass>0.198 && mass<0.2014) return 5;
  else if(mass>0.165 && mass<0.1678333) return 6;
  else return 0;

}

double KData::polynomialBestFit(vector<double>& x, vector<double>& y, vector<double>& coeff, int degree){
	if(degree>3){
		cout << "High order polynomials not supported with this function. Max=3" << endl;
		exit(1);
	}

	if(degree<2){
		cout << "Polynomials need at least two degrees. Min=2" << endl;
		exit(1);
	}

	int i,j,a;
	int n=(int)x.size();
	degree++;

	double sFactor=x[n/2];

	//set X matrix
	double** X = new double* [n];
	for(i=0;i<n;i++){
		X[i] = new double [degree];
		X[i][0] = 1.0;
		for(j=1;j<degree;j++) X[i][j]=X[i][j-1]*(x[i]-sFactor);
	}

	//make transpose of X
	double** Xt = new double* [degree];
	for(j=0;j<degree;j++) Xt[j] = new double [n];
	for(i=0;i<n;i++){
		for(j=0;j<degree;j++){
			Xt[j][i] = X[i][j];
		}
	}

	//matrix multiplication
	double** XtX = new double* [degree];
	for(i=0;i<degree;i++){
		XtX[i] = new double [degree];
		for(j=0;j<degree;j++){
			XtX[i][j]=0;
			for(a=0;a<n;a++) XtX[i][j]+=(Xt[i][a]*X[a][j]);
		}
	}

	//inverse using Gauss-Jordan Elimination
	double** XtXi = new double* [degree];
	for(i=0;i<degree;i++){
		XtXi[i] = new double [degree*2];
		for(j=0;j<degree*2;j++){
			if(j<degree) XtXi[i][j]=XtX[i][j];
			else if(j-degree==i) XtXi[i][j]=1;
			else XtXi[i][j]=0;
		}
	}
	double d;
	for(j=0;j<degree;j++){
		for(i=0;i<degree;i++){
			if(i==j) continue;
			if(XtXi[i][j]==0) continue;
			d=-XtXi[i][j]/XtXi[j][j];
			for(a=0;a<degree*2;a++) XtXi[i][a]+=(d*XtXi[j][a]);
		}
	}
	for(i=0;i<degree;i++){
		d=1/XtXi[i][i];
		for(j=0;j<degree*2;j++) XtXi[i][j]*=d;
	}

	//matrix multiplication
	double* Xty = new double [degree];
	for(i=0;i<degree;i++){
		Xty[i]=0;
		for(j=0;j<n;j++) Xty[i]+=Xt[i][j]*y[j];	
	}

	//matrix multiplication
	double* c = new double [degree];
	for(i=0;i<degree;i++){
		c[i]=0;
		for(j=0;j<degree;j++) c[i]+=XtXi[i][j+degree]*Xty[j];	
	}

	coeff.clear();
	for(i=0;i<degree;i++) {
		coeff.push_back(c[i]);
	}
	coeff.push_back(sFactor);

	vector<double> z;
	for(i=0;i<n;i++) z.push_back((x[i]-sFactor)*(x[i]-sFactor)*c[2]+(x[i]-sFactor)*c[1]+c[0]);

	//clean up memory
	delete [] c;
	delete [] Xty;
	for(i=0;i<degree;i++){
		delete [] XtXi[i];
		delete [] XtX[i];
		delete [] Xt[i];
	}
	for(i=0;i<n;i++) delete [] X[i];
	delete [] X;
	delete [] Xt;
	delete [] XtX;
	delete [] XtXi;

	double sxy=0;
  double sxx=0;
  double syy=0;
	double xavg=0;
	double yavg=0;
  for(i=0;i<n;i++){
		xavg+=z[i];
		yavg+=y[i];
	}
	xavg/=n;
	yavg/=n;
  for(i=0;i<n;i++){
    sxy += ((z[i]-xavg)*(y[i]-yavg));
    sxx += ((z[i]-xavg)*(z[i]-xavg));
    syy += ((y[i]-yavg)*(y[i]-yavg));
  }
	double r2 = (sxy*sxy)/(sxx*syy);

	return r2;

}

//Takes relative path and finds absolute path
bool KData::processPath(const char* in_path, char* out_path){
  char cwd[1024];
  if(getcwd(cwd,1024)==NULL) return false; //stop if failed to get CWD

  //if windows or unix in_path, just copy it to out_path
  if (strlen(in_path) > 0 && in_path[0] == '/'){ //unix
    strcpy(out_path,in_path);
    return true;
  }
  if (strlen(in_path) > 1 && in_path[1] == ':'){ //windows
    strcpy(out_path, in_path);
    return true;
  }

  //tokenize cwd
  char* tok;
  char str[1024];
  strcpy(str,cwd);
  string s;
  vector<string> v;

  tok = strtok(str, "\\/\n\r");
  while (tok != NULL){
    s=tok;
    v.push_back(s);
    tok = strtok(NULL, "\\/\n\r");
  }

  //tokenize in_path
  strcpy(str,in_path);

  tok = strtok(str, "\\/\n\r");
  while (tok != NULL){
    if (strcmp(tok, "..") == 0) {
      v.pop_back();
    } else if (strcmp(tok, ".") == 0){
      //do nothing
    } else {
      s=tok;
      v.push_back(s);
    }
    tok = strtok(NULL, "\\/\n\r");
  }

  //build absolute path
#ifdef _MSC_VER
  s.clear();
#else
  s.clear();
  s+=slashdir;
#endif
  for (size_t i = 0; i < v.size(); i++){
    s += v[i];
    s += slashdir;
  }
  s[s.size() - 1] = '\0';
  strcpy(out_path, &s[0]);
  return true;

}

string KData::processPeptide(kPeptide& pep, vector<kPepMod>* mod, KDatabase& db){
  char tmp[32];
  size_t j,k;
  string seq = "";
  string peptide;

  db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, peptide);

  if (pep.nTerm && aa.getFixedModMass('$') != 0) {
    sprintf(tmp, "n[%.2lf]", aa.getFixedModMass('$'));
    seq += tmp;
  }
  for (k = 0; k<mod->size(); k++){ //check for n-terminal peptide mod
    if (mod->at(k).pos == -1){
      sprintf(tmp, "n[%.2lf]", mod->at(k).mass);
      seq += tmp;
    }
  }
  for (j = 0; j<peptide.size(); j++) {
    seq += peptide[j];
    for (k = 0; k<mod->size(); k++){
      if(mod->at(k).pos<0) continue;
      if (j == (size_t)mod->at(k).pos){
        sprintf(tmp, "[%.2lf]", mod->at(k).mass);
        seq += tmp;
      }
    }
  }
  for (k = 0; k<mod->size(); k++){ //check for c-terminal peptide mod
    if (mod->at(k).pos ==-2){
      sprintf(tmp, "c[%.2lf]", mod->at(k).mass);
      seq += tmp;
    }
  }
  if (pep.cTerm && aa.getFixedModMass('%') != 0) {
    sprintf(tmp, "c[%.2lf]", aa.getFixedModMass('%'));
    seq += tmp;
  }

  return seq;
}

void KData::processProtein(int pepIndex, int site, char linkSite, string& prot, string& sites, bool& decoy, KDatabase& db){

  size_t j;
  kPeptide pep;
  char tmp[16];

  //automatically set decoyness to true; remains so until first non-decoy peptide is found
  decoy=true;

  //export protein
  pep = db.getPeptide(pepIndex);
  string peps;
  db.getPeptideSeq(pepIndex,peps);
  prot.clear();
  sites.clear();
  for (j = 0; j<pep.map->size(); j++){

    //for linkage to n- or c- termini, skip protein if not those things
    if(linkSite=='n' && pep.map->at(j).start>1) continue;
    if (linkSite == 'c' && pep.map->at(j).stop < db[pep.map->at(j).index].sequence.size()-1) continue;

    if(prot.size()>0) prot+=";"; //add spacer if appending a prior protein
    prot+=db[pep.map->at(j).index].name;

    if(site>-1){//add the protein site location
      if(sites.size()>0) sites+=";";
      sprintf(tmp, "%d", pep.map->at(j).start + site+1); 
      sites+=tmp;
    }

    //determine if target (if it is currently still decoy)
    if(decoy){
      if (db[pep.map->at(j).index].name.find(params->decoy) == string::npos) decoy=false;
    }

  }

}

void KData::writeMzIDDatabase(CMzIdentML& m, KDatabase& db){

  char outPath[1056];
  processPath(params->dbFile, outPath);
  string sDB = outPath;
  CSearchDatabase* m_db = m.dataCollection.inputs.addSearchDatabase(sDB);

  for(size_t a=0;a<db.getProteinDBSize();a++){
    CDBSequence m_dbs;
    sCvParam cv;

    string pName;
    string pDesc;
    if (db[a].name.find(' ') == string::npos){
      pName = db[a].name;
      pDesc.clear();
    } else {
      pName = db[a].name.substr(0, db[a].name.find(' '));
      pDesc = db[a].name.substr(db[a].name.find(' '), db[a].name.size());
    }

    m_dbs.accession = pName;
    m_dbs.searchDatabaseRef = m_db->id;
    sSeq ss;
    ss.text = db[a].sequence;
    m_dbs.seq.push_back(ss);
   
    if(!pDesc.empty()){
      cv.cvRef = "PSI-MS";
      cv.accession = "MS:1001088";
      cv.name = "protein description";
      cv.value = pDesc;
      m_dbs.cvParam.push_back(cv);
    }
    cv.cvRef = "PSI-MS";
    cv.accession = "MS:1001344";
    cv.name = "AA sequence";
    cv.value.clear();
    m_dbs.cvParam.push_back(cv);
    m.sequenceCollection.addDBSequence(m_dbs);
  }
}

bool KData::writeMzIDEnzyme(pxwBasicXMLTag t, CEnzymes& e){
  char str[256];
  char* tok;
  string valueA,valueB;

  strcpy(str,t.value.c_str());
  tok=strtok(str," \t\n\r");
  valueA=tok;
  tok=strtok(NULL," \t\n\r");
  if(tok!=NULL) valueB=tok;

  if(valueA.compare("[KR]|{P}")==0 || valueA.compare("[RK]|{P}")==0){
    CEnzyme ez;
    ez.id="SIP0_E0";
    ez.name="trypsin";
    ez.missedCleavages=params->miscleave;
    CEnzymeName en;
    sCvParam cv;
    cv.accession="MS:1001251";
    cv.cvRef="PSI-MS";
    cv.name="Trypsin";
    en.cvParam.push_back(cv);
    ez.enzymeName.push_back(en);
    e.enzyme.push_back(ez);
    return true;
  } else if (valueA.compare("[KR]") == 0 || valueA.compare("[RK]") == 0){
    CEnzyme ez;
    ez.id = "SIP0_E0";
    ez.name = "trypsin/p";
    ez.missedCleavages = params->miscleave;
    CEnzymeName en;
    sCvParam cv;
    cv.accession = "MS:1001313";
    cv.cvRef = "PSI-MS";
    cv.name = "Trypsin/P";
    en.cvParam.push_back(cv);
    ez.enzymeName.push_back(en);
    e.enzyme.push_back(ez);
    return true;
  }
  return false;
}

void KData::writeMzIDPE(CMzIdentML& m, CSpectrumIdentificationItem& m_sii, int pepID, KDatabase& db){
  //Add all proteins mapped by this peptide
  kPeptide pep = db.getPeptide(pepID);
  for (size_t i = 0; i<pep.map->size(); i++){
    if (pep.n15 && db[pep.map->at(i).index].name.find(params->n15Label) == string::npos) continue;
    if (!pep.n15 && strlen(params->n15Label)>0 && db[pep.map->at(i).index].name.find(params->n15Label) != string::npos) continue;

    CDBSequence m_dbs;
    if (db[pep.map->at(i).index].name.find(' ') == string::npos){
      m_dbs = m.getDBSequenceByAcc(db[pep.map->at(i).index].name);
    } else {
      string pName = db[pep.map->at(i).index].name.substr(0, db[pep.map->at(i).index].name.find(' '));
      m_dbs = m.getDBSequenceByAcc(pName);
    }
    char pre;
    char post;
    bool isDecoy;
    if (pep.map->at(i).start<1) pre = '-';
    else pre = db[pep.map->at(i).index].sequence[pep.map->at(i).start - 1];
    if ((size_t)pep.map->at(i).stop + 1 == db[pep.map->at(i).index].sequence.size()) post = '-';
    else post = db[pep.map->at(i).index].sequence[(size_t)pep.map->at(i).stop + 1];
    isDecoy = db[pep.map->at(i).index].decoy;

    m_sii.peptideEvidenceRef.push_back(m.addPeptideEvidence(m_dbs.id, m_sii.peptideRef, (int)pep.map->at(i).start + 1, (int)pep.map->at(i).stop + 1, pre, post, isDecoy));
  }
}

std::string KData::writeMzIDSIP(CMzIdentML& m, string& sRef, KParams& par){
  CSpectrumIdentificationProtocol* m_sip = m.analysisProtocolCollection.addSpectrumIdentificationProtocol(sRef);

  sCvParam cv;
  cv.accession="MS:1001494";
  cv.cvRef="PSI-MS";
  cv.name="no threshold";
  m_sip->threshold.cvParam.push_back(cv);

  //populate analysis software & protocol information if it is new
  cv.accession = "MS:1001083";
  cv.cvRef = "PSI-MS";
  cv.name = "ms-ms search";
  m_sip->searchType.cvParam=cv;

  //special case for cross-linking
  CAdditionalSearchParams m_asp;
  cv.accession = "MS:1002494";
  cv.cvRef = "PSI-MS";
  cv.name = "cross-linking search";
  m_asp.cvParam.push_back(cv);

  CModificationParams m_mp; //what if there are no modifications in the search?
  size_t i;
  string cStr;
  bool bTerm=false;
  char site;
  for(i=0;i<params->mods->size();i++){
    cStr.clear();
    site=(char)params->mods->at(i).index;
    if (site == '$') cStr += 'n';
    else if (site == '%') cStr += 'c';
    else cStr += site;
    if (site == 'n' || site == 'c') bTerm = true;
    m_mp.addSearchModification(false, params->mods->at(i).mass, cStr,bTerm);
  }
  for (i = 0; i<params->fMods->size(); i++){
    cStr.clear();
    site = (char)params->fMods->at(i).index;
    if (site == '$') cStr += 'n';
    else if (site == '%') cStr += 'c';
    else cStr += site;
    if (site == 'n' || site == 'c') bTerm = true;
    m_mp.addSearchModification(true, params->fMods->at(i).mass, cStr, bTerm);
  }
  
  vector<string> tokens;
  char str[1024];
  char* tok;
  bool nTerm;
  bool cTerm;
  for (i = 0; i<par.xmlParams.size(); i++){ //figure out how to write all parameters
    sUserParam u;
    u.name = par.xmlParams[i].name;
    u.value = par.xmlParams[i].value;
    m_asp.userParam.push_back(u);
    if (par.xmlParams[i].name.compare("cross_link") == 0){
      tokens.clear();
      strcpy(str, par.xmlParams[i].value.c_str());
      tok = strtok(str, " \t\n\r");
      while (tok != NULL){
        cStr = tok;
        tokens.push_back(cStr);
        tok = strtok(NULL, " \t\n\r");
      }
      m_mp.addSearchModificationXL(atof(tokens[2].c_str()), tokens[0], tokens[1]);
    } else if(par.xmlParams[i].name.compare("enzyme")==0){
      CEnzymes m_e;
      if(writeMzIDEnzyme(par.xmlParams[i],m_e)){
        m_sip->enzymes.push_back(m_e);
      }
    }
  }
  
  m_sip->modificationParams.push_back(m_mp);
  m_sip->additionalSearchParams.push_back(m_asp);

  //  m.consolidateSpectrumIdentificationProtocol();
  //  sip = m.getSpectrumIdentificationProtocol(si->spectrumIdentificationProtocolRef);
  //}

  return m_sip->id;
}



