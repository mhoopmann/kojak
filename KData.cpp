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

using namespace std;
using namespace MSToolkit;

Mutex KData::mutexMemoryPool;
bool* KData::memoryPool;
double** KData::tempRawData;
double** KData::tmpFastXcorrData;
float**  KData::fastXcorrData;
kPreprocessStruct** KData::preProcess;
kParams* KData::params;

KLog* KData::klog;

deque<Spectrum*> KData::dMS1;
vector<Spectrum*> KData::vMS1Buffer;
Mutex KData::mutexLockMS1;
CHardklor2** KData::h;
CHardklor** KData::hO;
CAveragine** KData::averagine;
CMercury8** KData::mercury;
CModelLibrary* KData::models;
Mutex* KData::mutexHardklor;
CHardklorSetting KData::hs;
CHardklorSetting KData::hs2;
CHardklorSetting KData::hs4;
bool* KData::bHardklor;

int KData::maxPrecursorMass;

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
  for(size_t a=0;a<spec.size();a++) delete spec[a];
}


/*============================
  Operators
============================*/
KSpectrum& KData::operator [](const int& i){
  return *spec[i];
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
      if (db[pep.map->at(j).index].name[i]==' ') break;
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

  ////Output for diagnostic purposes
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

void KData::convertToDiag(kScoreCard& sc, kDiag& d, KDatabase& db){
  size_t i, x;
  char strs[256];
  char st[32];
  string pep1, pep2, tmp;
  kPeptide pep;

  pep = db.getPeptide(sc.pep1);
  db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, strs);
  pep1.clear();
  if (pep.nTerm && aa.getFixedModMass('$') != 0) {
    sprintf(st, "[%.2lf]", aa.getFixedModMass('$'));
    pep1 += st;
  }
  for (i = 0; i < strlen(strs); i++) {
    pep1 += strs[i];
    for (x = 0; x < sc.mods1.size(); x++) {
      if (sc.mods1[x].pos == (char)i) {
        sprintf(st, "[%.2lf]", sc.mods1[x].mass);
        pep1 += st;
      }
    }
    if ((int)i == sc.k1 || (sc.pep2 < 0 && (int)i == sc.k2)) pep1 += "[x]";
  }
  if (pep.cTerm && aa.getFixedModMass('%') != 0) {
    sprintf(st, "[%.2lf]", aa.getFixedModMass('%'));
    pep1 += st;
  }

  pep2.clear();
  if (sc.pep2 > -1) {
    pep = db.getPeptide(sc.pep2);
    db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, strs);
    if (pep.nTerm && aa.getFixedModMass('$') != 0) {
      sprintf(st, "[%.2lf]", aa.getFixedModMass('$'));
      pep2 += st;
    }
    for (i = 0; i < strlen(strs); i++) {
      pep2 += strs[i];
      for (x = 0; x < sc.mods2.size(); x++) {
        if (sc.mods2[x].pos == (char)i) {
          sprintf(st, "[%.2lf]", sc.mods2[x].mass);
          pep2 += st;
        }
      }
      if ((int)i == sc.k2) pep2 += "[x]";
    }
    if (pep.cTerm && aa.getFixedModMass('%') != 0) {
      sprintf(st, "[%.2lf]", aa.getFixedModMass('%'));
      pep2 += st;
    }
  } else {
    //fprintf(f2, "\t(%d)", spec[i].getScoreCard(j).k2);
  }

  if (pep2.size() > 0 && pep1.compare(pep2) > 0) {
    tmp = pep1;
    pep1 = pep2;
    pep2 = tmp;
  }

  if (sc.link == -2) pep1 += "+";
  else if (pep2.size() > 0) pep1 += "--";
  pep1 += pep2;
  d.sequence=pep1;
  d.simpleScore= sc.simpleScore;
  d.cpScore= sc.cpScore1 + sc.cpScore2;
  d.evalue= sc.eVal;
  d.mass= sc.mass;
  if (sc.link < 0) d.xlMass=0;
  else d.xlMass=link[sc.link].mass;
}

//Returns closet mz value to desired point
int KData::findPeak(Spectrum* s, double mass){
  int sz = s->size();
  int lower = 0;
  int mid = sz / 2;
  int upper = sz;


  //binary search to closest mass
  while (s->at(mid).mz != mass){
    if (lower >= upper) break;
    if (mass<s->at(mid).mz){
      upper = mid - 1;
      mid = (lower + upper) / 2;
    } else {
      lower = mid + 1;
      mid = (lower + upper) / 2;
    }
    if (mid == sz) {
      mid--;
      break;
    }
  }

  //Check that mass is closest
  if (mid>0 && fabs(s->at(mid - 1).mz - mass)<fabs(s->at(mid).mz - mass)) return mid - 1;
  if (mid<s->size() - 1 && fabs(s->at(mid + 1).mz - mass)<fabs(s->at(mid).mz - mass)) return mid + 1;
  return mid;

}

//Returns point within precision or -1 if doesn't exist
int KData::findPeak(Spectrum* s, double mass, double prec){
  int sz = s->size();
  int lower = 0;
  int mid = sz / 2;
  int upper = sz;

  double minMass = mass - (mass / 1000000 * prec);
  double maxMass = mass + (mass / 1000000 * prec);

  //binary search to closest mass
  while (s->at(mid).mz<minMass || s->at(mid).mz>maxMass){
    if (lower >= upper) break;
    if (mass<s->at(mid).mz){
      upper = mid - 1;
      mid = (lower + upper) / 2;
    } else {
      lower = mid + 1;
      mid = (lower + upper) / 2;
    }
    if (mid == sz) {
      mid--;
      break;
    }
  }

  //Check that mass is correct
  if (s->at(mid).mz>minMass && s->at(mid).mz<maxMass) return mid;

  return -1;
}

void KData::formatMS2(MSToolkit::Spectrum* s, KSpectrum* pls){
  char nStr[256];
  string sStr;

  pls->setRTime(s->getRTime());
  pls->setScanNumber(s->getScanNumber());
  s->getNativeID(nStr,256);
  sStr=nStr;
  pls->setNativeID(sStr);

  bool doCentroid=false;
  switch(s->getCentroidStatus()){
  case 0:
    if(params->ms2Centroid) {
      char tmpStr[256];
      sprintf(tmpStr,"Kojak parameter indicates MS/MS data are centroid, but spectrum %d labeled as profile.",s->getScanNumber());
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
  int totalPeaks=0;
  if(doCentroid){
    centroid(s, pls, params->ms2Resolution, params->instrument);
    totalPeaks += pls->size();
  } else {
    kSpecPoint sp;
    for(int i=0;i<s->size();i++){
      sp.mass=s->at(i).mz;
      sp.intensity=s->at(i).intensity;
      pls->addPoint(sp);
    }
  }

  //remove precursor if requested
  if(params->removePrecursor>0){
    double pMin=s->getMZ()-params->removePrecursor;
    double pMax=s->getMZ()+params->removePrecursor;
    for(int i=0;i<pls->size();i++){
      if((*pls)[i].mass>pMin && (*pls)[i].mass<pMax) (*pls)[i].intensity=0;
    }
  }

  //Collapse the isotope peaks
  int collapsedPeaks=0;
  if (params->specProcess == 1 && pls->size()>1) {
    collapseSpectrum(*pls);
    collapsedPeaks += pls->size();
  }

  //If user limits number of peaks to analyze, sort by intensity and take top N
  if (params->maxPeaks>0 && pls->size()>params->maxPeaks){
    if (pls->size()>1) pls->sortIntensityRev();
    vector<kSpecPoint> v;
    for (int i = 0; i<params->maxPeaks; i++) v.push_back((*pls)[i]);
    pls->clear();
    for (int i = 0; i<params->maxPeaks; i++) pls->addPoint(v[i]);
    pls->sortMZ();
    pls->setMaxIntensity(v[0].intensity);
  } else {
    float max=0;
    for (int i = 0; i<pls->size(); i++){
      if((*pls)[i].intensity>max) max=(*pls)[i].intensity;
    }
    pls->setMaxIntensity(max);
  }

  //Get any additional information user requested
  pls->setCharge(s->getCharge());
  pls->setMZ(s->getMZ());
  if(params->preferPrecursor>0){
    if(s->getMonoMZ()>0 && s->getCharge()>0){
      kPrecursor pre;
      pre.monoMass=s->getMonoMZ()*s->getCharge()-s->getCharge()*1.007276466;
      pre.charge=s->getCharge();
      pre.corr=0;
      pls->addPrecursor(pre,params->topCount);
      for(int px=1;px<=params->isotopeError;px++){
        if(px==4) break;
        pre.monoMass -= 1.00335483;
        pre.corr -= 0.1;
        pls->addPrecursor(pre, params->topCount);
      }
      pls->setInstrumentPrecursor(true);
    }
  }

}

KSpectrum* KData::getSpectrum(const int& i){
  return spec[i];
}

char** KData::getXLTable(){
  return xlTable;
}

KSpectrum& KData::at(const int& i){
  return *spec[i];
}

//Places centroided peaks in bins, then finds average. Has potential for error when accuracy drifts into neighboring bins.
void KData::averageScansCentroid(vector<Spectrum*>& s, Spectrum& avg, double min, double max){
  unsigned int i;
  int j, k;
  double binWidth;
  double offset = -1.0;
  float* bin;
  float intensity;
  int binCount;
  int* pos;
  double  lowMZ = -1.0;
  kScanBin sb;
  vector<kScanBin> topList;

  avg.clear();

  //if vector is just one scan, then simply copy the peaks
  if (s.size() == 1){
    int index = findPeak(s[0], min);
    if (s[0]->at(index).mz<min) index++;
    while (index<s[0]->size() && s[0]->at(index).mz<max){
      avg.add(s[0]->at(index++));
    }
    return;
  }
  
  pos = new int[s.size()];

  //Set some really small bin width related to the mz region being summed.
  binWidth = 0.0001;

  binCount = (int)((max - min) / binWidth + 1);
  bin = new float[binCount];
  for (j = 0; j<binCount; j++) bin[j] = 0;

  //align all spectra to point closest to min and set offset
  for (i = 0; i<s.size(); i++)  {
    pos[i] = findPeak(s[i], min);
    if (s[i]->at(pos[i]).mz<min) pos[i]++;
    if (offset<0) offset = s[i]->at(pos[i]).mz;
    else if (s[i]->at(pos[i]).mz<offset) offset = s[i]->at(pos[i]).mz;
  }

  //Iterate all spectra and add peaks to bins
  for (i = 0; i<s.size(); i++) {
    while (pos[i]<s[i]->size() && s[i]->at(pos[i]).mz<max){
      j = (int)((s[i]->at(pos[i]).mz - offset) / binWidth);
      if (j<0){
        pos[i]++;
        continue;
      }
      if (j >= binCount) break;
      bin[j] += s[i]->at(pos[i]).intensity;
      pos[i]++;
    }
  }

  //Unsure of current efficiency. Finds bin of tallest peak, then combines with neighboring bins
  //to produce an average. Thus summing of neighboring bins allows flexibility when dealing with mass accuracy drift.
  //Using larger bins has the same effect, but perhaps fewer significant digits in the final result.
  for (j = 0; j<binCount; j++){
    if (bin[j]>0) {
      sb.index = j;
      sb.intensity = bin[j];
      topList.push_back(sb);
    }
  }
  if (topList.size()>0) {
    qsort(&topList[0], topList.size(), sizeof(kScanBin), compareScanBinRev2);
    for (i = 0; i<topList.size(); i++){
      if (bin[topList[i].index] == 0) continue;
      intensity = 0;
      j = topList[i].index - 50; //This will sum bins within 0.05 m/z... might be too wide...
      k = topList[i].index + 51;
      if (j<0) j = 0;
      if (k>binCount) k = binCount;
      for (j = j; j<k; j++){
        intensity += bin[j];
        bin[j] = 0;
      }
      intensity /= s.size();
      avg.add(offset + topList[i].index*binWidth, intensity);
    }
    avg.sortMZ();
  }

  //clean up memory
  delete[] bin;
  delete[] pos;

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
    for (int q = 0; q<spec[b]->sizePrecursor(); q++){
      if (spec[b]->getTopPeps(q)->singletList.size() == 0) continue;
      list<kSingletScoreCard>::iterator it = spec[b]->getTopPeps(q)->singletList.begin();
      if (it->simpleScore>spec[b]->getScoreCard(0).simpleScore){
        if (it->simpleScore>maxScore) {
          maxScore = it->simpleScore;
          unknownMass = spec[b]->getPrecursor(q).monoMass - it->mass;
        }
      }
    }
    if (maxScore>0) oddCount++;
    if (maxScore>3) bigCount++;
    if (maxScore>3 && spec[b]->getScoreCard(0).simpleScore>0 && maxScore>spec[b]->getScoreCard(0).simpleScore * 2) twoCount++;
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

void KData::initHardklor(){
  hs2.noBase = true;
  hs4.noBase = true;

  CHardklorVariant hv;
  hv.addAtom(107, 2);
  hv.addEnrich(107, 2, params->enrichment);
  hs2.variant->push_back(hv);

  hv.clear();
  hv.addAtom(107, 4);
  hv.addEnrich(107, 2, params->enrichment);
  hs4.variant->push_back(hv);

  hv.clear();
  hs.winSize = hs2.winSize = hs4.winSize = 10;
  hs.peptide = hs2.peptide = hs4.peptide = 4;
  hs.sn = hs2.sn = hs4.sn = 0;
  hs.depth = hs2.depth = hs4.depth = 3;
  hs.minCharge = hs2.minCharge = hs4.minCharge = 2;
  hs.maxCharge = hs2.maxCharge = hs4.maxCharge = 8;
  hs.algorithm = Version2;
  if (params->instrument == 1) hs.msType = hs2.msType = hs4.msType = FTICR;
  else hs.msType = hs2.msType = hs4.msType = OrbiTrap;
  hs.res400 = hs2.res400 = hs4.res400 = params->ms1Resolution;
  hs.corr = 0.875;
  hs2.corr = hs4.corr = 0.85;
  hs.centroid = hs2.centroid = hs4.centroid = true;

  strcpy(hs.inFile, "PLTmp.ms1");
  strcpy(hs2.inFile, "PLTmp.ms1");
  strcpy(hs4.inFile, "PLTmp.ms1");

  hs.fileFormat = hs2.fileFormat = hs4.fileFormat = ms1;

  CHardklorVariant hkv;
  vector<CHardklorVariant> pepVariants;
  pepVariants.clear();
  pepVariants.push_back(hkv);

  averagine = new CAveragine*[params->threads]();
  mercury = new CMercury8*[params->threads]();
  h = new CHardklor2*[params->threads]();
  hO = new CHardklor*[params->threads]();
  mutexHardklor = new Mutex[params->threads]();
  bHardklor = new bool[params->threads]();
  for (int a = 0; a<params->threads; a++){
    averagine[a] = new CAveragine(NULL, NULL);
    mercury[a] = new CMercury8(NULL);
    if (a == 0) models = new CModelLibrary(averagine[a], mercury[a]);
    h[a] = new CHardklor2(averagine[a], mercury[a], models);
    h[a]->Echo(false);
    h[a]->SetResultsToMemory(true);
    hO[a] = new CHardklor(averagine[a], mercury[a]);
    hO[a]->Echo(false);
    hO[a]->SetResultsToMemory(true);
    Threading::CreateMutex(&mutexHardklor[a]);
    bHardklor[a] = false;
  }
  models->eraseLibrary();
  models->buildLibrary(2, 8, pepVariants);

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

void KData::memoryAllocate(){
  //find largest possible array for a spectrum
  int threads=params->threads;
  double xlMass=0;
  for(size_t a=0;a<params->xLink->size();a++){
    if(params->xLink->at(a).mass>xlMass) xlMass=params->xLink->at(a).mass;
  }
  double td = params->maxPepMass * 2 + xlMass+1;
  maxPrecursorMass=(int)td;
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
  //kSingletScoreCard* sc;
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
    list<kSingletScoreCard>::iterator it=tp->singletList.begin();
    k=1;
    while (it != tp->singletList.end()){
      fprintf(f,"    <peptide rank=\"%d\" index=\"%d\", sequence=\"",k++,it->pep1);
      db.getPeptideSeq(db.getPeptideList()->at(it->pep1).map->at(0).index, db.getPeptideList()->at(it->pep1).map->at(0).start, db.getPeptideList()->at(it->pep1).map->at(0).stop, strs);
      for (i = 0; i<strlen(strs); i++){
        fprintf(f, "%c", strs[i]);
        for (x = 0; x<it->modLen; x++){
          if (it->mods[x].pos == char(i)) fprintf(f, "[%.2lf]", it->mods[x].mass);
        }
        if (char(i) == it->k1) fprintf(f, "[x]");
      }
      fprintf(f, "\" link_site=\"%d\" score=\"%.4f\" cpScore=\"%.4f\" matches=\"%d\" longest_run=\"%d\" mass=\"%.4lf\"/>\n", (int)it->k1+1, it->simpleScore, it->cpScore, it->matches, it->conFrag, it->mass);
      it++;
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
      for (x = 0; x<psm.mods1.size(); x++){
        if (psm.mods1[x].pos == (char)i) {
          sprintf(st, "[%.2lf]", psm.mods1[x].mass);
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
        for (x = 0; x<psm.mods2.size(); x++){
          if (psm.mods2[x].pos == (char)i) {
            sprintf(st, "[%.2lf]", psm.mods2[x].mass);
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
    fprintf(f, "sequence=\"%s\" score=\"%.4f\" cpScore=\"%.4f\" evalue=\"%.3e\" mass=\"%.4lf\"",&pep1[0], psm.simpleScore, psm.cpScore1+psm.cpScore2, psm.eVal, psm.mass);
    if (psm.link<0) fprintf(f, " crosslinker_mass=\"0\"/>\n"); 
    else fprintf(f, " crosslinker_mass=\"%.4lf\"/>\n", link[psm.link].mass);
  }
  fprintf(f, "  </results_list>\n");

  kDiag dsc;
  fprintf(f, "  <top_results>\n");
  for(size_t pr=0;pr<s.sizePrecursor();pr++){
    p = s.getPrecursor2(pr);
    fprintf(f, "   <precursor mass=\"%.4lf\" charge=\"%d\">\n", p->monoMass, p->charge);
    list<kScoreCard>::iterator it=p->topSingle.begin();
    while(it!=p->topSingle.end()){
      convertToDiag(*it, dsc, db);
      fprintf(f, "    <top_single sequence=\"%s\" score=\"%.4f\" cpScore=\"%.4f\" evalue=\"%.3e\" mass=\"%.4lf\" crosslinker_mass=\"%.4lf\"/>\n", dsc.sequence.c_str(), dsc.simpleScore, dsc.cpScore, dsc.evalue, dsc.mass, dsc.xlMass);
      for(size_t z=0;z<it->alternate.size();z++){
        convertToDiag(it->alternate[z], dsc, db);
        fprintf(f, "    <top_single sequence=\"%s\" score=\"%.4f\" cpScore=\"%.4f\" evalue=\"%.3e\" mass=\"%.4lf\" crosslinker_mass=\"%.4lf\"/>\n", dsc.sequence.c_str(), dsc.simpleScore, dsc.cpScore, dsc.evalue, dsc.mass, dsc.xlMass);
      }
      it++;
    }

    it = p->topLoop.begin();
    while (it != p->topLoop.end()) {
      convertToDiag(*it, dsc, db);
      fprintf(f, "    <top_loop sequence=\"%s\" score=\"%.4f\" cpScore=\"%.4f\" evalue=\"%.3e\" mass=\"%.4lf\" crosslinker_mass=\"%.4lf\"/>\n", dsc.sequence.c_str(), dsc.simpleScore, dsc.cpScore, dsc.evalue, dsc.mass, dsc.xlMass);
      for (size_t z = 0; z < it->alternate.size(); z++) {
        convertToDiag(it->alternate[z], dsc, db);
        fprintf(f, "    <top_loop sequence=\"%s\" score=\"%.4f\" cpScore=\"%.4f\" evalue=\"%.3e\" mass=\"%.4lf\" crosslinker_mass=\"%.4lf\"/>\n", dsc.sequence.c_str(), dsc.simpleScore, dsc.cpScore, dsc.evalue, dsc.mass, dsc.xlMass);
      }
      it++;
    }

    it=p->topXL.begin();
    while (it != p->topXL.end()) {
      convertToDiag(*it, dsc, db);
      fprintf(f, "    <top_xl sequence=\"%s\" score=\"%.4f\" cpScore=\"%.4f\" evalue=\"%.3e\" mass=\"%.4lf\" crosslinker_mass=\"%.4lf\"/>\n", dsc.sequence.c_str(), dsc.simpleScore, dsc.cpScore, dsc.evalue, dsc.mass, dsc.xlMass);
      for (size_t z = 0; z < it->alternate.size(); z++) {
        convertToDiag(it->alternate[z], dsc, db);
        fprintf(f, "    <top_xl sequence=\"%s\" score=\"%.4f\" cpScore=\"%.4f\" evalue=\"%.3e\" mass=\"%.4lf\" crosslinker_mass=\"%.4lf\"/>\n", dsc.sequence.c_str(), dsc.simpleScore, dsc.cpScore, dsc.evalue, dsc.mass, dsc.xlMass);
      }
      it++;
    }
    fprintf(f, "   </precursor>\n");
  }
  fprintf(f, "  </top_results>\n");

  /*
  fprintf(f, "  <histogramSinglet count=\"%d\" intercept=\"%.4f\" slope=\"%.4f\" rsq=\"%.4lf\" start=\"%.4f\" next=\"%.4f\" max=\"%d\">\n", s.histogramSingletCount, s.tmpSingletIntercept, s.tmpSingletSlope, s.tmpSingletRSquare, s.tmpSingletIStartCorr, s.tmpSingletINextCorr, s.tmpSingletIMaxCorr);
  for (j = 0; j<HISTOSZ; j++){
    fprintf(f, "   <bin id=\"%d\" value=\"%d\"/>\n", j, s.histogramSinglet[j]);
  }
  fprintf(f, "  </histogramSinglet>\n");
  */

  if(params->diagHistogram){
    fprintf(f, "  <histogram count=\"%d\" intercept=\"%.4f\" slope=\"%.4f\" rsq=\"%.4lf\" start=\"%.4f\" next=\"%.4f\" max=\"%d\">\n", s.histogramCount,s.tmpIntercept,s.tmpSlope,s.tmpRSquare,s.tmpIStartCorr,s.tmpINextCorr,s.tmpIMaxCorr);
    for (j = 0; j<HISTOSZ; j++){
      fprintf(f, "   <bin id=\"%d\" value=\"%d\" score=\"%.1lf\" count=\"%d\"/>\n", j, s.histogram[j],(double)j/10,s.histogramO[j]);
    }
    fprintf(f, "  </histogram>\n");
  }
  

  fprintf(f, " </scan>\n");
}

void KData::convertToResults(kScoreCard& sc, kResults& res, CnpxSpectrumQuery& sq, KDatabase& db, size_t index){
  int j;//, k, n, d;
  //char fName[1056];
  //char outPath[1056];
  char peptide[256];
  //char tmp[16];
  char specID[256];

  kPeptide pep;
  kPeptide pep2;
  kPrecursor precursor;
  //kScoreCard tmpSC2;

  //kEnzymeRules enzyme;
  //kResults res;
  res.scanNumber = spec[index]->getScanNumber();
  res.scanID = spec[index]->getNativeID();
  res.rTime = spec[index]->getRTime();

  sq.spectrum = res.baseName + "." + to_string(res.scanNumber) + "." + to_string(res.scanNumber) + "." + to_string(spec[index]->getPrecursor(0).charge);
  sq.start_scan = res.scanNumber;
  sq.end_scan = res.scanNumber;
  sq.precursor_neutral_mass = spec[index]->getPrecursor(0).monoMass;
  sq.assumed_charge = spec[index]->getPrecursor(0).charge;
  sq.index = 1;
  sq.retention_time_sec = spec[index]->getRTime() * 60;

  //Get precursor ion for the PSM
  precursor = spec[index]->getPrecursor((int)sc.precursor);
  res.obsMass = precursor.monoMass;
  res.charge = precursor.charge;
  res.ppm = (sc.mass - precursor.monoMass) / precursor.monoMass * 1e6;
  res.psmMass = sc.mass;
  res.hk = precursor.corr;

  sq.assumed_charge = res.charge;
  sq.precursor_neutral_mass = res.obsMass;
  //sq.spectrum is added later
  //sq.spectrum = specID;

  res.score = sc.simpleScore;
  res.scoreDelta = sc.dScore;
  res.eVal = sc.eVal;
  res.eVal1 = sc.eVal1;
  res.eVal2 = sc.eVal2;
  res.matches1 = sc.matches1;
  res.matches2 = sc.matches2;
  res.conFrag1 = sc.conFrag1;
  res.conFrag2 = sc.conFrag2;
  if (sc.score1 < sc.score2) res.scorePepDif = sc.score1;
  else res.scorePepDif = sc.score2;

  KTopPeps* tp = spec[index]->getTopPeps((int)sc.precursor);
  list<kSingletScoreCard>::iterator grr = tp->singletList.begin();
  int rank = 1;
  res.rankA = 0;
  res.rankB = 0;
  while (grr != tp->singletList.end()) {
    if (res.rankA == 0 && sc.pep1 == grr->pep1 && sc.score1 == grr->simpleScore) res.rankA = rank;
    if (res.rankB == 0 && sc.score2 > grr->simpleScore) res.rankB = rank;
    rank++;
    grr++;
  }
  if (res.rankB == 0) res.rankB = params->topCount;
  res.rank = res.rankA + res.rankB;
  res.scoreA = sc.score1;
  res.scoreB = sc.score2;
  res.massA = sc.mass1;
  res.massB = sc.mass2;

  //Get the peptide sequence(s)
  pep = db.getPeptide(sc.pep1);
  db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, peptide);
  res.peptide1 = peptide;
  res.mods1.clear();
  res.cTerm1 = pep.cTerm;
  res.nTerm1 = pep.nTerm;
  res.linkSite1 = sc.site1;
  if (sc.site2 > -1) res.linkSite2 = sc.site2; //loop-link
  res.n15Pep1 = pep.n15;
  for (j = 0; j < sc.mods1.size(); j++) res.mods1.push_back(sc.mods1[j]);
  res.peptide2 = "";
  if (sc.pep2 >= 0) {
    pep2 = db.getPeptide(sc.pep2);
    db.getPeptideSeq(pep2.map->at(0).index, pep2.map->at(0).start, pep2.map->at(0).stop, peptide);
    res.peptide2 = peptide;
    res.mods2.clear();
    res.cTerm2 = pep2.cTerm;
    res.nTerm2 = pep2.nTerm;
    res.linkSite2 = sc.site2;
    res.n15Pep2 = pep2.n15;
    for (j = 0; j < sc.mods2.size(); j++) res.mods2.push_back(sc.mods2[j]);
  }

  //Process the peptide
  res.modPeptide1 = processPeptide(pep, sc.mods1, db);
  res.modPeptide2 = "";
  if (res.peptide2.size() > 0) {
    res.modPeptide2 = processPeptide(pep2, sc.mods2, db);
  }

  //Get the link positions - relative to the peptide
  res.link1 = sc.k1;
  res.link2 = sc.k2;
  if (res.link1 >= 0) res.link1++;
  if (res.link2 >= 0) res.link2++;

  //set link type
  res.type = 0;
  if (sc.k1 >= 0 && sc.k2 >= 0) res.type = 1;
  if (sc.pep1 >= 0 && sc.pep2 >= 0) res.type = 2;
  if (res.type == 2 && sc.k1 == -1 && sc.k2 == -1) res.type = 3;

  if (res.type > 0 && res.type != 3) {
    res.xlMass = link[sc.link].mass;
    res.xlLabel = link[sc.link].label;
  }

  //Get the peptide indexes
  res.pep1 = sc.pep1;
  res.pep2 = sc.pep2;
  res.linkable1 = sc.linkable1;
  res.linkable2 = sc.linkable2;
  res.linkerID = sc.link;

  //Process the protein
  processProtein(res.pep1, res.link1 - 1, res.linkSite1, res.protein1, res.protPos1, res.decoy1, db);
  if (res.modPeptide2.size() > 1) {
    processProtein(res.pep2, res.link2 - 1, res.linkSite2, res.protein2, res.protPos2, res.decoy2, db);
    if (res.decoy1 || res.decoy2) res.decoy = true;
    else res.decoy = false;
  } else if (res.linkSite2 > -1) { //loop link special case.
    processProtein(res.pep1, res.link2 - 1, res.linkSite2, res.protein2, res.protPos2, res.decoy2, db);
    if (!res.decoy1 || !res.decoy2) res.decoy = false;
    else res.decoy = true;
  } else {
    res.decoy = res.decoy1;
  }

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
    sprintf(str, "%s_%d", m_sd->id.c_str(), (int)m_sil->spectrumIdentificationResult.size());
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
      if (db[pep.map->at(j).index].name[i]==' ') break;
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
        if(db[pep.map->at(j).index].name[i]==' ') break;
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

  if(count>1) return false; //currently only export first CSM

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
  if(r.type==2) fprintf(f,"\t%.4lf\t%.4lf",r.scoreA,r.scoreB);
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
  vector<string> vPr;
  pep = db.getPeptide(r.pep1);
  for(j=0;j<pep.map->size();j++){
    protein="";
    if(!r.decoy1){ //do not report decoy protein names if peptide also belongs to a target protein.
      if (db[pep.map->at(j).index].name.find(params->decoy)!=string::npos) continue;
    }
    for(i=0;i<db[pep.map->at(j).index].name.size();i++){
      if(params->truncate>0 && i==params->truncate) break;
      if(db[pep.map->at(j).index].name[i]==' ') protein+='_';
      else protein+=db[pep.map->at(j).index].name[i];
    }
    fprintf(f,"\t%s",protein.c_str());
    vPr.push_back(protein);
  }
  if(r.pep2>=0){
    pep = db.getPeptide(r.pep2);
    for(j=0;j<pep.map->size();j++){
      protein="";
      if (!r.decoy2){ //do not report decoy protein names if peptide also belongs to a target protein.
        if (db[pep.map->at(j).index].name.find(params->decoy) != string::npos) continue;
      }
      for(i=0;i<db[pep.map->at(j).index].name.size();i++){
        if(params->truncate>0 && i==params->truncate) break;
        if(db[pep.map->at(j).index].name[i]==' ') protein+='_';
        else protein+=db[pep.map->at(j).index].name[i];
      }
      size_t a;
      for(a=0;a<vPr.size();a++){
        if(vPr[a].compare(protein)==0) break;
      }
      if(a==vPr.size()){
        fprintf(f,"\t%s",protein.c_str());
        vPr.push_back(protein);
      }
    }
  }

  fprintf(f,"\n");

  return true;
}

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
  vector<kResults> vRes;

  NeoPepXMLParser p;
  NeoPepXMLParser pSingle;
  NeoPepXMLParser pLoop;
  NeoPepXMLParser pXL;
  string pepXML_fileName;
  string pepXML_fileName_single;
  string pepXML_fileName_loop;
  string pepXML_fileName_xl;

  CMzIdentML mzID;
  string analysisSoftware_ref;
  string sip_ref;

  bool bBadFiles;
  bool bInter;
  //bool bTarget1;
  bool bDupe;
  bool bDiag;

  int scoreIndex;
  int iDupe;

  double topScore;

  string tmpPep1;
  string tmpPep2;
  string outFile;
  //string dStr;

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
    sprintf(fName,"%s.perc.intra.pin",params->outFile);
    fIntra=fopen(fName,"wt");
    if(fIntra==NULL) bBadFiles=true;
    sprintf(fName,"%s.perc.inter.pin",params->outFile);
    fInter=fopen(fName,"wt");
    if(fInter==NULL) bBadFiles=true;
    sprintf(fName,"%s.perc.loop.pin",params->outFile);
    fLoop=fopen(fName,"wt");
    if(fLoop==NULL) bBadFiles=true;
    sprintf(fName,"%s.perc.single.pin",params->outFile);
    fSingle=fopen(fName,"wt");
    if(fSingle==NULL) bBadFiles=true;
    if (params->dimers){
      sprintf(fName, "%s.perc.dimer.pin", params->outFile);
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

    if(params->splitPepXML){
      sprintf(fName, "%s.single.pep.xml", params->outFile);
      pepXML_fileName_single = fName;
      ft = fopen(fName, "wt");
      if (ft == NULL) bBadFiles = true;
      else fclose(ft);
      pa.summary_xml = pepXML_fileName_single;
      pSingle.msms_pipeline_analysis.push_back(pa);

      sprintf(fName, "%s.loop.pep.xml", params->outFile);
      pepXML_fileName_loop = fName;
      ft = fopen(fName, "wt");
      if (ft == NULL) bBadFiles = true;
      else fclose(ft);
      pa.summary_xml = pepXML_fileName_loop;
      pLoop.msms_pipeline_analysis.push_back(pa);

      sprintf(fName, "%s.xl.pep.xml", params->outFile);
      pepXML_fileName_xl = fName;
      ft = fopen(fName, "wt");
      if (ft == NULL) bBadFiles = true;
      else fclose(ft);
      pa.summary_xml = pepXML_fileName_xl;
      pXL.msms_pipeline_analysis.push_back(pa);
    }
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
      fprintf(fIntra,"SpecId\tLabel\tscannr\tScore\tScoreA\tScoreB\tdScore\tnegLog10eVal\tnegLog10eValA\tnegLog10eValB\tIonMatch\tConIonMatch\tIonMatchA\tConIonMatchA\tIonMatchB\tConIonMatchB\t");
      fprintf(fInter,"SpecId\tLabel\tscannr\tScore\tScoreA\tScoreB\tdScore\tnegLog10eVal\tnegLog10eValA\tnegLog10eValB\tIonMatch\tConIonMatch\tIonMatchA\tConIonMatchA\tIonMatchB\tConIonMatchB\t");
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
  //dStr=params->decoy;
  for(i=0;i<spec.size();i++) {

    //update top hits so that a target result is always first among ties between targets and decoys
    //spec[i]->refreshScore(db, dStr);
    vRes.clear();

    //Check if we need to output diagnostic information
    bDiag=false;
    if(params->diag->size()>0){
      if(params->diag->at(0)==-1) {
        bDiag=true;
      } else {
        for (d = 0; d<params->diag->size(); d++){
          if (spec[i]->getScanNumber() == params->diag->at(d)){
            bDiag=true;
            break;
          }
        }
      }
    }
    if(bDiag) outputDiagnostics(fDiag,*spec[i],db);

    scoreIndex=0;
    tmpSC = spec[i]->getScoreCard(scoreIndex);
    res.scanNumber = spec[i]->getScanNumber();
    res.scanID = spec[i]->getNativeID();
    res.rTime = spec[i]->getRTime();
    
    CnpxSpectrumQuery sq;
    sq.spectrum = res.baseName + "." + to_string(res.scanNumber) + "." + to_string(res.scanNumber) + "." + to_string(spec[i]->getPrecursor(0).charge);
    sq.start_scan = res.scanNumber;
    sq.end_scan = res.scanNumber;
    sq.precursor_neutral_mass = spec[i]->getPrecursor(0).monoMass;
    sq.assumed_charge = spec[i]->getPrecursor(0).charge;
    sq.index = 1;
    sq.retention_time_sec = spec[i]->getRTime() * 60;

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
    //TODO: store all ties in a short vector so results can be sorted to always display in the same order.
    //      proper order should be Targets before Decoys, Sequence, then LinkPos
    topScore=tmpSC.simpleScore;
    int count=0;
    while(tmpSC.simpleScore==topScore){

      count++;

      //Get precursor ion for the PSM
      precursor = spec[i]->getPrecursor((int)tmpSC.precursor);
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
        tmpSC2 = spec[i]->getScoreCard(n++);
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
      
      KTopPeps* tp = spec[i]->getTopPeps((int)tmpSC.precursor);
      list<kSingletScoreCard>::iterator grr=tp->singletList.begin();
      int rank=1;
      res.rankA=0;
      res.rankB=0;
      while(grr!=tp->singletList.end()){
        if(res.rankA==0 && tmpSC.pep1==grr->pep1 && tmpSC.score1==grr->simpleScore) res.rankA=rank;
        if(res.rankB==0 && tmpSC.score2>grr->simpleScore) res.rankB=rank;
        rank++;
        grr++;
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
      for(j=0;j<tmpSC.mods1.size();j++) res.mods1.push_back(tmpSC.mods1[j]);
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
        for(j=0;j<tmpSC.mods2.size();j++) res.mods2.push_back(tmpSC.mods2[j]);
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
      res.linkerID = tmpSC.link;

      //Edge case where single peptide is shared between linked and non-linked peptide lists
      //This occurs when the peptide appears multiple times in a database: internally and on
      //the c-terminus for amine reactive cross-linkers, for example.
      bDupe=false;
      if(res.type==0){
        n=scoreIndex+1;
        iDupe=1;
        while(n<19){
          iDupe++;
          tmpSC2 = spec[i]->getScoreCard(n++);
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

      vRes.push_back(res);

      tmpPep1=res.peptide1;
      sprintf(tmp,"(%d)",res.link1);
      tmpPep1+=tmp;
      tmpPep2 = res.peptide2;
      sprintf(tmp, "(%d)", res.link2);
      tmpPep2 += tmp;

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
        if(params->splitPepXML && scoreIndex==0){
          if(spec[i]->topSingle.simpleScore>0) {
            kResults r;
            CnpxSpectrumQuery sq2;
            convertToResults(spec[i]->topSingle,r,sq2,db,i);
            //cout << sq2.start_scan << " .. " << res.scanNumber << endl;
            outputNeoPepXML(sq2,db,r);
            sq2.index = (int)pSingle.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.size() + 1;
            sq2.spectrum=sq.spectrum;
            pSingle.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.push_back(sq2);
            //if(sq2.start_scan<1) { cout << "single wtf " << sq2.start_scan << endl;}
          }
          if (spec[i]->topLoop.simpleScore > 0) {
            kResults r;
            CnpxSpectrumQuery sq2;
            convertToResults(spec[i]->topLoop, r, sq2, db, i);
            outputNeoPepXML(sq2, db, r);
            sq2.index = (int)pLoop.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.size() + 1;
            sq2.spectrum = sq.spectrum;
            pLoop.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.push_back(sq2);
            //if (sq2.start_scan < 1) { cout << "loop wtf " << sq2.start_scan << endl; }
          }
          if (spec[i]->topXL.simpleScore > 0) {
            kResults r;
            CnpxSpectrumQuery sq2;
            convertToResults(spec[i]->topXL, r, sq2, db, i);
            outputNeoPepXML(sq2, db, r);
            sq2.index = (int)pXL.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.size() + 1;
            sq2.spectrum = sq.spectrum;
            pXL.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.push_back(sq2);
            //if (sq2.start_scan < 1) { cout << "xl wtf " << sq2.start_scan << endl; }
          }
        }
      }

      //Get the next entry - it must also be exported if it has the same score
      if(bDupe) scoreIndex+=iDupe;
      else scoreIndex++;
      if(scoreIndex>=20) break;
      tmpSC = spec[i]->getScoreCard(scoreIndex);
    }

    outputTxt(fOut,db,vRes);

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
    if(params->splitPepXML){
      pSingle.write(pepXML_fileName_single.c_str(),true);
      pLoop.write(pepXML_fileName_loop.c_str(),true);
      pXL.write(pepXML_fileName_xl.c_str(),true);
    }
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

//TODO: decide if output should be styled this way. If so, fix tmpSC variables within.
bool KData::outputTxt(FILE* f, KDatabase& db, std::vector<kResults>& r){

  ////sort if needed
  //if(r.size()>1){
  //  for (size_t a = 0; a<r.size() - 1; a++){
  //    for (size_t b = a + 1; b<r.size(); b++){

  //      //sort according to decoy state, linker type (single, loop, cross), peptide sequence, link pos, second peptide, second link pos
  //      if (r[b].decoy<r[a].decoy) goto swap_hit;  //check decoy status first
  //      else if (r[a].decoy<r[b].decoy) continue;

  //      //check link type next
  //      if (r[b].type<r[a].type) goto swap_hit;
  //      else if (r[a].type<r[b].type) continue;

  //      //check peptide sequence (first only)
  //      int c = r[b].modPeptide1.compare(r[a].modPeptide1);
  //      if (c<0) goto swap_hit;
  //      else if (c>0) continue;

  //      //check link position (first peptide)
  //      if (r[b].link1<r[a].link1) goto swap_hit;
  //      else if (r[a].type==1 && r[b].type == 1 && r[a].link1 == r[b].link1){
  //        if (r[b].link2<r[a].link2) goto swap_hit; //check 2nd pos on loop links
  //        continue;
  //      }

  //      //check second peptide (cross-links only)
  //      if (r[a].type == 2 && r[b].type == 2){
  //        int c = r[b].modPeptide2.compare(r[a].modPeptide2);
  //        if (c<0) goto swap_hit;
  //        else if (c == 0){
  //          if (r[b].link2<r[a].link2) goto swap_hit;
  //        }
  //      }

  //      continue; //already in order

  //      //the swap code
  //      swap_hit:
  //      kResults tmp = r[a];
  //      r[a] = r[b];
  //      r[b] = tmp;
  //    }
  //  }
  //}

  for(size_t a=0;a<r.size();a++){
    //Export Results:
    fprintf(f, "%d", r[a].scanNumber);
    //fprintf(fOut, "\t%.4lf",res.hk); //this was for diagnostics of hardklor correlation results (or lack of)
    fprintf(f, "\t%.4f", r[a].rTime);
    fprintf(f, "\t%.4lf", r[a].obsMass);
    fprintf(f, "\t%d", r[a].charge);
    fprintf(f, "\t%.4lf", r[a].psmMass);
    fprintf(f, "\t%.4lf", r[a].ppm);
    fprintf(f, "\t%.4lf", r[a].score);
    fprintf(f, "\t%.4lf", r[a].scoreDelta);
    fprintf(f, "\t%.3e", r[a].eVal);
    //fprintf(fOut,"\t%.4lf",res.scorePepDif);
    if (r[a].scoreA == 0) fprintf(f, "\t%.4lf", r[a].score);
    else fprintf(f, "\t%.4lf", r[a].scoreA);
    fprintf(f, "\t%.3e", r[a].eVal1);
    fprintf(f, "\t%s", r[a].modPeptide1.c_str());
    if (r[a].n15Pep1) fprintf(f, "-15N");
    fprintf(f, "\t%d", r[a].link1);

    //export protein
    fprintf(f, "\t%s", r[a].protein1.c_str());

    if (r[a].link1>-1) fprintf(f, "\t%s", r[a].protPos1.c_str());
    else fprintf(f, "\t-");

    if (r[a].modPeptide2.size()>1) {
      fprintf(f, "\t%.4lf", r[a].scoreB);
      fprintf(f, "\t%.3e", r[a].eVal2);
      fprintf(f, "\t%s", r[a].modPeptide2.c_str());
      if (r[a].n15Pep2) fprintf(f, "-15N");
      fprintf(f, "\t%d", r[a].link2);
      fprintf(f, "\t%s", r[a].protein2.c_str());
      fprintf(f, "\t%s", r[a].protPos2.c_str());
      if (r[a].linkerID>-1)fprintf(f, "\t%.4lf", link[r[a].linkerID].mass);
      else fprintf(f, "\t0");
    } else if (r[a].link2>-1){
      fprintf(f, "\t0\t999\t-\t%d\t-\t%s", r[a].link2, r[a].protPos2.c_str());
      fprintf(f, "\t%.4lf", link[r[a].linkerID].mass);
    } else {
      fprintf(f, "\t0\t999\t-\t-1\t-\t-\t0");
    }

    fprintf(f, "\n");
  }

  return true;
}

void KData::processMS2(kMS2struct* s){

  int j;

  formatMS2(s->s, s->pls);

  if(s->pls->size()<=params->minPeaks){
    s->state=4;
    s->thread=false;
    return;
  }

  bool bAddHardklor = false;
  bool bAddEstimate = false;

  if (params->preferPrecursor == 1){
    if (s->pls->sizePrecursor() == 0){
      if (s->pls->getCharge()>0) bAddEstimate = true;
      else bAddHardklor = true;
    }
  } else if (params->preferPrecursor == 0){
    s->pls->clearPrecursors();
    bAddHardklor = true;
    if (s->pls->getCharge()) bAddEstimate = true;
  } else {
    bAddHardklor = true;
    if (s->pls->getCharge()) bAddEstimate = true;
  }

  int ret;
  if (bAddHardklor && params->precursorRefinement){
    //only do Hardklor analysis if data contain precursor scans
    //ret = pre.getSpecRange(*spec[i]);
    int tIndex;
    Threading::LockMutex(mutexLockMS1);
    for (tIndex = 0; tIndex<params->threads; tIndex++){
      if (!bHardklor[tIndex]){
        bHardklor[tIndex] = true;
        break;
      }
    }
    if (tIndex == params->threads) cout << "Thread overload" << endl;
    Threading::UnlockMutex(mutexLockMS1);

    Threading::LockMutex(mutexHardklor[tIndex]);
    ret = processPrecursor(s, tIndex);
    bHardklor[tIndex] = false;
    Threading::UnlockMutex(mutexHardklor[tIndex]);
    //Threading::LockMutex(mutexLockMS1); //is this necessary?
    //bHardklor[tIndex] = false;
    //Threading::UnlockMutex(mutexLockMS1);
  }

  if (bAddEstimate){
    kPrecursor pr;
    pr.monoMass = s->pls->getMZ()*s->pls->getCharge() - 1.007276466*s->pls->getCharge();
    pr.charge = s->pls->getCharge();
    pr.corr = -5;
    s->pls->setCharge(pr.charge);
    s->pls->addPrecursor(pr, params->topCount);
    for (int px = 1; px <= params->isotopeError; px++){
      if (px == 4) break;
      pr.monoMass -= 1.00335483;
      pr.corr -= 0.1;
      s->pls->addPrecursor(pr, params->topCount);
    }
  }

  //Now clean up any duplicate precursors. They should already be in order of priority. Use 5ppm as tolerance
  for (int k = 0; k<s->pls->sizePrecursor(); k++){
    for (int n = k + 1; n<s->pls->sizePrecursor(); n++){
      double m1 = s->pls->getPrecursor(k).monoMass;
      double m2 = s->pls->getPrecursor(n).monoMass;
      double m = (m1 - m2) / m1*1e6;
      if (fabs(m)<5){
        s->pls->erasePrecursor(n);
        n--;
      }
    }
  }

  if (s->pls->sizePrecursor()>0){
    //build singletList
    //s->pls->resetSingletList();
    s->pls->peakCounts = s->pls->size();

  } else {
    s->state = 4; //no precursors, so advance state past transform to delete.
    s->thread = false;
    return;
  }

  Threading::LockMutex(mutexMemoryPool);
  for (j = 0; j<params->threads; j++){
    if (!memoryPool[j]){
      memoryPool[j] = true;
      break;
    }
  }
  Threading::UnlockMutex(mutexMemoryPool);

  if (j == params->threads){
    cout << "Error in KData::processMS2::state==2" << endl;
    exit(-1);
  }
  s->pls->kojakXCorr(tempRawData[j], tmpFastXcorrData[j], fastXcorrData[j], preProcess[j]);
  memoryPool[j] = false;

  s->state=3;
  s->thread = false;
}

int KData::processPrecursor(kMS2struct* s, int tIndex){

  int j;
  float rt = s->pls->getRTime();
  double mz = s->pls->getMZ();
  float maxIntensity = 0;
  float maxRT = 0;
  int best = 0;
  int precursor = -1;
  int ret=0;

  //Find the MS1 scan that contains the precursor ion within 6 seconds of when it was acquired.
  for (int i = 0; i<dMS1.size(); i++){
    if (dMS1[i]->getRTime()<rt - 0.167) continue;
    if (dMS1[i]->getRTime()>rt + 0.167) break;
    int j = findPeak(dMS1[i], mz, 10);

    if (j>-1){
      if (dMS1[i]->at(j).intensity>maxIntensity){
        maxIntensity = dMS1[i]->at(j).intensity;
        maxRT = dMS1[i]->getRTime();
        best = dMS1[i]->getScanNumber();
        precursor = i;
      }
    }
  }

  if (precursor<0){
    //cout << "Warning: Precursor not found for " << s->pls->getScanNumber() << " " << mz << endl;
    return ret;
  }

  //Get up to +/-15 sec of spectra around the max precursor intensity
  //This is done by extending on both sides until a gap is found or time is reached.
  //Additionally, stop if 2 scans found flanking either side (maximum 5 scans per precursor).
  vector<Spectrum*> vs;
  float rtHigh;
  float rtLow;
  rtHigh = rtLow = dMS1[precursor]->getRTime();
  int k = 0;
  for (int i = precursor; i<dMS1.size(); i++){
    if (dMS1[i]->getRTime()>maxRT + 0.25) break;
    j = findPeak(dMS1[i], mz, 10);
    if (j<0) break;
    vs.push_back(dMS1[i]);
    rtHigh = dMS1[i]->getRTime();
    k++;
    if (k == 3) break;
  }

  k = 0;
  int i = precursor;
  while (i>0){
    i--;
    if (dMS1[i]->getRTime()<maxRT - 0.25) break;
    j = findPeak(dMS1[i], mz, 10);
    if (j<0) break;
    vs.push_back(dMS1[i]);
    rtLow = dMS1[i]->getRTime();
    k++;
    if (k == 2) break;
  }

  //Average points between mz-1.5 and mz+2
  Spectrum sp;
  if (params->enrichment>0) averageScansCentroid(vs, sp, mz - 1.75, mz + 1.75);
  else averageScansCentroid(vs, sp, mz - 1.0, mz + 1.5);
  if (sp.size() == 0) {
    cout << "\n   WARNING: Unexpected precursor scan data!";
    if (!params->ms1Centroid) cout << " Params are set to MS1 profile mode, but are MS1 scans centroided?" << endl;
    return ret;
  }
  sp.setScanNumber(dMS1[precursor]->getScanNumber());

  //Obtain the possible precursor charge states of the selected ion.
  //Find the index of the closest peak to the selected m/z.
  vector<int> preCharges;
  double tmz = fabs(mz - sp[0].mz);
  for (j = 1; j<sp.size(); j++){
    if (fabs(mz - sp[j].mz)<tmz) tmz = fabs(mz - sp[j].mz);
    else break;
  }
  j = j - 1;
  h[tIndex]->QuickCharge(sp, j, preCharges);

  //Clear corr
  double corr = 0;
  double monoMass = 0;
  int charge = 0;
  kPrecursor pre;

  //Perform 18O2 analysis with Hardklor. If enrichment is set to 0, store unenriched results in the *O2 variables.
  //This is done in a non-competitive to identify an 18O2 peptide without solving everything
  if (params->enrichment>0) {
    hO[tIndex]->GoHardklor(hs2, &sp);

    for (j = 0; j<hO[tIndex]->Size(); j++){
      if (hO[tIndex]->operator[](j).corr>corr){
        monoMass = hO[tIndex]->operator[](j).monoMass;
        charge = hO[tIndex]->operator[](j).charge;
        corr = hO[tIndex]->operator[](j).corr;
      }
    }
    if (corr>0){
      pre.monoMass = monoMass;
      pre.charge = charge;
      pre.corr = corr;
      if (params->enrichment>0) pre.label = 1;
      else pre.label = 0;
      s->pls->addPrecursor(pre, params->topCount);
    }

  } else {
    h[tIndex]->GoHardklor(hs, &sp);

    //If nothing was found, really narrow down the window and try again.
    if (h[tIndex]->Size() == 0){
      averageScansCentroid(vs, sp, mz - 0.6, mz + 1.2);
      sp.setScanNumber(dMS1[precursor]->getScanNumber());
      if (params->enrichment>0) h[tIndex]->GoHardklor(hs2, &sp);
      else h[tIndex]->GoHardklor(hs, &sp);
    }

    float intensity = 0;
    for (j = 0; j<h[tIndex]->Size(); j++){

      //Must have highest intensity and intersect isolated peak.
      if (h[tIndex]->operator[](j).intensity<intensity) continue;
      tmz = (h[tIndex]->operator[](j).monoMass + 1.007276466*h[tIndex]->operator[](j).charge) / h[tIndex]->operator[](j).charge;
      while (tmz<(s->pls->getMZ() + 0.01)){
        if (fabs(tmz - s->pls->getMZ())<0.01){
          monoMass = h[tIndex]->operator[](j).monoMass;
          charge = h[tIndex]->operator[](j).charge;
          corr = h[tIndex]->operator[](j).corr;
          intensity = h[tIndex]->operator[](j).intensity;
          ret = 1;
          break;
        }
        tmz += (1.00335483 / h[tIndex]->operator[](j).charge);
      }
    }

    //failing to match precursor peak, keep most intense precursor in presumed isolation window
    if (corr == 0){
      for (j = 0; j<h[tIndex]->Size(); j++){
        if (h[tIndex]->operator[](j).intensity>intensity){
          monoMass = h[tIndex]->operator[](j).monoMass;
          charge = h[tIndex]->operator[](j).charge;
          corr = h[tIndex]->operator[](j).corr;
          intensity = h[tIndex]->operator[](j).intensity;
          ret = 2;
        }
      }
    }

    if (corr>0){
      pre.monoMass = monoMass;
      pre.charge = charge;
      pre.corr = corr;
      if (params->enrichment>0) pre.label = 1;
      else pre.label = 0;
      s->pls->addPrecursor(pre, params->topCount);
      //also add isotope error
      if (params->isotopeError>0){
        pre.monoMass -= 1.00335483;
        pre.corr = -1;
        s->pls->addPrecursor(pre, params->topCount);
      }
      if (params->isotopeError>1){
        pre.monoMass -= 1.00335483;
        pre.corr = -2;
        s->pls->addPrecursor(pre, params->topCount);
      }
      if (params->isotopeError>2){
        pre.monoMass -= 1.00335483;
        pre.corr = -3;
        s->pls->addPrecursor(pre, params->topCount);
      }
    }

  }

  //Clear corr
  corr = 0;

  //Perform the 18O4 Hardklor analysis on the same data if needed
  if (params->enrichment>0) {
    hO[tIndex]->GoHardklor(hs4, &sp);
    for (j = 0; j<hO[tIndex]->Size(); j++){
      if (hO[tIndex]->operator[](j).corr>corr){
        monoMass = hO[tIndex]->operator[](j).monoMass;
        charge = hO[tIndex]->operator[](j).charge;
        corr = hO[tIndex]->operator[](j).corr;
        ret = 1;
      }
    }
  }
  if (corr>0){
    pre.monoMass = monoMass;
    pre.charge = charge;
    pre.corr = corr;
    pre.label = 2;
    s->pls->addPrecursor(pre, params->topCount);
  }

  //Assume two precursors with nearly identical mass (within precursor tolerance) are the same.
  //This can occur when checking multiple enrichment states.
  //Keep only the higher correlated precursor.
  if (s->pls->sizePrecursor()>1){
    bool bCheck = true;
    while (bCheck){
      bCheck = false;
      for (k = 0; k<s->pls->sizePrecursor() - 1; k++){
        for (j = k + 1; j<s->pls->sizePrecursor(); j++){
          if (fabs(s->pls->getPrecursor(k).monoMass - s->pls->getPrecursor(j).monoMass) / s->pls->getPrecursor(k).monoMass*1e6 < params->ppmPrecursor){
            if (s->pls->getPrecursor(k).corr>s->pls->getPrecursor(j).corr) s->pls->erasePrecursor(j);
            else s->pls->erasePrecursor(k);
            bCheck = true;
            break;
          }
        }
        if (bCheck) break;
      }
    }
  }

  return ret;
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
  Spectrum*   s;
  Spectrum   c;
  kPrecursor pre;

  int totalScans = 0;
  int totalPeaks = 0;
  int collapsedPeaks = 0;
  int iPercent = 0;
  int iTmp;

  //size_t index;

  deque<kMS2struct*> dMS2;
  vMS1Buffer.reserve(2000);
  memoryAllocate();
  initHardklor();

  Threading::CreateMutex(&mutexLockMS1);

  ThreadPool<kMS2struct*>* threadPool = new ThreadPool<kMS2struct*>(processMS2, params->threads, params->threads, 1);

  for(size_t a=0;a<spec.size();a++) delete spec[a];
  spec.clear();

  msr.setFilter(MS1);
  msr.addFilter(MS2);

  //Set progress meter
  printf("%2d%%", iPercent);
  fflush(stdout);

  s = new Spectrum;

  if (!msr.readFile(params->msFile, *s)) return false;

  //temporary
  int nextMS2=0;

  while (s->getScanNumber()>0){

    totalScans++;
    if (s->size()<1) {
      msr.readFile(NULL, *s);
      continue;
    }

    if (s->getMsLevel() == 1) {
      vMS1Buffer.emplace_back(s);
      if (vMS1Buffer.size() == 10){  //When buffer is full, transfer to MS1 memory pool
        for (int a = 0; a<params->threads; a++) Threading::LockMutex(mutexHardklor[a]);
        while (spec.size()>0 && dMS1.size()>0 && dMS1.front()->getRTime()<spec.back()->getRTime() - 1){ //clear old memory
          delete dMS1.front();
          dMS1.pop_front();
        }
        for (size_t a = 0; a<vMS1Buffer.size(); a++){
          dMS1.emplace_back(vMS1Buffer[a]);
          vMS1Buffer[a] = NULL;
        }
        vMS1Buffer.clear();
        for (int a = 0; a<params->threads; a++) Threading::UnlockMutex(mutexHardklor[a]);
      }

    } else {
      kMS2struct* ms = new kMS2struct(s, params->topCount, params->binSize, params->binOffset);
      dMS2.emplace_back(ms);
    }

    //Update progress meter
    iTmp = msr.getPercent();
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }

    while (dMS2.size()>0 && dMS2[0]->state >= 3){ //copy and/or clear finished MS2 spectra
      if (dMS2[0]->state == 3) spec.push_back(dMS2[0]->pls);
      else delete dMS2[0]->pls;
      dMS2.pop_front();
      nextMS2--;
    }
    //Launch next MS2
    while(nextMS2<dMS2.size()) {
      if (dMS1.size()>0 && (dMS2[nextMS2]->s->getRTime() + 2)<dMS1.back()->getRTime()){
        //only launch this if there are enough MS1 for precursor analysis
        dMS2[nextMS2]->thread = true;
        threadPool->Launch(dMS2[nextMS2]);
        nextMS2++;
      } else break; //we got here because we need more MS1 first.
    }

    s = new Spectrum;
    msr.readFile(NULL, *s);

  }

  //finish flushing buffer
  for (int a = 0; a<params->threads; a++) Threading::LockMutex(mutexHardklor[a]);
  while (spec.size()>0 && dMS1.size()>0 && dMS1.front()->getRTime()<spec.back()->getRTime() - 1){ //clear old memory
    delete dMS1.front();
    dMS1.pop_front();
  }
  for (size_t a = 0; a<vMS1Buffer.size(); a++){
    dMS1.emplace_back(vMS1Buffer[a]);
    vMS1Buffer[a] = NULL;
  }
  vMS1Buffer.clear();
  for (int a = 0; a<params->threads; a++)  Threading::UnlockMutex(mutexHardklor[a]);

  while (nextMS2<dMS2.size()) {
    dMS2[nextMS2]->thread = true;
    threadPool->Launch(dMS2[nextMS2]);
    nextMS2++;
  }

  threadPool->WaitForQueuedParams();
  threadPool->WaitForThreads();

  //finish processing last MS2 scans
  while (dMS2.size()>0){
    while (dMS2.size()>0 && dMS2[0]->state >= 3){ //copy and/or clear finished MS2 spectra
      if (dMS2[0]->state == 3) spec.push_back(dMS2[0]->pls);
      else delete dMS2[0]->pls;
      dMS2.pop_front();
      nextMS2--;
    }
  }

  //Finalize progress meter
  if (iPercent<100) printf("\b\b\b100%%");
  cout << endl;

  //clean up remaining memory
  while (dMS1.size()>0){
    delete dMS1.front();
    dMS1.pop_front();
  }

  delete threadPool;
  threadPool = NULL;
  memoryFree();
  releaseHardklor();
  Threading::DestroyMutex(mutexLockMS1);

  //Build mass list - this orders all precursor masses, with an index pointing to the actual
  //array position for the spectrum. This is because all spectra will have more than 1
  //precursor mass
  kMass m;
  massList.clear();
  for (size_t i = 0; i<spec.size(); i++){
    m.index = (int)i;
    for (int j = 0; j<spec[i]->sizePrecursor(); j++){
      m.mass = spec[i]->getPrecursor(j).monoMass;
      massList.push_back(m);
    }
  }

  //sort mass list from low to high
  qsort(&massList[0], massList.size(), sizeof(kMass), compareMassList);

  if (bScans != NULL) delete[] bScans;
  bScans = new bool[spec.size()];

	return true;
}

void KData::releaseHardklor(){
  for (int a = 0; a<params->threads; a++){
    delete h[a];
    delete hO[a];
    Threading::DestroyMutex(mutexHardklor[a]);
    delete averagine[a];
    delete mercury[a];
  }
  delete[] h;
  delete[] hO;
  delete[] mutexHardklor;
  delete[] bHardklor;
  delete[] averagine;
  delete[] mercury;
  delete models;
}

void KData::report(){
  char str[256];
  sprintf(str, "%d spectra with sufficient data points (%d peaks) will be analyzed.",(int)spec.size(),params->minPeaks);
  klog->addMessage(str, true);
  cout << "  " << str << endl;
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

/*============================
  Private Utilities
============================*/
//First derivative method, returns base peak intensity of the set
//TODO: rewrite this to not require reliance on resolution and expected peak shapes.
//TODO: add verbage in Kojak to warn the user to not rely on centroiding here, but from the authors of the raw data
void KData::centroid(Spectrum* s, KSpectrum* out, double resolution, int instrument){
  int i, j;
  float maxIntensity;
  int bestPeak;
  bool bLastPos;

  int nextBest;
  double FWHM;
  double maxMZ = (*s)[s->size() - 1].mz + 1.0;
  kSpecPoint centroid;

  vector<double> x;
  vector<double> y;
  vector<double> c;
  int left, right;
  bool bPoly;
  float lastIntensity;

  out->clear();

  bLastPos = false;
  for (i = 0; i<s->size() - 1; i++){

    if ((*s)[i].intensity<(*s)[i + 1].intensity) {
      bLastPos = true;
      continue;
    } else {
      if (bLastPos){
        bLastPos = false;

        //find max and add peak
        maxIntensity = 0;
        for (j = i; j<i + 1; j++){
          if ((*s)[j].intensity>maxIntensity){
            maxIntensity = (*s)[j].intensity;
            bestPeak = j;
          }
        }

        //walk left and right to find bounds above half max
        left = right = bestPeak;
        lastIntensity = maxIntensity;
        for (left = bestPeak - 1; left>0; left--){
          if ((*s)[left].intensity<(maxIntensity / 3) || (*s)[left].intensity>lastIntensity){
            left++;
            break;
          }
          lastIntensity = (*s)[left].intensity;
        }
        lastIntensity = maxIntensity;
        for (right = bestPeak + 1; right<s->size() - 1; right++){
          if ((*s)[right].intensity<(maxIntensity / 3) || (*s)[right].intensity>lastIntensity){
            right--;
            break;
          }
          lastIntensity = (*s)[right].intensity;
        }

        //if we have at least 5 data points, try polynomial fit
        double r2;
        bPoly = false;
        if ((right - left + 1)>4){
          x.clear();
          y.clear();
          for (j = left; j <= right; j++){
            x.push_back((*s)[j].mz);
            y.push_back(log((*s)[j].intensity));
          }
          r2 = polynomialBestFit(x, y, c);
          if (r2>0.95){
            bPoly = true;
            centroid.mass = -c[1] / (2 * c[2]) + c[3];
            centroid.intensity = (float)exp(c[0] - c[2] * (c[1] / (2 * c[2]))*(c[1] / (2 * c[2])));
          } else {

          }
        }

        if (!bPoly){
          //Best estimate of Gaussian centroid
          //Get 2nd highest point of peak
          if (bestPeak == s->size()) nextBest = bestPeak - 1;
          else if ((*s)[bestPeak - 1].intensity > (*s)[bestPeak + 1].intensity) nextBest = bestPeak - 1;
          else nextBest = bestPeak + 1;

          //Get FWHM
          switch (instrument){
          case 0: FWHM = (*s)[bestPeak].mz*sqrt((*s)[bestPeak].mz) / (20 * resolution); break;  //Orbitrap
          case 1: FWHM = (*s)[bestPeak].mz*(*s)[bestPeak].mz / (400 * resolution); break;				//FTICR
          default: break;
          }

          //Calc centroid MZ (in three lines for easy reading)
          centroid.mass = pow(FWHM, 2)*log((*s)[bestPeak].intensity / (*s)[nextBest].intensity);
          centroid.mass /= GAUSSCONST*((*s)[bestPeak].mz - (*s)[nextBest].mz);
          centroid.mass += ((*s)[bestPeak].mz + (*s)[nextBest].mz) / 2;

          //Calc centroid intensity
          centroid.intensity = (float)((*s)[bestPeak].intensity / exp(-pow(((*s)[bestPeak].mz - centroid.mass) / FWHM, 2)*GAUSSCONST));
        }

        //some peaks are funny shaped and have bad gaussian fit.
        //if error is more than 10%, keep existing intensity
        if (fabs(((*s)[bestPeak].intensity - centroid.intensity) / centroid.intensity * 100) > 10 ||
          //not a good check for infinity
          centroid.intensity>9999999999999.9 ||
          centroid.intensity < 0) {
          centroid.intensity = (*s)[bestPeak].intensity;
        }

        //Hack until I put in mass ranges
        if (centroid.mass<0 || centroid.mass>maxMZ) {
          //do nothing if invalid mz
        } else {
          out->addPoint(centroid);
        }

      }

    }
  }

}

//Function tries to remove isotopes of signals by stacking the intensities on the monoisotopic peak
//Also creates an equal n+1 peak in case wrong monoisotopic peak was identified.
void KData::collapseSpectrum(KSpectrum& s){
  int i, j, k, n;
  int charge, z;
  int maxIndex;
  float max;
  float cutoff;
  vector<int> dist;

  vector<kSpecPoint> s2;

  while (true){
    max = 0.1f;
    for (i = 0; i<s.size(); i++){
      if (s[i].intensity>max){
        max = s[i].intensity;
        maxIndex = i;
      }
    }

    //finish and exit function
    if (max<1) break;

    dist.clear();
    dist.push_back(maxIndex);

    //check right
    j = maxIndex + 1;
    while (j<s.size() && (s[j].mass - s[maxIndex].mass)<1.1){
      if (s[j].intensity<1) {
        j++;
        continue;
      }
      charge = getCharge(s, maxIndex, j);

      if (charge == 0){
        j++;
        continue;
      }

      //try stepping along at same charge state here out
      //note that if this doesn't work, it doesn't go back and look for a different charge state
      dist.push_back(j);
      k = j;
      n = j + 1;
      while (n<s.size() && (s[n].mass - s[k].mass)<1.1){
        if (s[n].intensity<1) {
          n++;
          continue;
        }
        z = getCharge(s, k, n);
        if (z>0 && z<charge) {
          break;
        } else if (z == charge && (s[n].mass - s[k].mass)>(0.99 / charge) && (s[n].mass - s[k].mass)<(1.0041 / charge)) {
          dist.push_back(n);
          k = n;
          n++;
        } else {
          n++;
        }
      }
      break;
    }

    //if nothing found to the right, quit here?
    if (dist.size() == 1){
      s2.push_back(s[dist[0]]);
      s[dist[0]].intensity = 0;
      continue;
    }

    //step to the left
    j = maxIndex - 1;
    while (j >= 0 && (s[maxIndex].mass - s[j].mass)<1.1){
      if (s[j].intensity<1) {
        j--;
        continue;
      }
      z = getCharge(s, j, maxIndex);
      if (z != charge){
        j--;
        continue;
      }

      //try stepping along at same charge state here out
      dist.push_back(j);
      k = j;
      n = j - 1;
      while (n >= 0 && (s[k].mass - s[n].mass)<1.1){
        if (s[n].intensity<1) {
          n--;
          continue;
        }
        z = getCharge(s, n, k);
        //printf("\tleft\t%.6lf\t%.6lf\t%d\n",s[n].mz,s[k].mz-s[n].mz,z);
        if (z>0 && z<charge) {
          break;
        } else if (z == charge && s[k].mass - s[n].mass > 0.99 / charge && s[k].mass - s[n].mass < 1.0041 / charge) {
          dist.push_back(n);
          k = n;
          n--;
        } else {
          n--;
        }
      }
      break;
    }


    //Only accept size of 2 if charge is 1 or 2
    if (dist.size() == 2){
      if (charge<3){
        max = s[dist[0]].intensity + s[dist[1]].intensity;
        kSpecPoint sp;
        sp.mass = s[dist[0]].mass;
        sp.intensity=max;
        s2.push_back(sp);
        s[dist[1]].intensity = 0;
        // s2.add(s[dist[1]].mz,max);
      } else {
        s2.push_back(s[dist[0]]);
        // s2.add(s[dist[1]]);
      }
      s[dist[0]].intensity = 0;
      //s[dist[1]].intensity=0;
    } else {
      cutoff = max / 20;
      max = 0;
      j = dist[0];
      k = dist[1];
      for (i = 0; i<(int)dist.size(); i++) {
        if (dist[i]<j && s[dist[i]].intensity>cutoff){
          k = j;
          j = dist[i];
        }
        if (s[dist[i]].intensity>cutoff){
          max += s[dist[i]].intensity;
          s[dist[i]].intensity = 0;
        }
      }
      kSpecPoint sp;
      sp.mass=s[j].mass;
      sp.intensity=max;
      s2.push_back(sp);
      //s2.add(s[k].mz,max);
    }

  }

  sort(s2.begin(),s2.end(),compareSpecPoint);
  s.clear();
  for (i = 0; i<s2.size(); i++) {
    if (i<s2.size() - 1 && s2[i].mass == s2[i + 1].mass){
      if (s2[i].intensity>s2[i + 1].intensity) s.addPoint(s2[i]);
      else s.addPoint(s2[i + 1]);
      i++;
    } else {
      s.addPoint(s2[i]);
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

int KData::getCharge(KSpectrum& s, int index, int next){
  double mass = s[next].mass - s[index].mass;
  if (mass>0.99 && mass<1.007) return 1;
  else if (mass>0.495 && mass<0.5035) return 2;
  else if (mass>0.33 && mass<0.335667) return 3;
  else if (mass>0.2475 && mass<0.25175) return 4;
  else if (mass>0.198 && mass<0.2014) return 5;
  else if (mass>0.165 && mass<0.1678333) return 6;
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

string KData::processPeptide(kPeptide& pep, vector<kPepMod>& mod, KDatabase& db){
  char tmp[32];
  size_t j,k;
  string seq = "";
  string peptide;

  db.getPeptideSeq(pep.map->at(0).index, pep.map->at(0).start, pep.map->at(0).stop, peptide);

  if (pep.nTerm && aa.getFixedModMass('$') != 0) {
    sprintf(tmp, "n[%.2lf]", aa.getFixedModMass('$'));
    seq += tmp;
  }
  for (k = 0; k<mod.size(); k++){ //check for n-terminal peptide mod
    if (mod[k].pos == -1){
      sprintf(tmp, "n[%.2lf]", mod[k].mass);
      seq += tmp;
    }
  }
  for (j = 0; j<peptide.size(); j++) {
    seq += peptide[j];
    for (k = 0; k<mod.size(); k++){
      if(mod[k].pos<0) continue;
      if (j == (size_t)mod[k].pos){
        sprintf(tmp, "[%.2lf]", mod[k].mass);
        seq += tmp;
      }
    }
  }
  for (k = 0; k<mod.size(); k++){ //check for c-terminal peptide mod
    if (mod[k].pos ==-2){
      sprintf(tmp, "c[%.2lf]", mod[k].mass);
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
    prot+=db[pep.map->at(j).index].short_name;

    if(site>-1){//add the protein site location
      if(sites.size()>0) sites+=";";
      sprintf(tmp, "%d", pep.map->at(j).start + site+1); 
      sites+=tmp;
    }

    //determine if target (if it is currently still decoy)
    if(decoy){
      if (db[pep.map->at(j).index].short_name.find(params->decoy) == string::npos) decoy=false;
    }

  }

}

void KData::writeMzIDDatabase(CMzIdentML& m, KDatabase& db){

  char outPath[1056];
  processPath(params->dbFile, outPath);
  string sDB = outPath;
  CSearchDatabase* m_db = m.dataCollection.inputs.addSearchDatabase(sDB);

  for(int a=0;a<db.getProteinDBSize();a++){
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
  //bool nTerm;
  //bool cTerm;
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


int KData::compareScanBinRev2(const void *p1, const void *p2){
  kScanBin d1 = *(kScanBin *)p1;
  kScanBin d2 = *(kScanBin *)p2;
  if (d1.intensity>d2.intensity) {
    return -1;
  } else if (d1.intensity<d2.intensity) {
    return 1;
  } else {
    return 0;
  }
}
