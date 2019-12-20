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

#include "KDB.h"

using namespace std;

//==============================
//  Constructors & Destructors
//==============================
KDatabase::KDatabase(){
  for(int i=0;i<128;i++) {
    AA[i]=0;
    AAn15[i]=0;
  }
  AA['A']=71.0371103;
  AA['C']=103.0091803;
  AA['D']=115.0269385;
  AA['E']=129.0425877;
  AA['F']=147.0684087;
  AA['G']=57.0214611;
  AA['H']=137.0589059;
  AA['I']=113.0840579;
  AA['K']=128.0949557;
  AA['L']=113.0840579;
  AA['M']=131.0404787;
  AA['N']=114.0429222;
  AA['P']=97.0527595;
  AA['Q']=128.0585714;
  AA['R']=156.1011021;
  AA['S']=87.0320244;
  AA['T']=101.0476736;
  AA['U']=150.9536303;
  AA['V']=99.0684087;
  AA['W']=186.0793065;
  AA['Y']=163.0633228;
  AA['n'] = 1.00782503;
  AA['c'] = 17.00273963;
  AA['$'] = 1.00782503;
  AA['%'] = 17.00273963;

  AAn15['A'] = 72.0341452;
  AAn15['C'] = 104.0062152;
  AAn15['D'] = 116.0239734;
  AAn15['E'] = 130.0396226;
  AAn15['F'] = 148.0654436;
  AAn15['G'] = 58.018496;
  AAn15['H'] = 140.0500106;
  AAn15['I'] = 114.0810928;
  AAn15['K'] = 130.0890255;
  AAn15['L'] = 114.0810928;
  AAn15['M'] = 132.0375136;
  AAn15['N'] = 116.036992;
  AAn15['P'] = 98.0497944;
  AAn15['Q'] = 130.0526412;
  AAn15['R'] = 160.0892417;
  AAn15['S'] = 88.0290593;
  AAn15['T'] = 102.0447085;
  AAn15['U'] = 151.9506652;
  AAn15['V'] = 100.0654436;
  AAn15['W'] = 188.0733763;
  AAn15['Y'] = 164.0603577;
  AAn15['n'] = 1.00782503;
  AAn15['c'] = 17.00273963;
  AAn15['$'] = 1.00782503;
  AAn15['%'] = 17.00273963;

  fixMassPepC=0;
  fixMassPepN=0;
  fixMassProtC=0;
  fixMassProtN=0;

  n15Label="NON15LABEL";
  klog=NULL;
  linkablePepCount=0;
}

KDatabase::~KDatabase(){
  klog=NULL;
}

//==============================
//  Operators
//==============================
kDB& KDatabase::operator[ ](const int& i){
  return vDB[i];
}

//==============================
//  User Functions
//==============================

//buildDB reads in a FASTA file and stores it in memory.
bool  KDatabase::buildDB(const char* fname, string decoyStr) {
  char    str[10240];
  char*   tok;
  FILE*   f;
  kDB    d;
  char   c;

  d.name="NIL";

  vDB.clear();
  f=fopen(fname,"rt");
  if(f==NULL) return false;
  while(!feof(f)){
    if(fgets(str,10240,f)==NULL) continue;
    if(strlen(str)>0){
      tok=strtok(str,"\r\n");
      if(tok==NULL) continue;
      strcpy(str,tok);
    }
    if(str[0]=='>') {
      if(d.name.compare("NIL")!=0) {
        if(d.sequence.length()>65000){
          if (klog != NULL) klog->addDBWarning(d.name + " has a sequence that is too long. It will be skipped.");
          else cout << "  WARNING: " << &d.name[0] << " has a sequence that is too long. It will be skipped." << endl;
        } else {
          if(d.name.find(decoyStr)!=string::npos) d.decoy=true;
          else d.decoy=false;
          vDB.push_back(d);
        }
      }
      d.name=&str[1];
      d.sequence="";
    } else {
      for(unsigned int i=0;i<strlen(str);i++){
        c=toupper(str[i]);
        if(AA[c]==0) {
          if (klog != NULL) klog->addDBWarning(d.name+" has an unexpected amino acid character or errant white space: '" + c + "'");
          else cout << "  WARNING: " << &d.name[0] << " has an unexpected amino acid character or errant white space: '" << c << "'" << endl;
        }
        if(c==' ' || c=='\t') continue;
        if (AA[c] == 0) {
          if (klog != NULL) {
            string tmpStr="Mass of '";
            klog->addDBWarning(tmpStr + c + "' is currently set to 0. Consider revising with the aa_mass parameter.");
          } else cout << "  WARNING: Mass of '" << c << "' is currently set to 0. Consider revising with the aa_mass parameter." << endl;
        }
        d.sequence+=c;
      }
    }
  }
  fclose(f);
  if(d.sequence.length()>65000){
    if (klog != NULL) klog->addDBWarning(d.name + " has a sequence that is too long. It will be skipped.");
    else cout << "  WARNING: " << &d.name[0] << " has a sequence that is too long. It will be skipped." << endl;
   } else {
    if (d.name.find(decoyStr) != string::npos) d.decoy = true;
    else d.decoy = false;
    vDB.push_back(d);
  }

  cout << "  Total Proteins: " << vDB.size() << endl;

  return true;
}

void KDatabase::buildDecoy(string decoy_label) {

  typedef struct clips {
    int start;
    int stop;
  } clips;

  size_t i, j,sz;
  vector<clips> cut;
  clips c;

  sz = vDB.size();
  for (i = 0; i < sz; i++) {

    cut.clear();
    c.start = -1;
    //if (vDB[i].sequence[0] == 'M')j = 1; //leave leading methionines in place
    //else j = 0;
    for (j = 1; j < vDB[i].sequence.size(); j++) {
      if (enzyme.cutN[vDB[i].sequence[j]] || enzyme.cutC[vDB[i].sequence[j]]) {
        if (c.start > -1) { //mark the space in between
          c.stop = (int)j - 1;
          cut.push_back(c);
          c.start = -1; //reset
        }
      } else {
        if (c.start < 0) c.start = (int)j;
      }
    }
    if (!enzyme.cutN[vDB[i].sequence[j - 1]] && !enzyme.cutC[vDB[i].sequence[j - 1]]) { //check last amino acid
      c.stop = (int)j - 1;
      cut.push_back(c);
    }

    //reverse the sequences
    string rev;
    kDB decoy = vDB[i];
    decoy.name = decoy_label + "_" + decoy.name;
    for (j = 0; j < cut.size(); j++) {
      rev.clear();

      //adjust ends for restrictive sites
      while (enzyme.exceptN[decoy.sequence[cut[j].start]] || enzyme.exceptC[decoy.sequence[cut[j].start]]) {
        cut[j].start++;
        if (cut[j].start == decoy.sequence.size()) break;
      }
      while (enzyme.exceptN[decoy.sequence[cut[j].stop]] || enzyme.exceptC[decoy.sequence[cut[j].stop]]) {
        cut[j].stop--;
        if (cut[j].stop == -1) break;
      }
      if (cut[j].start == decoy.sequence.size()) continue; //skip when out of bounds
      if (cut[j].stop == -1) continue; //skip when out of bounds
      if (cut[j].stop <= cut[j].start) continue; //skip if nothing will happen


      for (size_t k = cut[j].stop; k >= cut[j].start; k--) {
        rev += decoy.sequence[k];
        if (k == 0) break;
      }
      decoy.sequence.replace(cut[j].start, (size_t)cut[j].stop - (size_t)cut[j].start + 1, rev);
    }

    vDB.push_back(decoy);
  }

  cout << "  Adding Kojak-generated decoys. New Total Proteins: " << vDB.size() << endl;
}

//buildPeptides creates lists of peptides to search based on the user-defined enzyme rules
bool KDatabase::buildPeptides(double min, double max, int mis){

  double mass;
  bool bCutMarked;
  bool bNTerm;
  bool bCTerm;
  bool bN15;

  int mc;
  int next;

  char xlSites;

  kPepMap  pm;
  kPeptide p;
  
  size_t DBSize=vDB.size();
  size_t i;
  size_t k;
  size_t n;
  size_t seqSize;
  size_t start;

  vPep.clear();

  for(i=0;i<DBSize;i++){
    if (vDB[i].name.find(n15Label) == string::npos) bN15=false;
    else bN15=true;
    seqSize=vDB[i].sequence.size();
    start=0;
    n=0;
    k=0;
    mc=0;
    mass=18.0105633+fixMassPepN+fixMassProtN;
    if(vDB[i].sequence[0]=='M') next=0; //allow for next start site to be amino acid after initial M.
    else next = -1;

    pm.index=(int)i;
    bNTerm=false;
    bCTerm=false;
    xlSites=0;

    while(true){

      bCutMarked=false;

      //Check if we cut n-terminal to this AA
      if(n>0 && enzyme.cutN[vDB[i].sequence[start+n]] && !enzyme.exceptC[vDB[i].sequence[start+n-1]]){
        if(next==-1) next=(int)start+(int)n-1;
        if(!bCutMarked) mc++;
        bCutMarked=true;

        //Add the peptide now (if enough mass)
        if ((mass+fixMassPepC)>min) addPeptide((int)i, (int)start, (int)n - 1, mass+fixMassPepC, p, vPep, bNTerm, bCTerm, bN15, xlSites);

      }

      //Add the peptide mass
      if (bN15) mass += AAn15[vDB[i].sequence[start + n]];
      else mass+=AA[vDB[i].sequence[start+n]];

      //Check if we cut c-terminal to this AA
      if((start+n+1)<seqSize && enzyme.cutC[vDB[i].sequence[start+n]] && !enzyme.exceptN[vDB[i].sequence[start+n+1]]){
        if(next==-1) next=(int)(start+n);
        if(!bCutMarked) mc++;
        bCutMarked=true;

        //Add the peptide now (if enough mass)
        if((mass+fixMassPepC)>min && (mass+fixMassPepC)<max) addPeptide((int)i,(int)start,(int)n,mass+fixMassPepC,p,vPep,bNTerm,bCTerm, bN15, xlSites);

      }

      //Mark sites of cross-linker attachment (if searching for cross-links)
      if(checkAA(p,i,start,n,seqSize,bNTerm,bCTerm)) xlSites++;

      //Check if we are at the end of the sequence
      if((start+n+1)==seqSize) {

        //Add the peptide now (if enough mass)
        if ((mass+fixMassPepC+fixMassProtC)>min && (mass+fixMassPepC+fixMassProtC)<max) addPeptide((int)i, (int)start, (int)n, mass+fixMassPepC+fixMassProtC, p, vPep, bNTerm, bCTerm, bN15, xlSites);
        if(next>-1) {
          start=next+1;
          n=0;
          mc=0;
          mass=18.0105633+fixMassPepN;
          bNTerm = false;
          bCTerm = false;
          xlSites=0;
          next=-1;
          continue;
        } else {
          break;
        }

      }

      //Check if we exceeded peptide mass
      //Check if we exceeded the number of missed cleavages
      if((mass+fixMassPepC)>max || mc>mis ) {

        //if we know next cut site
        if(next>-1) {
          start=next+1;
          n=0;
          mc=0;
          mass=18.0105633+fixMassPepN;
          bNTerm = false;
          bCTerm = false;
          xlSites = 0;
          next=-1;

        //Otherwise, continue scanning until it is found
        } else {
          while((start+n)<seqSize-1){
            n++;
            if(n>0 && enzyme.cutN[vDB[i].sequence[start+n]] && !enzyme.exceptC[vDB[i].sequence[start+n-1]]){
              next=(int)(start+n);
              break;
            } else if((start+n+1)<seqSize && enzyme.cutC[vDB[i].sequence[start+n]] && !enzyme.exceptN[vDB[i].sequence[start+n+1]]){
              next=(int)(start+n);
              break;
            } 
          }
          if(next<0) break;

          start=next+1;
          n=0;
          mc=0;
          mass=18.0105633+fixMassPepN;
          bNTerm = false;
          bCTerm = false;
          xlSites = 0;
          next=-1;
        }  
      } else {
        n++;
        continue;
      }

    }

  }

  /* Diagnostic block - probably safe to give the boot
  cout << vPep.size() << endl;
  kPepSort ps2;
  for (i = 0; i < vPep.size(); i++) {
    getPeptideSeq(vPep[i].map->at(0).index, vPep[i].map->at(0).start, vPep[i].map->at(0).stop, ps2.sequence);
    cout << ps2.sequence << "\t" << vPep[i].map->at(0).start << endl;
  }
  */

  //merge duplicates
  kPepSort ps;
  ps.index=0;
  ps.sequence.clear();
  vector<kPepSort> vPS;
  for(i=0;i<vPep.size();i++){
    ps.index=(int)i;
    ps.n15=vPep[i].n15;
    getPeptideSeq(vPep[i].map->at(0).index,vPep[i].map->at(0).start,vPep[i].map->at(0).stop,ps.sequence);
    vPS.push_back(ps);
  }
  //qsort(&vPS[0],vPS.size(),sizeof(kPepSort),compareSequence);
  sort(vPS.begin(),vPS.end(),compareSequenceB);

  vector<kPeptide> vtp;
  for(i=vPS.size()-1;i>0;i--){
    if(vPS[i].sequence.compare(vPS[i-1].sequence)==0 && vPS[i].n15==vPS[i-1].n15){
      for(k=0;k<vPep[vPS[i].index].map->size();k++){
        vPep[vPS[i-1].index].map->push_back(vPep[vPS[i].index].map->at(k));
        if (vPep[vPS[i].index].cTerm) vPep[vPS[i - 1].index].cTerm = vPep[vPS[i].index].cTerm;
        if (vPep[vPS[i].index].nTerm) vPep[vPS[i - 1].index].nTerm = vPep[vPS[i].index].nTerm;
        if (vPep[vPS[i].index].xlSites>vPep[vPS[i - 1].index].xlSites) vPep[vPS[i - 1].index].xlSites = vPep[vPS[i].index].xlSites;
      }
      vPep[vPS[i].index].mass=-1;
    }
  }
  for(i=0;i<vPep.size();i++){
    if(vPep[i].mass>0) vtp.push_back(vPep[i]);
  }
  vPep.clear();
  n=0;
  for(i=0;i<vtp.size();i++) {
    if (vtp[i].xlSites>0)n++;
    vPep.push_back(vtp[i]);
  }

  cout << "  " << vPep.size() << " peptides to search (" << n << " linkable)." << endl;
  qsort(&vPep[0],vPep.size(),sizeof(kPeptide),compareMass);
  linkablePepCount=(int)n;

  //Reporting list
  /*
  char str[256];
  for(i=0;i<vPep.size();i++){
    getPeptideSeq(vPep[i].map->at(0).index,vPep[i].map->at(0).start,vPep[i].map->at(0).stop,str);
    cout << i << ", " << vPep[i].mass << "\t" << str;
    for(k=0;k<vPep[i].vA->size();k++) cout << "\t" << vPep[i].vA->at(k);
    cout << endl;
    //if(i==10) break;
  }
  for(i=0;i<vPepK.size();i++){
    getPeptideSeq(vPepK[i].map->at(0).index,vPepK[i].map->at(0).start,vPepK[i].map->at(0).stop,str);
    cout << i << ", " << vPepK[i].mass << "\t" << str;
    for(k=0;k<vPepK[i].vA->size();k++) cout << "\t" << vPepK[i].vA->at(k);
    cout << endl;
    //if(i==10) break;
  }

  exit(0);
  */
  return true;

}

void KDatabase::exportDB(string fName) {
  size_t i;
  FILE* f = fopen(fName.c_str(), "wt");
  for (i = 0; i < vDB.size(); i++) {
    fprintf(f,">%s\n",vDB[i].name.c_str());
    fprintf(f,"%s\n",vDB[i].sequence.c_str());
  }
  fclose(f);
}

//==============================
//  Accessors & Modifiers
//==============================
void KDatabase::addFixedMod(char mod, double mass){
  if (mod == 'n') fixMassPepN=mass;
  else if(mod=='c') fixMassPepC=mass;
  else if(mod=='$') fixMassProtN=mass;
  else if(mod=='%') fixMassProtC=mass;
  else {
    AA[mod]+=mass;
    AAn15[mod]+=mass;
  }
}

kDB& KDatabase::at(const int& i){
  return vDB[i];
}

double KDatabase::getAAMass(char aa, bool n15){
  if (n15) return AAn15[aa];
  return AA[aa];
}

kEnzymeRules& KDatabase::getEnzymeRules(){
  return enzyme;
}

kPeptide& KDatabase::getPeptide(int index){
  return vPep[index];
}

vector<kPeptide>* KDatabase::getPeptideList(){
  return &vPep;
}

int KDatabase::getPeptideListSize(){
  return (int)vPep.size();
}

bool KDatabase::getPeptideSeq(int index, int start, int stop, char* str){
  if ((size_t)index>vDB.size()) return false;
  string str1=vDB[index].sequence.substr(start,stop-start+1);
  strcpy(str,&str1[0]);
  return true;
}

bool KDatabase::getPeptideSeq(int index, int start, int stop, string& str){
  if((size_t)index>vDB.size()) return false;
  str=vDB[index].sequence.substr(start,stop-start+1);
  return true;
}

bool KDatabase::getPeptideSeq(kPeptide& p, string& str){
  str=vDB[p.map->at(0).index].sequence.substr(p.map->at(0).start,p.map->at(0).stop-p.map->at(0).start+1);
  return true;
}

bool KDatabase::getPeptideSeq(int pepIndex, string& str){
  if((size_t)pepIndex>vPep.size()) return false;
  kPeptide p = vPep[(size_t)pepIndex];
  str = vDB[p.map->at(0).index].sequence.substr(p.map->at(0).start, p.map->at(0).stop - p.map->at(0).start + 1);
  return true;
}

int KDatabase::getProteinDBSize(){
  return (int)vDB.size();
}

void KDatabase::setAAMass(char aa, double mass, bool n15){
  if (n15) AAn15[aa] = mass;
  else AA[aa] = mass;
}

bool KDatabase::setEnzyme(char* str){
  bool stateNTerm=false;
  int stateRule=0;

  unsigned int i;
  for(i=0;i<128;i++){
    enzyme.cutC[i]=enzyme.cutN[i]=enzyme.exceptC[i]=enzyme.exceptN[i]=false;
  }

  for(i=0;i<strlen(str);i++){
    switch(str[i]){
      case '[':
        if(stateRule>0){
          cout << "Error in enzyme string: " << str << endl;
          return false;
        }
        stateRule=1;
        break;
      case ']':
        if(stateRule!=1){
          cout << "Error in enzyme string: " << str << endl;
          return false;
        }
        stateRule=0;
        break;
      case '{':
        if(stateRule>0){
          cout << "Error in enzyme string: " << str << endl;
          return false;
        }
        stateRule=2;
        break;
      case '}':
        if(stateRule!=2){
          cout << "Error in enzyme string: " << str << endl;
          return false;
        }
        stateRule=0;
        break;
      case '|':
        if(stateNTerm) {
          cout << "Error in enzyme string: " << str << endl;
          return false;
        }
        stateNTerm=true;
        break;
      default:
        if(stateRule==0){
          cout << "Error in enzyme string: " << str << endl;
          return false;
        }
        if(stateNTerm){
          if(stateRule==1) enzyme.cutN[str[i]]=true;
          else enzyme.exceptN[str[i]]=true;
        } else {
          if(stateRule==1) enzyme.cutC[str[i]]=true;
          else enzyme.exceptC[str[i]]=true;
        }
        break;
    }
  }

  //for(i=65;i<90;i++){
  //  cout << (char)i << "\t" << enzyme.cutC[i] << "\t" << enzyme.cutN[i] << "\t" << enzyme.exceptC[i] << "\t" << enzyme.exceptN[i] << endl;
  //}

  return true;

}

void KDatabase::setLog(KLog* c){
  klog=c;
}

void KDatabase::setN15Label(char* str){
  n15Label=str;
}

void KDatabase::setXLTable(char** arr, int szA, int szB){
  for (int i = 0; i < szA; i++){
    for (int j = 0; j < szB; j++){
      xlTable[i][j] = arr[i][j];
    }
  }
}

//==============================
//  Private Functions
//==============================
void KDatabase::addPeptide(int index, int start, int len, double mass, kPeptide& p, vector<kPeptide>& vP, bool bN, bool bC, bool bN15, char xlSites){
  kPepMap  pm;
                           
  pm.index=index;
  pm.start=start;
  pm.stop=start+len;
  
  p.nTerm=bN;
  p.cTerm=bC;
  p.xlSites=xlSites;
  p.mass=mass;
  p.n15=bN15;
  p.map->clear();
  p.map->push_back(pm);
  vP.push_back(p);

  //char str[256];
  //getPeptideSeq(p.map->at(0).index,p.map->at(0).start,p.map->at(0).stop,str);
  //cout << "Adding: " << str << "\t" << p.vA->size()+p.vB->size() << endl;

}

bool KDatabase::checkAA(kPeptide& p, size_t i, size_t start, size_t n, size_t seqSize, bool& bN, bool& bC){
  if (start + n == 0){
    bN=true;
    if (xlTable['n'][0]>-1) return true;
  }
  if (start + n == 1 && start == 1) {
    bN=true;
    if (xlTable['n'][0]>-1) return true;
  }
  if (start + n == seqSize-1) {
    bC=true;
    if (xlTable['c'][0]>-1) return true;
  }
  if (xlTable[vDB[i].sequence[start + n]][0] > -1) return true;
  return false;
}

//==============================
//  Utility Functions
//==============================
int KDatabase::compareMass(const void *p1, const void *p2){ //sort high to low
  const kPeptide d1 = *(kPeptide *)p1;
  const kPeptide d2 = *(kPeptide *)p2;
  if(d1.mass<d2.mass) {
		return 1;
	} else if(d1.mass>d2.mass) {
  	return -1;
  } else {
	  return 0;
  }
}

int KDatabase::compareSequence(const void *p1, const void *p2){
  const kPepSort d1 = *(kPepSort *)p1;
  const kPepSort d2 = *(kPepSort *)p2;
  return d1.sequence.compare(d2.sequence);
}

bool KDatabase::compareSequenceB(const kPepSort& p1, const kPepSort& p2){
  int i=p1.sequence.compare(p2.sequence);
  if(i==0) return p1.n15<p2.n15;
  else return i<0;
}


