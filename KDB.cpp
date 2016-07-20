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

//==============================
//  Constructors & Destructors
//==============================
KDatabase::KDatabase(){
  for(int i=0;i<128;i++) AA[i]=0;
  AA['a']=AA['A']=71.0371103;
  AA['c']=AA['C']=103.0091803;
  AA['d']=AA['D']=115.0269385;
  AA['e']=AA['E']=129.0425877;
  AA['f']=AA['F']=147.0684087;
  AA['g']=AA['G']=57.0214611;
  AA['h']=AA['H']=137.0589059;
  AA['i']=AA['I']=113.0840579;
  AA['k']=AA['K']=128.0949557;
  AA['l']=AA['L']=113.0840579;
  AA['m']=AA['M']=131.0404787;
  AA['n']=AA['N']=114.0429222;
  AA['p']=AA['P']=97.0527595;
  AA['q']=AA['Q']=128.0585714;
  AA['r']=AA['R']=156.1011021;
  AA['s']=AA['S']=87.0320244;
  AA['t']=AA['T']=101.0476736;
  AA['v']=AA['V']=99.0684087;
  AA['w']=AA['W']=186.0793065;
  AA['y']=AA['Y']=163.0633228;

  setA=0;
  setB=0;
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
bool  KDatabase::buildDB(char* fname) {
  char    str[10240];
  char*   tok;
  FILE*   f;
  kDB    d;

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
          cout << "  WARNING: " << &d.name[0] << " has a sequence that is too long. It will be skipped." << endl;
        } else {
          vDB.push_back(d);
        }
      }
      d.name=&str[1];
      d.sequence="";
    } else {
      for(unsigned int i=0;i<strlen(str);i++){
        if(AA[str[i]]==0) cout << "  WARNING: " << &d.name[0] << " has an unexpected amino acid character or errant white space." << endl;
        if(str[i]==' ' || str[i]=='\t') continue;
        d.sequence+=toupper(str[i]);
      }
    }
  }
  fclose(f);
  if(d.sequence.length()>65000){
    cout << "  WARNING: " << &d.name[0] << " has a sequence that is too long. It will be skipped." << endl;
   } else {
    vDB.push_back(d);
  }

  cout << "  Total Proteins: " << vDB.size() << endl;
  return true;
}


//buildPeptides creates lists of peptides to search based on the user-defined enzyme rules
bool  KDatabase::buildPeptides(double min, double max, int mis){

  double mass;
  bool bCutMarked;
  bool bNTerm;
  bool bCTerm;

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
    seqSize=vDB[i].sequence.size();
    start=0;
    n=0;
    k=0;
    mc=0;
    mass=18.0105633;
    next=0; //allow for next start site to be amino acid after initial M.

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
        if (mass>min) addPeptide((int)i, (int)start, (int)n - 1, mass, p, vPep, bNTerm, bCTerm, xlSites);

      }

      //Add the peptide mass
      mass+=AA[vDB[i].sequence[start+n]];

      //Check if we cut c-terminal to this AA
      if((start+n+1)<seqSize && enzyme.cutC[vDB[i].sequence[start+n]] && !enzyme.exceptN[vDB[i].sequence[start+n+1]]){
        if(next==-1) next=(int)(start+n);
        if(!bCutMarked) mc++;
        bCutMarked=true;

        //Add the peptide now (if enough mass)
        if(mass>min && mass<max) addPeptide((int)i,(int)start,(int)n,mass,p,vPep,bNTerm,bCTerm,xlSites);

      }

      //Mark sites of cross-linker attachment (if searching for cross-links)
      if(checkAA(p,i,start,n,seqSize,bNTerm,bCTerm)) xlSites++;

      //Check if we are at the end of the sequence
      if((start+n+1)==seqSize) {

        //Add the peptide now (if enough mass)
        if (mass>min && mass<max) addPeptide((int)i, (int)start, (int)n, mass, p, vPep, bNTerm, bCTerm, xlSites);
        if(next>-1) {
          start=next+1;
          n=0;
          mc=0;
          mass=18.0105633;
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
      if(mass>max || mc>mis ) {

        //if we know next cut site
        if(next>-1) {
          start=next+1;
          n=0;
          mc=0;
          mass=18.0105633;
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
          mass=18.0105633;
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

  //merge duplicates
  kPepSort ps;
  vector<kPepSort> vPS;
  for(i=0;i<vPep.size();i++){
    ps.index=(int)i;
    getPeptideSeq(vPep[i].map->at(0).index,vPep[i].map->at(0).start,vPep[i].map->at(0).stop,ps.sequence);
    vPS.push_back(ps);
  }
  qsort(&vPS[0],vPS.size(),sizeof(kPepSort),compareSequence);

  vector<kPeptide> vtp;
  for(i=vPS.size()-1;i>0;i--){
    if(vPS[i].sequence.compare(vPS[i-1].sequence)==0){
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

//==============================
//  Accessors & Modifiers
//==============================
void KDatabase::addFixedMod(char mod, double mass){
  AA[mod]+=mass;
}

kDB& KDatabase::at(const int& i){
  return vDB[i];
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

  string str1=vDB[index].sequence.substr(start,stop-start+1);
  strcpy(str,&str1[0]);
  return true;
}

bool KDatabase::getPeptideSeq(int index, int start, int stop, string& str){

  str=vDB[index].sequence.substr(start,stop-start+1);
  return true;
}

bool KDatabase::getPeptideSeq(kPeptide& p, string& str){

  str=vDB[p.map->at(0).index].sequence.substr(p.map->at(0).start,p.map->at(0).stop-p.map->at(0).start+1);
  return true;
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

void KDatabase::setXLTable(char** arr, int szA, int szB){
  for (int i = 0; i < szA; i++){
    for (int j = 0; j < szB; j++){
      xlTable[i][j] = arr[i][j];
    }
  }
}

bool KDatabase::setXLType(int a, int b){
  setA=a;
  setB=b;
  return true;
}

//==============================
//  Private Functions
//==============================
void KDatabase::addPeptide(int index, int start, int len, double mass, kPeptide& p, vector<kPeptide>& vP, bool bN, bool bC, char xlSites){
  kPepMap  pm;
                           
  pm.index=index;
  pm.start=start;
  pm.stop=start+len;
  
  p.nTerm=bN;
  p.cTerm=bC;
  p.xlSites=xlSites;
  p.mass=mass;
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
int KDatabase::compareMass(const void *p1, const void *p2){
  const kPeptide d1 = *(kPeptide *)p1;
  const kPeptide d2 = *(kPeptide *)p2;
  if(d1.mass<d2.mass) {
		return -1;
	} else if(d1.mass>d2.mass) {
  	return 1;
  } else {
	  return 0;
  }
}

int KDatabase::compareSequence(const void *p1, const void *p2){
  const kPepSort d1 = *(kPepSort *)p1;
  const kPepSort d2 = *(kPepSort *)p2;
  return d1.sequence.compare(d2.sequence);
}

