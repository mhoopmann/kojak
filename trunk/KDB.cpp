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

KDatabase::KDatabase(){
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
}

kDB& KDatabase::operator[ ](const int& i){
  return vDB[i];
}

bool  KDatabase::buildDB(char* fname) {
  char    str[1024];
  FILE*   f;
  kDB    d;

  d.name="NIL";
  cout << "Reading DB...";

  vDB.clear();
  f=fopen(fname,"rt");
  if(f==NULL) return false;
  while(!feof(f)){
    fgets(str,1024,f);
    if(feof(f)) continue;
    while(strlen(str)>0){
      if(str[strlen(str)-1]=='\n') {
        str[strlen(str)-1]='\0';
        continue;
      }
      if(str[strlen(str)-1]==' ') {
        str[strlen(str)-1]='\0';
        continue;
      }
      break;
    }
    if(str[0]=='>') {
      if(d.name.compare("NIL")!=0) {
        if(d.sequence.length()>65000){
          cout << "WARNING: " << &d.name[0] << " has a sequence that is too long. It will be skipped." << endl;
        } else {
          vDB.push_back(d);
        }
      }
      d.name=str;
      d.sequence="";
    } else {
      for(unsigned int i=0;i<strlen(str);i++){
        d.sequence+=toupper(str[i]);
      }
    }
  }
  fclose(f);
  if(d.sequence.length()>65000){
    cout << "WARNING: " << &d.name[0] << " has a sequence that is too long. It will be skipped." << endl;
   } else {
    vDB.push_back(d);
  }

  cout << "Done!  Total Proteins: " << vDB.size() << endl;
  return true;
}

bool  KDatabase::buildPeptides(double min, double max, int mis){
  //char str[256];

  double mass;

  int mc;
  int next;
  int tmp;

  kPepMap  pm;
  kPeptide p;
  
  unsigned int DBSize=vDB.size();
  unsigned int i;
  unsigned int k;
  unsigned int n;
  unsigned int seqSize;
  unsigned int start;

  vPep.clear();
  vPepK.clear();

  for(i=0;i<DBSize;i++){

    seqSize=vDB[i].sequence.size();
    start=0;
    n=0;
    k=0;
    mc=0;
    mass=18.0105633;
    next=0; //allow for next start site to be amino acid after initial M.

    pm.index=i;
    p.vK->clear();

    while(true){

      //Check if we are at the end of the sequence
      if((start+n)==seqSize) {
        //Add to database
        if(mass>min && mass<max) {
        p.mass=mass;
          pm.start=start;
          pm.stop=start+n;
          p.map->clear();
          p.map->push_back(pm);
          if(p.vK->size()>0) vPepK.push_back(p);
          else vPep.push_back(p);
        }
        if(next==-1) break;
        start=next+1;
        n=0;
        mc=0;
        mass=18.0105633;
        p.vK->clear();
        next=-1;
        continue;
      }

      //add amino acid mass
      mass+=AA[vDB[i].sequence[start+n]];

      //mark position if a lysine is the current aa
      if(vDB[i].sequence[start+n]=='K') p.vK->push_back(start+n);
      else if(start+n==0) p.vK->push_back(start+n); //mark N-terminus, too
      else if(start+n==1 && start==1) p.vK->push_back(start+n); //mark Met-removed N-terminus, too
				
			//skip non-tryptic peptides (ignoring proline rule)
      if(vDB[i].sequence[start+n]!='K' && vDB[i].sequence[start+n]!='R') {
        n++;
        continue;
      } else {
        if(next==-1) next=start+n;
        mc++;
      }
 
      if(mass>min && mass<max && (mc-1)<=mis){
        //Add to database
        p.mass=mass;
        pm.start=start;
        pm.stop=start+n;
        p.map->clear();
        p.map->push_back(pm);
        if(vDB[i].sequence[start+n]=='K' && start+n!=seqSize-1) {
          //unmark last lysine if it is the peptide cut site, unless it is the C-terminus
          tmp=p.vK->at(p.vK->size()-1);
          p.vK->pop_back();
          if(p.vK->size()>0) vPepK.push_back(p);
          else vPep.push_back(p);
          p.vK->push_back(tmp); //add it back for when peptide continues to build
        } else if(start+n!=seqSize-1) {
          //add c-terminal peptides to the database on the next iteration
          //otherwise, add this non-c-terminal peptide now.
          if(p.vK->size()>0) vPepK.push_back(p);
          else vPep.push_back(p);
        }

        n++;
        continue;
      }

      //When we reach the max peptide
      if(mass>max || mc>mis ) {

        //if we know next K or R
        if(next>-1) {
          start=next+1;
          n=0;
          mc=0;
          mass=18.0105633;
          p.vK->clear();
          next=-1;

        //Otherwise, continue scanning until it is found
        } else {
          while((start+n)<seqSize-1){
            n++;
            if(vDB[i].sequence[start+n]=='K' && vDB[i].sequence[start+n]=='R') {
              next=start+n;
              break;
            }
          }
          if(next<0) break;

          start=next+1;
          n=0;
          mc=0;
          mass=18.0105633;
          p.vK->clear();
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
  vector<kPepSort> vPSK;
  for(i=0;i<vPep.size();i++){
    ps.index=i;
    getPeptideSeq(vPep[i].map->at(0).index,vPep[i].map->at(0).start,vPep[i].map->at(0).stop,ps.sequence);
    vPS.push_back(ps);
  }
  for(i=0;i<vPepK.size();i++){
    ps.index=i;
    getPeptideSeq(vPepK[i].map->at(0).index,vPepK[i].map->at(0).start,vPepK[i].map->at(0).stop,ps.sequence);
    vPSK.push_back(ps);
  }
  qsort(&vPS[0],vPS.size(),sizeof(kPepSort),compareSequence);
  qsort(&vPSK[0],vPSK.size(),sizeof(kPepSort),compareSequence);

  vector<kPeptide> vtp;
  for(i=vPS.size()-1;i>0;i--){
    if(vPS[i].sequence.compare(vPS[i-1].sequence)==0){
      for(k=0;k<vPep[vPS[i].index].map->size();k++){
        vPep[vPS[i-1].index].map->push_back(vPep[vPS[i].index].map->at(k));
      }
      vPep[vPS[i].index].mass=-1;
    }
  }
  for(i=0;i<vPep.size();i++){
    if(vPep[i].mass>0) vtp.push_back(vPep[i]);
  }
  vPep.clear();
  for(i=0;i<vtp.size();i++) vPep.push_back(vtp[i]);

  for(i=vPSK.size()-1;i>0;i--){
    if(vPSK[i].sequence.compare(vPSK[i-1].sequence)==0){
      for(k=0;k<vPepK[vPSK[i].index].map->size();k++){
        vPepK[vPSK[i-1].index].map->push_back(vPepK[vPSK[i].index].map->at(k));
      }
      vPepK[vPSK[i].index].mass=-1;
    }
  }
  vtp.clear();
  for(i=0;i<vPepK.size();i++){
    if(vPepK[i].mass>0) vtp.push_back(vPepK[i]);
  }
  vPepK.clear();
  for(i=0;i<vtp.size();i++) vPepK.push_back(vtp[i]);

  cout << vPep.size() << " peptides to search." << endl;
  cout << vPepK.size() << " linkable peptides to search." << endl;

  qsort(&vPep[0],vPep.size(),sizeof(kPeptide),compareMass);
  qsort(&vPepK[0],vPepK.size(),sizeof(kPeptide),compareMass);

  //Reporting list
  /*
  char str[256];
  for(i=0;i<vPep.size();i++){
    getPeptideSeq(vPep[i].map->at(0).index,vPep[i].map->at(0).start,vPep[i].map->at(0).stop,str);
    cout << i << ", " << vPep[i].mass << "\t" << str;
    for(k=0;k<vPep[i].vK->size();k++) cout << "\t" << vPep[i].vK->at(k);
    cout << endl;
    if(i==10) break;
  }
  for(i=0;i<vPepK.size();i++){
    getPeptideSeq(vPepK[i].map->at(0).index,vPepK[i].map->at(0).start,vPepK[i].map->at(0).stop,str);
    cout << i << ", " << vPepK[i].mass << "\t" << str;
    for(k=0;k<vPepK[i].vK->size();k++) cout << "\t" << vPepK[i].vK->at(k);
    cout << endl;
    //if(i==10) break;
  }

  exit(0);
  */
  return true;

}

void KDatabase::combo(){
  qsort(&vPep[0],vPep.size(),sizeof(kPeptide),compareMass);
  unsigned int i,j,k;

  k=0;
  for(i=0;i<vPep.size();i++){
    if(vPep[i].mass<8000) k++;
  }

  for(i=0;i<vPepK.size();i++){
    //check mono links
    if(vPepK[i].mass<8000) k+=3;

    //check crosslinks
    for(j=i+1;j<vPepK.size();j++){
      if(vPepK[i].mass+vPepK[j].mass<8000) k++;
      else break;
    }
  }
  cout << "Combos to check: " << k << endl;


  /*
  for(i=0;i<vPep.size();i++){
    getPeptideSeq(vPep[i].index,vPep[i].start,vPep[i].stop,str);
    cout << vPep[i].mass << "\t" << str << endl;
  }
  */
}

void KDatabase::addFixedMod(char mod, double mass){
  AA[mod]+=mass;
}

kPeptide& KDatabase::getPeptide(int index, bool linkable){
  if(linkable) return vPepK[index];
  else return vPep[index];
}

vector<kPeptide>& KDatabase::getPeptideList(bool linkable){
  if(linkable) return vPepK;
  else return vPep;
}

int KDatabase::getPeptideListSize(bool linkable){
  if(linkable) return vPepK.size();
  else return vPep.size();
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

