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

/*============================
  Constructors
============================*/
KData::KData(){
  bScans=NULL;
  params=NULL;
  xlTable = new char*[128];
  for (int i = 0; i<128; i++) xlTable[i] = new char[20];
}

KData::KData(kParams* p){
  bScans=NULL;
  params=p;
  size_t i;
  for(i=0;i<p->fMods->size();i++) aa.addFixedMod((char)p->fMods->at(i).index,p->fMods->at(i).mass);
  xlTable = new char*[128];
  for (int j = 0; j<128; j++) xlTable[j] = new char[20];
}

KData::~KData(){
  params=NULL;
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


/*============================
  Functions
============================*/
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
  /*
  for (i = 0; i < motifCount; i++){
    cout << "Motif " << i << ": " << &motifs[i].motif[0] << endl;
    for(k=0;k<10;k++){
      if(motifs[i].counterMotif[k]<0) break;
      cout << "  XLink index: " << (int)motifs[i].xlIndex[k] << "\tCounter Motif: " << (int)motifs[i].counterMotif[k] << endl;
    }
  }
  */
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

//This function tries to assign best possible 18O2 and 18O4 precursor ion mass values
//for all MS2 spectra
bool KData::mapPrecursors(){
  
  int iPercent=0;
  int iTmp;
  
  unsigned int i;
  int j,k;

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

    //If instrument determined precursors are preferred, only compute precursors if none supplied
    if(params->preferPrecursor==1 && spec[i].sizePrecursor()>0){
      prePre++;
      specCounts++;
      peakCounts+=spec[i].size();
      continue;
    }

    //if (spec[i].getScanNumber() == 25028) {
    //  cout << "Current precursors: " << spec[i].sizePrecursor() << endl;
    //  for (j = 0; j<spec[i].sizePrecursor(); j++) cout << spec[i].getPrecursor(j).monoMass << " " << spec[i].getPrecursor(j).charge << endl;
    //}

    //Find precursor using object function. Take results and copy them to spectra
    ret=pre.getSpecRange(spec[i]);

    //if (spec[i].getScanNumber() == 25028) {
    //  cout << "Current precursors: " << spec[i].sizePrecursor() << "\tret = " << ret << endl;
    //  for (j = 0; j<spec[i].sizePrecursor(); j++) cout << spec[i].getPrecursor(j).monoMass << " " << spec[i].getPrecursor(j).charge << endl;
    //}

    if(ret>0){

      //if supplementing instrument predicted precursor, then chance for precursor
      //to be seen twice (first by instrument, then by Hardklor). Keep Hardklor result.
      if(spec[i].getInstrumentPrecursor() && spec[i].sizePrecursor()>1){
        if(spec[i].getPrecursor(0).charge==spec[i].getPrecursor(1).charge && 
           fabs((spec[i].getPrecursor(0).monoMass-spec[i].getPrecursor(1).monoMass)/spec[i].getPrecursor(0).monoMass*1e6)<10.0 ){
          spec[i].erasePrecursor(0);
          spec[i].setInstrumentPrecursor(false);
        }
      }

      //if precursor prediction doesn't overlap selected ion, predict additional
      //precursors using presumed charge states and the selected ion.
      if(ret==2) {
        pre.estimatePrecursor(spec[i]);

        //if (spec[i].getScanNumber() == 25028) {
        //  cout << "Current precursors: " << spec[i].sizePrecursor() << endl;
        //  for (j = 0; j<spec[i].sizePrecursor(); j++) cout << spec[i].getPrecursor(j).monoMass << " " << spec[i].getPrecursor(j).charge << endl;
        //}

        //if supplementing instrument predicted precursor, then chance for precursor
        //to be seen twice (first by instrument, then by charge prediction). Keep instrument.
        if(spec[i].getInstrumentPrecursor() && spec[i].sizePrecursor()>1){
          for(k=1;k<spec[i].sizePrecursor();k++){
            if(spec[i].getPrecursor(0).charge==spec[i].getPrecursor(k).charge && 
               fabs((spec[i].getPrecursor(0).monoMass-spec[i].getPrecursor(k).monoMass)/spec[i].getPrecursor(0).monoMass*1e6)<10.0 ){
              spec[i].erasePrecursor(k);
              k--;
            }
          }
        }

      }

    } else { 
      
      //If no precursors found, estimate using mercury and selected mass
      pre.estimatePrecursor(spec[i]);
      
      //if supplementing instrument predicted precursor, then chance for precursor
      //to be seen twice (first by instrument, then by charge prediction). Keep instrument.
      if(spec[i].getInstrumentPrecursor() && spec[i].sizePrecursor()>1){
        for(k=1;k<spec[i].sizePrecursor();k++){
          if(spec[i].getPrecursor(0).charge==spec[i].getPrecursor(k).charge && 
             fabs(spec[i].getPrecursor(0).monoMass-spec[i].getPrecursor(k).monoMass)<0.01 ){
            spec[i].erasePrecursor(k);
            k--;
          }
        }
      }

    }

    if(spec[i].sizePrecursor()>0){
      foundPre++;
      specCounts++;
      peakCounts+=spec[i].size();
    }    

  }
 

  //Finalize the progress
  printf("\b\b\b100%%");
  cout << endl;

  cout << "  " << specCounts << " spectra with " << peakCounts << " peaks will be analyzed." << endl;

  //Build mass list - this orders all precursor masses, with an index pointing to the actual
  //array position for the spectrum. This is because all spectra will have more than 1
  //precursor mass
  massList.clear();
  for(i=0;i<spec.size();i++){
    m.index=i;
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
    sh.modInfo.addMod((int)r.mods1[i].pos+1,r.mods1[i].mass+aa.getAAMass(r.peptide1[r.mods1[i].pos]));
  }
  for(i=0;i<r.peptide1.size();i++){
    if(aa.getFixedModMass(r.peptide1[i])>0) {
      sh.modInfo.addMod(i+1,aa.getAAMass(r.peptide1[i]));
    }
  }

  if(r.type>1){
    shB.hit_rank=1;
    shB.peptide=r.peptide2;
    shB.calc_neutral_pep_mass=r.massB;
    shB.massdiff=r.obsMass-r.massB;
    
    for(i=0;i<r.mods2.size();i++){
      shB.modInfo.addMod((int)r.mods2[i].pos+1,r.mods2[i].mass+aa.getAAMass(r.peptide2[r.mods2[i].pos]));
    }
    for(i=0;i<r.peptide2.size();i++){
      if(aa.getFixedModMass(r.peptide2[i])>0) {
        shB.modInfo.addMod(i+1,aa.getAAMass(r.peptide2[i]));
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

      shB.xlink_type="xl";
      sprintf(score,"%.4lf",r.scoreB);
      shB.addXLScore("score",score);
      sprintf(score,"%d",r.rankB);
      shB.addXLScore("rank",score);
      sprintf(score,"%d",r.link2);
      shB.addXLScore("link",score);
    }
    
  }

  sprintf(score,"%.4lf",r.score);
  sh.addScore("kojak_score",score);
  sprintf(score,"%.4lf",r.scoreDelta);
  sh.addScore("delta_score",score);
  sprintf(score,"%.4lf",r.ppm);
  sh.addScore("ppm_error",score);

  //Get proteins
  pep = db.getPeptide(r.pep1);
  sh.num_tot_proteins=(int)pep.map->size();
  for(j=0;j<pep.map->size();j++){
    protein="";
    for(i=0;i<db[pep.map->at(j).index].name.size();i++){
      if(params->truncate>0 && i==params->truncate) break;
      protein+=db[pep.map->at(j).index].name[i];
    }
    if(pep.map->at(j).start<1) n='-';
    else n=db[pep.map->at(j).index].sequence[pep.map->at(j).start-1];
    if(pep.map->at(j).stop+1==db[pep.map->at(j).index].sequence.size()) c='-';
    else c=db[pep.map->at(j).index].sequence[pep.map->at(j).stop+1];
    sh.addProtein(protein,c,n);
  }

  if(r.type>1){
    pep = db.getPeptide(r.pep2);
    shB.num_tot_proteins=(int)pep.map->size();
    for(j=0;j<pep.map->size();j++){
      protein="";
      for(i=0;i<db[pep.map->at(j).index].name.size();i++){
        if(params->truncate>0 && i==params->truncate) break;
        protein+=db[pep.map->at(j).index].name[i];
      }
      if(pep.map->at(j).start<1) n='-';
      else n=db[pep.map->at(j).index].sequence[pep.map->at(j).start-1];
      if(pep.map->at(j).stop+1==db[pep.map->at(j).index].sequence.size()) c='-';
      else c=db[pep.map->at(j).index].sequence[pep.map->at(j).stop+1];
      shB.addProtein(protein,c,n);
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

  kPeptide pep;
  kScoreCard sc;
  kScoreCard sc2;

  //Export Results:
  if(r.decoy) fprintf(f,"D-");
  else fprintf(f,"T-");
  fprintf(f,"%d-%.2f",r.scanNumber,r.rTime);
  if(count>1) fprintf(f,"-%d",count);
  if(r.decoy) fprintf(f,"\t-1");
  else fprintf(f,"\t1");
  if(params->percVersion>2.04) fprintf(f,"\t%d",r.scanNumber);
  fprintf(f,"\t%.4lf",r.score);
  fprintf(f,"\t%.4lf",r.scoreDelta);
  if(r.type==2 || r.type==3) fprintf(f,"\t%d\t%.4lf",r.rank,r.scorePepDif);
  //if(r.type==1) fprintf(f,"\t1\t0");
  //else if(r.type==2) fprintf(f,"\t0\t1");
  //else if(r.type==3) fprintf(f,"\t0\t0");
  //else fprintf(f,"\t0\t0");
  fprintf(f,"\t%d",r.charge);
  fprintf(f,"\t%.4lf",r.psmMass);
  fprintf(f,"\t%.4lf",r.ppm);
  if(r.pep2>=0 && r.peptide2.size()>0){
    if(r.peptide2.size()<r.peptide1.size()) {
      if(r.type==3) fprintf(f,"\t%d\t%d\t%d\t-.%s+%s.-",r.peptide2.size(),r.peptide1.size(),r.peptide1.size()+r.peptide2.size(),&r.modPeptide1[0],&r.modPeptide2[0]);
      else fprintf(f,"\t%d\t%d\t%d\t-.%s(%d)--%s(%d).-",r.peptide2.size(),r.peptide1.size(),r.peptide1.size()+r.peptide2.size(),&r.modPeptide1[0],r.link1,&r.modPeptide2[0],r.link2);
    } else {
      if(r.type==3) fprintf(f,"\t%d\t%d\t%d\t-.%s+%s.-",r.peptide1.size(),r.peptide2.size(),r.peptide1.size()+r.peptide2.size(),&r.modPeptide2[0],&r.modPeptide1[0]);
      else fprintf(f,"\t%d\t%d\t%d\t-.%s(%d)--%s(%d).-",r.peptide1.size(),r.peptide2.size(),r.peptide1.size()+r.peptide2.size(),&r.modPeptide2[0],r.link2,&r.modPeptide1[0],r.link1);
    }
  } else {
    fprintf(f,"\t%d\t-.%s",r.peptide1.size(),&r.modPeptide1[0]);
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
  int j,k,n,x,d;
  char fPath[1024];
  char fName[256];
  char peptide[256];
  char tmp[16];
  char specID[256];

  kPeptide pep;
  kPeptide pep2;
  kScoreCard tmpSC;
  kScoreCard tmpSC2;

  kResults res;

  PepXMLWriter p;
  pxwMSMSRunSummary rs;
  PXWSearchSummary ss;
  PXWSpectrumQuery sq;

  bool b;
  bool bBadFiles;
  bool bInter;
  bool bTarget1;
  bool bDupe;

  int preIndex;
  int scoreIndex;
  int iDupe;

  double dbl;
  double topScore;
  double ppm1;
  double ppm2;

  string tmpPep1;
  string tmpPep2;
  string outFile;

  FILE* fOut    = NULL;
  FILE* fIntra  = NULL;
  FILE* fInter  = NULL;
  FILE* fLoop   = NULL;
  FILE* fSingle = NULL;

  //Open all the required output files.
  bBadFiles=false;
  sprintf(fName,"%s.kojak.txt",params->outFile);
  fOut=fopen(fName,"wt");
  if(fOut==NULL) bBadFiles=true;
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
  }
  if(params->exportPepXML) {
    //cout << "export PepXML" << endl;
    getcwd(fPath,1024);
    //sprintf(fName,"%s%c%s.pep.xml",fPath,slashdir,params->outFile);
    //cout << fName << endl;
    outFile=params->outFile;
    i = outFile.find_last_of("/\\");
    outFile = outFile.substr(i+1);
    sprintf(fName, "%s.pep.xml", params->outFile);
    rs.base_name=fPath;
    rs.base_name+=slashdir;
    rs.base_name+=outFile;
    rs.raw_data=params->ext;
    rs.raw_data_type="raw";
    rs.search_engine="Kojak";
    ss.base_name=rs.base_name;
    ss.search_engine="Kojak";
    ss.search_engine_version=version;
    ss.search_database=params->dbFile;
    i=ss.search_database.find_last_of("/\\");
    ss.search_database=ss.search_database.substr(i+1);
    for(i=0;i<par.xmlParams.size();i++){
      ss.parameters->push_back(par.xmlParams[i]);
    }
    if(!p.createPepXML(fName,rs,&ss)) bBadFiles=true;
  }
  if(bBadFiles){
    cout << "Error exporting results. Please make sure drive is writable." << endl;
    if(fOut!=NULL)    fclose(fOut);
    if(fIntra!=NULL)  fclose(fIntra);
    if(fInter!=NULL)  fclose(fInter);
    if(fLoop!=NULL)   fclose(fLoop);
    if(fSingle!=NULL) fclose(fSingle);
    return false;
  }

  //Put the headers on all the files
  fprintf(fOut,"Kojak version %s\n",version);
  fprintf(fOut,"Scan Number\tRet Time\tObs Mass\tCharge\tPSM Mass\tPPM Error\tScore\tdScore\tPep. Diff.\tPeptide #1\tLink #1\tProtein #1\tPeptide #2\tLink #2\tProtein #2\tLinker Mass\n");
  if(params->exportPercolator){
    if(params->percVersion>2.04) {
      fprintf(fIntra,"SpecId\tLabel\tscannr\tScore\tdScore\t");
      fprintf(fInter,"SpecId\tLabel\tscannr\tScore\tdScore\t");
      fprintf(fLoop,"SpecId\tLabel\tscannr\tScore\tdScore\t");
      fprintf(fSingle,"SpecId\tLabel\tscannr\tScore\tdScore\t");
    } else {
      fprintf(fIntra,"SpecId\tLabel\tScore\tdScore\t");
      fprintf(fInter,"SpecId\tLabel\tScore\tdScore\t");
      fprintf(fLoop,"SpecId\tLabel\tScore\tdScore\t");
      fprintf(fSingle,"SpecId\tLabel\tScore\tdScore\t");
    }
    fprintf(fIntra,"NormRank\tPPScoreDiff\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum\tPeptide\tProteins\n");
    fprintf(fInter,"NormRank\tPPScoreDiff\tCharge\tMass\tPPM\tLenShort\tLenLong\tLenSum\tPeptide\tProteins\n");
    fprintf(fLoop,"Charge\tMass\tPPM\tLen\tPeptide\tProteins\n");
    fprintf(fSingle,"Charge\tMass\tPPM\tLen\tPeptide\tProteins\n");
  }
  
  //Output top score for each spectrum
  //Must iterate through all possible precursors for that spectrum
  for(i=0;i<spec.size();i++) {

    //This block is strictly for exporting the diagnostics
    //<--------- Start Diagnostics
    for(d=0;d<params->diag->size();d++){
      if(spec[i].getScanNumber()==params->diag->at(d)){
        char diagStr[256];
        sprintf(diagStr,"diagnostic_%d.txt",params->diag->at(d));

        FILE* f2=fopen(diagStr,"at");
        fprintf(f2,"\n");
        fprintf(f2,"# precursors: %d\n",spec[i].sizePrecursor());
        for(j=0;j<spec[i].sizePrecursor();j++) fprintf(f2,"Precursor: %.4lf\t%d\t%.4lf\t%d\n",spec[i].getPrecursor(j).monoMass,spec[i].getPrecursor(j).charge,spec[i].getPrecursor(j).corr,(int)spec[i].getPrecursor(j).label);
        fprintf(f2,"\n");

        for(j=0;j<20;j++){
          tmpSC = spec[i].getScoreCard(j);
          pep = db.getPeptide(tmpSC.pep1);
          db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
          for(n=0;n<strlen(peptide);n++){
            fprintf(f2,"%c",peptide[n]);
            for(x=0;x<tmpSC.mods1->size();x++){
              if(tmpSC.mods1->at(x).pos==n) fprintf(f2,"[%.2lf]",tmpSC.mods1->at(x).mass);
            }
            if(n==tmpSC.k1 || (tmpSC.pep2<0 && n==tmpSC.k2)) fprintf(f2,"[x]");
          }
          fprintf(f2,"(%d)",spec[i].getScoreCard(j).k1);
          if(tmpSC.k1>-1 && link[tmpSC.link].mono==0){
            if(tmpSC.pep2>-1){
              fprintf(f2,"\t");
              pep = db.getPeptide(tmpSC.pep2);
              db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
              for(n=0;n<strlen(peptide);n++){
                fprintf(f2,"%c",peptide[n]);
                for(x=0;x<tmpSC.mods2->size();x++){
                  if(tmpSC.mods2->at(x).pos==n) fprintf(f2,"[%.2lf]",tmpSC.mods2->at(x).mass);
                }
                if(n==tmpSC.k2) fprintf(f2,"[x]");
              }
              fprintf(f2,"(%d)",spec[i].getScoreCard(j).k2);
            } else {
              fprintf(f2,"\t(%d)",spec[i].getScoreCard(j).k2);
            }
          }
          if(tmpSC.link==-2){
            fprintf(f2,"+");
            pep = db.getPeptide(tmpSC.pep2);
            db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
            for(n=0;n<strlen(peptide);n++){
              fprintf(f2,"%c",peptide[n]);
              for(x=0;x<tmpSC.mods2->size();x++){
                if(tmpSC.mods2->at(x).pos==n) fprintf(f2,"[%.2lf]",tmpSC.mods2->at(x).mass);
              }
            }
            fprintf(f2,"(%d)",spec[i].getScoreCard(j).k2);
          }
          if(tmpSC.link>-1) fprintf(f2,"\t%.4lf\t%.4lf\t%.4lf\n",tmpSC.simpleScore,tmpSC.mass,link[tmpSC.link].mass);
          else fprintf(f2,"\t%.4lf\t%.4lf\t--\n",tmpSC.simpleScore,tmpSC.mass);
        }
        fclose(f2);
        break;
      }
    }
    //<--------- End Diagnostics

    scoreIndex=0;
    tmpSC=spec[i].getScoreCard(scoreIndex);
    res.scanNumber=spec[i].getScanNumber();
    res.rTime=spec[i].getRTime();

    //if there are no matches to the spectrum, return null result and continue
    if(tmpSC.simpleScore==0){
      fprintf(fOut,"%d\t%.4f\t0\t0\t0\t0\t0\t0\t0\t-\t-\t-\t-\t-\t-\t0\n",res.scanNumber,res.rTime);
      continue;
    }

    if(params->exportPepXML){
      sq.clear();
      sq.end_scan=res.scanNumber;
      sq.retention_time_sec=res.rTime;
      sq.start_scan=res.scanNumber;
    }

    //Export top scoring peptide, plus any ties that occur after it.
    topScore=tmpSC.simpleScore;
    int count=0;
    while(tmpSC.simpleScore==topScore){

      count++;

      //Get the best precursor ion for the PSM
      ppm1=(tmpSC.mass-spec[i].getPrecursor(0).monoMass)/spec[i].getPrecursor(0).monoMass*1e6;
      preIndex=0;
      for(j=1;j<spec[i].sizePrecursor();j++){
        ppm2=(tmpSC.mass-spec[i].getPrecursor(j).monoMass)/spec[i].getPrecursor(j).monoMass*1e6;
        if(fabs(ppm1)>fabs(ppm2)){
          preIndex=j;
          ppm1=ppm2;
        } else if(fabs(ppm1)==fabs(ppm2)){
          if(spec[i].getPrecursor(j).corr>spec[i].getPrecursor(preIndex).corr){
            preIndex=j;
            ppm1=ppm2;
          }
        }
      }
      res.obsMass = spec[i].getPrecursor(preIndex).monoMass;
      res.charge  = spec[i].getPrecursor(preIndex).charge;
      res.ppm     = ppm1;
      res.psmMass = tmpSC.mass;
      res.hk = spec[i].getPrecursor(preIndex).corr;

      if(params->exportPepXML){
        sq.assumed_charge=res.charge;
        sq.precursor_neutral_mass=res.obsMass;
        sprintf(specID,"%s.%d.%d.%d",&outFile[0],res.scanNumber,res.scanNumber,res.charge);
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
        ppm1=(tmpSC2.mass-spec[i].getPrecursor(preIndex).monoMass)/spec[i].getPrecursor(preIndex).monoMass*1e6;
        if(fabs(ppm1)>params->ppmPrecursor) continue;

        //if peptides and link sites are the same, go to the next one
        if(tmpSC2.pep1==tmpSC.pep1 && tmpSC2.pep2==tmpSC.pep2 && tmpSC2.k1==tmpSC.k1 && tmpSC2.k2==tmpSC.k2){
          continue;
        }
        break;
      }
      res.score       = tmpSC.simpleScore;
      res.scoreDelta  = tmpSC.simpleScore-tmpSC2.simpleScore;
      if(tmpSC.score1<tmpSC.score2) res.scorePepDif = tmpSC.score1;
      else res.scorePepDif = tmpSC.score2;
      res.rank        = tmpSC.rank1+tmpSC.rank2;
      res.rankA       = tmpSC.rank1+1;
      res.rankB       = tmpSC.rank2+1;
      res.scoreA      = tmpSC.score1;
      res.scoreB      = tmpSC.score2;
      res.massA       = tmpSC.mass1;
      res.massB       = tmpSC.mass2;

      //Get the peptide sequence(s)
      pep = db.getPeptide(tmpSC.pep1);
      db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
      res.peptide1 = peptide;
      res.mods1.clear();
      for(j=0;j<tmpSC.mods1->size();j++) res.mods1.push_back(tmpSC.mods1->at(j));
      res.peptide2 = "";
      if(tmpSC.pep2>=0){
        pep = db.getPeptide(tmpSC.pep2);
        db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
        res.peptide2 = peptide;
        res.mods2.clear();
        for(j=0;j<tmpSC.mods2->size();j++) res.mods2.push_back(tmpSC.mods2->at(j));
      }

      res.modPeptide1 = "";
      res.modPeptide2 = "";
      for(j=0;j<res.peptide1.size();j++) {
        res.modPeptide1 += res.peptide1[j];
        for(k=0;k<tmpSC.mods1->size();k++){
          if(j==(unsigned int)tmpSC.mods1->at(k).pos){
            sprintf(tmp,"[%.2lf]",tmpSC.mods1->at(k).mass);
            res.modPeptide1 += tmp;
          }
        }
      }
      for(j=0;j<res.peptide2.size();j++) {
        res.modPeptide2+=res.peptide2[j];
        for(k=0;k<tmpSC.mods2->size();k++){
          if(j==(unsigned int)tmpSC.mods2->at(k).pos){
            sprintf(tmp,"[%.2lf]",tmpSC.mods2->at(k).mass);
            res.modPeptide2 += tmp;
          }
        }
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

      if(res.type>0) {
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
              bDupe=true;
              break;
            }
          }
        }
      }

      //Determine if target or decoy
      bTarget1=false;
      pep = db.getPeptide(res.pep1);
      for(j=0;j<pep.map->size();j++) if(db[pep.map->at(j).index].name.find(params->decoy)==string::npos) bTarget1=true;
      if(bTarget1 && res.pep2>=0){
        if(!bDupe) bTarget1=false;
        pep = db.getPeptide(res.pep2);
        for(j=0;j<pep.map->size();j++) if(db[pep.map->at(j).index].name.find(params->decoy)==string::npos) bTarget1=true;
      }
      if(bTarget1) res.decoy=false;
      else res.decoy=true;

      //reorder peptides alphabetically
      if(res.type>1 && res.peptide1.compare(res.peptide2)>0){
        tmpPep1=res.peptide1;
        res.peptide1=res.peptide2;
        res.peptide2=tmpPep1;
        x=res.link1;
        res.link1=res.link2;
        res.link2=x;
        b=res.linkable1;
        res.linkable1=res.linkable2;
        res.linkable2=b;
        tmpPep1=res.modPeptide1;
        res.modPeptide1=res.modPeptide2;
        res.modPeptide2=tmpPep1;
        x=res.pep1;
        res.pep1=res.pep2;
        res.pep2=x;
        x=res.rankA;
        res.rankA=res.rankB;
        res.rankB=x;
        dbl=res.scoreA;
        res.scoreA=res.scoreB;
        res.scoreB=dbl;
        dbl=res.massA;
        res.massA=res.massB;
        res.massB=dbl;
        vector<kPepMod> vpm;
        for(j=0;j<res.mods1.size();j++) vpm.push_back(res.mods1[j]);
        res.mods1.clear();
        for(j=0;j<res.mods2.size();j++) res.mods1.push_back(res.mods2[j]);
        res.mods2.clear();
        for(j=0;j<vpm.size();j++) res.mods2.push_back(vpm[j]);
      }

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
      fprintf(fOut,"\t%.4lf",res.scorePepDif);
      fprintf(fOut,"\t%s",&res.modPeptide1[0]);
      fprintf(fOut,"\t%d",res.link1);

      fprintf(fOut,"\t");
      pep = db.getPeptide(res.pep1);
      for(j=0;j<pep.map->size();j++){
        fprintf(fOut,"%s",&db[pep.map->at(j).index].name[0]);
        if(res.link1>=0) fprintf(fOut,"(%d);",pep.map->at(j).start+res.link1); //put position from start of protein
        else fprintf(fOut,";");
      }
      if(bDupe){
        pep = db.getPeptide(res.pep2);
        for(j=0;j<pep.map->size();j++){
          fprintf(fOut,"%s;",&db[pep.map->at(j).index].name[0]);
          //if(res.link1>=0) fprintf(fOut,"(%d);",pep.map->at(j).start+res.link1); //only non-linked peptides
        }
      }

      if(res.modPeptide2.size()>1) {
        fprintf(fOut,"\t%s",&res.modPeptide2[0]);
        fprintf(fOut,"\t%d",res.link2);
        fprintf(fOut,"\t");
        pep = db.getPeptide(res.pep2);
        for(j=0;j<pep.map->size();j++){
          fprintf(fOut,"%s",&db[pep.map->at(j).index].name[0]);
          if(res.link1>=0) fprintf(fOut,"(%d);",pep.map->at(j).start+res.link2); //put position from start of protein
        }
        if(tmpSC.link>-1)fprintf(fOut,"\t%.4lf",link[tmpSC.link].mass);
        else fprintf(fOut,"\t-");
      } else if(res.link2>-1){
        fprintf(fOut,"\t-\t%d\t-",res.link2);
        fprintf(fOut,"\t%.4lf",link[tmpSC.link].mass);
      } else {
        fprintf(fOut,"\t-\t-1\t-\t0");
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
      
      if(params->exportPercolator) {
        switch(res.type){
          case 1:   outputPercolator(fLoop,db,res,count);   break;
          case 2:
            if(bInter)  outputPercolator(fInter,db,res,count);
            else        outputPercolator(fIntra,db,res,count);
            break;
          default:  outputPercolator(fSingle,db,res,count); break;
        }
      }

      if(params->exportPepXML){
        outputPepXML(sq,db,res);
      }

      //Get the next entry - it must also be exported if it has the same score
      if(bDupe) scoreIndex+=iDupe;
      else scoreIndex++;
      if(scoreIndex>=20) break;
      tmpSC=spec[i].getScoreCard(scoreIndex);
    }

    if(params->exportPepXML) {
      p.writeSpectrumQuery(sq);
    }

  }

  fclose(fOut);
  if(params->exportPercolator) {
    fclose(fIntra);
    fclose(fInter);
    fclose(fLoop);
    fclose(fSingle);
  }
  if(params->exportPepXML){
    p.closePepXML();
  }

  return true;

}

void KData::readLinkers(char* fn){
  FILE* f;
  kLinker k;

  cout << "Reading linker file...";

  f=fopen(fn,"rt");
  while(!feof(f)){
    fscanf(f,"%lf\t%d\n",&k.mass,&k.mono);
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
  int finalPeaks=0;
  int iPercent=0;
  int iTmp;

  int i;
  int j;

  spec.clear();
  msr.setFilter(MS2);

  //Set progress meter
  printf("%2d%%", iPercent);
  fflush(stdout);

  if(!msr.readFile(params->msFile,s)) return false;
  while(s.getScanNumber()>0){

    totalScans++;
    if(s.size()<1) {
      msr.readFile(NULL,s);
      continue;
    }

    //This is for the methods used in old grad school data
    /*
    if(s.getRTime()<15){
      msr.readFile(NULL,s);
      continue;
    }
    */

    pls.clear();
    pls.setRTime(s.getRTime());
    pls.setScanNumber(s.getScanNumber());
    max=0;

    //Check whether scans are centroided or not (info supplied by user in params)
    if(!params->ms2Centroid) {

      //If not centroided, do so now.
      centroid(s,c,params->ms2Resolution,params->instrument);

      totalPeaks+=c.size();
      
      //Collapse the isotope peaks
      if(params->specProcess==1 && c.size()>1) {
        collapseSpectrum(c);
        collapsedPeaks+=c.size();
      }

      //If user limits number of peaks to analyze, sort by intensity and take top N
      if(params->maxPeaks>0){
        if(c.size()>1) c.sortIntensityRev();
        if(c.size()<params->maxPeaks) j=c.size();
        else j=params->maxPeaks;
      } else {
        j=c.size();
      }
      for(i=0;i<j;i++){
        sp.mass=c[i].mz;
        sp.intensity=c[i].intensity;
        pls.addPoint(sp);
        if(sp.intensity>max) max=sp.intensity;
      }
      pls.setMaxIntensity(max);

      //Sort again by MZ, if needed
      if(pls.size()>1 && params->maxPeaks>0) pls.sortMZ();

      finalPeaks+=pls.size();

    } else {

      //Collapse the isotope peaks
      if(params->specProcess==1 && c.size()>1) collapseSpectrum(c);

      //If user limits number of peaks to analyze, sort by intensity and take top N
      if(params->maxPeaks>0){
        if(s.size()>1) s.sortIntensityRev();
        if(s.size()<params->maxPeaks) j=s.size();
        else j=params->maxPeaks;
      } else {
        j=s.size();
      }
      for(i=0;i<j;i++){
        sp.mass=s[i].mz;
        sp.intensity=s[i].intensity;
        pls.addPoint(sp);
        if(sp.intensity>max) max=sp.intensity;
      }
      pls.setMaxIntensity(max);
      
      //Sort again by MZ, if needed
      if(pls.size()>1 && params->maxPeaks>0) pls.sortMZ();
    }

    //Get any additional information user requested
    pls.setCharge(s.getCharge());
    pls.setMZ(s.getMZ());
    if(params->preferPrecursor>0){
      if(s.getMonoMZ()>0 && s.getCharge()>0){
        pre.monoMass=s.getMonoMZ()*s.getCharge()-s.getCharge()*1.007276466;
        pre.charge=s.getCharge();
        pre.corr=0;
        pls.addPrecursor(pre);
        pls.setInstrumentPrecursor(true);
      } else if(s.sizeZ()>0){
        for(j=0; j<s.sizeZ(); j++){
          pre.monoMass=s.atZ(j).mh-1.007276466;
          pre.charge=s.atZ(j).z;
          pre.corr=-5;
          pls.setCharge(pre.charge);
          pls.addPrecursor(pre);
          pls.setInstrumentPrecursor(true);
        }
      }
    }

    //Add spectrum (if it has enough data points) to data object and read next file
    if(pls.size()>12) spec.push_back(pls);

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

    //Update progress meter
    iTmp = msr.getPercent();
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }

    msr.readFile(NULL,s);
  }

  //Finalize progress meter
  if(iPercent<100) printf("\b\b\b100%%");
  cout << endl;

  cout << "  " << spec.size() << " total spectra have enough data points for searching." << endl;
  //cout << totalScans << " total scans were loaded." <<  endl;
  //cout << totalPeaks << " total peaks in original data." << endl;
  //cout << collapsedPeaks << " peaks after collapsing." << endl;
  //cout << finalPeaks << " peaks after top N." << endl;
	return true;
}

void KData::setLinker(kLinker x){
  if(x.mono==0) link.push_back(x);
}

void KData::setVersion(char* v){
  strcpy(version,v);
}

int KData::size(){
  return (int)spec.size();
}

int KData::sizeLink(){
  return (int)link.size();
}

void KData::xCorr(bool b){
  if(b) cout << "  Using XCorr scores." << endl;
  else  cout << "  Using Kojak modified XCorr scores." << endl;

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



