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
  params=NULL;
}

KData::KData(kParams* p){
  params=p;
}

KData::~KData(){
  params=NULL;
}


/*============================
  Operators
============================*/
KSpectrum& KData::operator [](const int &i){
  return spec[i];
}


/*============================
  Functions
============================*/
bool KData::getBoundaries(double mass1, double mass2, vector<int>& index){
  int sz=massList.size();
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

}

//Get the list of spectrum array indexes to search based on desired mass
bool KData::getBoundaries2(double mass, double prec, vector<int>& index){
  int sz=massList.size();
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

//This function tries to assign best possible 18O2 and 18O4 precursor ion mass values
//for all MS2 spectra
bool KData::mapPrecursors(){
  
  int iPercent=0;
  int iTmp;
  
  unsigned int i;
  int j;

  KPrecursor pre(params);
  kMass      m;

  int peakCounts=0;
  int specCounts=0;

  int prePre=0;
  int foundPre=0;
  int noPre=0;

  //Open the data file in the precursor mapping object
  //if(!pre.setFile(&p)) return false;

  //Print progress
  printf("Mapping precursors: %2d%%",iPercent);
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

    //if precursors already existed (user request at spectrum load time),
    //don't compute new precursors
    if(spec[i].sizePrecursor()>0){
      prePre++;
      specCounts++;
      peakCounts+=spec[i].size();
      continue;
    }

    //Find precursor using object function. Take results and copy them to spectra
    if(pre.getSpecRange(spec[i])){
    //if(pre.getSpecRange(spec[i].getScanNumber(),monoO2,chargeO2,corrO2,monoO4,chargeO4,corrO4)){
    //  if(corrO2>0){
    //    spec[i].setCharge(chargeO2,true);
    //    spec[i].setMass(monoO2,true);
    //    spec[i].setCorr(corrO2,true);
    //  }
    //  if(corrO4>0){
    //    spec[i].setCharge(chargeO4,false);
    //    spec[i].setMass(monoO4,false);
    //    spec[i].setCorr(corrO4,false);
    //  }    
    //  if(corrO2>0 || corrO4>0){
      foundPre++;
        specCounts++;
        peakCounts+=spec[i].size();
    //  }
    } else { //If no precursors found, estimate using mercury and selected mass
      if(params->enrichment==0 && pre.estimatePrecursor(spec[i])){
        noPre++;
        specCounts++;
        peakCounts+=spec[i].size();
      }
    } 

  }
 

  //Finalize the progress
  printf("\b\b\b100%%");
  cout << endl;

  cout << specCounts << " spectra with " << peakCounts << " peaks will be analyzed." << endl;
  //cout << prePre << " spectra had precomputed precursor mass." << endl;
  //cout << foundPre << " spectra had PepLynx processed precursor mass." << endl;
  //cout << noPre << " spectra used low-res isolation mass to predict precursor mass." << endl;

  /*
  FILE* f=fopen("data.ms2","wt");
  fprintf(f,"H\tExtractor\tPepLynx\n");
  fprintf(f,"H\tExtractorVersion\t2.52a\n");
  for(i=0;i<spec.size();i++){
    fprintf(f,"S\t%d\t%d\n",spec[i].getScanNumber(),spec[i].getScanNumber());
    for(j=0;j<spec[i].sizePrecursor();j++){
      fprintf(f,"Z\t%d\t%.4lf\n",spec[i].getPrecursor(j).charge,spec[i].getPrecursor(j).monoMass+1.007276466);
    }
    for(j=0;j<spec[i].size();j++){
      fprintf(f,"%.4lf\t%.1f\n",spec[i][j].mass,spec[i][j].intensity);
    }
  }
  fclose(f);
  exit(1);
  */

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

  return true;
}

bool KData::outputPercolator(char* out, char* tag, KDatabase& db, int pLen){

  bool bTarget1;
  bool bTarget2;

  double dScore;
  double ppm;
  double ppm2;
  double mass1;

  int cross;
  int dimer;
  int label;
  unsigned int len1;
  unsigned int len2;
  unsigned int len3;
  int loop;
  int mono;
  int pos;
  int preIndex;
  //int preIndex2;

  unsigned int i;
  unsigned int j;
  unsigned int k;

  char ID[16];
  char tmp[16];

  string peptide;
  string protein;
  string sequence;
  string tStr;

  kPeptide pep;
  kScoreCard sc;
  kScoreCard sc2;

  FILE* f=fopen(out,"wt");
  if(f==NULL) return false;

  if(params->percVersion>2.04) fprintf(f,"SpecId\tLabel\tScanNum\tdScore\t");
  else fprintf(f,"SpecId\tLabel\tScore\tdScore\t");
  if(params->relaxedAnalysis) fprintf(f,"NormRank\t");
  fprintf(f,"Dimer\tMono\tLoop\tCross\tCharge\tMass\tPPM\tLenShort\tLenLong\tPeptide\tProteins\n");
  for(i=0;i<spec.size();i++){

    bTarget1=false;
    bTarget2=false;
    cross=0;
    dimer=0;
    mono=0;
    loop=0;
    len1=0;
    len2=0;

    sc=spec[i].getScoreCard(0);
    if(sc.simpleScore==0) continue;
  
    preIndex=0;
    ppm=(sc.mass-spec[i].getPrecursor(0).monoMass)/spec[i].getPrecursor(0).monoMass*1e6;
    for(j=1;j<(unsigned int)spec[i].sizePrecursor();j++){
      if(fabs(ppm)>fabs((sc.mass-spec[i].getPrecursor(j).monoMass)/spec[i].getPrecursor(j).monoMass*1e6)) {
        ppm = (sc.mass-spec[i].getPrecursor(j).monoMass)/spec[i].getPrecursor(j).monoMass*1e6;
        preIndex=j;
      }
    }

    k=1;
    while(k<19){
      sc2=spec[i].getScoreCard(k++);
      ppm2=(sc2.mass-spec[i].getPrecursor(preIndex).monoMass)/spec[i].getPrecursor(preIndex).monoMass*1e6;
      if(fabs(ppm2)>params->ppmPrecursor) continue;

      //if peptides and link sites are the same, go to the next one
      if(sc2.pep1==sc.pep1 && sc2.pep2==sc.pep2 && sc2.k1==sc.k1 && sc2.k2==sc.k2 && sc2.linkable1==sc.linkable1 && sc2.linkable2==sc.linkable2){
        continue;
      }
      break;
    }
    dScore=sc.simpleScore-sc2.simpleScore;

    peptide="";
    pep = db.getPeptide(sc.pep1,sc.linkable1);
    db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,tStr);
    len1=tStr.length();
    mass1=pep.mass;
    for(j=0;j<len1;j++){
      peptide+=tStr[j];
      for(k=0;k<sc.mods->size();k++){
        if(j==(unsigned int)sc.mods->at(k).pos){
          sprintf(tmp,"[%.2lf]",sc.mods->at(k).mass);
          peptide+=tmp;
        }
      }
    }
    if(sc.k1>-1){
      sprintf(tmp,"(%d)",sc.k1-pep.map->at(0).start+1);
      peptide+=tmp;
    }

    protein="";
    for(j=0;j<pep.map->size();j++){
      if(db[pep.map->at(j).index].name.find(tag)==string::npos) bTarget1=true;
      if(j>0) protein+="\t";
      if(pLen>0) protein+=db[pep.map->at(j).index].name.substr(0,pLen);
      else protein+=db[pep.map->at(j).index].name;
      if(sc.k1>-1){
        pos=sc.k1-pep.map->at(0).start+1;
        sprintf(tmp,"(%d)",pep.map->at(j).start+pos);
        protein+=tmp;
      }
    }

    if(sc.link==-2){ //two non-covalently linked peptides
      dimer=1;
      pep = db.getPeptide(sc.pep2,sc.linkable2);
      db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,tStr);
      sequence="";
      for(j=0;j<tStr.length();j++){
        sequence+=tStr[j];
        for(k=0;k<sc.mods->size();k++){
          if(j==(unsigned int)sc.mods->at(k).pos){
            sprintf(tmp,"[%.2lf]",sc.mods->at(k).mass);
            sequence+=tmp;
          }
        }
      }
      if((int)tStr.length()>len1 || ((int)tStr.length()==len1 && pep.mass>mass1)){
        peptide=sequence+"+"+peptide;
        len2=len1;
        len1=tStr.length();
      } else {
        peptide+="+"+sequence;
        len2=tStr.length();
      }

      bTarget2=false;
      for(j=0;j<pep.map->size();j++){
        if(db[pep.map->at(j).index].name.find(tag)==string::npos) bTarget2=true;
        protein+="\t";
        if(pLen>0) protein+=db[pep.map->at(j).index].name.substr(0,pLen);
        else protein+=db[pep.map->at(j).index].name;
      }

    } else if(sc.link==-1) { //single peptides
      bTarget2=true;
    } else if(link[sc.link].mono>0) { //mono-links
      bTarget2=true;
      mono=1;
      peptide+="-MONO";
    } else { //Cross- and loop-links
      if(sc.pep2==-1) { //Loop
        bTarget2=true;
        loop=1;
        sprintf(tmp,",%d)-LOOP",sc.k2-pep.map->at(0).start+1);
        peptide.erase(peptide.length()-1,1);
        peptide+=tmp;
      } else {
        cross=1;
        pep = db.getPeptide(sc.pep2,sc.linkable2);
        db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,sequence);
        len3=sequence.length();
        sprintf(tmp,"(%d)",sc.k2-pep.map->at(0).start+1);
        sequence+=tmp;
        if(len3>len1 || (len3==len1 && pep.mass>mass1)){
          peptide=sequence+"--"+peptide;
          len2=len1;
          len1=len3;
        } else {
          peptide+="--"+sequence;
          len2=len3;
        }
      }

      for(j=0;j<pep.map->size();j++){
        if(db[pep.map->at(j).index].name.find(tag)==string::npos) bTarget2=true;
        protein+="\t";
        if(pLen>0) protein+=db[pep.map->at(j).index].name.substr(0,pLen);
        else protein+=db[pep.map->at(j).index].name;
        if(sc.k2>-1){
          pos=sc.k2-pep.map->at(0).start+1;
          sprintf(tmp,"(%d)",pep.map->at(j).start+pos);
          protein+=tmp;
        }
      }
    }


    if(bTarget1 && bTarget2) {
      sprintf(ID,"T-%d",spec[i].getScanNumber());
      label=1;
    } else {
      sprintf(ID,"D-%d",spec[i].getScanNumber());
      label=-1;
    }

    for(j=0;j<protein.size();j++){
      if(protein[j]==' ') protein[j]='_';
    }

    if(params->percVersion>2.04) fprintf(f,"%s\t%d\t%d\t%.4lf\t%.4lf\t",ID,label,spec[i].getScanNumber(),sc.simpleScore,dScore);
    else fprintf(f,"%s\t%d\t%.4lf\t%.4lf\t",ID,label,sc.simpleScore,dScore);
    if(params->relaxedAnalysis) fprintf(f,"%d\t",sc.rank);
    fprintf(f,"%d\t%d\t%d\t%d\t%d\t%.4lf\t%.2lf\t%d\t%d\t-.%s.-\t%s\n",dimer,mono,loop,cross,
      spec[i].getPrecursor(preIndex).charge,sc.mass,ppm,len1,len2,&peptide[0],&protein[0]);
                                                           
  }

  fclose(f);
  return true;
}

bool KData::outputResults(char* out, KDatabase& db, bool bEnrich){

  unsigned int i,j,n;
  char peptide[256];
  char tmp[16];

  kPeptide pep;
  kScoreCard best[5];
  kScoreCard tmpSC;
  kScoreCard tmpSC2;

  int charge;
  int link1;
  int link2;
  int pos1;
  int pos2;
  int preIndex;
  //int preIndex2;
  double ppm1;
  double ppm2;
  string pep1;
  string pep2;
  string pepB;
  string prot1;
  string prot2;

  FILE* f=fopen(out,"wt");
  if(f==NULL) return false;

  fprintf(f,"PepLynx version %s\n",version);

  fprintf(f,"Scan Number\tObs Mass\tCorr\t");
  if(bEnrich) fprintf(f,"Label\t");
  fprintf(f,"Peptide Mass\tCharge\tPPM Error\tProtein #1\tProtein #2\tPeptide #1\tPeptide #2\tLink #1\tLink #2\tScore\tdScore\t");
  if(params->relaxedAnalysis>0) fprintf(f,"Norm Rank\t");
  fprintf(f,"Link Mass\tMod Mass\n");

  //Output top n (currently 1) scores for each spectrum
  //Must iterate through all possible precursors for that spectrum
  for(i=0;i<spec.size();i++) {

    
    if(spec[i].getScanNumber()==params->diagnostic){
      FILE* f2=fopen("diagnostic.txt","at");
      fprintf(f2,"\n");
      for(j=0;j<(unsigned int)spec[i].sizePrecursor();j++) fprintf(f2,"Precursor: %.4lf\t%d\t%.4lf\n",spec[i].getPrecursor(j).monoMass,spec[i].getPrecursor(j).charge,spec[i].getPrecursor(j).corr);
      fprintf(f2,"\n");
      for(j=0;j<20;j++){
        pep = db.getPeptide(spec[i].getScoreCard(j).pep1,spec[i].getScoreCard(j).linkable1);
        db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
        fprintf(f2,"%s(%d)",peptide,spec[i].getScoreCard(j).k1);
        if(spec[i].getScoreCard(j).k1>-1 && link[spec[i].getScoreCard(j).link].mono==0){
          pep = db.getPeptide(spec[i].getScoreCard(j).pep2,spec[i].getScoreCard(j).linkable2);
          db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
          fprintf(f2,"\t%s(%d)",peptide,spec[i].getScoreCard(j).k2);
        }
        fprintf(f2,"\t%.4lf\n",spec[i].getScoreCard(j).simpleScore);
      }
      fclose(f2);
    }

    fprintf(f,"%d\t",spec[i].getScanNumber());
    tmpSC=spec[i].getScoreCard(0);

    if(tmpSC.simpleScore>0){  
      ppm1=(tmpSC.mass-spec[i].getPrecursor(0).monoMass)/spec[i].getPrecursor(0).monoMass*1e6;
      preIndex=0;
      for(j=1;j<(unsigned int)spec[i].sizePrecursor();j++){
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
      charge=spec[i].getPrecursor(preIndex).charge;
      fprintf(f,"%.4lf\t%.4lf\t",spec[i].getPrecursor(preIndex).monoMass,spec[i].getPrecursor(preIndex).corr);
      if(bEnrich){
        if(spec[i].getPrecursor(preIndex).label==0) fprintf(f,"0\t");
        else if(spec[i].getPrecursor(preIndex).label==1) fprintf(f,"O2\t");
        else fprintf(f,"O4\t");
      }
      fprintf(f,"%.4lf\t%d\t%.4lf\t",tmpSC.mass,spec[i].getPrecursor(preIndex).charge,ppm1);
      pep = db.getPeptide(tmpSC.pep1,tmpSC.linkable1);
      db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
      
      pep1=peptide;
      prot1="";

      for(j=0;j<pep.map->size();j++){
        prot1+=db[pep.map->at(j).index].name;
        pos1=tmpSC.k1-pep.map->at(0).start+1;           //get position from start of peptide
        sprintf(tmp,"(%d);",pep.map->at(j).start+pos1); //put position from start of protein
        prot1+=tmp;
      }
      
      if(tmpSC.link==-2){
        link1=-1;
        link2=-1;
        pep = db.getPeptide(tmpSC.pep2,tmpSC.linkable2);
        db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
        pep2=peptide;
        prot2="";
        for(j=0;j<pep.map->size();j++){
          prot2+=db[pep.map->at(j).index].name;
          pos2=tmpSC.k2-pep.map->at(0).start+1;           //get position from start of peptide
          sprintf(tmp,"(%d);",pep.map->at(j).start+pos2); //put position from start of protein
          prot2+=tmp;
        }
      } else if(tmpSC.link==-1){
        link1=tmpSC.k1;
        pep2="-";
        prot2="-";
        link2=-1;
      } else if(link[tmpSC.link].mono>0){
        link1=tmpSC.k1-pep.map->at(0).start;
        pep2="Mono";
        prot2="-";
        link2=-1;
      } else {
        link1=tmpSC.k1-pep.map->at(0).start;
        if(tmpSC.pep2>-1){
          pep = db.getPeptide(tmpSC.pep2,tmpSC.linkable2);
          db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
          pep2=peptide;
          prot2="";
          for(j=0;j<pep.map->size();j++){
            prot2+=db[pep.map->at(j).index].name;
            pos2=tmpSC.k2-pep.map->at(0).start+1;           //get position from start of peptide
            sprintf(tmp,"(%d);",pep.map->at(j).start+pos2); //put position from start of protein
            prot2+=tmp;
          }
        } else {
          pep2="Loop";
          prot2="-";
        }
        link2=tmpSC.k2-pep.map->at(0).start;
      }
      fprintf(f,"%s\t%s\t%s\t%s\t%d\t%d\t",&prot1[0],&prot2[0],&pep1[0],&pep2[0],link1,link2);
      fprintf(f,"%.4lf\t",tmpSC.simpleScore);
      
      //grab the next highest score that matches to the same precursor ion for the delta score
      //if no other match has the same precursor, just take the lowest score in the list
      n=1;
      while(n<19){
        tmpSC2=spec[i].getScoreCard(n++);
        if(tmpSC2.simpleScore==0) break;
        ppm1=(tmpSC2.mass-spec[i].getPrecursor(preIndex).monoMass)/spec[i].getPrecursor(preIndex).monoMass*1e6;
        if(fabs(ppm1)>params->ppmPrecursor) continue;

        //if peptides and link sites are the same, go to the next one
        if(tmpSC2.pep1==tmpSC.pep1 && tmpSC2.pep2==tmpSC.pep2 && tmpSC2.k1==tmpSC.k1 && tmpSC2.k2==tmpSC.k2 && tmpSC2.linkable1==tmpSC.linkable1 && tmpSC2.linkable2==tmpSC.linkable2){
          continue;
        }
        break;
      }
      fprintf(f,"%.4lf\t",tmpSC.simpleScore-tmpSC2.simpleScore);
      
      if(params->relaxedAnalysis) fprintf(f,"%d\t",tmpSC.rank);
      if(tmpSC.link<0) fprintf(f,"0\t");
      else fprintf(f,"%.2lf\t",link[tmpSC.link].mass);

      fprintf(f,"0\n"); //tmpSC.modMass);

    } else {
      
      fprintf(f,"0\t0\t");
      if(bEnrich) fprintf(f,"-\t");
      fprintf(f,"0\t0\t-100000\t-\t-\t-\t-\t-\t-\t0\t0\t");
      if(params->relaxedAnalysis) fprintf(f,"0\t");
      fprintf(f,"0\t0\n");

    }
  }

  fclose(f);
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

//Reads in Hardklor files. May change this to some other format.
/*
bool KData::readSpectra(char* fname){
  FILE*       hkr;
  KSpectrum  pls();
  specPoint   sp;

	bool        firstScan;

	char        line[256];
  char        tag;
	char*       tok;

	int count=0;

  //Read in the Hardklor results
  firstScan=true;
	hkr = fopen(fname,"rt");
	if(hkr==NULL) {
		cout << "Problem reading " << fname << endl;
		return false;
	}

	while(!feof(hkr)){
    if(fgets(line,256,hkr)==NULL) continue;
		tag=line[0];

    if(tag=='S') {
			if(firstScan) firstScan=false;
      else {
        if(count>12) spec.push_back(pls);	
        count=0;
        pls.clear();
      }

			tok=strtok(line,"\t\n");
			tok=strtok(NULL,"\t\n");
      pls.setScanNumber(atoi(tok));

			tok=strtok(NULL,"\t\n");
			pls.setRTime((float)atof(tok));

		} else {
			count++;
      tok=strtok(line,"\t\n");
      tok=strtok(NULL,"\t\n");
      sp.mass=atof(tok);

      tok=strtok(NULL,"\t\n");
      tok=strtok(NULL,"\t\n");
      sp.intensity=(float)atof(tok);

      pls.addPoint(sp);
		}
	}
  if(count>12) spec.push_back(pls);	//Only accept spectra with this minimum number of data points
	fclose(hkr);

  cout << spec.size() << " total spectra have enough data points for searching." << endl;
	return true;
}
*/

//Reads in raw files.
bool KData::readSpectra2(){
  MSReader    msr;
  Spectrum    s;
  Spectrum    c;
  KSpectrum  pls(params->topCount);
  specPoint   sp;
  float       max;
  kPrecursor pre;

  int totalScans=0;
  int totalPeaks=0;
  int collapsedPeaks=0;
  int finalPeaks=0;

  int i;
  int j;

  cout << "Reading spectra ..." << endl;

  msr.setFilter(MS2);
  msr.readFile(params->msFile,s);
  while(s.getScanNumber()>0){

    totalScans++;

    //This is for the methods used in old grad school data
    /*
    if(s.getRTime()<15){
      msr.readFile(NULL,s);
      continue;
    }
    */

    //printf("Scan Number:  %d\n",s.getScanNumber());

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

    /*
    if(s.getScanNumber()==2661){
      for(int xx=0;xx<pls.size();xx++){
        printf("%.4lf\t0\n%.4lf\t%.0f\n%.4lf\t0\n",pls[xx].mass-0.001,pls[xx].mass,pls[xx].intensity,pls[xx].mass+0.001);
      }
    }
    */
    //Get any additional information user requested
    pls.setCharge(s.getCharge());
    pls.setMZ(s.getMZ());
    if(params->preferPrecursor){
      if(s.getMonoMZ()>0 && s.getCharge()>0){
        pre.monoMass=s.getMonoMZ()*s.getCharge()-s.getCharge()*1.007276466;
        pre.charge=s.getCharge();
        pls.addPrecursor(pre);
      }
    }

    //Add spectrum (if it has enough data points) to data object and read next file
    if(pls.size()>12) spec.push_back(pls);

    if(pls.getScanNumber()==params->diagnostic){
      FILE* f=fopen("diagnostic_spectrum.txt","wt");
      fprintf(f,"Scan: %d\t%d\n",pls.getScanNumber(),pls.size());
      for(int k=0;k<pls.size();k++) fprintf(f,"%.6lf\t%.0f\n",pls[k].mass,pls[k].intensity);
      fclose(f);
    }

    msr.readFile(NULL,s);
  }

  cout << spec.size() << " total spectra have enough data points for searching." << endl;
  //cout << totalScans << " total scans were loaded." <<  endl;
  //cout << totalPeaks << " total peaks in original data." << endl;
  //cout << collapsedPeaks << " peaks after collapsing." << endl;
  //cout << finalPeaks << " peaks after top N." << endl;
	return true;
}

void KData::setLinker(double mass, int type){
  kLinker p;
  p.mass=mass;
  p.mono=type;
  link.push_back(p);
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

void KData::xCorr(){
  cout << "Using XCorr scores." << endl;
  for(unsigned int i=0;i<spec.size();i++) spec[i].xCorrScore();
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

    //printf("Checking %.6lf %.0f\n",s[maxIndex].mz,s[maxIndex].intensity);

    //check right
    j=maxIndex+1;
    while(j<s.size() && (s[j].mz-s[maxIndex].mz)<1.1){
      if(s[j].intensity<1) {
        j++;
        continue;
      }
      charge=getCharge(s,maxIndex,j);

        //printf("\t%.6lf\t%.6lf\t%d\n",s[j].mz,s[j].mz-s[maxIndex].mz,charge);

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
        //printf("\tright\t%.6lf\t%.6lf\t%d\n",s[n].mz,s[n].mz-s[k].mz,z);
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
      //printf("\t%.6lf\t%.6lf\t%d\n",s[j].mz,s[maxIndex].mz-s[j].mz,charge);
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
        //printf("Adding %.6lf\t%d+\t%d\n",s[dist[0]].mz,charge,dist.size());
        s[dist[1]].intensity=0;
       // s2.add(s[dist[1]].mz,max);
      } else {
        s2.add(s[dist[0]]);
       // s2.add(s[dist[1]]);
      }
      s[dist[0]].intensity=0;
      //s[dist[1]].intensity=0;
    } else {
      //printf("collapsing\n");
      cutoff=max/20;
      max=0;
      j=dist[0];
      k=dist[1];
      for(i=0;i<(int)dist.size();i++) {
        if(dist[i]<j && s[dist[i]].intensity>cutoff){
          k=j;
          j=dist[i];
          //printf("new index: %.6lf\n",s[j].mz);
        }
        if(s[dist[i]].intensity>cutoff){
          //printf("removing: %.6lf\n",s[dist[i]].mz);
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

/*
int KData::compareSpecMass(const void *p1, const void *p2){
  KSpectrum d1 = *(KSpectrum *)p1;
  KSpectrum d2 = *(KSpectrum *)p2;
  if(d1.getMass()<d2.getMass()) {
		return -1;
	} else if(d1.getMass()>d2.getMass()) {
  	return 1;
  } else {
	  return 0;
  }
}
*/

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

