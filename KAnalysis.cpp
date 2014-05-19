#include "KAnalysis.h"

/*============================
  Constructors & Destructors
============================*/
KAnalysis::KAnalysis(kParams& p){
  unsigned int i;
  params=p;
  for(i=0;i<params.fMods->size();i++) ions.addFixedMod((char)params.fMods->at(i).index,params.fMods->at(i).mass);
  for(i=0;i<params.mods->size();i++) ions.addMod((char)params.mods->at(i).index,params.mods->at(i).mass);
  ions.setMaxModCount(params.maxMods);
  //xCorrCount=0;
}

KAnalysis::~KAnalysis(){
}

/*============================
  Public Functions
============================*/
bool KAnalysis::analyzePeptides(KDatabase& db, KData& d, bool crossLink){
  unsigned int i,j,k,k2;
  int n;
  int iPercent;
  int iTmp;
  vector<kPeptide>* p;
  vector<int> index;
  vector<kPepMod> mods;
  double totalMass;
  double highLinkMass=0;
  double lowLinkMass=0;
  double maxMass = d.getMaxMass()+0.25;
  double minMass = d.getMinMass()-0.25;
  kScoreCard sc;
  double mm;
  bool bt;

  //if(!crossLink) return true;

  //Set progress meter
  iPercent=0;
  printf("Progress: %2d%%",iPercent);
  fflush(stdout);

  //set cross-linker mass boundaries.
  for(n=0;n<d.sizeLink();n++){
    if(d.getLink(n).mono==0){
      if(lowLinkMass==0) lowLinkMass=d.getLink(n).mass;
      if(highLinkMass==0) highLinkMass=d.getLink(n).mass;
      if(d.getLink(n).mass<lowLinkMass) lowLinkMass=d.getLink(n).mass;
      if(d.getLink(n).mass>highLinkMass) highLinkMass=d.getLink(n).mass;
    }
  }

  //Set which list of peptides to search (with and without internal lysine)
  p=&db.getPeptideList(crossLink);

  //Iterate entire peptide list
  for(i=0;i<p->size();i++){

    //Update progress meter
    iTmp=(int)(i*100.0/p->size());
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }

    //Peptides are sorted by mass. If greater than max mass, stop checking peptides
    if(p->at(i).mass>maxMass) break;
      
    //Get a set of spectra to analyze. Spectra were sorted by precursor mass
    if(p->at(i).mass<minMass) bt=false;
    else bt=d.getBoundaries2(p->at(i).mass,params.ppmPrecursor,index);

    //Set the peptide, calc the ions, and score it against the spectra
    ions.setPeptide(true,&db[p->at(i).map->at(0).index].sequence[p->at(i).map->at(0).start],p->at(i).map->at(0).stop-p->at(i).map->at(0).start+1,p->at(i).mass);
    if(bt){
      ions.buildIons(0,0);
      scoreSpectra(index,d,p->at(i).mass,0,crossLink,i,-1,-1,-1,-1);
    }
    
    //Intercept here to check for links using singlet analysis
    if(params.relaxedAnalysis>0){

      //Must search singlets for all peptides if searching for dimers.
      if(params.dimers>0) {
        if(!bt) ions.buildIons(0,0);
        analyzeSingletsNoLysine(db,d,p->at(i),i,crossLink);
      }

      //Crosslinked peptides must also search singlets with reciprocol mass on each lysine
      if(crossLink) analyzeSinglets(db,d,p->at(i),i,lowLinkMass,highLinkMass);

    }

    //Iterate across all differential mods
    if(params.maxMods>0){
      mm=ions.buildModIons(0,0);
      while(mm>0){
        if(mm>maxMass || mm<minMass){
          mm=ions.buildModIons(0,0);
          continue;
        }
        bt=d.getBoundaries2(mm,params.ppmPrecursor,index);
        if(bt) {
          ions.getPeptideMods(mods);
          scoreSpectra(index,d,mm,&mods,crossLink,i,-1,-1,-1,-1);
        }
        mm=ions.buildModIons(0,0);
      }
    }

    //if peptides cannot be crosslinked, move on to next peptide now
    if(!crossLink) continue;
    
    //check monolinks & looplinks by iterating through each cross-linker mass
    for(n=0;n<d.sizeLink();n++){

      totalMass = p->at(i).mass+d.getLink(n).mass;

      //skip looplinks on peptides with fewer than 2 linkable lysines
      if(d.getLink(n).mono==0 && p->at(i).vK->size()<2) continue;
      if(totalMass>maxMass) continue;

      //Get a set of spectra to analyze. Spectra were sorted by precursor mass
      if(totalMass<minMass) bt=false;
      else bt=d.getBoundaries2(totalMass,params.ppmPrecursor,index);

      if(bt){

        //Handle looplinks differently
        if(d.getLink(n).mono==0){

          //Iterate all possible loops across the same peptide
          for(k=0;k<p->at(i).vK->size()-1;k++){
            for(k2=k+1;k2<p->at(i).vK->size();k2++){               
              ions.buildLoopIons(d.getLink(n).mass,p->at(i).vK->at(k)-p->at(i).map->at(0).start,p->at(i).vK->at(k2)-p->at(i).map->at(0).start);
              scoreSpectra(index,d,totalMass,NULL,crossLink,i,-1,p->at(i).vK->at(k),p->at(i).vK->at(k2),n);
            }
          }

        } else {

          //Iterate every lysine as a site for the monolink
          for(k=0;k<p->at(i).vK->size();k++){
            ions.buildIons(d.getLink(n).mass,p->at(i).vK->at(k) - p->at(i).map->at(0).start);
            scoreSpectra(index,d,totalMass,NULL,crossLink,i,-1,p->at(i).vK->at(k),-1,n);
          }
        }
      }
    }

    //Intercept here to continue if using singlet analysis
    if(params.relaxedAnalysis>0) continue;

    //check crosslinks (including self links). cross-links are to every next peptide in the list until mass is greater than max
    for(j=i;j<p->size();j++){

      //Iterate all cross-link masses.
      for(n=0;n<d.sizeLink();n++){

        //skip mono-link masses
        if(d.getLink(n).mono==1) continue;
        totalMass=p->at(i).mass+p->at(j).mass+d.getLink(n).mass;
        if(totalMass>maxMass) continue;

        //Get a set of spectra to analyze. Spectra were sorted by precursor mass
        if(totalMass<minMass) bt=false;
        else bt=d.getBoundaries2(totalMass,params.ppmPrecursor,index);

        if(bt){

          //set second peptide
          ions.setPeptide(false,&db[p->at(j).map->at(0).index].sequence[p->at(j).map->at(0).start],p->at(j).map->at(0).stop-p->at(j).map->at(0).start+1,p->at(j).mass);

          //Iterate every combination of lysine linkages
          for(k=0;k<p->at(i).vK->size();k++){
            for(k2=0;k2<p->at(j).vK->size();k2++){

              //Calc the ions, and score it against the spectra
              ions.buildXIons(d.getLink(n).mass,p->at(i).vK->at(k) - p->at(i).map->at(0).start,p->at(j).vK->at(k2) - p->at(j).map->at(0).start);
              scoreSpectra(index,d,totalMass,NULL,crossLink,i,j,p->at(i).vK->at(k),p->at(j).vK->at(k2),n);

            } //k2
          }//k

        }//if(bt)

      }
    }
  }

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;

  //clean up memory
  p=NULL;
  return true;
}

bool KAnalysis::analyzePeptidesNonCovalent(KDatabase& db, KData& d){
  unsigned int i,j;
  int k;
  int iPercent;
  int iTmp;
  kPeptide pep;
  kPeptideB pepB;
  vector<kPeptideB> p;
  vector<int> index;
  double totalMass;
  double maxMass = d.getMaxMass()+0.25;
  double minMass = d.getMinMass()-0.25;
  bool bt;

  //Combine list of linkable and non-linkable peptides
  pepB.linkable=false;
  for(k=0;k<db.getPeptideListSize(false);k++){
    pepB.mass=db.getPeptide(k,false).mass;
    pepB.index=k;
    p.push_back(pepB);
  }
  pepB.linkable=true;
  for(k=0;k<db.getPeptideListSize(true);k++){
    pepB.mass=db.getPeptide(k,true).mass;
    pepB.index=k;
    p.push_back(pepB);
  }

  //sort list by mass
  qsort(&p[0],p.size(),sizeof(kPeptideB),comparePeptideBMass);

  /*
  char peptide[256];
  for(i=0;i<p.size();i++){
    pep=db.getPeptide(p[i].index,p[i].linkable);
    db.getPeptideSeq( pep.map->at(0).index,pep.map->at(0).start,pep.map->at(0).stop,peptide);
    cout << peptide << "\t" << p[i].index << "\t" << p[i].linkable << "\t" << p[i].mass << endl;
  }
  exit(0);
  */

  //Set progress meter
  iPercent=0;
  printf("Progress: %2d%%",iPercent);
  fflush(stdout);

  //Iterate entire peptide list
  for(i=0;i<p.size();i++){

    //Update progress meter
    iTmp=(int)(i*100.0/p.size());
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }

    //Peptides are sorted by mass. If greater than max mass, stop checking peptides
    if(p[i].mass>maxMass) break;

    //Set the peptide, calc the ions, and score it against the spectra
    pep=db.getPeptide(p[i].index,p[i].linkable);
    ions.setPeptide(true,&db[pep.map->at(0).index].sequence[pep.map->at(0).start],pep.map->at(0).stop-pep.map->at(0).start+1,p[i].mass);

    //check dimerization with another peptide (including itself)
    for(j=i;j<p.size();j++){

      totalMass=p[i].mass+p[j].mass;
      if(totalMass>maxMass) break;
      //cout << totalMass << endl;

      //Get a set of spectra to analyze. Spectra were sorted by precursor mass
      if(totalMass<minMass) bt=false;
      else bt=d.getBoundaries2(totalMass,params.ppmPrecursor,index);

      if(bt){
        //cout << "In bounds" << endl;

        //set second peptide
        pep=db.getPeptide(p[j].index,p[j].linkable);
        ions.setPeptide(false,&db[pep.map->at(0).index].sequence[pep.map->at(0).start],pep.map->at(0).stop-pep.map->at(0).start+1,pep.mass);
              
        //Calc the ions, and score it against the spectra   
        ions.buildNCIons();
        //cout << "Build ions ok" << endl;
        scoreNCSpectra(index,d,totalMass,p[i].linkable,p[j].linkable,p[i].index,p[j].index);
        //cout << "Score ok" << endl;

      }//if(bt)

    }
  }

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;

  return true;
}

bool KAnalysis::analyzeRelaxed(KDatabase& db, KData& d){
  int i,j,k,n,m;
  int index;
  int iPercent;
  int iTmp;

  double ppm;
  double totalMass;

  kSingletScoreCardPlus* s=new kSingletScoreCardPlus[params.topCount];
  kSingletScoreCard sc1;
  kScoreCard sc;

  //Set progress meter
  iPercent=0;
  printf("Progress: %2d%%",iPercent);
  fflush(stdout);

  for(i=0;i<d.size();i++){

    //cout << d[i].getScanNumber() << endl;
   
    if(d[i].getScanNumber()==params.diagnostic){
      FILE* f=fopen("diagnostic.txt","wt");
      fprintf(f,"Scan: %d\n",d[i].getScanNumber());
      char strs[256];
      for(k=0;k<params.topCount;k++){
        sc1=d[i].getSingletScoreCard(k);
        db.getPeptideSeq( db.getPeptideList(sc1.linkable)[sc1.pep1].map->at(0).index,db.getPeptideList(sc1.linkable)[sc1.pep1].map->at(0).start,db.getPeptideList(sc1.linkable)[sc1.pep1].map->at(0).stop,strs);
        fprintf(f,"%s\t%d\t%.6lf\t%.4lf\t%.4lf\n",strs,sc1.k1,db.getPeptideList(sc1.linkable)[sc1.pep1].mass,sc1.simpleScore,sc1.simpleScore*sc1.len);
      }
      fclose(f);
    }

    //Update progress meter
    iTmp=(int)(i*100.0/d.size());
    if(iTmp>iPercent){
      iPercent=iTmp;
      printf("\b\b\b%2d%%",iPercent);
      fflush(stdout);
    }

    //Make a sortable list of top hits to compare
    for(j=0;j<params.topCount;j++){
      sc1=d[i].getSingletScoreCard(j);
      s[j].len=sc1.len;
      s[j].k1=sc1.k1;
      s[j].linkable=sc1.linkable;
      s[j].pep1=sc1.pep1;
      s[j].rank=j;
      s[j].simpleScore=sc1.simpleScore;
      if(sc1.simpleScore>0)s[j].mass=db.getPeptideList(sc1.linkable)[sc1.pep1].mass;
      else s[j].mass=0;
    }
    qsort(s,params.topCount,sizeof(kSingletScoreCardPlus),compareSSCPlus);

    //Check true cross-links
    for(k=0;k<d.sizeLink();k++){
      if(d.getLink(k).mono>0) continue;

      for(j=0;j<params.topCount;j++){
        
        if(s[j].simpleScore<=0) continue; //skips anything we've already considered
        if(!s[j].linkable) continue;
        if(s[j].k1<0) continue;

        for(m=0;m<d[i].sizePrecursor();m++){
          index=findMass(s,d[i].getPrecursor(m).monoMass-s[j].mass-d.getLink(k).mass);
          n=index;
          while(n<params.topCount){
            if(s[n].simpleScore<0 || s[n].k1<0){
              n++;
              continue;
            }
            totalMass=s[j].mass+s[n].mass+d.getLink(k).mass;
            ppm = (totalMass-d[i].getPrecursor(m).monoMass)/d[i].getPrecursor(m).monoMass*1e6;
            //printf("%.6lf\t%.6lf\t%.4lf\t%.2lf\n",s[j].mass,s[n].mass,ppm,s[j].simpleScore+s[n].simpleScore);
            if(fabs( ppm )<=params.ppmPrecursor) {
              sc.simpleScore=s[j].simpleScore*s[j].len+s[n].simpleScore*s[n].len;
              sc.k1=s[j].k1;
              sc.k2=s[n].k1;
              sc.mass=totalMass;
              //sc.modMass=0;
              sc.linkable1=s[j].linkable;
              sc.linkable2=s[n].linkable;
              sc.pep1=s[j].pep1;
              sc.pep2=s[n].pep1;
              sc.link=k;
              sc.rank=s[j].rank+s[n].rank;
              d[i].checkScore(sc);
            } else if(ppm>params.ppmPrecursor) {
              //cout << "break" << endl;
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
            totalMass=s[j].mass+s[n].mass+d.getLink(k).mass;
            ppm = (totalMass-d[i].getPrecursor(m).monoMass)/d[i].getPrecursor(m).monoMass*1e6;
            //printf("%.6lf\t%.6lf\t%.4lf\t%.2lf\n",s[j].mass,s[n].mass,ppm,s[j].simpleScore+s[n].simpleScore);
            if(fabs( ppm )<=params.ppmPrecursor) {
              sc.simpleScore=s[j].simpleScore*s[j].len+s[n].simpleScore*s[n].len;
              sc.k1=s[j].k1;
              sc.k2=s[n].k1;
              sc.mass=totalMass;
              //sc.modMass=0;
              sc.linkable1=s[j].linkable;
              sc.linkable2=s[n].linkable;
              sc.pep1=s[j].pep1;
              sc.pep2=s[n].pep1;
              sc.link=k;
              sc.rank=s[j].rank+s[n].rank;
              d[i].checkScore(sc);
            } else if(ppm<-params.ppmPrecursor) {
              break;
              //cout << "break" << endl;
            }
            n--;
          }
        }
        s[j].simpleScore=-s[j].simpleScore;
      }
    }

    //reset scores
    for(j=0;j<params.topCount;j++){
      if(s[j].simpleScore<0) s[j].simpleScore=-s[j].simpleScore;
    }

    //Check non-covalent dimers
    if(params.dimers==0) continue;
    for(j=0;j<params.topCount;j++){
      if(s[j].simpleScore<=0) continue;
      if(s[j].k1>-1) continue;

      for(m=0;m<d[i].sizePrecursor();m++){
        index=findMass(s,d[i].getPrecursor(m).monoMass-s[j].mass);
        n=index;
        while(n<params.topCount){
          if(s[n].simpleScore<0){
            n++;
            continue;
          }
          totalMass=s[j].mass+s[n].mass;
          ppm = (totalMass-d[i].getPrecursor(m).monoMass)/d[i].getPrecursor(m).monoMass*1e6;
          if(fabs( ppm )<=params.ppmPrecursor) {
            sc.simpleScore=s[j].simpleScore*s[j].len+s[n].simpleScore*s[n].len;
            sc.k1=-1;
            sc.k2=-1;
            sc.mass=totalMass;
            //sc.modMass=0;
            sc.linkable1=s[j].linkable;
            sc.linkable2=s[n].linkable;
            sc.pep1=s[j].pep1;
            sc.pep2=s[n].pep1;
            sc.link=-2;
            sc.rank=s[j].rank+s[n].rank;
            d[i].checkScore(sc);
          } else if(ppm>params.ppmPrecursor) {
            break;
          }
          n++;
        }
        n=index-1;
        while(n>-1){
          if(s[n].simpleScore<0){
            n--;
            continue;
          }
          totalMass=s[j].mass+s[n].mass;
          ppm = (totalMass-d[i].getPrecursor(m).monoMass)/d[i].getPrecursor(m).monoMass*1e6;
          if(fabs( ppm )<=params.ppmPrecursor) {
            sc.simpleScore=s[j].simpleScore*s[j].len+s[n].simpleScore*s[n].len;
            sc.k1=-1;
            sc.k2=-1;
            sc.mass=totalMass;
            //sc.modMass=0;
            sc.linkable1=s[j].linkable;
            sc.linkable2=s[n].linkable;
            sc.pep1=s[j].pep1;
            sc.pep2=s[n].pep1;
            sc.link=-2;
            sc.rank=s[j].rank+s[n].rank;
            d[i].checkScore(sc);
          } else if(ppm<-params.ppmPrecursor) {
            break;
          }
          n--;
        }
      }
      s[j].simpleScore=-s[j].simpleScore;
    }

  }

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;

  delete [] s;
  return true;
}

bool KAnalysis::analyzeSinglets(KDatabase& db, KData& d, kPeptide& pep, int index, double lowLinkMass, double highLinkMass){
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
  ions.setPeptide(true,&db[pep.map->at(0).index].sequence[pep.map->at(0).start],len,pep.mass);

  //Iterate all spectra from (peptide mass + low linker + minimum mass) to (peptide mass + high linker + maximum mass)
  if(!d.getBoundaries(minMass,maxMass,scanIndex)) return false;

  //Iterate every lysine as a site for the monolink
  for(k=0;k<pep.vK->size();k++){

    //build fragment ions and score against all potential spectra
    ions.buildSingletIons(pep.vK->at(k)-pep.map->at(0).start);

    for(j=0;j<scanIndex.size();j++){
      scoreSingletSpectra(scanIndex[j],d,pep.mass,len,index,pep.vK->at(k),true,minMass);
    }

  }

  return true;
}

bool KAnalysis::analyzeSingletsNoLysine(KDatabase& db, KData& d, kPeptide& pep, int index, bool linkable){
  unsigned int j;
  double minMass;
  double maxMass;
  vector<int> scanIndex;

  //Set Mass boundaries
  minMass=pep.mass+params.minPepMass;
  maxMass=pep.mass+params.maxPepMass;
  minMass-=(minMass/1000000*params.ppmPrecursor);
  maxMass+=(maxMass/1000000*params.ppmPrecursor);

  //Iterate all spectra from (peptide mass + low linker + minimum mass) to (peptide mass + high linker + maximum mass)
  if(!d.getBoundaries(minMass,maxMass,scanIndex)) return false;
  for(j=0;j<scanIndex.size();j++){
    scoreSingletSpectra(scanIndex[j],d,pep.mass,pep.map->at(0).stop-pep.map->at(0).start+1,index,-1,linkable,minMass);
  }
  return true;

}


/*============================
  Private Functions
============================*/
int KAnalysis::findMass(kSingletScoreCardPlus* s, double mass){
  int sz=params.topCount;
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

void KAnalysis::scoreNCSpectra(vector<int>& index, KData& d, double mass, bool linkable1, bool linkable2, int pep1, int pep2){
  unsigned int a;
  kScoreCard sc;
 
  //score spectra
  for(a=0;a<index.size();a++){
    if(params.xcorr==1) sc.simpleScore=xCorrScoring(d[index[a]]);
    else sc.simpleScore=simpleScoring(d[index[a]]);
    sc.k1=-1;
    sc.k2=-1;
    sc.mass=mass;
    //sc.modMass=0;
    sc.linkable1=linkable1;
    sc.linkable2=linkable2;
    sc.pep1=pep1;
    sc.pep2=pep2;
    sc.link=-2;
    d[index[a]].checkScore(sc);
  }
}

void KAnalysis::scoreSingletSpectra(int index, KData& d, double mass, int len, int pep, int k, bool linkable, double minMass){
  //cout << "scoreSingletSpectra " << d[index].getScanNumber() << "\t" << ions.getIonCount() << "\t" << pep << "\t" << k << "\t" << linkable << "\t" << mass << endl;
  kSingletScoreCard sc;
  double score1=0;
  double score2=0;
  int i;

  for(i=0;i<d[index].sizePrecursor();i++){
    if(d[index].getPrecursor(i).monoMass>minMass){
      if(params.xcorr==1) {
        if(k>-1) score1=xCorrScoring2(d[index],d[index].getPrecursor(i).monoMass-mass);
        else score1=xCorrScoring(d[index]);
      } else {
        score1=simpleScoring(d[index]);
      }
      if(score1>score2) score2=score1;
    }
  }

  //Check the highest score
  sc.len=len;
  sc.simpleScore=(float)score2/len;
  sc.k1=k;
  sc.linkable=linkable;
  sc.pep1=pep;
  if(sc.simpleScore>0) d[index].checkSingletScore(sc);

}

void KAnalysis::scoreSpectra(vector<int>& index, KData& d, double mass, vector<kPepMod>* v, bool linkable, int pep1, int pep2, int k1, int k2, int link){
  unsigned int a;
  kScoreCard sc;
 
  //scoreSpectra(index,d,p->at(i).mass,0,crossLink,i,-1,-1,-1,-1);

  //score spectra
  for(a=0;a<index.size();a++){
    sc.mods->clear();
    if(params.xcorr==1) sc.simpleScore=xCorrScoring(d[index[a]]);
    else sc.simpleScore=simpleScoring(d[index[a]]);
    sc.k1=k1;
    sc.k2=k2;
    sc.mass=mass;
    sc.linkable1=sc.linkable2=linkable;
    sc.pep1=pep1;
    sc.pep2=pep2;
    sc.link=link;
    if(v!=NULL){
      for(unsigned int i=0;i<v->size();i++) sc.mods->push_back(v->at(i));
    }
    d[index[a]].checkScore(sc);
  }
}

//This is the PepLynx specific scoring algorithm. It uses number of matched ions weighted by intensity (with some extra details).
float KAnalysis::simpleScoring(KSpectrum& s){
  int k;
  int maxCharge;
  int ionCount=ions.getIonCount();

  double massXplus;
  double massXminus;
  double mass2Xminus;
  double mass2Xplus;
  double tol=params.ppmFragment/1000000.0;
  double tol2X=tol*2;
  
  float seriesScoresB[6];
  float seriesScoresY[6];
  float simpleScore=0;

	unsigned int i;
  int j[6];
  int n[6];
  unsigned int SpecSize=s.size();

  //Get the maximum charge state for fragment ions.
  //This is PrecursorCharge - 1
  maxCharge=s.getCharge();
  if(maxCharge>6) maxCharge=6;
  for(k=0;k<6;k++){
    j[k]=0;
    n[k]=0;
    seriesScoresB[k]=0;
    seriesScoresY[k]=0;
  }

  i=0;

  //Iterate through the spectrum, matching to any b- and y-ions in a ratchet fashion.
  //Scores are weighted by the percentile of the peak being matched. There are 10 percentile groups with
  //weights ranging from 0.1 to 1.0
  while(i<SpecSize){
    massXminus=s[i].mass-tol*s[i].mass;
    massXplus=s[i].mass+tol*s[i].mass;
    mass2Xminus=s[i].mass-tol2X*s[i].mass;
    mass2Xplus=s[i].mass+tol2X*s[i].mass;
    for(k=1;k<maxCharge;k++){

      //Check b-ions
      while(j[k]<ionCount && ions.bIons[k][j[k]]<mass2Xminus) j[k]++;
      
      if(j[k]<ionCount && ions.bIons[k][j[k]]<mass2Xplus){
        if(ions.bIons[k][j[k]] > massXminus && ions.bIons[k][j[k]] < massXplus) seriesScoresB[k]+=s[i].intensity;
        else if(ions.bIons[k][j[k]] > mass2Xminus && ions.bIons[k][j[k]] < mass2Xplus) seriesScoresB[k]+=s[i].intensity*0.25f;
      }

      //Check y-ions
      while(n[k]<ionCount && ions.yIons[k][n[k]]<mass2Xminus) n[k]++;
      
      if(n[k]<ionCount && ions.yIons[k][n[k]]<mass2Xplus){
        if(ions.yIons[k][n[k]] > massXminus && ions.yIons[k][n[k]] < massXplus) seriesScoresY[k]+=s[i].intensity;
        else if(ions.yIons[k][n[k]] > mass2Xminus && ions.yIons[k][n[k]] < mass2Xplus) seriesScoresY[k]+=s[i].intensity*0.25f;
      }

    }
    i++;
  }
  
  //Sum the b- and y-ion scores at each charge state;
  //In the future, this can be changed to use only the top N charge series scores.
  for(i=0;i<6;i++){
    simpleScore+=seriesScoresB[i];
    simpleScore+=seriesScoresY[i];
  }

	return simpleScore;
}

//An alternative score uses the XCorr metric from the Comet algorithm
float KAnalysis::xCorrScoring(KSpectrum& s) { 

  //xCorrCount++;

  double dXcorr;
  double invBinSize=s.getInvBinSize();

  int ionCount=ions.getIonCount();
  int k;
  int maxCharge;
  int xx;

  int i;
  unsigned int j;

  unsigned int SpecSize=s.size();

  dXcorr=0.0;
  
  //Get the number of fragment ion series to analyze
  //The number is PrecursorCharge-1
  maxCharge=s.getCharge();
  if(maxCharge>6) maxCharge=6;

  //Iterate all series
  for(k=1;k<maxCharge;k++){

    //Ratchet through pfFastXcorrData
    xx=0;
    for(i=0;i<ionCount;i++){
      j = (unsigned int)(ions.bIons[k][i]*invBinSize+1.0);
      while( j >=  (unsigned) s.xCorrSparseArray[xx].bin) xx++;
      dXcorr += s.xCorrSparseArray[xx-1].fIntensity;
    }

    xx=0;
    for(i=0;i<ionCount;i++){
      j = (unsigned int)(ions.yIons[k][i]*invBinSize+1.0);
      while( j >=  (unsigned) s.xCorrSparseArray[xx].bin) xx++;
      dXcorr += s.xCorrSparseArray[xx-1].fIntensity;
    }
  }

  //Scale score appropriately
  if(dXcorr <= 0.0) dXcorr=0.0;
  else dXcorr *= 0.005;

  return float(dXcorr);
}

//An alternative score uses the XCorr metric from the Comet algorithm
float KAnalysis::xCorrScoring2(KSpectrum& s, double modMass) { 

  //xCorrCount++;

  double dXcorr;
  double invBinSize=s.getInvBinSize();

  int ionCount=ions.getIonCount();
  int k;
  int maxCharge;
  int xx;

  int i;
  unsigned int j;

  unsigned int SpecSize=s.size();

  dXcorr=0.0;
  
  //Get the number of fragment ion series to analyze
  //The number is PrecursorCharge-1
  maxCharge=s.getCharge();
  if(maxCharge>6) maxCharge=6;

  //Iterate all series
  for(k=1;k<maxCharge;k++){

    //Ratchet through pfFastXcorrData
    xx=0;
    for(i=0;i<ionCount;i++){
      if(ions.bIons[k][i]<0) j = (unsigned int)((modMass/k-ions.bIons[k][i])*invBinSize+1.0);
      else j = (unsigned int)(ions.bIons[k][i]*invBinSize+1.0);
      while( j >=  (unsigned) s.xCorrSparseArray[xx].bin) xx++;
      dXcorr += s.xCorrSparseArray[xx-1].fIntensity;
    }

    xx=0;
    for(i=0;i<ionCount;i++){
      if(ions.yIons[k][i]<0) j = (unsigned int)((modMass/k-ions.yIons[k][i])*invBinSize+1.0);
      else j = (unsigned int)(ions.yIons[k][i]*invBinSize+1.0);
      while( j >=  (unsigned) s.xCorrSparseArray[xx].bin) xx++;
      dXcorr += s.xCorrSparseArray[xx-1].fIntensity;
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
