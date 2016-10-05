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

#include "KPrecursor.h"

/*============================
  Constructors & Destructors
============================*/
KPrecursor::KPrecursor(kParams* p){
  params=p;

  hs2.noBase=true;
  hs4.noBase=true;

  hv.addAtom(107,2);
  hv.addEnrich(107,2,params->enrichment);
  hs2.variant->push_back(hv);

  hv.clear();
  hv.addAtom(107,4);
  hv.addEnrich(107,2,params->enrichment);
  hs4.variant->push_back(hv);

  hv.clear();
  hs.winSize=hs2.winSize=hs4.winSize=10;
  hs.peptide=hs2.peptide=hs4.peptide=4;
  hs.sn=hs2.sn=hs4.sn=0;
  hs.depth=hs2.depth=hs4.depth=3;
  hs.minCharge=hs2.minCharge=hs4.minCharge=2;
  hs.maxCharge=hs2.maxCharge=hs4.maxCharge=8;
  hs.algorithm=Version2;
  if(params->instrument==1) hs.msType=hs2.msType=hs4.msType=FTICR;
  else hs.msType=hs2.msType=hs4.msType=OrbiTrap;
  hs.res400=hs2.res400=hs4.res400=params->ms1Resolution;
  hs.corr=0.875;
  hs2.corr=hs4.corr=0.85;
  hs.centroid=hs2.centroid=hs4.centroid=true;

  strcpy(hs.inFile,"PLTmp.ms1");
  strcpy(hs2.inFile,"PLTmp.ms1");
  strcpy(hs4.inFile,"PLTmp.ms1");

  hs.fileFormat=hs2.fileFormat=hs4.fileFormat=ms1;

  averagine = new CAveragine(NULL,NULL);
  mercury = new CMercury8(NULL);
  models = new CModelLibrary(averagine,mercury);

  h = new CHardklor2(averagine,mercury,models);
  hO = new CHardklor(averagine,mercury);

  CHardklorVariant hkv;
  vector<CHardklorVariant> pepVariants;
  pepVariants.clear();
  pepVariants.push_back(hkv);

  models->eraseLibrary();
  models->buildLibrary(2,8,pepVariants);

  h->Echo(false);
  h->SetResultsToMemory(true);

  hO->Echo(false);
  hO->SetResultsToMemory(true);

  centBuf = new deque<Spectrum>;

  setFile();

}

KPrecursor::~KPrecursor(){
  delete h;
  delete hO;
  delete averagine;
  delete mercury;
  delete models;
  delete centBuf;
  params = NULL;
}

/*============================
  Public Functions
============================*/

bool KPrecursor::estimatePrecursor(KSpectrum& s){
  char formula[256];
  double mass;
  double offset;
  double dif;
  unsigned int i,j;
  kPrecursor pre;

  //only estimate if precursor charge is known or assumed
  if(s.getCharge()<1 && preCharges.size()==0) return false;  

  //Check if precursor charge is already known. If so, use only it.
  if(s.getCharge()>0){
    preCharges.clear();
    preCharges.push_back(s.getCharge());
  }

  for(j=0;j<preCharges.size();j++){

    mass = s.getMZ()*preCharges[j]-preCharges[j]*1.007276466;
    averagine->clear();
    averagine->calcAveragine(mass,hv);
    averagine->getAveragine(formula);

    mercury->GoMercury(formula);
    for(i=0;i<mercury->FixedData.size();i++){
      if(mercury->FixedData[i].data>99.99) break;
    }
    offset=mercury->FixedData[i].mass-mass;

    pre.charge=preCharges[j];
    pre.monoMass=mass;
    pre.corr=-0.5;
    s.addPrecursor(pre);

    //if base peak is not monoisotopic, try one peak to the left.
    //Could try all peaks until monoisotopic is found, but might not make any difference,
    //only add time.
    if(i>0){
      dif=mercury->FixedData[i].mass-mercury->FixedData[i-1].mass;
      pre.monoMass=mercury->FixedData[i-1].mass-offset;
      s.addPrecursor(pre);
      if(pre.monoMass>3000 && i>1){
        pre.monoMass=mercury->FixedData[i-2].mass-offset;
        s.addPrecursor(pre);
      }
      if(pre.monoMass>5000 && i>2){
        pre.monoMass=mercury->FixedData[i-3].mass-offset;
        s.addPrecursor(pre);
      }
    }

    //also add the next peak?
    //pre.monoMass=mercury->FixedData[i+1].mass-offset;
    //s.addPrecursor(pre);
  }

  return true;
}

//Finds the 18O2 and 18O4 precursor monoisotopic masses and charge states for a given MS/MS scan.
//Return values:
//0 = No precursor found
//1 = Found expected precursor
//2 = Found alternate precursor(s)
int KPrecursor::getSpecRange(KSpectrum& pls){

  unsigned int i;

  int charge;
  int j;
  int k;
  int precursor;
  int best;
  int scanNum=pls.getScanNumber();
  int ret=0;

  float maxIntensity;
  float maxRT;
  float rt;
  float rtHigh;
  float rtLow;
  
  double corr;
  double monoMass;
  double mz;
  double tmz;

  kPrecursor pre;
  Spectrum s;
  vector<Spectrum*> vs;

  //Look up original MS2
  msr.setFilter(MS2);
  msr.readFile(NULL,spec,scanNum);
  rt=spec.getRTime();
  mz=spec.getMZ();

  //clear scan buffer of unneeded spectra
  while(centBuf->size()>0 && centBuf->at(0).getRTime()<rt-1.0){
    centBuf->pop_front();
  }

  //if whole buffer was emptied, read in files until we find the scan range needed.  
  if(centBuf->size()==0){
    msr.setFilter(MS1);
    msr.readFile(NULL,spec,lastScan);
    msr.readFile(NULL,spec);
    while(spec.getScanNumber()>0 && spec.getRTime()<rt-1.0){
      msr.readFile(NULL,spec);
    }
    if(spec.getScanNumber()>0) {
      if(!params->ms1Centroid){
        centroid(spec,cent,params->ms1Resolution,params->instrument);
        cent.setScanNumber(spec.getScanNumber());
        cent.setRTime(spec.getRTime());
        centBuf->push_back(cent);
      } else {
        centBuf->push_back(spec);
      }
      lastScan=spec.getScanNumber();
    }
  }

  //Sanity check. Exit if you can't find precursor scans within range of MS/MS retention time
  if(centBuf->size()==0){
    //cout << "Reached end of file buffer and no spectra remain." << endl;
    return ret;
  }

  if(centBuf->at(centBuf->size()-1).getRTime()<rt+0.6){
    msr.setFilter(MS1);
    msr.readFile(NULL,spec,lastScan);
    msr.readFile(NULL,spec);
    while(spec.getScanNumber()>0 && spec.getRTime()<rt+2.0){
      if(!params->ms1Centroid){
        centroid(spec,cent,params->ms1Resolution,params->instrument);
        cent.setScanNumber(spec.getScanNumber());
        cent.setRTime(spec.getRTime());
        centBuf->push_back(cent);
      } else {
        centBuf->push_back(spec);
      }
      lastScan=spec.getScanNumber();
      msr.readFile(NULL,spec);
    }
    if(spec.getScanNumber()>0){
      if(!params->ms1Centroid){
        centroid(spec,cent,params->ms1Resolution,params->instrument);
        cent.setScanNumber(spec.getScanNumber());
        cent.setRTime(spec.getRTime());
        centBuf->push_back(cent);
      } else {
        centBuf->push_back(spec);
      }
      lastScan=spec.getScanNumber();
    }
  }

  //Find a precursor within 10 seconds for a given scan 
  maxIntensity=0;
  precursor=-1;
  for(i=0;i<centBuf->size();i++){
    if(centBuf->at(i).getRTime()<rt-0.167) continue;
    if(centBuf->at(i).getRTime()>rt+0.167) break;
    j=findPeak(centBuf->at(i),mz,10);
    
    if(j>-1){
      if(centBuf->at(i)[j].intensity>maxIntensity){
        maxIntensity=centBuf->at(i)[j].intensity;
        maxRT=centBuf->at(i).getRTime();
        best=centBuf->at(i).getScanNumber();
        precursor=i;
      }
    }
  }

  if(precursor<0){
    //cout << "Warning: Precursor not found for " << scanNum << " " << mz << endl;
    return ret;
  }
  
  //Get up to +/-15 sec of spectra around the max precursor intensity
  //This is done by extending on both sides until a gap is found or time is reached.
  //Additionally, stop if 2 scans found flanking either side (maximum 5 scans per precursor).
  vs.clear();
  rtHigh=rtLow=centBuf->at(precursor).getRTime();
  k=0;
  for(i=precursor;i<centBuf->size();i++){
    if(centBuf->at(i).getRTime()>maxRT+0.25) break;
    j=findPeak(centBuf->at(i),mz,10);
    if(j<0) break;
    vs.push_back(&centBuf->at(i));
    rtHigh=centBuf->at(i).getRTime();
    k++;
    if(k==3) break;
  }

  k=0;
  i=precursor;
  while(i>0){
    i--;
    if(centBuf->at(i).getRTime()<maxRT-0.25) break;
    j=findPeak(centBuf->at(i),mz,10);
    if(j<0) break;
    vs.push_back(&centBuf->at(i));
    rtLow=centBuf->at(i).getRTime();
    k++;
    if(k==2) break;
  }

  //If total spectra to be combined is small (5 scan events), try 
  //grabbing a +/- 5 sec window around max, regardless of precursor observation
  /*
  if(vs.size()<5 && (rtHigh-rtLow)<0.167){
    vs.clear();
    for(i=precursor;i<scanBuf->size();i++){
      if(scanBuf->at(i).getRTime()<maxRT+0.083) vs.push_back(&scanBuf->at(i));
      else break;
    }
    for(i=precursor-1;i>-1;i--){
      if(scanBuf->at(i).getRTime()<maxRT-0.083) vs.push_back(&scanBuf->at(i));
      else break;
    }
  }
  */

  //Average points between mz-1.5 and mz+2
  if(params->enrichment>0)averageScansCentroid(vs,s,mz-1.75,mz+1.75);
  else averageScansCentroid(vs,s,mz-1.0,mz+1.5);
  s.setScanNumber(centBuf->at(precursor).getScanNumber());

  //Obtain the possible precursor charge states of the selected ion.
  //Find the index of the closest peak to the selected m/z.
  preCharges.clear();
  tmz=fabs(mz-s[0].mz);
  for(j=1;j<s.size();j++){
    if(fabs(mz-s[j].mz)<tmz) tmz=fabs(mz-s[j].mz);
    else break;
  }
  j=j-1;
  h->QuickCharge(s,j,preCharges);

  //Clear corr
  corr=0;

  //Perform 18O2 analysis with Hardklor. If enrichment is set to 0, store unenriched results in the *O2 variables.
  //This is done in a non-competitive to identify an 18O2 peptide without solving everything
  if(params->enrichment>0) {
    hO->GoHardklor(hs2,&s);

    for(j=0;j<hO->Size();j++){
      if(hO->operator[](j).corr>corr){
        monoMass=hO->operator[](j).monoMass;
        charge=hO->operator[](j).charge;
        corr=hO->operator[](j).corr;
      }
    }
    if(corr>0){
      pre.monoMass=monoMass;
      pre.charge=charge;
      pre.corr=corr;
      if(params->enrichment>0) pre.label=1;
      else pre.label=0;
      pls.addPrecursor(pre);
    }

  } else {
    
    h->GoHardklor(hs,&s);

    //If nothing was found, really narrow down the window and try again.
    if(h->Size()==0){
      averageScansCentroid(vs,s,mz-0.6,mz+1.2);
      s.setScanNumber(centBuf->at(precursor).getScanNumber());
      if(params->enrichment>0) h->GoHardklor(hs2,&s);
      else h->GoHardklor(hs,&s);
    }

    float intensity=0;
    for(j=0;j<h->Size();j++){

      //Must have highest intensity and intersect isolated peak.
      if(h->operator[](j).intensity<intensity) continue;
      tmz=(h->operator[](j).monoMass+1.007276466*h->operator[](j).charge)/h->operator[](j).charge;
      while(tmz<(pls.getMZ()+0.01)){
        if(fabs(tmz-pls.getMZ())<0.01){
          monoMass=h->operator[](j).monoMass;
          charge=h->operator[](j).charge;
          corr=h->operator[](j).corr;
          intensity=h->operator[](j).intensity;
          ret=1;
          break;
        }
        tmz+=(1.00335483/h->operator[](j).charge);
      }
    }

    //failing to match precursor peak, keep most intense precursor in presumed isolation window
    if(corr==0){
      for(j=0;j<h->Size();j++){
        if(h->operator[](j).intensity>intensity){
          monoMass=h->operator[](j).monoMass;
          charge=h->operator[](j).charge;
          corr=h->operator[](j).corr;
          intensity=h->operator[](j).intensity;
          ret=2;
        }
      }
    }

    if(corr>0){
      pre.monoMass=monoMass;
      pre.charge=charge;
      pre.corr=corr;
      if(params->enrichment>0) pre.label=1;
      else pre.label=0;
      pls.addPrecursor(pre);
      //also add isotope error
      if (params->isotopeError>0){
        pre.monoMass-=1.00335483;
        pre.corr=-1;
        pls.addPrecursor(pre);
      }
      if (params->isotopeError>1){
        pre.monoMass -= 1.00335483;
        pre.corr = -2;
        pls.addPrecursor(pre);
      }
    }

  }

  //Clear corr
  corr=0;

  //Perform the 18O4 Hardklor analysis on the same data if needed
  if(params->enrichment>0) {
    hO->GoHardklor(hs4,&s);
    for(j=0;j<hO->Size();j++){
      if(hO->operator[](j).corr>corr){
        monoMass=hO->operator[](j).monoMass;
        charge=hO->operator[](j).charge;
        corr=hO->operator[](j).corr;
        ret=1;
      }
    }
  }
  if(corr>0){
    pre.monoMass=monoMass;
    pre.charge=charge;
    pre.corr=corr;
    pre.label=2;
    pls.addPrecursor(pre);
  }

  //Assume two precursors with nearly identical mass (within precursor tolerance) are the same.
  //This can occur when checking multiple enrichment states.
  //Keep only the higher correlated precursor.
  if(pls.sizePrecursor()>1){
    bool bCheck=true;
    while(bCheck){
      bCheck=false;
      for(k=0;k<pls.sizePrecursor()-1;k++){
        for(j=k+1;j<pls.sizePrecursor();j++){
          if(fabs(pls.getPrecursor(k).monoMass-pls.getPrecursor(j).monoMass)/pls.getPrecursor(k).monoMass*1e6 < params->ppmPrecursor){
            if(pls.getPrecursor(k).corr>pls.getPrecursor(j).corr) pls.erasePrecursor(j);
            else pls.erasePrecursor(k);
            bCheck=true;
            break;
          }
        }
        if(bCheck) break;
      }
    }
  }

  //If we have something, return it
  return ret;
  //if(pls.sizePrecursor()>0)  return true;
  //else return false;

}

double KPrecursor::polynomialBestFit(vector<double>& x, vector<double>& y, vector<double>& coeff, int degree){
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

//This function starts the file buffer
bool KPrecursor::setFile() {

  msr.setFilter(MS1);
  msr.readFile(params->msFile,spec);
  strcpy(fileName,params->msFile);

  if(spec.getScanNumber()==0) return false;
  lastScan=spec.getScanNumber();

  centBuf->clear();
  if(!params->ms1Centroid){
    centroid(spec,cent,params->ms1Resolution,params->instrument);
    cent.setScanNumber(spec.getScanNumber());
    centBuf->push_back(cent);
  } else {
    centBuf->push_back(spec);
  }

  return true;
}

/*============================
  Private Functions
============================*/
//Places centroided peaks in bins, then finds average. Has potential for error when accuracy drifts into neighboring bins.
void KPrecursor::averageScansCentroid(vector<Spectrum*>& s, Spectrum& avg, double min, double max){

  unsigned int i;
  int j,k;
  double binWidth;
  double offset=-1.0;
  float* bin;
  float intensity;
  int binCount;
  int* pos;
  double  lowMZ=-1.0;
  kScanBin sb;
  vector<kScanBin> topList;

  pos=new int[s.size()];

  avg.clear();

  //Set some really small bin width related to the mz region being summed.
  binWidth=0.0001;
  
  binCount=(int)((max-min)/binWidth+1);
  bin = new float[binCount];
  for(j=0;j<binCount;j++) bin[j]=0;

  //align all spectra to point closest to min and set offset
  for(i=0;i<s.size();i++)  {
    pos[i]=findPeak(*s[i],min);
    if(s[i]->at(pos[i]).mz<min) pos[i]++;
    if(offset<0) offset=s[i]->at(pos[i]).mz;
    else if(s[i]->at(pos[i]).mz<offset) offset=s[i]->at(pos[i]).mz;
    //cout << s[i]->getScanNumber() << "\t" << pos[i] << "\t" << s[i]->at(pos[i]).mz << endl;
  }
  //printf("Offset: %.5lf\n",offset);

  //Iterate all spectra and add peaks to bins
  for(i=0;i<s.size();i++) {
    while(pos[i]<s[i]->size() && s[i]->at(pos[i]).mz<max){ 
      j=(int)((s[i]->at(pos[i]).mz-offset)/binWidth);
      //cout << i << "\t" << pos[i] << "\t" << j << endl;
      if(j<0){
        pos[i]++;
        continue;
      }
      if(j>=binCount) break;
      bin[j]+=s[i]->at(pos[i]).intensity;
      pos[i]++;
    }
  }

  //Unsure of current efficiency. Finds bin of tallest peak, then combines with neighboring bins
  //to produce an average. Thus summing of neighboring bins allows flexibility when dealing with mass accuracy drift.
  //Using larger bins has the same effect, but perhaps fewer significant digits in the final result.
  for(j=0;j<binCount;j++){
    if(bin[j]>0) {
      sb.index=j;
      sb.intensity=bin[j];
      topList.push_back(sb);
    }
  }
  qsort(&topList[0],topList.size(),sizeof(kScanBin),compareScanBinRev);
  for(i=0;i<topList.size();i++){
    if(bin[topList[i].index]==0) continue;
    intensity=0;
    j=topList[i].index-50;
    k=topList[i].index+51;
    if(j<0) j=0;
    if(k>binCount) k=binCount;
    for(j=j;j<k;j++){
      intensity+=bin[j];
      bin[j]=0;
    }
    intensity/=s.size();
    avg.add(offset+topList[i].index*binWidth,intensity);
  }
  avg.sortMZ();

  //clean up memory
  delete [] bin;
  delete [] pos;

}


//First derivative method, returns base peak intensity of the set
void KPrecursor::centroid(Spectrum& s, Spectrum& out, double resolution, int instrument){
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

//Returns closet mz value to desired point
int KPrecursor::findPeak(Spectrum& s, double mass){
  int sz=s.size();
  int lower=0;
  int mid=sz/2;
  int upper=sz;


  //binary search to closest mass
  while(s[mid].mz!=mass){
		if(lower>=upper) break;
    if(mass<s[mid].mz){
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

	//Check that mass is closest
  if(mid>0 && fabs(s[mid-1].mz-mass)<fabs(s[mid].mz-mass)) return mid-1;
  if(mid<s.size()-1 && fabs(s[mid+1].mz-mass)<fabs(s[mid].mz-mass)) return mid+1;
	return mid;

}

//Returns point within precision or -1 if doesn't exist
int KPrecursor::findPeak(Spectrum& s, double mass, double prec){
  int sz=s.size();
  int lower=0;
  int mid=sz/2;
  int upper=sz;

  double minMass = mass - (mass/1000000*prec);
  double maxMass = mass + (mass/1000000*prec);

  //binary search to closest mass
  while(s[mid].mz<minMass || s[mid].mz>maxMass){
		if(lower>=upper) break;
    if(mass<s[mid].mz){
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
	if(s[mid].mz>minMass && s[mid].mz<maxMass) return mid;

  return -1;
}

/*============================
  Private Utilities
============================*/
int KPrecursor::comparePLPScanNum(const void *p1, const void *p2){
  kRTProfile d1 = *(kRTProfile *)p1;
  kRTProfile d2 = *(kRTProfile *)p2;
  if(d1.scan<d2.scan) {
		return -1;
	} else if(d1.scan>d2.scan) {
  	return 1;
  } else {
	  return 0;
  }
}

int KPrecursor::compareScanBinRev(const void *p1, const void *p2){
  kScanBin d1 = *(kScanBin *)p1;
  kScanBin d2 = *(kScanBin *)p2;
  if(d1.intensity>d2.intensity) {
		return -1;
  } else if(d1.intensity<d2.intensity) {
  	return 1;
  } else {
	  return 0;
  }
}

