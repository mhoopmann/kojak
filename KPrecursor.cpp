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
  hs.depth=hs2.depth=hs4.depth=1;
  hs.minCharge=hs2.minCharge=hs4.minCharge=2;
  hs.maxCharge=hs2.maxCharge=hs4.maxCharge=8;
  if(params->instrument==1) hs.msType=hs2.msType=hs4.msType=FTICR;
  else hs.msType=hs2.msType=hs4.msType=OrbiTrap;
  hs.res400=hs2.res400=hs4.res400=params->ms1Resolution;
  hs.corr=hs2.corr=hs4.corr=0.85;
  if(params->ms1Centroid) hs.centroid=hs2.centroid=hs4.centroid=true;
  else hs.centroid=hs2.centroid=hs4.centroid=false;

  strcpy(hs.inFile,"PLTmp.ms1");
  strcpy(hs2.inFile,"PLTmp.ms1");
  strcpy(hs4.inFile,"PLTmp.ms1");

  hs.fileFormat=hs2.fileFormat=hs4.fileFormat=ms1;

  averagine = new CAveragine(NULL,NULL);
	mercury = new CMercury8(NULL);

  h = new CHardklor(averagine,mercury);
  h->Echo(false);
  h->SetResultsToMemory(true);

  scanBuf = new deque<Spectrum>;
  centBuf = new deque<Spectrum>;

  setFile();

  //params = NULL;
}

KPrecursor::~KPrecursor(){
  delete h;
  delete averagine;
  delete mercury;
  delete scanBuf;
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
  unsigned int i;
  kPrecursor pre;

  if(s.getCharge()<1) return false;  //only estimate if precursor charge is known

  mass = s.getMZ()*s.getCharge()-s.getCharge()*1.007276466;
  averagine->clear();
  averagine->calcAveragine(mass,hv);
  averagine->getAveragine(formula);

  mercury->GoMercury(formula);
  for(i=0;i<mercury->FixedData.size();i++){
    if(mercury->FixedData[i].data>99.99) break;
  }
  offset=mercury->FixedData[i].mass-mass;

  pre.charge=s.getCharge();
  pre.monoMass=mass;
  s.addPrecursor(pre);

  //if base peak is not monoisotopic, try one peak to the left.
  //Could try all peaks until monoisotopic is found, but might not make any difference,
  //only add time.
  if(i>0){
    dif=mercury->FixedData[i].mass-mercury->FixedData[i-1].mass;
    pre.monoMass=mercury->FixedData[i-1].mass-offset;
    s.addPrecursor(pre);
    if(pre.monoMass>3000){
      pre.monoMass=mercury->FixedData[i-1].mass-dif-offset;
      s.addPrecursor(pre);
    }
  }

  //also add the next peak
  pre.monoMass=mercury->FixedData[i+1].mass-offset;
  s.addPrecursor(pre);

  /*
  if(s.getScanNumber()==18211){
    cout << "\n" << i << endl;
    for(i=0;i<s.sizePrecursor();i++) cout << s.getPrecursor(i).monoMass << "\t" << s.getPrecursor(i).charge << "\t" << s.getPrecursor(i).corr << endl;
    exit(0);
  }
  */

  return true;
}

//Finds the 18O2 and 18O4 precursor monoisotopic masses and charge states for a given MS/MS scan.
bool KPrecursor::getSpecRange(KSpectrum& pls){

  unsigned int i;

  int charge;
  int j;
  int k;
  int precursor;
  int best;
  int lastScan;
  int scanNum=pls.getScanNumber();

  float maxIntensity;
  float maxRT;
  float rt;
  float rtHigh;
  float rtLow;
  
  double corr;
  double monoMass;
  double mz;

  kPrecursor pre;
  kRTProfile p;
  Spectrum s;
  vector<kRTProfile> v;
  vector<Spectrum*> vs;

  //Look up original MS2
  msr.setFilter(MS2);
  msr.readFile(NULL,spec,scanNum);
  rt=spec.getRTime();
  mz=spec.getMZ();

  //clear scan buffer of unneeded spectra
  while(scanBuf->size()>0 && scanBuf->at(0).getRTime()<rt-1.0){
    lastScan=scanBuf->at(0).getScanNumber();
    scanBuf->pop_front();
    centBuf->pop_front();
  }

  //if whole buffer was emptied, read in files until we find the scan range needed.  
  if(scanBuf->size()==0){
    msr.setFilter(MS1);
    msr.readFile(NULL,spec,lastScan);
    msr.readFile(NULL,spec);
    while(spec.getScanNumber()>0 && spec.getRTime()<rt-1.0){
      msr.readFile(NULL,spec);
    }
    if(spec.getScanNumber()>0) {
      scanBuf->push_back(spec);
      if(!params->ms1Centroid){
        centroid(spec,cent,params->ms1Resolution,params->instrument);
        centBuf->push_back(cent);
      } else {
        centBuf->push_back(spec);
      }
    }
  }

  //Sanity check. Exit if you can't find precursor scans within range of MS/MS retention time
  if(scanBuf->size()==0){
    //cout << "Reached end of file buffer and no spectra remain." << endl;
    return false;
  }

  //add spectra that need to be buffered
  lastScan=scanBuf->at(scanBuf->size()-1).getScanNumber();
  if(scanBuf->at(scanBuf->size()-1).getRTime()<rt+0.6){
    msr.setFilter(MS1);
    msr.readFile(NULL,spec,lastScan);
    msr.readFile(NULL,spec);
    while(spec.getScanNumber()>0 && spec.getRTime()<rt+2.0){
      scanBuf->push_back(spec);
      if(!params->ms1Centroid){
        centroid(spec,cent,params->ms1Resolution,params->instrument);
        centBuf->push_back(cent);
      } else {
        centBuf->push_back(spec);
      }
      msr.readFile(NULL,spec);
    }
    if(spec.getScanNumber()>0){
      scanBuf->push_back(spec);
      if(!params->ms1Centroid){
        centroid(spec,cent,params->ms1Resolution,params->instrument);
        centBuf->push_back(cent);
      } else {
        centBuf->push_back(spec);
      }
    }
  }

  //Find a precursor within 10 seconds for a given scan 
  v.clear();
  maxIntensity=0;
  for(i=0;i<scanBuf->size();i++){
    if(scanBuf->at(i).getRTime()<rt-0.167) continue;
    if(scanBuf->at(i).getRTime()>rt+0.167) break;
    j=findPeak(centBuf->at(i),mz,10);
    if(j>-1){
      p.rt=scanBuf->at(i).getRTime();
      p.scan=scanBuf->at(i).getScanNumber();
      v.push_back(p);
      if(centBuf->at(i)[j].intensity>maxIntensity){
        maxIntensity=centBuf->at(i)[j].intensity;
        maxRT=p.rt;
        best=p.scan;
        precursor=i;
      }
    }
  }

  //If no precursor was found within 10 seconds of MS/MS scan event, exit.
  //In the future, perhaps add verbosity parameter to capture this message.
  if(v.size()<1){
    //cout << "Warning: Precursor not found for " << scanNum << " " << mz << endl;
    return false;
  }
  
  //Get up to +/-15 sec of spectra around the max precursor intensity
  //This is done by extending on both sides until a gap is found or time is reached.
  vs.clear();
  rtHigh=rtLow=scanBuf->at(precursor).getRTime();
  for(i=precursor;i<scanBuf->size();i++){
    if(scanBuf->at(i).getRTime()>maxRT+0.25) break;
    j=findPeak(centBuf->at(i),mz,10);
    if(j<0) break;
    vs.push_back(&scanBuf->at(i));
    rtHigh=scanBuf->at(i).getRTime();
  }

  for(k=precursor-1;k>-1;k--){
    if(scanBuf->at(k).getRTime()<maxRT-0.25) break;
    j=findPeak(centBuf->at(k),mz,10);
    if(j<0) break;
    vs.push_back(&scanBuf->at(k));
    rtLow=scanBuf->at(k).getRTime();
  }

  //If total spectra to be combined is small (5 scan events), try 
  //grabbing a +/- 5 sec window around max, regardless of precursor observation
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

  //Average points between mz-1.5 and mz+2
  if(params->ms1Centroid) averageScansCentroid(vs,s,mz-1.5,mz+2);
  else averageScansBin(vs,s,mz-1.5,mz+2);
  s.setScanNumber(scanBuf->at(precursor).getScanNumber());

  //Clear corr
  corr=0;

  //Perform 18O2 analysis with Hardklor. If enrichment is set to 0, store unenriched results in the *O2 variables.
  //This is done in a non-competitive to identify an 18O2 peptide without solving everything
  if(params->enrichment>0) h->GoHardklor(hs2,&s);
  else h->GoHardklor(hs,&s);

  for(j=0;j<h->Size();j++){
    if(h->operator[](j).corr>corr){
      monoMass=h->operator[](j).monoMass;
      charge=h->operator[](j).charge;
      corr=h->operator[](j).corr;
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

  //Clear corr
  corr=0;

  //Perform the 18O4 Hardklor analysis on the same data if needed
  if(params->enrichment>0) {
    h->GoHardklor(hs4,&s);
    for(j=0;j<h->Size();j++){
      if(h->operator[](j).corr>corr){
        monoMass=h->operator[](j).monoMass;
        charge=h->operator[](j).charge;
        corr=h->operator[](j).corr;
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

  //If we have something, return it
  if(pls.sizePrecursor()>0)  return true;
  else return false;

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

  centBuf->clear();
  scanBuf->clear();
  scanBuf->push_back(spec);
  if(!params->ms1Centroid){
    centroid(spec,cent,params->ms1Resolution,params->instrument);
    centBuf->push_back(cent);
  } else {
    centBuf->push_back(spec);
  }

  return true;
}

/*============================
  Private Functions
============================*/
void KPrecursor::averageScansBin(vector<Spectrum*>& s, Spectrum& avg, double min, double max){

  unsigned int i;
  int j;
  double binWidth;
  double offset;
  float* bin;
  float intensity;
  int binCount;
  int* pos;

  pos=new int[s.size()];

  avg.clear();

  //Get some point near peak of interest to determine bin size. Make sure it is not a 0 intensity point.
  pos[0]=findPeak(*s[0],(min+max)/2);
  while(pos[0]<s[0]->size()-2 && (s[0]->at(pos[0]).intensity<0.1 || s[0]->at(pos[0]-1).intensity<0.1 || s[0]->at(pos[0]+1).intensity<0.1)){
    pos[0]++;
  }
  
  //binWidth is the distance between two points in profile data divided by 2
  binWidth=s[0]->at(pos[0]+1).mz-s[0]->at(pos[0]).mz;
  binWidth/=2;
  
  binCount=(int)((max-min)/binWidth+1);
  bin = new float[binCount];
  for(j=0;j<binCount;j++) bin[j]=0;

  //align all spectra to point closest to min and set offset
  for(i=0;i<s.size();i++)  {
    pos[i]=findPeak(*s[i],min);
    if(s[i]->at(pos[i]).mz<min) pos[i]++;
  }
  offset=s[0]->at(pos[0]).mz;

  //Iterate all spectra, adding intensities to the bins
  for(i=0;i<s.size();i++) {
    while(pos[i]<s[i]->size() && s[i]->at(pos[i]).mz<max){
      j=(int)((s[i]->at(pos[i]).mz-offset)/binWidth);
      if(j<0){
        pos[i]++;
        continue;
      }
      if(j>=binCount) break;
      bin[j]+=s[i]->at(pos[i]).intensity;
      pos[i]++;
    }
  }

  //average scan is sum of two bins at a time
  for(j=1;j<binCount-1;j+=2){
    intensity=bin[j]+bin[j+1];
    intensity/=s.size();
    if(intensity>0) avg.add(offset+j*binWidth+binWidth/2,intensity);
  }

  //clear memory
  delete [] bin;
  delete [] pos;
}

//Almost identical to averageScansBin, except centroided spectra need to be treated differently
void KPrecursor::averageScansCentroid(vector<Spectrum*>& s, Spectrum& avg, double min, double max){

  unsigned int i;
  int j;
  double binWidth;
  double offset;
  float* bin;
  float intensity;
  int binCount;
  int* pos;

  pos=new int[s.size()];

  avg.clear();

  //Set some really small bin width related to the mz region being summed.
  binWidth=(min+max)/2*0.00001;
  binWidth/=2;
  
  binCount=(int)((max-min)/binWidth+1);
  bin = new float[binCount];
  for(j=0;j<binCount;j++) bin[j]=0;

  //align all spectra to point closest to min and set offset
  for(i=0;i<s.size();i++)  {
    pos[i]=findPeak(*s[i],min);
    if(s[i]->at(pos[i]).mz<min) pos[i]++;
  }
  offset=s[0]->at(pos[0]).mz;

  //Iterate all spectra and add peaks to bins
  for(i=0;i<s.size();i++) {
    while(pos[i]<s[i]->size() && s[i]->at(pos[i]).mz<max){
      j=(int)((s[i]->at(pos[i]).mz-offset)/binWidth);
      if(j<0){
        pos[i]++;
        continue;
      }
      if(j>=binCount) break;
      bin[j]+=s[i]->at(pos[i]).intensity;
      pos[i]++;
    }
  }

  //average scan is sum of two bins at a time
  for(j=1;j<binCount-1;j+=2){
    intensity=bin[j]+bin[j+1];
    intensity/=s.size();
    if(intensity>0) avg.add(offset+j*binWidth+binWidth/2,intensity);
  }

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

