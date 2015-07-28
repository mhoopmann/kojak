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

#include "KParams.h"

//==============================
//  Constructors & Destructors
//==============================
KParams::KParams(){
  params = NULL;
}

KParams::KParams(kParams* p){
  params = p;
}

KParams::~KParams(){
  params = NULL;
}

//==============================
//  Public Functions
//==============================
bool KParams::parseConfig(char* fname){
  FILE* f;
  char str[512];

  if(params==NULL){
    printf("kParams structure not set in KParams object. No parameters read.\n");
    return false;
  }

  cout << "Parsing: " << fname << endl;
  f=fopen(fname,"rt");
  if(f==NULL){
    printf("Cannot open config file!\n");
    return false;
  }

  while(!feof(f)) {
    if(!fgets(str,512,f)) continue;
    if(str[0]=='#') continue;
    parse(str);
  }

  fclose(f);
  return true;
}

//==============================
//  Private Functions
//==============================
bool KParams::checkMod(kMass m){
  unsigned int i;
  for(i=0;i<params->mods->size();i++){
    if(params->mods->at(i).index == m.index && params->mods->at(i).mass == m.mass && params->mods->at(i).xl == m.xl) return true;
  }
  return false;
}

void KParams::parse(char* cmd) {

  char *tok;

  char param[32];
  char tmpstr[256];

  int i;

  kMass m;
  kLinker x;

  string tstr;
  vector<string> values;

  //Pre-process entire line to remove characters that should not be read
	//Replace first # with a terminator
	tok=strstr(cmd,"#");
	if(tok!=NULL) strncpy(tok,"\0",1);

	//if we have only white space, exit here
	strcpy(tmpstr,cmd);
	tok=strtok(tmpstr," \t\n\r");
	if(tok==NULL) return;

	//Check if we have a parameter (has '=' in it) or lots of random text.
	tok=strstr(cmd,"=");
  if(tok==NULL) {
    printf("Unknown parameter line in config file: %s\n",cmd);
    return;
  }

  //Process parameters
	//Read parameter into param name (before = sign) and value (after = sign)
	tok=strtok(cmd," \t=\n\r");
	if(tok==NULL) return;
	strcpy(param,tok);
	tok=strtok(NULL," \t=\n\r");
	if(tok==NULL) {
		warn(param,0);
		return;
  } else {
    while(tok!=NULL){
      tstr=tok;
      values.push_back(tstr);
      tok=strtok(NULL," \t=\n\r");
    }
  }

	//Look up parameter and assign the value
  if(strcmp(param,"cross_link")==0){
    //check first site
    if(values.size()!=3){
      printf("Error in cross_link parameter(s)\n");
      exit(-5);
    }
    i=atoi(&values[0][0]);
    if(i<1 || i>5) {
      printf("Error in cross_link parameter(s)\n");
      exit(-5);
    }
    if(params->setA==0) {
      params->setA=i;
      x.siteA=1;
    } else if(params->setA==i) {
      x.siteA=1;
    } else if(params->setB==0) {
      params->setB=i;
      x.siteA=2;
    } else if(params->setB==i) {
      x.siteA=2;
    } else {
      printf("Error in cross_link parameter(s)\n");
      exit(-5);
    }
    //check second site
    i=atoi(&values[1][0]);
    if(i<1 || i>5) {
      printf("Error in cross_link parameter(s)\n");
      exit(-5);
    }
    if(params->setA==0) {
      params->setA=i;
      x.siteB=1;
    } else if(params->setA==i) {
      x.siteB=1;
    } else if(params->setB==0) {
      params->setB=i;
      x.siteB=2;
    } else if(params->setB==i) {
      x.siteB=2;
    } else {
      printf("Error in cross_link parameter(s)\n");
      exit(-5);
    }
    x.mass=atof(&values[2][0]);
    x.mono=0;
    params->xLink->push_back(x);

	} else if(strcmp(param,"database")==0){
    strcpy(params->dbFile,&values[0][0]);

  } else if(strcmp(param,"diagnostic")==0){
    //params->diagnostic=atoi(&values[0][0]);
    params->diag->push_back(atoi(&values[0][0]));

  } else if(strcmp(param,"diff_mods_on_xl")==0){
    if(atoi(&values[0][0])!=0) params->diffModsOnXL=true;
    else params->diffModsOnXL=false;

	} else if(strcmp(param,"decoy_filter")==0){
    strcpy(params->decoy,&values[0][0]);
 
  } else if(strcmp(param,"enrichment")==0){
    params->enrichment=atof(&values[0][0]);

  } else if(strcmp(param,"enzyme")==0){
    strcpy(params->enzyme,&values[0][0]);

  } else if(strcmp(param,"export_percolator")==0){
    if(atoi(&values[0][0])!=0) params->exportPercolator=true;
    else params->exportPercolator=false;

  } else if(strcmp(param,"fixed_modification")==0){
    m.index=(int)values[0][0];
    m.mass=atof(&values[1][0]);
    params->fMods->push_back(m);

	} else if(strcmp(param,"fragment_bin_offset")==0){
    params->binOffset=1.0-atof(&values[0][0]);

	} else if(strcmp(param,"fragment_bin_size")==0){
    params->binSize=atof(&values[0][0]);
    if(params->binSize<=0){
      printf("Invalid value for fragment_bin_size parameter\n");
      exit(-5);
    }

	} else if(strcmp(param,"instrument")==0){
    params->instrument=atoi(&values[0][0]);
    if(params->instrument<0 || params->instrument>1){
      warn("Value out of range for instrument. Defaulting to 1=Orbitrap",3);
      params->instrument=1;
    }

  } else if(strcmp(param,"ion_series_A")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[0]=false;
    else params->ionSeries[0]=true;

  } else if(strcmp(param,"ion_series_B")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[1]=false;
    else params->ionSeries[1]=true;

  } else if(strcmp(param,"ion_series_C")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[2]=false;
    else params->ionSeries[2]=true;

  } else if(strcmp(param,"ion_series_X")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[3]=false;
    else params->ionSeries[3]=true;

  } else if(strcmp(param,"ion_series_Y")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[4]=false;
    else params->ionSeries[4]=true;

  } else if(strcmp(param,"ion_series_Z")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[5]=false;
    else params->ionSeries[5]=true;

	} else if(strcmp(param,"max_miscleavages")==0){
    params->miscleave=atoi(&values[0][0]);

  } else if(strcmp(param,"max_mods_per_peptide")==0){
    params->maxMods=atoi(&values[0][0]);

	} else if(strcmp(param,"max_peptide_mass")==0){
    params->maxPepMass=atof(&values[0][0]);

  } else if(strcmp(param,"max_spectrum_peaks")==0){
    params->maxPeaks=atoi(&values[0][0]);

	} else if(strcmp(param,"min_peptide_mass")==0){
    params->minPepMass=atof(&values[0][0]);

  } else if(strcmp(param,"modification")==0){
    m.xl=false;
    m.index=(int)values[0][0];
    m.mass=atof(&values[1][0]);
    if(!checkMod(m)) params->mods->push_back(m);

  } else if(strcmp(param,"mono_link")==0){
    m.xl=true;
    m.mass=atof(&values[1][0]);
    i=atoi(&values[0][0]);
    if(i==1){
      m.index=75;
      if(!checkMod(m)) params->mods->push_back(m);
      m.index=64;
      if(!checkMod(m)) params->mods->push_back(m);
    } else if(i==2){
      m.index=68;
      if(!checkMod(m)) params->mods->push_back(m);
      m.index=69;
      if(!checkMod(m)) params->mods->push_back(m);
      m.index=36;
      if(!checkMod(m)) params->mods->push_back(m);
    } else if(i==3){
      m.index=67;
      if(!checkMod(m)) params->mods->push_back(m);
    } else if(i==4){
      m.index=81;
      if(!checkMod(m)) params->mods->push_back(m);
    } else if(i==5){
      m.index=75;
      if(!checkMod(m)) params->mods->push_back(m);
      m.index=64;
      if(!checkMod(m)) params->mods->push_back(m);
      m.index=83;
      if(!checkMod(m)) params->mods->push_back(m);
      m.index=84;
      if(!checkMod(m)) params->mods->push_back(m);
      m.index=89;
      if(!checkMod(m)) params->mods->push_back(m);
    }

  } else if(strcmp(param,"mono_links_on_xl")==0){
    if(atoi(&values[0][0])!=0) params->monoLinksOnXL=true;
    else params->monoLinksOnXL=false;

	} else if(strcmp(param,"MS_data_file")==0){
    strcpy(params->msFile,&values[0][0]);

  } else if(strcmp(param,"MS1_centroid")==0){
    params->ms1Centroid=atoi(&values[0][0]);

  } else if(strcmp(param,"MS2_centroid")==0){
    params->ms2Centroid=atoi(&values[0][0]);

	} else if(strcmp(param,"MS1_resolution")==0){
    params->ms1Resolution=atoi(&values[0][0]);

	} else if(strcmp(param,"MS2_resolution")==0){
		params->ms2Resolution=atoi(&values[0][0]);

	} else if(strcmp(param,"output_file")==0){
    //strcpy(params->outFile,&values[0][0]);
    warn(param,2);

	} else if(strcmp(param,"percolator_file")==0){
    //strcpy(params->percolator,&values[0][0]);
    warn(param,2);

  } else if(strcmp(param,"percolator_version")==0){
    params->percVersion=atof(&values[0][0]);

	} else if(strcmp(param,"ppm_tolerance_pre")==0){
    params->ppmPrecursor=atof(&values[0][0]);

  } else if(strcmp(param,"prefer_precursor_pred")==0){
    params->preferPrecursor=atoi(&values[0][0]);

  } else if(strcmp(param,"search_dimers")==0){
    params->dimers=atoi(&values[0][0]);

  } else if(strcmp(param,"spectrum_processing")==0) {
    params->specProcess=atoi(&values[0][0]);

  } else if(strcmp(param,"threads")==0) {
    params->threads=atoi(&values[0][0]);
    if(params->threads<1) {
      printf("Invalid threads parameter. Setting to 1\n");
      params->threads=1;
    }

  } else if(strcmp(param,"top_count")==0) {
    params->topCount=atoi(&values[0][0]);

  } else if(strcmp(param,"truncate_prot_names")==0) {
    params->truncate=atoi(&values[0][0]);

  } else if(strcmp(param,"use_comet_xcorr")==0){
    if(atoi(&values[0][0])!=0) params->xcorr=true;
    else params->xcorr=false;

	} else {
		warn(param,1);
	}

}

void KParams::warn(char* c, int i){
	switch(i){
		case 0:
			printf("Parameter %s has no value.",c);
			break;
		case 1:
			printf("Unknown parameter: %s\n",c);
			break;
		case 2:
      printf("Parameter %s has been deprecated and will be ignored.\n",c);
      break;
    case 3:
		default:
			printf("%s\n",c);
			break;
	}
}
