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
void KParams::parse(char* cmd) {

  char *tok;

  char param[32];
  char tmpstr[256];
  char value[256];
  char value2[256];

  kMass m;

  string tstr;
  vector<string> vs;

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
    strcpy(value,tok);
  }
  
  //get second parameter for modification
	if(strcmp(param,"modification")==0 || strcmp(param,"fixed_modification")==0) {
		tok=strtok(NULL," \t=\n\r");
    if(tok==NULL) {
		  warn(param,0);
		  return;
    } else {
      strcpy(value2,tok);
    }
  }

	//Look up parameter and assign the value
	if(strcmp(param,"cross_link_mass")==0){
    params->xLink->push_back(atof(value));

	} else if(strcmp(param,"database")==0){
    strcpy(params->dbFile,value);

  } else if(strcmp(param,"diagnostic")==0){
    params->diagnostic=atoi(value);

	} else if(strcmp(param,"decoy_filter")==0){
    strcpy(params->decoy,value);
 
  } else if(strcmp(param,"enrichment")==0){
    params->enrichment=atof(value);

  } else if(strcmp(param,"fixed_modification")==0){
    m.index=(int)value[0];
    m.mass=atof(value2);
    params->fMods->push_back(m);

	} else if(strcmp(param,"instrument")==0){
    params->instrument=atoi(value);
    if(params->instrument<0 || params->instrument>1){
      warn("Value out of range for instrument. Defaulting to 1=Orbitrap",2);
      params->instrument=1;
    }

	} else if(strcmp(param,"max_miscleavages")==0){
    params->miscleave=atoi(value);

  } else if(strcmp(param,"max_mods_per_peptide")==0){
    params->maxMods=atoi(value);

	} else if(strcmp(param,"max_peptide_mass")==0){
    params->maxPepMass=atof(value);

  } else if(strcmp(param,"max_spectrum_peaks")==0){
    params->maxPeaks=atoi(value);

	} else if(strcmp(param,"min_peptide_mass")==0){
    params->minPepMass=atof(value);

  } else if(strcmp(param,"modification")==0){
    m.index=(int)value[0];
    m.mass=atof(value2);
    params->mods->push_back(m);

	} else if(strcmp(param,"mono_link_mass")==0){
		params->mLink->push_back(atof(value));

	} else if(strcmp(param,"MS_data_file")==0){
    strcpy(params->msFile,value);

  } else if(strcmp(param,"MS1_centroid")==0){
    params->ms1Centroid=atoi(value);

  } else if(strcmp(param,"MS2_centroid")==0){
    params->ms2Centroid=atoi(value);

	} else if(strcmp(param,"MS1_resolution")==0){
    params->ms1Resolution=atoi(value);

	} else if(strcmp(param,"MS2_resolution")==0){
		params->ms2Resolution=atoi(value);

	} else if(strcmp(param,"output_file")==0){
    strcpy(params->outFile,value);

	} else if(strcmp(param,"percolator_file")==0){
    strcpy(params->percolator,value);

  } else if(strcmp(param,"percolator_version")==0){
    params->percVersion=atof(value);

	} else if(strcmp(param,"ppm_tolerance_frag")==0){
    params->ppmFragment=atof(value);

	} else if(strcmp(param,"ppm_tolerance_pre")==0){
    params->ppmPrecursor=atof(value);

  } else if(strcmp(param,"prefer_precursor_pred")==0){
    params->preferPrecursor=atoi(value);

	} else if(strcmp(param,"relaxed_analysis")==0){
    params->relaxedAnalysis=atoi(value);

  } else if(strcmp(param,"search_dimers")==0){
    params->dimers=atoi(value);

  } else if(strcmp(param,"spectrum_processing")==0) {
    params->specProcess=atoi(value);

  } else if(strcmp(param,"top_count")==0) {
    params->topCount=atoi(value);

  } else if(strcmp(param,"truncate_prot_names")==0) {
    params->truncate=atoi(value);

  } else if(strcmp(param,"use_comet_xcorr")==0){
    params->xcorr=atoi(value);

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
		default:
			printf("%s\n",c);
			break;
	}
}
