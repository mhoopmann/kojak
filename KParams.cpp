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

  pxwBasicXMLTag xml;

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
    //Check number of parameters
    if (values.size() != 4){
      warn("ERROR: bad cross_link parameter(s).",3);
      exit(-5);
    }
    i = atoi(&values[0][0]);
    if (i > 0) {
      warn("ERROR: bad cross_link parameter(s). Suspected use of deprecated format.",3);
      exit(-5);
    }
    x.motifA = values[0];
    i = atoi(&values[1][0]);
    if (i > 0) {
      warn("ERROR: bad cross_link parameter(s). Suspected use of deprecated format.",3);
      exit(-5);
    }
    x.motifB = values[1];
    x.mass = atof(&values[2][0]);
    x.mono = 0;
    x.label = values[3];
    params->xLink->push_back(x);
    xml.name="cross_link";
    xml.value=values[0]+" "+values[1]+" "+values[2]+" "+values[3];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"database")==0){
    strcpy(params->dbFile,&values[0][0]);
    xml.name = "database";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"diagnostic")==0){
    params->diag->push_back(atoi(&values[0][0]));
    xml.name = "diagnostic";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"diff_mods_on_xl")==0){
    if(atoi(&values[0][0])!=0) params->diffModsOnXL=true;
    else params->diffModsOnXL=false;
    xml.name = "diff_mods_on_xl";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"decoy_filter")==0){
    strcpy(params->decoy,&values[0][0]);
    xml.name = "decoy_filter";
    xml.value = values[0];
    xmlParams.push_back(xml);
 
  } else if(strcmp(param,"enrichment")==0){
    params->enrichment=atof(&values[0][0]);
    xml.name = "enrichment";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"enzyme")==0){
    strcpy(params->enzyme,&values[0][0]);
    xml.name = "enzyme";
    xml.value = values[0];
    if (values.size() > 1) strcpy(params->enzymeName, &values[1][0]);
    else strcpy(params->enzymeName, "Unnamed");
    xmlParams.push_back(xml);

  } else if(strcmp(param, "export_pepXML")==0 || strcmp(param, "export_pepxml")==0){
    if(atoi(&values[0][0])!=0) params->exportPepXML=true;
    else params->exportPepXML=false;
    xml.name = "export_pepXML";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"export_percolator")==0){
    if(atoi(&values[0][0])!=0) params->exportPercolator=true;
    else params->exportPercolator=false;
    xml.name = "export_percolator";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"fixed_modification")==0){
    m.index=(int)values[0][0];
    m.mass=atof(&values[1][0]);
    if(m.mass!=0) params->fMods->push_back(m);
    xml.name = "fixed_modification";
    xml.value = values[0]+" "+values[1];
    xmlParams.push_back(xml);

  } else if (strcmp(param, "fixed_modification_protC") == 0){
    m.index = (int)'%';
    m.mass = atof(&values[0][0]);
    if (m.mass != 0) params->fMods->push_back(m);
    xml.name = "fixed_modification_protC";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if (strcmp(param, "fixed_modification_protN") == 0){
    m.index = (int)'$';
    m.mass = atof(&values[0][0]);
    if (m.mass != 0) params->fMods->push_back(m);
    xml.name = "fixed_modification_protN";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"fragment_bin_offset")==0){
    params->binOffset=1.0-atof(&values[0][0]);
    xml.name = "fragment_bin_offset";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"fragment_bin_size")==0){
    params->binSize=atof(&values[0][0]);
    if(params->binSize<=0){
      warn("ERROR: Invalid value for fragment_bin_size parameter. Stopping analysis.",3);
      exit(-5);
    }
    xml.name = "fragment_bin_size";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"instrument")==0){
    params->instrument=atoi(&values[0][0]);
    if(params->instrument<0 || params->instrument>1){
      warn("ERROR: instrument value out of range for instrument. Stopping analysis.",3);
      exit(-5);
    }
    xml.name = "instrument";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if (strcmp(param, "intermediate") == 0){
    params->intermediate = atoi(&values[0][0]);
    if (params->intermediate < 0){
      warn("ERROR: intermediate value out of range. Stopping analysis.", 3);
      exit(-5);
    }

  } else if(strcmp(param,"ion_series_A")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[0]=false;
    else params->ionSeries[0]=true;
    xml.name = "ion_series_A";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"ion_series_B")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[1]=false;
    else params->ionSeries[1]=true;
    xml.name = "ion_series_B";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"ion_series_C")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[2]=false;
    else params->ionSeries[2]=true;
    xml.name = "ion_series_C";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"ion_series_X")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[3]=false;
    else params->ionSeries[3]=true;
    xml.name = "ion_series_X";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"ion_series_Y")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[4]=false;
    else params->ionSeries[4]=true;
    xml.name = "ion_series_Y";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"ion_series_Z")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[5]=false;
    else params->ionSeries[5]=true;
    xml.name = "ion_series_Z";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if (strcmp(param, "isotope_error")==0){
    params->isotopeError = atoi(&values[0][0]);
    if (params->isotopeError < 0 || params->isotopeError>2){
      warn("ERROR: isotope_error has invalid value. Stopping analysis.",3);
      exit(-5);
    }
    xml.name = "isotope_error";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"max_miscleavages")==0){
    params->miscleave=atoi(&values[0][0]);
    xml.name = "max_miscleavages";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"max_mods_per_peptide")==0){
    params->maxMods=atoi(&values[0][0]);
    xml.name = "max_mods_per_peptide";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"max_peptide_mass")==0){
    params->maxPepMass=atof(&values[0][0]);
    xml.name = "max_peptide_mass";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"max_spectrum_peaks")==0){
    params->maxPeaks=atoi(&values[0][0]);
    xml.name = "max_spectrum_peaks";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"min_peptide_mass")==0){
    params->minPepMass=atof(&values[0][0]);
    xml.name = "min_peptide_mass";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"modification")==0){
    m.xl=false;
    m.index=(int)values[0][0];
    m.mass=atof(&values[1][0]);
    if(m.mass!=0 && !checkMod(m)) params->mods->push_back(m);
    xml.name = "modification";
    xml.value = values[0]+" "+values[1];
    xmlParams.push_back(xml);

  } else if (strcmp(param, "modification_protC") == 0){
    m.xl = false;
    m.index = (int)'%';
    m.mass = atof(&values[0][0]);
    if (m.mass!=0) { //acceptable to use placeholder; don't load it as a parameter
      if (!checkMod(m)) params->mods->push_back(m);
    }
    xml.name = "modification_protC";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if (strcmp(param, "modification_protN") == 0){
    m.xl = false;
    m.index = (int)'$';
    m.mass = atof(&values[0][0]);
    if(m.mass != 0) { //acceptable to use placeholder; don't load it as a parameter
      if (!checkMod(m)) params->mods->push_back(m);
    }
    xml.name = "modification_protN";
    xml.value = values[0];
    xmlParams.push_back(xml);
    
  } else if(strcmp(param,"mono_link")==0){
    //Check number of parameters
    if (values.size() != 2){
      warn("ERROR: bad mono_link parameter(s)",3);
      exit(-5);
    }
    i = atoi(&values[0][0]);
    if (i > 0) {
      warn("ERROR: bad mono_link parameter(s). Suspected use of deprecated format.",3);
      exit(-5);
    }
    m.xl = true;
    m.mass = atof(&values[1][0]);
    for (i = 0; i < values[0].size(); i++){
      if (values[0][i] == 'c') m.index = (int)'%';
      else if (values[0][i] == 'n') m.index = (int)'$';
      else m.index = (int)values[0][i];
      if (!checkMod(m)) params->mods->push_back(m);
    }
    xml.name = "mono_link";
    xml.value = values[0] + " " + values[1];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"mono_links_on_xl")==0){
    if(atoi(&values[0][0])!=0) params->monoLinksOnXL=true;
    else params->monoLinksOnXL=false;
    xml.name = "mono_links_on_xl";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"MS_data_file")==0){
    strcpy(params->msFile,&values[0][0]);
    xml.name = "MS_data_file";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"MS1_centroid")==0){
    params->ms1Centroid=atoi(&values[0][0]);
    xml.name = "MS1_centroid";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"MS2_centroid")==0){
    params->ms2Centroid=atoi(&values[0][0]);
    xml.name = "MS2_centroid";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"MS1_resolution")==0){
    params->ms1Resolution=atoi(&values[0][0]);
    xml.name = "MS1_resolution";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"MS2_resolution")==0){
		params->ms2Resolution=atoi(&values[0][0]);
    xml.name = "MS2_resolution";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"output_file")==0){
    warn(param,2);

	} else if(strcmp(param,"percolator_file")==0){
    warn(param,2);

  } else if(strcmp(param,"percolator_version")==0){
    params->percVersion=atof(&values[0][0]);
    xml.name = "percolator_version";
    xml.value = values[0];
    xmlParams.push_back(xml);

	} else if(strcmp(param,"ppm_tolerance_pre")==0){
    params->ppmPrecursor=atof(&values[0][0]);
    xml.name = "ppm_tolerance_pre";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if (strcmp(param, "precursor_refinement") == 0){
    if (atoi(&values[0][0]) == 0) params->precursorRefinement = false;
    else params->precursorRefinement = true;
    xml.name = "precursor_refinement";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"prefer_precursor_pred")==0){
    params->preferPrecursor=atoi(&values[0][0]);
    xml.name = "prefer_precursor_pred";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"search_dimers")==0){
    warn(param,2);

  } else if(strcmp(param,"search_dimers_xl")==0){
    if(atoi(&values[0][0])==0) params->dimersXL=false;
    else params->dimersXL=true;
    xml.name = "search_dimers_xl";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"spectrum_processing")==0) {
    params->specProcess=atoi(&values[0][0]);
    xml.name = "spectrum_processing";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"threads")==0) {
    params->threads=atoi(&values[0][0]);
    if(params->threads<1) {
      warn("WARNING: Invalid threads parameter. Setting to 1.",3);
      params->threads=1;
    }
    xml.name = "threads";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"top_count")==0) {
    params->topCount=atoi(&values[0][0]);
    xml.name = "top_count";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"truncate_prot_names")==0) {
    params->truncate=atoi(&values[0][0]);
    xml.name = "truncate_prot_names";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if (strcmp(param, "turbo_button") == 0){
    if (atoi(&values[0][0]) != 0) params->turbo = true;
    else params->turbo = false;
    xml.name = "turbo_button";
    xml.value = values[0];
    xmlParams.push_back(xml);

  } else if(strcmp(param,"use_comet_xcorr")==0){
    warn(param, 2);
    //if(atoi(&values[0][0])!=0) params->xcorr=true;
    //else params->xcorr=false;
    //xml.name = "use_comet_xcorr";
    //xml.value = values[0];
    //xmlParams.push_back(xml);

	} else {
		warn(param,1);
	}

}

void KParams::warn(const char* c, int i){
	switch(i){
		case 0:
			printf("  WARNING: Parameter %s has no value.",c);
			break;
		case 1:
			printf("  WARNING: Unknown parameter: %s\n",c);
			break;
		case 2:
      printf("  WARNING: Parameter %s has been deprecated and will be ignored.\n",c);
      break;
    case 3:
		default:
			printf("  %s\n",c);
			break;
	}
}
