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

using namespace std;

//==============================
//  Constructors & Destructors
//==============================
KParams::KParams(){
  log=NULL;
  params = NULL;
}

KParams::KParams(kParams* p){
  params = p;
}

KParams::~KParams(){
  log=NULL;
  params = NULL;
}

//==============================
//  Public Functions
//==============================
bool KParams::parseConfig(const char* fname){
  FILE* f;
  char str[4096];

  if(params==NULL){
    printf("kParams structure not set in KParams object. No parameters read.\n");
    return false;
  }

  f=fopen(fname,"rt");
  if(f==NULL){
    if (log != NULL) log->addError("Cannot open config file!");
    else printf("Cannot open config file!\n");
    return false;
  }

  vector<string> vParamLines;
  bool bReadXLTable=false;;
  while(!feof(f)) {
    if(!fgets(str,4096,f)) continue;
    if(strlen(str)==0) continue;
    if(str[0]=='#') continue;
    //if(strlen(str)>4096) {
    //  if(log!=NULL) log->addError("Parameter line too long!");
    //  else printf("Parameter line too long!\n");
    //  return false;
    //}
    string s=str;

    if(s.find("[XL_PARAMS]")==0){
      bReadXLTable=true;
      continue;
    }
    if (s.find("[END_XL_PARAMS]") == 0){
      bReadXLTable = false;
      continue;
    }

    if(bReadXLTable){
      char* tok=strtok(str," \t\n\r");
      kPXLTable x;
      x.id=atoi(tok);
      size_t a;
      for(a=0;a<pXLTable.size();a++) if(pXLTable[a].id==x.id) break;
      if(a<pXLTable.size()){
        if (log != NULL) log->addError("Duplicate ID number in [XL_PARAMS]");
        else printf("Duplicate ID number in [XL_PARAMS]\n");
        return false;
      }

      tok = strtok(NULL, " \t\n\r");
      if (checkToken(tok)) x.xlName = tok;
      else return false;

      tok = strtok(NULL, " \t\n\r");
      if (checkToken(tok)) x.xlQuench = tok;
      else return false;

      tok = strtok(NULL, " \t\n\r");
      if(checkToken(tok)) x.targetA=tok;
      else return false;

      tok = strtok(NULL, " \t\n\r");
      if (checkToken(tok)) x.targetB = tok;
      else return false;

      tok = strtok(NULL, " \t\n\r");
      if (checkToken(tok)) x.linkMass = atof(tok);
      else return false;

      tok = strtok(NULL, " \t\n\r");
      if (checkToken(tok)) splitMasses(tok,x.monoA);
      else return false;

      tok = strtok(NULL, " \t\n\r");
      if (checkToken(tok)) splitMasses(tok, x.monoB);
      else return false;

      tok = strtok(NULL, " \t\n\r");
      if (checkToken(tok)) splitMasses(tok, x.cleavageMass);
      else return false;

      pXLTable.push_back(x);
      continue;
    }

    vParamLines.push_back(s);
    parse(str);
  }

  fclose(f);

  cout << pXLTable.size() << " preconfigured crosslinkers." << endl;
  return true;
}

//==============================
//  Private Functions
//==============================
bool KParams::buildOutput(char* in, char* base, char* ext){
  char cwd[1024];
  char str[1056];
  char outPath[1024];
  char* ret;
  string tmp;
  string outFile;
  size_t i;
  FILE* f;

  //get current working directory and process input file
  strcpy(params->msFile, in);
  strcpy(params->ext, ext);
  ret=getcwd(cwd, 1024);
  processPath(cwd, in, str);
  strcpy(params->inFile, str);
  
  //process output file
  //if user did not specify output path, use cwd as full path for output
  if (strlen(params->resPath)==0){
    processPath(cwd, base, outPath);  
  } else {
    if(params->resPath[0]=='.'){
      processPath(cwd, params->resPath, outPath);
      tmp = outPath;
    } else {
      tmp = params->resPath;
    }
    params->fullPath=tmp;
    outFile = base;
    i = outFile.find_last_of("/\\");
    if (i != string::npos) outFile = outFile.substr(i + 1);
    processPath(&tmp[0], &outFile[0], outPath);
  }
  strcpy(params->outFile, outPath);

  /*
  cout << "\nMS_data_file: " << params->msFile << endl;
  cout << "results_path: " << params->resPath << endl;
  cout << "input:        " << params->inFile << endl;
  cout << "output:       " << params->outFile << endl;
  */

  //Create Kojak output file to test output paths
  sprintf(str, "%s.kojak.txt", params->outFile);
  f = fopen(str, "wt");
  if (f == NULL) {
    if (log != NULL) log->addError("\nERROR: cannot open " + string(str) + " for output. Please ensure path and write permissions are correct.");
    else cout << "\nERROR: cannot open " << str << " for output. Please ensure path and write permissions are correct." << endl;
    return false;
  } else {
    fclose(f);
  }

  //Create Kojak log file
  sprintf(str,"%s.log",params->outFile);
  f = fopen(str, "wt");
  if (f == NULL) {
    if (log != NULL) log->addError("\nERROR: cannot open " + string(str) + " to log Kojak. Please ensure path and write permissions are correct.");
    else cout << "\nERROR: cannot open " << str << " to log Kojak. Please ensure path and write permissions are correct." << endl;
    return false;
  } else {
    fclose(f);
  }
  logFile=str;

  return true;
}

bool KParams::checkMod(kMass m){
  unsigned int i;
  for(i=0;i<params->mods->size();i++){
    if(params->mods->at(i).index == m.index && params->mods->at(i).mass == m.mass && params->mods->at(i).xl == m.xl) return true;
  }
  return false;
}

void KParams::exportDefault(string ver){
  kParams def;
  FILE* f = fopen("kojak_default_params.conf", "wt");

  fprintf(f, "# Kojak %s parameter file\n\n", ver.c_str());
  fprintf(f, "# This file was autogenerated. Please alter and rename before use.\n");
  fprintf(f, "# All parameters are separated from their values by an equals sign ('=')\n");
  fprintf(f, "# Anything after a '#' will be ignored for the remainder of the line.\n");
  fprintf(f, "# Remove the '#' at the beginning of a parameter line to activate that parameter.\n");
  fprintf(f, "# Please see online documentation at:\n");
  fprintf(f, "# http://www.kojak-ms.org/param");

  fprintf(f, "\n\n#\n");
  fprintf(f, "# XL_PARAMS defines pre - configured crosslinker settings.The table can be expanded to suit new crosslinkers.\n");
  fprintf(f, "# This table MUST precede the predefined_crosslink parameter.\n");
  fprintf(f, "# Columns are : [ID Num] [Name] [Quenching Reagent] [Target A] [Target B] [XL Mass] [MonoMasses A] [MonoMasses B] [Cleavage_Product_Masses]\n");
  fprintf(f, "#\n");
  fprintf(f, "[XL_PARAMS]\n");
  fprintf(f, "1   BS3/DSS  NH2       nK  nK  138.068074  155.094629             x	 x\n");
  fprintf(f, "2   BS3/DSS  NH2+H2O   nK  nK  138.068074  155.094629,156.078644  x  x\n");
  fprintf(f, "3   BS3/DSS  Tris      nK  nK  138.068074  259.141973             x  x\n");
  fprintf(f, "4   BS3/DSS  Tris+H2O  nK  nK  138.068074  259.141973,156.078644  x  x\n");
  fprintf(f, "5   DSSO     NH2       nK  nK  158.003765  175.030314             x  54.010565,85.982635,103.993200\n");
  fprintf(f, "6   DSSO     NH2+H2O   nK  nK  158.003765  175.030314,176.014330  x  54.010565,85.982635,103.993200\n");
  fprintf(f, "7   DSSO     Tris      nK  nK  158.003765  279.077658             x  54.010565,85.982635,103.993200\n");
  fprintf(f, "8   DSSO     Tris+H2O  nK  nK  158.003765  279.077658,176.014330  x  54.010565,85.982635,103.993200\n");
  fprintf(f, "9   PhoX     NH2       nK  nK  209.972     226.998                x  x\n");
  fprintf(f, "10  PhoX     NH2+H2O   nK  nK  209.972     226.998,227.982        x  x\n");
  fprintf(f, "[END_XL_PARAMS]\n");

  fprintf(f, "\n\n#\n# Computational Settings\n#\n");
  fprintf(f, "threads = %d\n", def.threads);
  
  fprintf(f, "\n\n#\n# Input and Output files - specify full path for input files if not in current working directory\n#\n");
  fprintf(f, "MS_data_file = yourData.mzML            #users specify their data here.\n");
  fprintf(f, "database = SearchDatabase.fasta         #users specify their proteins here.\n");
  fprintf(f, "results_path = .                        #path must exist. Use '.' for current working directory.\n");
  fprintf(f, "export_pepXML = %d                       #0=no, 1=yes\n", (int)def.exportPepXML);
  fprintf(f, "export_percolator = %d                   #0=no, 1=yes\n", (int)def.exportPercolator);
  fprintf(f, "export_mzID = %d                         #0=no, 1=yes\n", (int)def.exportMzID);
  fprintf(f, "percolator_version = %.2lf\n", def.percVersion);
  fprintf(f, "split_pepXML = %d                        #0=no, 1=yes\n", (int)def.splitPepXML);
  fprintf(f, "truncate_prot_names = %d                 #shorten protein names in output to this number of characters, 0 = off\n", def.truncate);
  
  fprintf(f, "\n\n#\n# Parameters describing data analyzed by Kojak\n#\n");
  fprintf(f, "enrichment = 0            #Values between 0 and 1 to describe 18O APE. For example, 0.25 equals 25 APE.\n");
  fprintf(f, "instrument = %d            #0=Orbitrap, 1=FTICR\n", def.instrument);
  fprintf(f, "MS1_centroid = %d          #0=no, 1=yes\n", def.ms1Centroid);
  fprintf(f, "MS2_centroid = %d          #0=no, 1=yes\n", def.ms2Centroid);
  fprintf(f, "MS1_resolution = %d    #resolution at 400 m/z, value ignored if data are centroided\n", def.ms1Resolution);
  fprintf(f, "MS2_resolution = %d    #resolution at 400 m/z, value ignored if data are centroided\n", def.ms2Resolution);
  
  fprintf(f, "\n\n#\n# Crosslinker information\n#\n");
  fprintf(f, "predefined_crosslink = 1      #Value chosen from the XL_PARAMS table above.\n");
  fprintf(f, "mono_links_on_xl = %d          #allow monolinks to be combined with a crosslink. 0=no, 1=yes\n",(int)def.monoLinksOnXL);

  fprintf(f, "\n# Alternatively, a custom crosslinker can be specified using the parameters below.\n");
  fprintf(f, "# Format for cross_link is[amino acids][amino acids][mass mod][identifier]\n");
  fprintf(f, "# Format for mono_link is[amino acids][mass mod]\n");
  fprintf(f, "# One or more amino acids(uppercase only!!) can be specified for each linkage moiety\n");
  fprintf(f, "# Use lowercase 'n' or 'c' to indicate protein N - terminus or C - terminus\n");
  fprintf(f, "# Each parameter can be used multiple times for multiple crosslinkers or chemistries in the analysis\n");
  fprintf(f, "# Below is an example of how to replicate DSSO quenched with AmBic and H2O\n#\n");
  fprintf(f, "#cross_link = nK nK 158.003765 DSSO     #Typical DSSO crosslinker settings\n");
  fprintf(f, "#mono_link = nK 176.014330              #DSSO_H20_monolink\n");
  fprintf(f, "#mono_link = nK 175.030314              #DSSO_NH2_monolink\n");
  fprintf(f, "#xl_cleavage_product_mass = 54.010565   #DSSO cleavage product mass #1\n"); 
  fprintf(f, "#xl_cleavage_product_mass = 85.982635   #DSSO cleavage product mass #2\n");
  fprintf(f, "#xl_cleavage_product_mass = 103.993200  #DSSO cleavage product mass #3\n");
  
  fprintf(f, "\n\n#\n# Amino acid search modification. Use uppercase amino acid letters. n=peptide N-terminus, c=peptide C-terminus\n");
  fprintf(f, "# These parameters can be used multiple times to accomodate many possible modifications.\n#\n");
  fprintf(f, "fixed_modification = C 57.02146       #fixed modifications are applied to all amino acid instances.\n");
  fprintf(f, "#fixed_modification_protC = 15.994915  #fixed modifications to protC are applied to all protein C-termini.\n");
  fprintf(f, "#fixed_modification_protN = 15.994915  #fixed modifications to protN are applied to all protein N-termini.\n");
  fprintf(f, "#aa_mass = X 101.074328                #overwrite the default mass of amino acid. Also specify non-canonical amino acids (e.g. B,J,X,Z) to your search.\n");
  fprintf(f, "\nmodification = M 15.994915            #modifications are possible (differential) mass alterations.\n");
  fprintf(f, "#modification = N 0.984016             #use multiple modification lines for each possible modification mass.\n");
  fprintf(f, "#modification_protN = 42.010565        #modifications to protN are differential on protein N-termini.\n");
  fprintf(f, "#modification_protC = 15.994915        #modifications to protC are differential on protein C-termini.\n");
  fprintf(f, "\nmax_mods_per_peptide = %d              #limit the number of modifications on a peptide (does not include fixed_modifications)\n", def.maxMods);
  fprintf(f, "diff_mods_on_xl = %d                   #allow differential mods on crosslinks. 0=no, 1=yes\n",(int)def.diffModsOnXL);

  fprintf(f, "\n\n#\n# Digestion enzyme rules: see http://kojak-ms.org/param/enzyme.html\n#\n");
  fprintf(f, "enzyme = [KR]|{P} Trypsin     #must specify rule, followed by enzyme name.\n");
  fprintf(f, "max_miscleavages = %d          #number of missed enzyme cleavages to allow.\n", def.miscleave);
  fprintf(f, "                              #must be >=1 to find crosslinks on cleavage sites, and >=2 to find looplinks on cleavage sites.\n");

  fprintf(f, "\n\n#\n# Scoring Algorithm Parameters: specifies resolution and fragment ions to search.\n#\n");
  fprintf(f, "# fragment_bin_offset and fragment_bin_size influence algorithm precision and memory usage.\n");
  fprintf(f, "# They should be set appropriately for the data analyzed.\n");
  fprintf(f, "# For ion trap ms/ms:  1.0005 size, 0.4 offset\n");
  fprintf(f, "# For high res ms/ms:    0.01 size, 0.0 offset\n");
  fprintf(f, "fragment_bin_size = %.2lf      #in m/z units\n", def.binSize);
  fprintf(f, "fragment_bin_offset = %.1lf     #between 0.0 and 1.0\n", def.binOffset);
  fprintf(f, "ion_series_A = %d              #0=do not search, #1=search\n", (int)def.ionSeries[0]);
  fprintf(f, "ion_series_B = %d\n", (int)def.ionSeries[1]);
  fprintf(f, "ion_series_C = %d\n", (int)def.ionSeries[2]);
  fprintf(f, "ion_series_X = %d\n", (int)def.ionSeries[3]);
  fprintf(f, "ion_series_Y = %d\n", (int)def.ionSeries[4]);
  fprintf(f, "ion_series_Z = %d              #Z-dot values are used\n", (int)def.ionSeries[5]);

  fprintf(f, "\n\n#\n# Kojak Search Space Prameters: specifies breadth of data analysis.\n#\n");
  fprintf(f, "decoy_filter = %s %d     #identifier for all decoy protein sequences. 0=database has decoys, 1=must generate decoys in Kojak\n", def.decoy,(int)def.buildDecoy);
  fprintf(f, "e_value_depth = %d       #robustness of e-value histogram. Larger number improves e-value estimates, but increases computation time.\n", def.decoySize);
  fprintf(f, "\nppm_tolerance_pre = %.1lf   #mass tolerance on precursor when searching (in ppm)\n", def.ppmPrecursor);
  fprintf(f, "precursor_refinement = %d   #0 = off, 1 = attempt to correct precursor prediction errors from MS1 scan.\n", (int)def.precursorRefinement);
  fprintf(f, "isotope_error = %d          #search isotope peak offsets. 0=off, 1=one offset, 2=two offsets, 3=three offsets\n", def.isotopeError);
  fprintf(f, "prefer_precursor_pred = %d  #prefer precursor mono mass predicted by instrument software.\n", def.preferPrecursor);
  fprintf(f, "                           #  0 = ignore previous predictions\n");
  fprintf(f, "                           #  1 = use only previous predictions\n");
  fprintf(f, "                           #  2 = supplement predictions with additional analysis\n");
  fprintf(f, "\nspectrum_processing = %d    #0 = no, 1 = collapse MS2 isotope distributions to a single monoisotopic peak.\n", def.specProcess);
  fprintf(f, "min_spectrum_peaks = %d    #minimum peaks in a MS2 scan to be searched.\n", def.minPeaks);
  fprintf(f, "max_spectrum_peaks = %d     #maximum number of MS2 peaks to use during analysis. 0 uses all peaks.\n", def.maxPeaks);
  //TODO
  //fprintf(f, "\nmin_peptide_length = %d     #minimum number of amino acids per peptide searched.\n", def.minPepLen);
  //fprintf(f, "max_peptide_length = %d    #maximum number of amino acids per peptide searched.\n", def.maxPepLen);
  fprintf(f, "\nmin_peptide_mass = %.1lf   #minimum allowed peptide mass in Daltons.\n", def.minPepMass);
  fprintf(f, "max_peptide_mass = %.1lf  #maximum allowed peptide mass in Daltons.\n", def.maxPepMass);
  fprintf(f, "min_peptide_score = %.1lf    #lowest Xcorr allowed for alpha peptide consideration in the first pass.\n", def.minPepScore);
  fprintf(f, "top_count = %d             #number of candidate alpha peptides from first pass to attempt to crosslink in second pass.\n", def.topCount);
  
  fprintf(f, "\n\n#\n# Diagnostics: Only recommended for advanced users. If enabled, a diagnostic XML file is output with the Kojak results.\n#\n");
  fprintf(f, "# diagnostic = 502          #diagnose intermediate peptide calculations for any scan number.\n");
  fprintf(f, "# diagnostic = 503          #list multiple MS2 scan numbers, one per line.\n");
  fprintf(f, "# diagnostic = -1           #or specify -1 to diagnose all MS2 scans.\n");
  fprintf(f, "# diagnostic_histogram = 1  #include score histogram with the diagnostics.\n");

  fclose(f);
}

void KParams::logParam(string name, string value){
  pxwBasicXMLTag t;
  t.name=name;
  t.value=value;
  if (log != NULL) log->addParameter(t.name + " = " + t.value);
  xmlParams.push_back(t);
}

void KParams::logParam(pxwBasicXMLTag& t){
  if(log!=NULL) log->addParameter(t.name + " = " + t.value);
  xmlParams.push_back(t);
}

void KParams::parse(const char* cmd) {

  char *tok;
  char c_cmd[4096];
  char param[32];
  char tmpstr[256];

  int i;

  kMass m;
  kLinker x;

  string tstr;
  vector<string> values;

  pxwBasicXMLTag xml;

  strcpy(c_cmd,cmd);

  //Pre-process entire line to remove characters that should not be read
	//Replace first # with a terminator
	tok=strstr(c_cmd,"#");
	if(tok!=NULL) strncpy(tok,"\0",1);

	//if we have only white space, exit here
	strcpy(tmpstr,c_cmd);
	tok=strtok(tmpstr," \t\n\r");
	if(tok==NULL) return;

	//Check if we have a parameter (has '=' in it) or lots of random text.
	tok=strstr(c_cmd,"=");
  if(tok==NULL) {
    if (log != NULL) log->addParameterWarning("Unknown parameter line in config file: " + string(cmd));
    else printf("Unknown parameter line in config file: %s\n",cmd);
    return;
  }

  //Process parameters
	//Read parameter into param name (before = sign) and value (after = sign)
	tok=strtok(c_cmd," \t=\n\r");
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
  if (strcmp(param, "aa_mass") == 0){
    m.index = (int)values[0][0];
    m.mass = atof(&values[1][0]);
    if (values.size() > 2 && values[2][0]!='0') m.xl=true;
    else m.xl=false;
    params->aaMass->push_back(m);
    logParam("aa_mass",values[0] + " " + values[1]);

  } else if(strcmp(param,"cross_link")==0){
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
    logParam("cross_link",values[0]+" "+values[1]+" "+values[2]+" "+values[3]);

	} else if(strcmp(param,"database")==0){
    strcpy(params->dbFile,&values[0][0]);
    logParam("database",values[0]);

  } else if(strcmp(param,"diagnostic")==0){  //a value of -1 means diagnose all spectra, overriding any existing or following spectrum specifications
    if (atoi(&values[0][0])==-1) params->diag->clear();
    params->diag->push_back(atoi(&values[0][0]));
    logParam("diagnostic",values[0]);

  } else if (strcmp(param, "diagnostic_histogram") == 0) {  //any non-zero value means include the histogram in the diagnostics
    if (atoi(values[1].c_str()) != 0) params->diagHistogram = true;
    else params->diagHistogram = false;
    logParam("diagnostic_histogram", values[0]);

  } else if(strcmp(param,"diff_mods_on_xl")==0){
    if(atoi(&values[0][0])!=0) params->diffModsOnXL=true;
    else params->diffModsOnXL=false;
    logParam("diff_mods_on_xl",values[0]);

	} else if(strcmp(param,"decoy_filter")==0){
    if (values.size()!=2) {
      warn("ERROR: bad decoy_filter parameter. Suspected use of deprecated format.", 3);
      exit(-5);
    }
    strcpy(params->decoy,values[0].c_str());
    if (atoi(values[1].c_str()) == 0) params->buildDecoy = false;
    else params->buildDecoy = true;
    logParam("decoy_filter",values[0] + " " + values[1]);
 
  } else if(strcmp(param,"enrichment")==0){
    params->enrichment=atof(&values[0][0]);
    logParam("enrichment",values[0]);

  } else if(strcmp(param,"enzyme")==0){
    strcpy(params->enzyme,&values[0][0]);
    xml.name = "enzyme";
    xml.value = values[0];
    if (values.size() > 1) {
      strcpy(params->enzymeName, &values[1][0]);
      xml.value+= " " + values[1];
    } else {
      strcpy(params->enzymeName, "Unnamed");
      warn("Using deprecated enzyme paramter format. Enzyme name will be set to \"Unnamed\"",4);
    }
    logParam(xml);

  } else if (strcmp(param, "e_value_depth") == 0){
    params->decoySize = atoi(&values[0][0]);
    xml.name = "e_value_depth";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "export_mzID") == 0 || strcmp(param, "export_mzid") == 0){
    if (atoi(&values[0][0]) != 0) params->exportMzID = true;
    else params->exportMzID = false;
    xml.name = "export_mzID";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param, "export_pepXML")==0 || strcmp(param, "export_pepxml")==0){
    if(atoi(&values[0][0])!=0) params->exportPepXML=true;
    else params->exportPepXML=false;
    xml.name = "export_pepXML";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"export_percolator")==0){
    if(atoi(&values[0][0])!=0) params->exportPercolator=true;
    else params->exportPercolator=false;
    xml.name = "export_percolator";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"fixed_modification")==0){
    if(values.size()!=2){
      warn("ERROR: Invalid value(s) for fixed_modifcation parameter. Stopping analysis.", 3);
      exit(-5);
    } else {
      m.index=(int)values[0][0];
      m.mass=atof(&values[1][0]);
      if(m.mass!=0) params->fMods->push_back(m);
      xml.name = "fixed_modification";
      xml.value = values[0]+" "+values[1];
      logParam(xml);
    }

  } else if (strcmp(param, "fixed_modification_protC") == 0){
    m.index = (int)'%';
    m.mass = atof(&values[0][0]);
    if (m.mass != 0) params->fMods->push_back(m);
    xml.name = "fixed_modification_protC";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "fixed_modification_protN") == 0){
    m.index = (int)'$';
    m.mass = atof(&values[0][0]);
    if (m.mass != 0) params->fMods->push_back(m);
    xml.name = "fixed_modification_protN";
    xml.value = values[0];
    logParam(xml);

	} else if(strcmp(param,"fragment_bin_offset")==0){
    params->binOffset=1.0-atof(&values[0][0]);
    logParam("fragment_bin_offset",values[0]);

	} else if(strcmp(param,"fragment_bin_size")==0){
    params->binSize=atof(&values[0][0]);
    if(params->binSize<=0){
      warn("ERROR: Invalid value for fragment_bin_size parameter. Stopping analysis.",3);
      exit(-5);
    }
    logParam("fragment_bin_size",values[0]);

	} else if(strcmp(param,"instrument")==0){
    params->instrument=atoi(&values[0][0]);
    if(params->instrument<0 || params->instrument>1){
      warn("ERROR: instrument value out of range for instrument. Stopping analysis.",3);
      exit(-5);
    }
    logParam("instrument",values[0]);

  } else if (strcmp(param, "intermediate") == 0){
    params->intermediate = atoi(&values[0][0]);
    if(params->intermediate>0) params->turbo=false; //turn off turbo mode when using intermediate analysis.
    if (params->intermediate < 0){
      warn("ERROR: intermediate value out of range. Stopping analysis.", 3);
      exit(-5);
    }

  } else if(strcmp(param,"ion_series_A")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[0]=false;
    else params->ionSeries[0]=true;
    xml.name = "ion_series_A";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"ion_series_B")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[1]=false;
    else params->ionSeries[1]=true;
    xml.name = "ion_series_B";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"ion_series_C")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[2]=false;
    else params->ionSeries[2]=true;
    xml.name = "ion_series_C";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"ion_series_X")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[3]=false;
    else params->ionSeries[3]=true;
    xml.name = "ion_series_X";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"ion_series_Y")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[4]=false;
    else params->ionSeries[4]=true;
    xml.name = "ion_series_Y";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"ion_series_Z")==0){
    if(atoi(&values[0][0])==0) params->ionSeries[5]=false;
    else params->ionSeries[5]=true;
    xml.name = "ion_series_Z";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "isotope_error")==0){
    params->isotopeError = atoi(&values[0][0]);
    if (params->isotopeError < 0 || params->isotopeError>3){
      warn("ERROR: isotope_error has invalid value. Stopping analysis.",3);
      exit(-5);
    }
    xml.name = "isotope_error";
    xml.value = values[0];
    logParam(xml);

	} else if(strcmp(param,"max_miscleavages")==0){
    params->miscleave=atoi(&values[0][0]);
    xml.name = "max_miscleavages";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"max_mods_per_peptide")==0){
    params->maxMods=atoi(&values[0][0]);
    xml.name = "max_mods_per_peptide";
    xml.value = values[0];
    logParam(xml);

	} else if(strcmp(param,"max_peptide_mass")==0){
    params->maxPepMass=atof(&values[0][0]);
    xml.name = "max_peptide_mass";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"max_spectrum_peaks")==0){
    params->maxPeaks=atoi(&values[0][0]);
    xml.name = "max_spectrum_peaks";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "min_peptide_contribution") == 0){
    params->minPepUnique = atof(&values[0][0]);
    xml.name = "min_peptide_contribution";
    xml.value = values[0];
    logParam(xml);

	} else if(strcmp(param,"min_peptide_mass")==0){
    params->minPepMass=atof(&values[0][0]);
    xml.name = "min_peptide_mass";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "min_peptide_score") == 0){
    params->minPepScore = atof(&values[0][0]);
    xml.name = "min_peptide_score";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "min_spectrum_peaks") == 0) {
    params->minPeaks = atoi(&values[0][0]);
    xml.name = "min_spectrum_peaks";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"modification")==0){
    if (values.size() != 2){
      warn("ERROR: Invalid value(s) for modifcation parameter. Stopping analysis.", 3);
      exit(-5);
    } else {
      m.xl=false;
      m.index=(int)values[0][0];
      m.mass=atof(&values[1][0]);
      if(m.mass!=0 && !checkMod(m)) params->mods->push_back(m);
      xml.name = "modification";
      xml.value = values[0]+" "+values[1];
      logParam(xml);
    }

  } else if (strcmp(param, "modification_protC") == 0){
    m.xl = false;
    m.index = (int)'%';
    m.mass = atof(&values[0][0]);
    if (m.mass!=0) { //acceptable to use placeholder; don't load it as a parameter
      if (!checkMod(m)) params->mods->push_back(m);
    }
    xml.name = "modification_protC";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "modification_protN") == 0){
    m.xl = false;
    m.index = (int)'$';
    m.mass = atof(&values[0][0]);
    if(m.mass != 0) { //acceptable to use placeholder; don't load it as a parameter
      if (!checkMod(m)) params->mods->push_back(m);
    }
    xml.name = "modification_protN";
    xml.value = values[0];
    logParam(xml);
    
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
    logParam(xml);

  } else if(strcmp(param,"mono_links_on_xl")==0){
    if(atoi(&values[0][0])!=0) params->monoLinksOnXL=true;
    else params->monoLinksOnXL=false;
    xml.name = "mono_links_on_xl";
    xml.value = values[0];
    logParam(xml);

	} else if(strcmp(param,"MS_data_file")==0){
    strcpy(params->msFile,&values[0][0]);
    xml.name = "MS_data_file";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"MS1_centroid")==0){
    params->ms1Centroid=atoi(&values[0][0]);
    xml.name = "MS1_centroid";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"MS2_centroid")==0){
    params->ms2Centroid=atoi(&values[0][0]);
    xml.name = "MS2_centroid";
    xml.value = values[0];
    logParam(xml);

	} else if(strcmp(param,"MS1_resolution")==0){
    params->ms1Resolution=atoi(&values[0][0]);
    xml.name = "MS1_resolution";
    xml.value = values[0];
    logParam(xml);

	} else if(strcmp(param,"MS2_resolution")==0){
		params->ms2Resolution=atoi(&values[0][0]);
    xml.name = "MS2_resolution";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "15N_filter") == 0){
    strcpy(params->n15Label, &values[0][0]);
    xml.name = "15N_filter";
    xml.value = values[0];
    logParam(xml);

	} else if(strcmp(param,"output_file")==0){
    warn(param,2);

	} else if(strcmp(param,"percolator_file")==0){
    warn(param,2);

  } else if(strcmp(param,"percolator_version")==0){
    params->percVersion=atof(&values[0][0]);
    xml.name = "percolator_version";
    xml.value = values[0];
    logParam(xml);

	} else if(strcmp(param,"ppm_tolerance_pre")==0){
    params->ppmPrecursor=atof(&values[0][0]);
    xml.name = "ppm_tolerance_pre";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "precursor_refinement") == 0){
    if (atoi(&values[0][0]) == 0) params->precursorRefinement = false;
    else params->precursorRefinement = true;
    xml.name = "precursor_refinement";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "predefined_crosslink") == 0){
    int iID=atoi(values[0].c_str());
    size_t a;
    for(a=0;a<pXLTable.size();a++) if(pXLTable[a].id==iID) break;
    if(a==pXLTable.size()){
      warn("ERROR: predifined_crosslink ID number not found in XL_PARAMS table. Please make sure XL_PARAMS precedes the predifined_crosslink parameter and contains desired ID number.", 3);
      exit(-5);
    }
    xml.name = "predefined_crosslink";
    xml.value = values[0];
    logParam(xml);

    //expand predefinition to add a crosslink
    x.motifA = pXLTable[a].targetA;
    x.motifB = pXLTable[a].targetB;
    x.mass = pXLTable[a].linkMass;
    x.mono = 0;
    x.label = pXLTable[a].xlName;
    params->xLink->push_back(x);
    logParam("predefined_crosslink:cross_link", x.motifA + " " + x.motifB + " " + to_string(x.mass) + " " + x.label);

    //monolinks from the first half of the crosslinker
    for(size_t b=0;b<pXLTable[a].monoA.size();b++){
      m.xl = true;
      m.mass = pXLTable[a].monoA[b];
      for (i = 0; i < pXLTable[a].targetA.size(); i++){
        if (pXLTable[a].targetA[i] == 'c') m.index = (int)'%';
        else if (pXLTable[a].targetA[i] == 'n') m.index = (int)'$';
        else m.index = (int)pXLTable[a].targetA[i];
        if (!checkMod(m)) params->mods->push_back(m);
      }
      xml.name = "predefined_crosslink:mono_link";
      xml.value = pXLTable[a].targetA + " " + to_string(m.mass);
      logParam(xml);
    }

    //monolinks from the second half of the crosslinker (heterobifunctional)
    for (size_t b = 0; b<pXLTable[a].monoB.size(); b++){
      m.xl = true;
      m.mass = pXLTable[a].monoB[b];
      for (i = 0; i < pXLTable[a].targetB.size(); i++){
        if (pXLTable[a].targetB[i] == 'c') m.index = (int)'%';
        else if (pXLTable[a].targetB[i] == 'n') m.index = (int)'$';
        else m.index = (int)pXLTable[a].targetB[i];
        if (!checkMod(m)) params->mods->push_back(m);
      }
      xml.name = "predefined_crosslink:mono_link";
      xml.value = pXLTable[a].targetB + " " + to_string(m.mass);
      logParam(xml);
    }

    //cleavage products if cleavable crosslinker
    for(size_t b=0;b<pXLTable[a].cleavageMass.size();b++){
      params->cleavageProducts->push_back(pXLTable[a].cleavageMass[b]);
      xml.name = "predefined_crosslink:xl_cleavage_product_mass";
      xml.value = to_string(pXLTable[a].cleavageMass[b]);
      logParam(xml);
    }
    

  } else if(strcmp(param,"prefer_precursor_pred")==0){
    params->preferPrecursor=atoi(&values[0][0]);
    xml.name = "prefer_precursor_pred";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "remove_precursor") == 0){
    params->removePrecursor = atof(&values[0][0]);
    xml.name = "remove_precursor";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "results_path") == 0){
    strcpy(params->resPath, &values[0][0]);
    xml.name = "results_path";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"search_dimers")==0){
    if (atoi(&values[0][0]) == 0) params->dimers = false;
    else params->dimers = true;
    xml.name = "search_dimers";
    xml.value = values[0];
    logParam(xml);
    //warn(param,2);

  } else if(strcmp(param,"search_dimers_xl")==0){
    if(atoi(&values[0][0])==0) params->dimersXL=false;
    else params->dimersXL=true;
    xml.name = "search_dimers_xl";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"spectrum_processing")==0) {
    params->specProcess=atoi(&values[0][0]);
    xml.name = "spectrum_processing";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "split_pepXML") == 0 || strcmp(param, "split_pepxml") == 0) {
    if (atoi(&values[0][0]) != 0) params->splitPepXML = true;
    else params->splitPepXML = false;
    xml.name = "split_pepXML";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"threads")==0) {
    params->threads=atoi(&values[0][0]);
    int iCores;
#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    iCores = sysinfo.dwNumberOfProcessors;
#else
    iCores = sysconf(_SC_NPROCESSORS_ONLN);
#endif
    if(params->threads==0){
      params->threads=iCores;
      char tst[8];
      sprintf(tst,"%d",iCores);
      string ts="Setting number of threads from 0 to maximum of ";
      ts+=tst;
      warn(ts,4);
    } else if(params->threads<0) {
      params->threads+=iCores;
      if(params->threads<1) {
        warn("Net number of threads is less than 1. Forcing thread value of 1.", 4);
        params->threads=1;
      } else {
        char tst[8];
        sprintf(tst, "%d", iCores);
        string ts = "Adjusting number of threads to " + values[0] + " below to maximum of ";
        sprintf(tst, "%d", iCores);
        ts+=tst;
        warn(ts,4);
      }
    } else if (params->threads>iCores) {
      warn("Requested threads exceeds number of CPU cores. Performance may be degraded.",4);
    }
    xml.name = "threads";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"top_count")==0) {
    params->topCount=atoi(&values[0][0]);
    if(params->topCount>50) {
      warn("WARNING: top_count is larger than expected. Are you using a pre-version 1.6.0 value? An appropriate value is likely between 5 and 50.",4);
    }
    if (params->topCount<1){
      warn("ERROR: top_count must be greater than zero. Stopping Kojak.", 3);
      params->topCount=1;
    }
    xml.name = "top_count";
    xml.value = values[0];
    logParam(xml);

  } else if(strcmp(param,"truncate_prot_names")==0) {
    params->truncate=atoi(&values[0][0]);
    xml.name = "truncate_prot_names";
    xml.value = values[0];
    logParam(xml);

  } else if (strcmp(param, "turbo_button") == 0){
    warn(param,2);

  } else if(strcmp(param,"use_comet_xcorr")==0){
    warn(param, 2);

  } else if (strcmp(param, "xl_cleavage_product_mass") == 0){
    params->cleavageProducts->push_back(atof(&values[0][0]));
    xml.name = "xl_cleavage_product_mass";
    xml.value = values[0];
    logParam(xml);

	} else {
		warn(param,1);
	}

}

void KParams::splitMasses(char*& c, std::vector<double>& v){
  string s;
  v.clear();
  if(c[0]=='x') return;
  for(size_t a=0;a<strlen(c);a++){
    if(c[a]==','){
      v.push_back(atof(s.c_str()));
      s.clear();
    } else s+=c[a];
  }
  v.push_back(atof(s.c_str()));
}

void KParams::setLog(KLog* c){
  log = c;
}

void KParams::setParams(kParams* p){
  params=p;
}

void KParams::warn(const char* c, int i){
	switch(i){
		case 0:
      if(log!=NULL) log->addParameterWarning("Parameter " + string(c) + " has no value.");
      else printf("  WARNING: Parameter %s has no value.",c);
			break;
		case 1:
      if (log != NULL) log->addParameterWarning("Unknown parameter: " + string(c));
			else printf("  WARNING:  %s\n",c);
			break;
		case 2:
      if (log != NULL) log->addParameterWarning("Parameter " + string(c) + "has been deprecated and will be ignored.");
      else printf("  WARNING: Parameter %s has been deprecated and will be ignored.\n",c);
      break;
    case 3:
      if(log!=NULL) log->addError(string(c));
      else printf("  %s\n", c);
      break;
    case 4:
      if(log!=NULL) log->addParameterWarning(string(c));
      else printf("  WARNING: %s\n", c);
      break;
		default:
      if (log != NULL) log->addParameterWarning(string(c));
      else printf("  %s\n", c);
			break;
	}
}

void KParams::warn(string c, int i){
  warn(c.c_str(),i);
}

//Takes relative path and finds absolute path
bool KParams::processPath(const char* cwd, const char* in_path, char* out_path){
  //if windows or unix in_path, just copy it to out_path
  if (strlen(in_path) > 0 && in_path[0] == '/'){ //unix
    strcpy(out_path, in_path);
    return true;
  }
  if (strlen(in_path) > 1 && in_path[1] == ':'){ //windows
    strcpy(out_path, in_path);
    return true;
  }

  //tokenize cwd
  char* tok;
  char str[1024];
  strcpy(str, cwd);
  string s;
  vector<string> v;

  tok = strtok(str, "\\/\n\r");
  while (tok != NULL){
    s = tok;
    v.push_back(s);
    tok = strtok(NULL, "\\/\n\r");
  }

  //tokenize in_path
  strcpy(str, in_path);

  tok = strtok(str, "\\/\n\r");
  while (tok != NULL){
    if (strcmp(tok, "..") == 0) {
      v.pop_back();
    } else if (strcmp(tok, ".") == 0){
      //do nothing
    } else {
      s = tok;
      v.push_back(s);
    }
    tok = strtok(NULL, "\\/\n\r");
  }

  //build absolute path
#ifdef _MSC_VER
  s.clear();
#else
  s.clear();
  s += slashdir;
#endif
  for (size_t i = 0; i < v.size(); i++){
    s += v[i];
    s += slashdir;
  }
  s[s.size() - 1] = '\0';
  strcpy(out_path, &s[0]);
  return true;

}
