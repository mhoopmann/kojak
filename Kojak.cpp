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

#include "KAnalysis.h"
#include "KData.h"
#include "KDB.h"
#include "KIons.h"
#include "KParams.h"

#define VERSION "2.0.0-dev"
#define BDATE "Oct 16 2018"

bool getBaseFileName(string& base, char* fName, string& extP);

int main(int argc, char* argv[]){

  cout << "\nKojak version " << VERSION << ", " << BDATE << endl;
  cout << "Copyright Michael Hoopmann, Institute for Systems Biology" << endl;
  cout << "Visit http://kojak-ms.org for full documentation." << endl;
  if(argc<2){
    cout << "Usage: Kojak <Config File> [<Data File>...]" << endl;
    return 1;
  }

  cout << "\n****** Begin Kojak Analysis ******" << endl;

  time_t timeNow;
  time(&timeNow);
  cout << " Time at start: " << ctime(&timeNow) << endl;

  unsigned int i;
  int k;
  kFile f;
  vector<kFile> files;

  //Step #1: Read in parameters & batch files for analysis
  cout << " Parameter file: " << argv[1] << endl;
  kParams params;
  KParams p(&params);
  p.parseConfig(argv[1]);

  if(argc==2){
    f.input=params.msFile;
    if(!getBaseFileName(f.base,params.msFile,f.ext)){
      cout << "  Error with input file parameter: " << params.msFile << " - unknown file or extension." << endl;
      return -1;
    }
    files.push_back(f);
  } else {
    for(k=2;k<argc;k++){
      f.input=argv[k];
      if(!getBaseFileName(f.base,argv[k],f.ext)){
        cout << "  Error with batch file parameter:  " << argv[k] << " - unknown file or extension." << endl;
        return -1;
      }
      files.push_back(f);
    }
  }
  KData spec(&params);
  spec.setVersion(VERSION);
  for(i=0;i<params.xLink->size();i++) spec.setLinker(params.xLink->at(i));
  spec.buildXLTable();

  //Step #2: Read in database and generate peptide lists
  KDatabase db;
  for(i=0;i<params.fMods->size();i++) db.addFixedMod(params.fMods->at(i).index,params.fMods->at(i).mass);
  for(i=0;i<params.aaMass->size();i++) db.setAAMass((char)params.aaMass->at(i).index,params.aaMass->at(i).mass,params.aaMass->at(i).xl);
  if (strlen(params.n15Label)>0) db.setN15Label(params.n15Label);
  if(!db.setEnzyme(params.enzyme)) exit(-3);
  db.setXLTable(spec.getXLTable(),128,20);
  cout << "\n Reading FASTA database: " << params.dbFile << endl;
  if(!db.buildDB(params.dbFile)){
    cout << "  Error opening database file: " << params.dbFile << endl;
    return -1;
  }
  db.buildPeptides(params.minPepMass,params.maxPepMass,params.miscleave);


  //Step #3: Read in spectra and map precursors
  //Iterate over all input files
  for(i=0;i<files.size();i++){
    p.buildOutput(&files[i].input[0], &files[i].base[0], &files[i].ext[0]);
 
    //strcpy(params.msFile,&files[i].input[0]);
    //strcpy(params.outFile,&files[i].base[0]);
    //strcpy(params.ext,&files[i].ext[0]);

    if (strcmp(params.ext, ".mgf") == 0 && params.precursorRefinement){
      cout << "\n ERROR: Cannot perform precursor refinement using MGF files. Please disable by setting precursor_refinement=0" << endl;
      exit(-10);
    }
    cout << "\n Reading spectra data file: " << &files[i].input[0] << " ... ";
    if(!spec.readSpectra()){
      cout << "  Error reading MS_data_file. Exiting." << endl;
      return -2;
    }
    spec.mapPrecursors();
    spec.xCorr(params.xcorr);

    //Step #4: Analyze single peptides, monolinks, and crosslinks
    KAnalysis anal(params, &db, &spec);

    time(&timeNow);
    cout << "\n Start spectral search: " << ctime(&timeNow);
    cout << "  Scoring peptides ... ";
    anal.doPeptideAnalysis();

    //if(params.intermediate>0) spec.outputIntermediate(db);

    cout << "  Calculating e-values (" << DECOY_SIZE << ")... ";
    anal.doEValueAnalysis();

    time(&timeNow);
    cout << " Finished spectral search: " << ctime(&timeNow) << endl;

    //Future diagnostics
    //spec.diagSinglet();

    //Step #5: Output results
    cout << " Exporting Results." << endl;
    spec.outputResults(db,p);

  }

  cout << "\n****** Finished Kojak Analysis ******" << endl;

  return 0;
}

bool getBaseFileName(string& base, char* fName, string& extP) {
  char file[256];
	char ext[256];
	char *tok;
	char preExt[256];
	unsigned int i;

	strcpy(ext,"");

	strcpy(file,fName);
	tok=strtok(file,".\n");
	while(tok!=NULL){
		strcpy(preExt,ext);
		strcpy(ext,tok);
		tok=strtok(NULL,".\n");
	}

	for(i=0;i<strlen(ext);i++) ext[i]=toupper(ext[i]);
	for(i=0;i<strlen(preExt);i++) preExt[i]=toupper(preExt[i]);

  base=fName;
  if(strcmp(ext,"MZML")==0) {
    base[base.size()-5]='\0';
    extP=".mzML";
    return true;
  }
  if(strcmp(ext,"MZXML")==0) {
    base[base.size()-6]='\0';
    extP=".mzXML";
    return true;
  }
	if(strcmp(ext,"GZ")==0) {
    if(strcmp(preExt,"MZML")==0){
      base[base.size()-8]='\0';
      extP=".mzML.gz";
      return true;
    }
    if(strcmp(preExt,"MZXML")==0) {
      base[base.size()-9]='\0';
      extP=".mzXML.gz";
      return true;
    }
	}
  if(strcmp(ext,"RAW")==0) {
    base[base.size()-4]='\0';
    extP=".raw";
    return true;
  }
  if(strcmp(ext,"MGF")==0) {
    base[base.size()-4]='\0';
    extP=".mgf";
    return true;
  }
  if(strcmp(ext, "MS2")==0) {
    base[base.size()-4]='\0';
    extP=".ms2";
    return true;
  }
	return false;
}

