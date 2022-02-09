#include "KojakManager.h"
#include "KAnalysis.h"
#include "KData.h"
#include "KDB.h"
#include "KIons.h"

using namespace std;

KojakManager::KojakManager(){
  param_obj.setParams(&params);
  param_obj.setLog(&log);
}

void KojakManager::clearFiles(){
  files.clear();
}

int KojakManager::setFile(const char* fn){
  kFile f;
  f.input = fn;
  if (!getBaseFileName(f.base, fn, f.ext)){
    cout << "  Error with batch file parameter:  " << fn << " - unknown file or extension." << endl;
    return -1;
  }
  files.push_back(f);
  return (int)files.size();
}

int KojakManager::setFile(string& s){
  return setFile(s.c_str());
}

void KojakManager::setParam(const char* p){
  param_obj.parse(p);
}

void KojakManager::setParam(string& s){
  setParam(s.c_str());
}

bool KojakManager::setParams(const char* fn){
  kFile f;
  paramFile=fn; //hang onto this

  cout << " Parameter file: " << fn << endl;
  param_obj.parseConfig(fn);
  f.input = params.msFile;
  if (!getBaseFileName(f.base, params.msFile, f.ext)){
    cout << "  Error with input file parameter: " << params.msFile << " - unknown file or extension." << endl;
    return false;
  }
  files.push_back(f);
  return true;

}

bool KojakManager::setParams(string& s){
  return setParams(s.c_str());
}

int KojakManager::run(){
  time_t timeNow;
  size_t i;

  //Step 1: Prepare from settings
  KData spec(&params);
  spec.setLog(&log);
  spec.setVersion(VERSION);
  for (i = 0; i<params.xLink->size(); i++) spec.setLinker(params.xLink->at(i));
  spec.buildXLTable();

  //Step #2: Read in database and generate peptide lists
  KDatabase db;
  db.setLog(&log);
  for (i = 0; i<params.fMods->size(); i++) db.addFixedMod(params.fMods->at(i).index, params.fMods->at(i).mass);
  for (i = 0; i<params.aaMass->size(); i++) db.setAAMass((char)params.aaMass->at(i).index, params.aaMass->at(i).mass, params.aaMass->at(i).xl);
  if (strlen(params.n15Label)>0) db.setN15Label(params.n15Label);
  if (!db.setEnzyme(params.enzyme)) exit(-3);
  db.setXLTable(spec.getXLTable(), 128, 20);
  cout << "\n Reading FASTA database: " << params.dbFile << endl;
  string str = params.decoy;
  if (!db.buildDB(params.dbFile,str)){
    cout << "  Error opening database file: " << params.dbFile << endl;
    return -1;
  }
  if (params.buildDecoy) db.buildDecoy(str);
  db.buildPeptides(params.minPepMass, params.maxPepMass, params.miscleave);
  log.setDBinfo(string(params.dbFile),db.getProteinDBSize(),db.getPeptideListSize(),db.linkablePepCount);

  //Step #3: Read in spectra and map precursors
  //Iterate over all input files
  for (i = 0; i<files.size(); i++){

    //set up our log
    log.clear();
    param_obj.buildOutput(&files[i].input[0], &files[i].base[0], &files[i].ext[0]);
    log.setLog(param_obj.logFile);
    log.addMessage("Kojak version: " + string(VERSION), true);
    log.addMessage("Parameter file: " + paramFile,true);

    if (strcmp(params.ext, ".mgf") == 0 && params.precursorRefinement){
      log.addError("Cannot perform precursor refinement using MGF files. Please disable by setting precursor_refinement=0");
      return -10;
    }
    log.addMessage("Reading spectra data file: " + files[i].input,true);
    cout << "\n Reading spectra data file: " << files[i].input.c_str() << " ... ";
    if (!spec.readSpectra()){
      log.addError("Error reading MS_data_file: " + files[i].input);
      return -2;
    }
    spec.mapPrecursors();
    spec.doXCorr(params);
    //spec.xCorr(params.xcorr);

    //Step #4: Analyze single peptides, monolinks, and crosslinks
    KAnalysis anal(params, &db, &spec);
    anal.setLog(&log);

    log.addMessage("Start spectral search.",true);
    time(&timeNow);
    cout << "\n Start spectral search: " << ctime(&timeNow);
    log.addMessage("Scoring peptides (first pass).",true);
    cout << "  Scoring peptides ... ";
    anal.doPeptideAnalysis();

    //if(params.intermediate>0) spec.outputIntermediate(db);

    char ts[16];
    sprintf(ts,"%d",params.decoySize);
    log.addMessage("Calculating e-values (" + string(ts) + ")",true);
    cout << "  Calculating e-values (" << params.decoySize << ")... ";
    anal.doEValueAnalysis();

    log.addMessage("Finish spectral search.",true);
    time(&timeNow);
    cout << " Finished spectral search: " << ctime(&timeNow) << endl;

    //Future diagnostics
    //spec.diagSinglet();

    //Step #5: Output results
    log.addMessage("Exporting results.",true);
    cout << " Exporting Results." << endl;
    spec.outputResults(db, param_obj);

    log.addMessage("Finished Kojak analysis.",true);
    log.exportLog();

  }

  return 0;
}

bool KojakManager::getBaseFileName(string& base, const char* fName, string& extP) {
  char file[256];
  char ext[256];
  char *tok;
  char preExt[256];
  unsigned int i;

  strcpy(ext, "");

  strcpy(file, fName);
  tok = strtok(file, ".\n");
  while (tok != NULL){
    strcpy(preExt, ext);
    strcpy(ext, tok);
    tok = strtok(NULL, ".\n");
  }

  for (i = 0; i<strlen(ext); i++) ext[i] = toupper(ext[i]);
  for (i = 0; i<strlen(preExt); i++) preExt[i] = toupper(preExt[i]);

  base = fName;
  if (strcmp(ext, "MZML") == 0) {
    base[base.size() - 5] = '\0';
    extP = ".mzML";
    return true;
  }
  if (strcmp(ext, "MZXML") == 0) {
    base[base.size() - 6] = '\0';
    extP = ".mzXML";
    return true;
  }
  if (strcmp(ext, "GZ") == 0) {
    if (strcmp(preExt, "MZML") == 0){
      base[base.size() - 8] = '\0';
      extP = ".mzML.gz";
      return true;
    }
    if (strcmp(preExt, "MZXML") == 0) {
      base[base.size() - 9] = '\0';
      extP = ".mzXML.gz";
      return true;
    }
  }
  if (strcmp(ext, "RAW") == 0) {
    base[base.size() - 4] = '\0';
    extP = ".raw";
    return true;
  }
  if (strcmp(ext, "MGF") == 0) {
    base[base.size() - 4] = '\0';
    extP = ".mgf";
    return true;
  }
  if (strcmp(ext, "MS2") == 0) {
    base[base.size() - 4] = '\0';
    extP = ".ms2";
    return true;
  }
  return false;
}
