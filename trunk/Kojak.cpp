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

int main(int argc, char* argv[]){

  cout << "Kojak version 1.0, May 19 2014" << endl;
  cout << "Copyright Michael Hoopmann, Institute for Systems Biology" << endl;
  if(argc<2){
    cout << "Usage: Kojak <Config File>" << endl;
    return 1;
  }

  time_t timeNow;
  time(&timeNow);
  cout << "Time at start: " << ctime(&timeNow) << endl;
  
  unsigned int i;

  //Step #1: Read in parameters
  kParams params;
  KParams p(&params);
  p.parseConfig(argv[1]);

  //Step #2: Read in database and generate peptide lists
  KDatabase db;
  for(i=0;i<params.fMods->size();i++) db.addFixedMod(params.fMods->at(i).index,params.fMods->at(i).mass);
  db.buildDB(params.dbFile);
  db.buildPeptides(params.minPepMass,params.maxPepMass,params.miscleave);
  
  //Step #3: Read in spectra and map precursors
  KData spec(&params);
  spec.setVersion("1.0");
  spec.readSpectra2();
  spec.mapPrecursors();
  if(params.xcorr==1) spec.xCorr();
  for(i=0;i<params.mLink->size();i++) spec.setLinker(params.mLink->at(i),1);
  for(i=0;i<params.xLink->size();i++) spec.setLinker(params.xLink->at(i),0);

  time(&timeNow);
  cout << "Time at analysis: " << ctime(&timeNow) << endl;

  //Step #4: Analyze single peptides, monolinks, and crosslinks
  KAnalysis anal(params);
  cout << "Scoring non-linked peptides. ";
  anal.analyzePeptides(db,spec,false);
  if(params.dimers>0 && params.relaxedAnalysis==0){
    cout << "Scoring non-covalent dimerizations. ";
    anal.analyzePeptidesNonCovalent(db,spec);
  }
  cout << "Scoring linked peptides. ";
  anal.analyzePeptides(db,spec,true);
  if(params.relaxedAnalysis>0){
    cout << "Finalizing relaxed analysis. ";
    anal.analyzeRelaxed(db,spec);
  }

  //Step #5: Output results
  if(params.percolator[0]!='\0') spec.outputPercolator(params.percolator,params.decoy,db,params.truncate);
  if(params.enrichment>0)  spec.outputResults(params.outFile,db,true);
  else spec.outputResults(params.outFile,db,false);

  cout << "Finished." << endl;

  time(&timeNow);
  cout << "Time at finish: " << ctime(&timeNow) << endl;

  //printf("XCorr was called %lld times.\n",anal.xCorrCount);

  return 0;
}
