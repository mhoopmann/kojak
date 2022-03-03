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
#include "KojakManager.h"

using namespace std;

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

  int i;
  int fc=1;
  KojakManager manager;
  if(!manager.setParams(argv[1])) return -3;
  if(argc>2){
    manager.clearFiles();
    for(i=2;i<argc;i++) fc=manager.setFile(argv[i]);
  }
  manager.run();

  time(&timeNow);
  cout << " Time at finish: " << ctime(&timeNow) << endl;
  cout << "\n****** Finished Kojak Analysis ******" << endl;
  return 0;
}
