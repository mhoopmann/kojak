#include "KLog.h"

using namespace std;

KLog::KLog(){
  clear();
}

void KLog::addDBWarning(std::string msg){
  string myMsg = "WARNING: " + msg;
  cout << myMsg << endl;
  vDBWarnings.push_back(myMsg);
}

void KLog::addError(std::string msg){
  strError="\n";
  time_t timeNow;
  time(&timeNow);
  strError += ctime(&timeNow);
  strError.pop_back();
  strError += "\tERROR: " + msg;
  cout << strError << endl;
  exportLog();
  exit(-4004);
}

void KLog::addMessage(string msg, bool silent){
  if(!silent) cout << msg << endl;

  time_t timeNow;
  time(&timeNow);
  string myMsg=ctime(&timeNow);
  myMsg.pop_back();
  myMsg+="\t"+msg;

  vMsg.push_back(myMsg);
}

void KLog::addParameter(std::string msg){
  vParams.push_back(msg);
}

void KLog::addParameterWarning(std::string msg){
  string myMsg = "WARNING: " + msg;
  cout << myMsg << endl;
  vParamWarnings.push_back(myMsg);
}

void KLog::addWarning(size_t id, std::string msg){
  if(idIndex[id]!=SIZE_MAX){
    vWarnings[idIndex[id]].count++;
  } else {
    string myMsg = "WARNING: " + msg;
    kWarning w;
    w.count=1;
    w.msg=myMsg;
    idIndex[id]=vWarnings.size();
    vWarnings.push_back(w);
    cout << myMsg << endl;
  }
}

void KLog::clear(){
  for(size_t i=0;i<LOGSZ;i++) idIndex[i]=SIZE_MAX;
  logFile.clear();
  vMsg.clear();
  vWarnings.clear();
  strError.clear();
}

void KLog::exportLog(){
  if(logFile.size()==0) return;
  size_t i;
  FILE* f=fopen(logFile.c_str(),"wt");

  fprintf(f,"***** PARAMETERS *****\n");
  for(i=0;i<vParams.size();i++) fprintf(f,"%s\n",vParams[i].c_str());
  fprintf(f,"\n");
  for(i=0;i<vParamWarnings.size();i++) fprintf(f,"%s\n",vParamWarnings[i].c_str());

  fprintf(f, "\n\n***** DATABASE *****\n");
  fprintf(f, "%s",dbInfo.c_str());
  fprintf(f, "\n");
  for (i = 0; i<vDBWarnings.size(); i++) fprintf(f, "%s\n", vDBWarnings[i].c_str());

  fprintf(f,"\n\n***** KOJAK LOG *****\n");
  for(i=0;i<vMsg.size();i++) fprintf(f,"%s\n",vMsg[i].c_str());

  if(strError.size()>0) {
    fprintf(f,"\n\n***** KOJAK EXECUTION HALTED DUE TO ERROR!!! *****");
    fprintf(f,"%s\n", strError.c_str());
  }

  if(vWarnings.size()>0 || vParamWarnings.size()>0){
    fprintf(f, "\n\n***** WARNINGS *****\n");
    for (i = 0; i<vWarnings.size(); i++) fprintf(f, "(%d instances) %s\n", vWarnings[i].count,vWarnings[i].msg.c_str());
    if(vParamWarnings.size()>0) fprintf(f,"Check Parameter log above for warnings originating in configuration file.\n");
    if (vDBWarnings.size()>0) fprintf(f, "Check Database log above for warnings originating in the FASTA database file.\n");
  }

  fclose(f);
}

void KLog::setDBinfo(std::string fn, int prot, int pep, int linkPep){
  dbInfo="FASTA database: ";
  dbInfo+=fn;
  char tStr[512];
  sprintf(tStr,"\nTotal proteins: %d\nTotal peptides to search: %d (%d linkable).",prot,pep,linkPep);
  dbInfo+=tStr;
}

void KLog::setLog(char* fn){
  string f=fn;
  setLog(f);
}

void KLog::setLog(std::string fn){
  logFile=fn;
}