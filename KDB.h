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

#ifndef _KDB_H
#define _KDB_H

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "KStructs.h"

class KDatabase{
public:

  //Constructors & Destructors
  KDatabase();

  //Operators
  kDB& operator[ ](const int& i);

  //User Functions
  bool  buildDB       (char* fname);                     //Reads FASTA file and populates vDB
  bool  buildPeptides (double min, double max, int mis); //Make peptide list within mass boundaries and miscleavages.

  //Accessors & Modifiers
  void                addFixedMod         (char mod, double mass);
  kDB&                at                  (const int& i);
  double              getAAMass           (char aa, bool n15=false);
  kEnzymeRules&       getEnzymeRules      ();
  kPeptide&           getPeptide          (int index);
  std::vector<kPeptide>*   getPeptideList();
  int                 getPeptideListSize  ();
  bool                getPeptideSeq       (int index, int start, int stop, char* str);
  bool                getPeptideSeq       (int index, int start, int stop, std::string& str);
  bool                getPeptideSeq       (kPeptide& p, std::string& str);
  bool                getPeptideSeq       (int pepIndex, std::string& str);
  kDB&                getProtein          (int index);
  int                 getProteinDBSize    ();
  void                setAAMass           (char aa, double mass, bool n15 = false);
  bool                setEnzyme           (char* str);
  void                setN15Label         (char* str);
  void                setXLTable          (char** arr, int szA, int szB);

private:
  
  //Data Members
  double        AA[128];   //Amino acid masses
  double        AAn15[128];
  double        fixMassPepC;
  double        fixMassPepN;
  double        fixMassProtC;
  double        fixMassProtN;
  char          xlTable[128][20];
  kEnzymeRules  enzyme;    //Where to cut to generate peptides
  std::string   n15Label;

  std::vector<kDB>      vDB;    //Entire FASTA database stored in memory
  std::vector<kPeptide> vPep;   //List of all peptides

  void addPeptide(int index, int start, int len, double mass, kPeptide& p, std::vector<kPeptide>& vP, bool bN, bool bC, bool bN15, char xlSites);
  bool checkAA(kPeptide& p, size_t i, size_t start, size_t n, size_t seqSize, bool& bN, bool& bC);

  //Utility functions (for sorting)
  static int compareMass      (const void *p1, const void *p2);
  static int compareSequence  (const void *p1, const void *p2);
  static bool compareSequenceB(const kPepSort& p1, const kPepSort& p2);

};

#endif
