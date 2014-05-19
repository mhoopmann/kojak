#ifndef _KDB_H
#define _KDB_H

#include <iostream>
#include <string>
#include <vector>

#include "KStructs.h"

using namespace std;


class KDatabase{
public:

  //Constructors & Destructors
  KDatabase();

  kDB& operator[ ](const int& i);

  //Accessors & Modifiers
  void                addFixedMod         (char mod, double mass);
  kPeptide&           getPeptide          (int index, bool linkable);
  vector<kPeptide>&   getPeptideList      (bool linkable);
  int                 getPeptideListSize  (bool linkable);
  bool                getPeptideSeq       (int index, int start, int stop, char* str);
  bool                getPeptideSeq       (int index, int start, int stop, string& str);
  bool                getPeptideSeq       (kPeptide& p, string& str);
  kDB&                getProtein          (int index);
  int                 getProteinDBSize    ();

  //User Functions
  bool  buildDB(char* fname);                           //Reads FASTA file and populates vDB
  bool  buildPeptides(double min, double max, int mis); //Make peptide list within mass boundaries and miscleavages.
  void  combo();

private:
  
  //Data Members
  double  AA[128];            //Amino acid masses

  vector<kDB>      vDB;    //Entire FASTA database stored in memory
  vector<kPeptide> vPep;   //List of all peptides not linkable
  vector<kPeptide> vPepK;  //List of all peptides with a linkable K

  //Utility functions (for sorting)
  static int compareMass      (const void *p1, const void *p2);
  static int compareSequence  (const void *p1, const void *p2);

};

#endif
