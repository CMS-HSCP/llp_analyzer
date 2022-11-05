//Class for analyzing ntuples produced by the llp_ntupler framework
//
//Author: Cristian Pena & Si Xie

#ifndef RazorAnalyzer_h
#define RazorAnalyzer_h

#include "HscpCandidates.h" //This is a MakeClass of the llp tree in the ntuple to be analyzed


//ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include "TRandom3.h"

//C++ includes
#include <map>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

class RazorAnalyzer: public HscpCandidates {
    public :
        RazorAnalyzer(TTree *tree=0);
        virtual ~RazorAnalyzer();

        //------ LIST OF ANALYSES ------//
        virtual void Analyze(bool isData, int option, string outputFileName, string label);



        enum RazorBox { //boxes for razor inclusive analysis
	  MuEle = 0,
	  MuMu = 1,
	  EleEle = 2,
	  MuSixJet = 3,
	  MuFourJet = 4,
	  MuJet = 5,
	  EleSixJet = 6,
	  EleFourJet = 7,
	  EleJet = 8,
	  LooseLeptonSixJet = 9,
	  LooseLeptonFourJet = 10,
	  SixJet = 11,
	  FourJet = 12,
	  LooseLeptonDiJet = 13,
	  DiJet = 14,
	  TwoBJet = 15,
	  OneBJet = 16,
	  ZeroBJet = 17,
	  MuMultiJet = 18,
	  EleMultiJet = 19,
	  LooseLeptonMultiJet = 20,
	  MultiJet = 21,
	  NONE = 999
        };
};

#endif
