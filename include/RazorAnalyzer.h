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

        double deltaPhi(double phi1, double phi2);
      	double deltaR(double eta1, double phi1, double eta2, double phi2);
        TLorentzVector makeTLorentzVector(double pt, double eta, double phi, double energy);
      	TLorentzVector makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass);
};

#endif
