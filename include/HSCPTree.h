// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef HSCPTree_H
#define HSCPTree_H

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define N_MAX_CSC 200
#define N_MAX_CSCRECHITS 5000
#define N_MAX_DTRECHITS 20000
#define NTriggersMAX 1201 // Number of trigger in the .dat file
#define N_CSC_CUT 20
#define JET_PT_CUT 10
#define MUON_PT_CUT 20
#define N_MAX_GPARTICLES 5000
#define N_MAX_HSCP 10

#include <iostream>
#include <string>
#include <sys/stat.h>
#include "assert.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TTree.h"
#include "RazorAnalyzer.h"
#include "RazorHelper.h"

class HSCPTree
{

public:
  HSCPTree();
  ~HSCPTree();
  // HSCPTree::HSCPTree()
  // {
  //   InitVariables();
  // };
  // HSCPTree::~HSCPTree()
  // {
  //   if (f_) f_->Close();
  // };
  TTree *tree_;
  TFile *f_;

  UInt_t  runNum, lumiSec, evtNum;

  float weight, pileupWeight;
  bool HLT_Mu50, HLT_PFMET120_PFMHT120_IDTight, HLT_PFHT500_PFMET100_PFMHT100_IDTight, HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, HLT_MET105_IsoTrk50;

  float met, metPhi;

  // HSCP
  int nHSCP;
  float HSCP_pt[N_MAX_HSCP];
  float HSCP_eta[N_MAX_HSCP];
  UInt_t HSCP_nPixelHit[N_MAX_HSCP];
  UInt_t HSCP_nHits[N_MAX_HSCP];
  bool HSCP_isHighPurity[N_MAX_HSCP];
  float HSCP_fracValidHits[N_MAX_HSCP];
  float HSCP_EOverP[N_MAX_HSCP];
  float HSCP_chi2[N_MAX_HSCP];
  UInt_t HSCP_nDof[N_MAX_HSCP];
  float HSCP_track_genTrackMiniIsoSumPt[N_MAX_HSCP];
  float HSCP_pfMiniIso_relative[N_MAX_HSCP];
  float HSCP_probQ[N_MAX_HSCP];
  float HSCP_ptErr[N_MAX_HSCP];
  float HSCP_ias_StripOnly[N_MAX_HSCP];
  float HSCP_mass[N_MAX_HSCP];
  float HSCP_dZ[N_MAX_HSCP];
  float HSCP_dXY[N_MAX_HSCP];


  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
