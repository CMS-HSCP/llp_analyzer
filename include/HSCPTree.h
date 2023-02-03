// Class to manage files for b-tag scale factors, lepton scale factors, pileup weights, and other information

#ifndef HSCPTree_H
#define HSCPTree_H

#define N_MAX_LEPTONS 20
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

  TTree *tree_;
  TFile *f_;

  UInt_t  runNum, lumiSec, evtNum;

  float weight, pileupWeight;
  bool HLT_Mu50, HLT_PFMET120_PFMHT120_IDTight, HLT_PFHT500_PFMET100_PFMHT100_IDTight, HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, HLT_MET105_IsoTrk50;

  float met, metPhi;
  int nLeptons;
  float lepE[N_MAX_LEPTONS];
  float lepPt[N_MAX_LEPTONS];
  float lepEta[N_MAX_LEPTONS];
  float lepPhi[N_MAX_LEPTONS];
  int lepPdgId[N_MAX_LEPTONS];
  float lepDZ[N_MAX_LEPTONS];
  bool lepLooseId[N_MAX_LEPTONS];
  bool lepTightId[N_MAX_LEPTONS];
  bool lepPassLooseIso[N_MAX_LEPTONS];
  bool lepPassTightIso[N_MAX_LEPTONS];
  //bool lepPassVTightIso[N_MAX_LEPTONS];
  //bool lepPassVVTightIso[N_MAX_LEPTONS];

  float ZMass;
  float ZPt;
  float ZEta;
  float ZPhi;
  int ZleptonIndex1;
  int ZleptonIndex2;


  float genZMass;
  float genZPt;
  float genZEta;
  float genZPhi;
  int genZleptonIndex1;
  int genZleptonIndex2;

  int nGenHSCP;
  int genHSCP_pdgid[N_MAX_HSCP];
  int genHSCP_type[N_MAX_HSCP];
  float genHSCP_pt[N_MAX_HSCP];
  float genHSCP_eta[N_MAX_HSCP];
  float genHSCP_phi[N_MAX_HSCP];
  float genHSCP_e[N_MAX_HSCP];

  // HSCP
  int nHSCP;
  float HSCP_pt[N_MAX_HSCP];
  float HSCP_eta[N_MAX_HSCP];
  float HSCP_phi[N_MAX_HSCP];
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
  float HSCP_dZ_pv[N_MAX_HSCP];
  float HSCP_dXY_pv[N_MAX_HSCP];
  float HSCP_Ih[N_MAX_HSCP];
  int HSCP_ErrorHisto_bin[N_MAX_HSCP];
  int HSCP_type[N_MAX_HSCP];
  bool HSCP_match_genHSCP[N_MAX_HSCP];
  float HSCP_match_genHSCP_minDeltaR[N_MAX_HSCP];
  int HSCP_match_genHSCP_index[N_MAX_HSCP];
  bool HLTDecision[NTriggersMAX];
  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
