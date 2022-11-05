#include "RazorHelper.h"
#include "HSCPTree.h"
#include "assert.h"
#include "TTree.h"

// Constructor
HSCPTree::HSCPTree()
{
  InitVariables();
};
HSCPTree::~HSCPTree()
{
  if (f_) f_->Close();
};
void HSCPTree::InitVariables()
{
  runNum=0; lumiSec=0; evtNum=0; weight = 0;
  HLT_Mu50 = false; HLT_PFMET120_PFMHT120_IDTight = false; HLT_PFHT500_PFMET100_PFMHT100_IDTight = false;
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = false; HLT_MET105_IsoTrk50 = false;
  met = 0; metPhi = 0;
  nHSCP = 0;
  pileupWeight = 0;
  for( int i = 0; i < N_MAX_HSCP; i++ )
  {
   HSCP_pt[i] = -999.;
   HSCP_eta[i] = -999.;
   HSCP_nPixelHit[i] = -999;
   HSCP_nHits[i] = -999;
   HSCP_isHighPurity[i] =false;
   HSCP_fracValidHits[i] = -999.;
   HSCP_EOverP[i] = -999.;
   HSCP_chi2[i] = -999.;
   HSCP_nDof[i] = -999;
   HSCP_track_genTrackMiniIsoSumPt[i] = -999.;
   HSCP_pfMiniIso_relative[i] = -999.;
   HSCP_probQ[i] = -999.;
   HSCP_ptErr[i] = -999.;
   HSCP_ias_StripOnly[i] = -999.;
   HSCP_mass[i] = -999.;
   HSCP_dZ[i] = -999.;
   HSCP_dXY[i] = -999.;
  }


};

void HSCPTree::InitTree()
{
  assert(tree_);
  InitVariables();

  tree_->SetBranchAddress("runNum",      &runNum);
  tree_->SetBranchAddress("lumiSec",     &lumiSec);
  tree_->SetBranchAddress("evtNum",      &evtNum);
  tree_->SetBranchAddress("weight",      &weight);
  tree_->SetBranchAddress("pileupWeight",      &pileupWeight);
  tree_->SetBranchAddress("met",      &met);
  tree_->SetBranchAddress("metPhi",      &metPhi);

  tree_->SetBranchAddress("HLT_Mu50",      &HLT_Mu50);
  tree_->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight",      &HLT_PFMET120_PFMHT120_IDTight);
  tree_->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight",      &HLT_PFHT500_PFMET100_PFMHT100_IDTight);
  tree_->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",      &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
  tree_->SetBranchAddress("HLT_MET105_IsoTrk50",      &HLT_MET105_IsoTrk50);

  tree_->SetBranchAddress("nHSCP",      &nHSCP);
  tree_->SetBranchAddress("HSCP_pt",      HSCP_pt);
  tree_->SetBranchAddress("HSCP_eta",      HSCP_eta);
  tree_->SetBranchAddress("HSCP_nPixelHit",      HSCP_nPixelHit);
  tree_->SetBranchAddress("HSCP_nHits",      HSCP_nHits);
  tree_->SetBranchAddress("HSCP_isHighPurity",      HSCP_isHighPurity);
  tree_->SetBranchAddress("HSCP_fracValidHits",      HSCP_fracValidHits);
  tree_->SetBranchAddress("HSCP_EOverP",      HSCP_EOverP);
  tree_->SetBranchAddress("HSCP_chi2",      HSCP_chi2);
  tree_->SetBranchAddress("HSCP_nDof",      HSCP_nDof);
  tree_->SetBranchAddress("HSCP_track_genTrackMiniIsoSumPt",      HSCP_track_genTrackMiniIsoSumPt);
  tree_->SetBranchAddress("HSCP_pfMiniIso_relative",      HSCP_pfMiniIso_relative);
  tree_->SetBranchAddress("HSCP_probQ",      HSCP_probQ);
  tree_->SetBranchAddress("HSCP_ptErr",      HSCP_ptErr);
  tree_->SetBranchAddress("HSCP_ias_StripOnly",      HSCP_ias_StripOnly);
  tree_->SetBranchAddress("HSCP_mass",      HSCP_mass);
  tree_->SetBranchAddress("HSCP_dZ",      HSCP_dZ);
  tree_->SetBranchAddress("HSCP_dXY",      HSCP_dXY);

};
void HSCPTree::LoadTree(const char* file)
{
  f_ = TFile::Open(file);
  assert(f_);
  tree_ = dynamic_cast<TTree*>(f_->Get("tree"));
  InitTree();
  assert(tree_);
};

void HSCPTree::CreateTree()
{
  tree_ = new TTree("tree","tree");
  f_ = 0;

  tree_->Branch("runNum",      &runNum,     "runNum/i");
  tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");
  tree_->Branch("evtNum",      &evtNum,     "evtNum/i");
  tree_->Branch("weight",      &weight,     "weight/F");
  tree_->Branch("pileupWeight",      &pileupWeight,     "pileupWeight/F");


  tree_->Branch("met",      &met,     "met/F");
  tree_->Branch("metPhi",      &metPhi,     "metPhi/F");
  tree_->Branch("HLT_Mu50",      &HLT_Mu50,     "HLT_Mu50/O");
  tree_->Branch("HLT_PFMET120_PFMHT120_IDTight",      &HLT_PFMET120_PFMHT120_IDTight,     "HLT_PFMET120_PFMHT120_IDTight/O");      // event number
  tree_->Branch("HLT_PFHT500_PFMET100_PFMHT100_IDTight",      &HLT_PFHT500_PFMET100_PFMHT100_IDTight,     "HLT_PFHT500_PFMET100_PFMHT100_IDTight/O");      // event number
  tree_->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",      &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,     "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60/O");      // event number
  tree_->Branch("HLT_MET105_IsoTrk50",      &HLT_MET105_IsoTrk50,     "HLT_MET105_IsoTrk50/O");      // event number


  tree_->Branch("nHSCP",      &nHSCP, "nHSCP/I");
  tree_->Branch("HSCP_pt",      HSCP_pt, "HSCP_pt[nHSCP]/F");
  tree_->Branch("HSCP_eta",      HSCP_eta, "HSCP_eta[nHSCP]/F");
  tree_->Branch("HSCP_nPixelHit",      HSCP_nPixelHit, "HSCP_nPixelHit[nHSCP]/i");
  tree_->Branch("HSCP_nHits",      HSCP_nHits, "HSCP_nHits[nHSCP]/i");
  tree_->Branch("HSCP_isHighPurity",      HSCP_isHighPurity, "HSCP_isHighPurity[nHSCP]/O");
  tree_->Branch("HSCP_fracValidHits",      HSCP_fracValidHits, "HSCP_fracValidHits[nHSCP]/F");
  tree_->Branch("HSCP_EOverP",      HSCP_EOverP, "HSCP_EOverP[nHSCP]/F");
  tree_->Branch("HSCP_chi2",      HSCP_chi2, "HSCP_chi2[nHSCP]/F");
  tree_->Branch("HSCP_nDof",      HSCP_nDof, "HSCP_nDof[nHSCP]/i");
  tree_->Branch("HSCP_track_genTrackMiniIsoSumPt",      HSCP_track_genTrackMiniIsoSumPt, "HSCP_track_genTrackMiniIsoSumPt[nHSCP]/F");
  tree_->Branch("HSCP_pfMiniIso_relative",      HSCP_pfMiniIso_relative, "HSCP_pfMiniIso_relative[nHSCP]/F");
  tree_->Branch("HSCP_probQ",      HSCP_probQ, "HSCP_probQ[nHSCP]/F");
  tree_->Branch("HSCP_ptErr",      HSCP_ptErr, "HSCP_ptErr[nHSCP]/F");
  tree_->Branch("HSCP_ias_StripOnly",      HSCP_ias_StripOnly, "HSCP_ias_StripOnly[nHSCP]/F");
  tree_->Branch("HSCP_mass",      HSCP_mass, "HSCP_mass[nHSCP]/F");
  tree_->Branch("HSCP_dZ",      HSCP_dZ, "HSCP_dZ[nHSCP]/F");
  tree_->Branch("HSCP_dXY",      HSCP_dXY, "HSCP_dXY[nHSCP]/F");


};
