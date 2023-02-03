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
  // HLT_Mu50 = false; HLT_PFMET120_PFMHT120_IDTight = false; HLT_PFHT500_PFMET100_PFMHT100_IDTight = false;
  // HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = false; HLT_MET105_IsoTrk50 = false;
  met = 0; metPhi = 0;
  pileupWeight = 0;
  nLeptons = 0;

  ZMass = -999.;
  ZPt = -999.;
  ZEta = -999.;
  ZPhi = -999.;
  ZleptonIndex1= 999;
  ZleptonIndex2 = 999;

  genZMass = -999.;
  genZPt = -999.;
  genZEta = -999.;
  genZPhi = -999.;
  genZleptonIndex1= 999;
  genZleptonIndex2 = 999;

  nHSCP = 0;
  nGenHSCP = 0;
  for( int i = 0; i < N_MAX_HSCP; i++ )
  {
   HSCP_pt[i] = -999.;
   HSCP_eta[i] = -999.;
   HSCP_phi[i] = -999.;
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
   HSCP_dZ_pv[i] = -999.;
   HSCP_dXY_pv[i] = -999.;
   HSCP_Ih[i] = -999.;
   HSCP_ErrorHisto_bin[i] = 0;
   HSCP_type[i] = 999;
   HSCP_match_genHSCP[i] = false;
   HSCP_match_genHSCP_minDeltaR[i] = -999.;
   HSCP_match_genHSCP_index[i] = -999;

   genHSCP_pt[i] = -999.;
   genHSCP_eta[i] = -999.;
   genHSCP_phi[i] = -999.;
   genHSCP_e[i] = -999.;
   genHSCP_type[i] = -999;
   genHSCP_pdgid[i] = -999;

  }
  for( int i = 0; i < N_MAX_LEPTONS; i++ )
  {
    lepE[i] = -999.;
    lepPt[i] = -999.;
    lepEta[i] = -999.;
    lepPhi[i] = -999.;
    lepPdgId[i] = -999;
    lepDZ[i] = -999.;
    lepLooseId[i] = false;
    lepTightId[i]= false;
    lepPassLooseIso[i]= false;
    lepPassTightIso[i]= false;
    //lepPassVTightIso[i]= false;
    //lepPassVVTightIso[i]= false;
  }
  for( int i = 0; i < NTriggersMAX; i++ )
  {
    HLTDecision[i] = false;
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



  tree_->SetBranchAddress("ZMass",      &ZMass);
  tree_->SetBranchAddress("ZPt",      &ZPt);
  tree_->SetBranchAddress("ZEta",      &ZEta);
  tree_->SetBranchAddress("ZPhi",      &ZPhi);
  tree_->SetBranchAddress("ZleptonIndex1",      &ZleptonIndex1);
  tree_->SetBranchAddress("ZleptonIndex2",      &ZleptonIndex2);



    tree_->SetBranchAddress("genZMass",      &genZMass);
    tree_->SetBranchAddress("genZPt",      &genZPt);
    tree_->SetBranchAddress("genZEta",      &genZEta);
    tree_->SetBranchAddress("genZPhi",      &genZPhi);
    tree_->SetBranchAddress("genZleptonIndex1",      &genZleptonIndex1);
    tree_->SetBranchAddress("genZleptonIndex2",      &genZleptonIndex2);
  // tree_->SetBranchAddress("HLT_Mu50",      &HLT_Mu50);
  // tree_->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight",      &HLT_PFMET120_PFMHT120_IDTight);
  // tree_->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight",      &HLT_PFHT500_PFMET100_PFMHT100_IDTight);
  // tree_->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",      &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
  // tree_->SetBranchAddress("HLT_MET105_IsoTrk50",      &HLT_MET105_IsoTrk50);
  tree_->SetBranchAddress("HLTDecision",      HLTDecision);


    tree_->SetBranchAddress("nLeptons",    &nLeptons);
    tree_->SetBranchAddress("lepE",        lepE);
    tree_->SetBranchAddress("lepPt",       lepPt);
    tree_->SetBranchAddress("lepEta",      lepEta);
    tree_->SetBranchAddress("lepPhi",      lepPhi);
    tree_->SetBranchAddress("lepPdgId",  lepPdgId);
    tree_->SetBranchAddress("lepDZ",     lepDZ);
    tree_->SetBranchAddress("lepLooseId", lepLooseId);
    tree_->SetBranchAddress("lepTightId", lepTightId);
    tree_->SetBranchAddress("lepPassLooseIso", lepPassLooseIso);
    tree_->SetBranchAddress("lepPassTightIso", lepPassTightIso);
    //tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);
    //tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);





    tree_->SetBranchAddress("nGenHSCP",      &nGenHSCP);
    tree_->SetBranchAddress("genHSCP_pt",      genHSCP_pt);
    tree_->SetBranchAddress("genHSCP_eta",      genHSCP_eta);
    tree_->SetBranchAddress("genHSCP_phi",      genHSCP_phi);
    tree_->SetBranchAddress("genHSCP_e",      genHSCP_e);
    tree_->SetBranchAddress("genHSCP_type",      genHSCP_type);
    tree_->SetBranchAddress("genHSCP_pdgid",      genHSCP_pdgid);


  tree_->SetBranchAddress("nHSCP",      &nHSCP);
  tree_->SetBranchAddress("HSCP_pt",      HSCP_pt);
  tree_->SetBranchAddress("HSCP_eta",      HSCP_eta);
  tree_->SetBranchAddress("HSCP_phi",      HSCP_phi);
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
  tree_->SetBranchAddress("HSCP_dZ_pv",      HSCP_dZ_pv);
  tree_->SetBranchAddress("HSCP_dXY_pv",      HSCP_dXY_pv);
  tree_->SetBranchAddress("HSCP_Ih",      HSCP_Ih);
  tree_->SetBranchAddress("HSCP_ErrorHisto_bin",      HSCP_ErrorHisto_bin);
  tree_->SetBranchAddress("HSCP_type",      HSCP_type);
  tree_->SetBranchAddress("HSCP_match_genHSCP",      HSCP_match_genHSCP);
  tree_->SetBranchAddress("HSCP_match_genHSCP_index",      HSCP_match_genHSCP_index);
  tree_->SetBranchAddress("HSCP_match_genHSCP_minDeltaR",      HSCP_match_genHSCP_minDeltaR);


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
  // tree_->Branch("HLT_Mu50",      &HLT_Mu50,     "HLT_Mu50/O");
  // tree_->Branch("HLT_PFMET120_PFMHT120_IDTight",      &HLT_PFMET120_PFMHT120_IDTight,     "HLT_PFMET120_PFMHT120_IDTight/O");      // event number
  // tree_->Branch("HLT_PFHT500_PFMET100_PFMHT100_IDTight",      &HLT_PFHT500_PFMET100_PFMHT100_IDTight,     "HLT_PFHT500_PFMET100_PFMHT100_IDTight/O");      // event number
  // tree_->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",      &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,     "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60/O");      // event number
  // tree_->Branch("HLT_MET105_IsoTrk50",      &HLT_MET105_IsoTrk50,     "HLT_MET105_IsoTrk50/O");      // event number

  tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[1201]/O"); //hardcoded

  tree_->Branch("ZMass",  &ZMass, "ZMass/F");
  tree_->Branch("ZPt",  &ZPt, "ZPt/F");
  tree_->Branch("ZEta",  &ZEta, "ZEta/F");
  tree_->Branch("ZPhi",  &ZPhi, "ZPhi/F");
  tree_->Branch("ZleptonIndex1",  &ZleptonIndex1, "ZleptonIndex1/I");
  tree_->Branch("ZleptonIndex2",  &ZleptonIndex2, "ZleptonIndex2/I");


    tree_->Branch("genZMass",  &genZMass, "genZMass/F");
    tree_->Branch("genZPt",  &genZPt, "genZPt/F");
    tree_->Branch("genZEta",  &genZEta, "genZEta/F");
    tree_->Branch("genZPhi",  &genZPhi, "genZPhi/F");
    tree_->Branch("genZleptonIndex1",  &genZleptonIndex1, "genZleptonIndex1/I");
    tree_->Branch("genZleptonIndex2",  &genZleptonIndex2, "genZleptonIndex2/I");

  //leptons
  tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
  tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
  tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
  tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
  tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
  tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");
  tree_->Branch("lepDZ",     lepDZ,     "lepDZ[nLeptons]/F");
  tree_->Branch("lepLooseId", lepLooseId, "lepLooseId[nLeptons]/O");
  tree_->Branch("lepTightId", lepTightId, "lepTightId[nLeptons]/O");
  tree_->Branch("lepPassLooseIso", lepPassLooseIso, "lepPassLooseIso[nLeptons]/O");
  tree_->Branch("lepPassTightIso", lepPassTightIso, "lepPassTightIso[nLeptons]/O");
  //tree_->Branch("lepPassVTightIso", lepPassVTightIso, "lepPassVTightIso[nLeptons]/O");
  //tree_->Branch("lepPassVVTightIso", lepPassVVTightIso, "lepPassVVTightIso[nLeptons]/O");




  tree_->Branch("nGenHSCP",      &nGenHSCP, "nGenHSCP/I");
  tree_->Branch("genHSCP_pt",      genHSCP_pt, "genHSCP_pt[nGenHSCP]/F");
  tree_->Branch("genHSCP_eta",      genHSCP_eta, "genHSCP_eta[nGenHSCP]/F");
  tree_->Branch("genHSCP_phi",      genHSCP_phi, "genHSCP_phi[nGenHSCP]/F");
  tree_->Branch("genHSCP_e",      genHSCP_e, "genHSCP_e[nGenHSCP]/F");
  tree_->Branch("genHSCP_type",      genHSCP_type, "genHSCP_type[nGenHSCP]/I");
  tree_->Branch("genHSCP_pdgid",      genHSCP_pdgid, "genHSCP_pdgid[nGenHSCP]/I");

  tree_->Branch("nHSCP",      &nHSCP, "nHSCP/I");
  tree_->Branch("HSCP_pt",      HSCP_pt, "HSCP_pt[nHSCP]/F");
  tree_->Branch("HSCP_eta",      HSCP_eta, "HSCP_eta[nHSCP]/F");
  tree_->Branch("HSCP_phi",      HSCP_phi, "HSCP_phi[nHSCP]/F");
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
  tree_->Branch("HSCP_dZ_pv",      HSCP_dZ_pv, "HSCP_dZ_pv[nHSCP]/F");
  tree_->Branch("HSCP_dXY_pv",      HSCP_dXY_pv, "HSCP_dXY_pv[nHSCP]/F");
  tree_->Branch("HSCP_Ih",      HSCP_Ih, "HSCP_Ih[nHSCP]/F");
  tree_->Branch("HSCP_ErrorHisto_bin",      HSCP_ErrorHisto_bin, "HSCP_ErrorHisto_bin[nHSCP]/I");
  tree_->Branch("HSCP_type",      HSCP_type, "HSCP_type[nHSCP]/I");
  tree_->Branch("HSCP_match_genHSCP",      HSCP_match_genHSCP, "HSCP_match_genHSCP[nHSCP]/O");
  tree_->Branch("HSCP_match_genHSCP_index",      HSCP_match_genHSCP_index, "HSCP_match_genHSCP_index[nHSCP]/I");
  tree_->Branch("HSCP_match_genHSCP_minDeltaR",      HSCP_match_genHSCP_minDeltaR, "HSCP_match_genHSCP_minDeltaR[nHSCP]/F");



};
