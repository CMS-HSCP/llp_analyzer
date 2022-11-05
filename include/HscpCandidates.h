//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov  3 11:50:40 2022 by ROOT version 6.14/09
// from TTree HscpCandidates/HscpCandidates
// found on file: root://cmsxrootd.fnal.gov//store/group/lpchscp/ntuples/V1p0/MC_UL18/HSCPstopOnlyNeutral_M-100_TuneCP5_13TeV-pythia8/HSCPNtupler_2018_HSCPstopOnlyNeutral_M-100_wProbQ_v1p0_v1/221015_125013/0000/Histos_5.root
//////////////////////////////////////////////////////////

#ifndef HscpCandidates_h
#define HscpCandidates_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
using std::vector;

class HscpCandidates {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          Trig;
   UInt_t          Run;
   ULong64_t       Event;
   UInt_t          Lumi;
   UInt_t          PileUp;
   UInt_t          nofVtx;
   UInt_t          Hscp;
   UInt_t          nmuons;
   UInt_t          njets;
   Float_t         Weight;
   Float_t         GeneratorWeight;
   Float_t         GeneratorBinningValues;
   Bool_t          HLT_Mu50;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight;
   Bool_t          HLT_PFHT500_PFMET100_PFMHT100_IDTight;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
   Bool_t          HLT_MET105_IsoTrk50;
   Float_t         RecoCaloMET;
   Float_t         RecoCaloMET_phi;
   Float_t         RecoCaloMET_sigf;
   Float_t         RecoPFMET;
   Float_t         RecoPFMET_phi;
   Float_t         RecoPFMET_sigf;
   Float_t         RecoPFMHT;
   Float_t         HLTCaloMET;
   Float_t         HLTCaloMET_phi;
   Float_t         HLTCaloMET_sigf;
   Float_t         HLTCaloMETClean;
   Float_t         HLTCaloMETClean_phi;
   Float_t         HLTCaloMETClean_sigf;
   Float_t         HLTCaloMHT;
   Float_t         HLTCaloMHT_phi;
   Float_t         HLTCaloMHT_sigf;
   Float_t         HLTPFMET;
   Float_t         HLTPFMET_phi;
   Float_t         HLTPFMET_sigf;
   Float_t         HLTPFMHT;
   Float_t         HLTPFMHT_phi;
   Float_t         HLTPFMHT_sigf;
   Bool_t          matchedMuonWasFound;
   Float_t         Muon1_Pt;
   Float_t         Muon1_eta;
   Float_t         Muon1_phi;
   Float_t         Muon2_Pt;
   Float_t         Muon2_eta;
   Float_t         Muon2_phi;
   vector<float>   *Jet_pt;
   vector<float>   *Jet_eta;
   vector<float>   *Jet_phi;
   vector<float>   *Jet_mass;
   vector<float>   *Jet_energy;
   vector<float>   *Jet_pdgId;
   vector<float>   *Jet_et;
   vector<float>   *Jet_chargedEmEnergyFraction;
   vector<float>   *Jet_neutralEmEnergyFraction;
   vector<float>   *mT;
   vector<bool>    *passCutPt55;
   vector<bool>    *passPreselection;
   vector<bool>    *passSelection;
   vector<float>   *Charge;
   vector<float>   *Pt;
   vector<float>   *PtErr;
   vector<float>   *Ias;
   vector<float>   *Ias_noTIBnoTIDno3TEC;
   vector<float>   *Ias_PixelOnly;
   vector<float>   *Ias_StripOnly;
   vector<float>   *Ias_PixelOnly_noL1;
   vector<float>   *Ih;
   vector<float>   *Ick;
   vector<float>   *Fmip;
   vector<float>   *ProbXY;
   vector<float>   *ProbXY_noL1;
   vector<float>   *ProbQ;
   vector<float>   *ProbQ_noL1;
   vector<float>   *Ndof;
   vector<float>   *Chi2;
   vector<int>     *QualityMask;
   vector<bool>    *isHighPurity;
   vector<float>   *EoverP;
   vector<bool>    *isMuon;
   vector<bool>    *isPhoton;
   vector<bool>    *isElectron;
   vector<bool>    *isChHadron;
   vector<bool>    *isNeutHadron;
   vector<bool>    *isPfTrack;
   vector<bool>    *isUndefined;
   vector<float>   *ECAL_energy;
   vector<float>   *HCAL_energy;
   vector<float>   *TOF;
   vector<float>   *TOFErr;
   vector<unsigned int> *TOF_ndof;
   vector<float>   *DTTOF;
   vector<float>   *DTTOFErr;
   vector<unsigned int> *DTTOF_ndof;
   vector<float>   *CSCTOF;
   vector<float>   *CSCTOFErr;
   vector<unsigned int> *CSCTOF_ndof;
   vector<float>   *Mass;
   vector<float>   *MassErr;
   vector<float>   *dZ;
   vector<float>   *dXY;
   vector<float>   *dR;
   vector<float>   *p;
   vector<float>   *eta;
   vector<float>   *phi;
   vector<unsigned int> *NOH;
   vector<unsigned int> *NOPH;
   vector<float>   *FOVH;
   vector<unsigned int> *NOMH;
   vector<float>   *FOVHD;
   vector<unsigned int> *NOM;
   vector<float>   *matchTrigMuon_minDeltaR;
   vector<float>   *matchTrigMuon_pT;
   vector<float>   *iso_TK;
   vector<float>   *iso_ECAL;
   vector<float>   *iso_HCAL;
   vector<float>   *track_genTrackMiniIsoSumPt;
   vector<float>   *PFMiniIso_relative;
   vector<float>   *PFMiniIso_wMuon_relative;
   vector<float>   *TrackPFIsolationR005_sumChargedHadronPt;
   vector<float>   *TrackPFIsolationR005_sumNeutralHadronPt;
   vector<float>   *TrackPFIsolationR005_sumPhotonPt;
   vector<float>   *TrackPFIsolationR005_sumPUPt;
   vector<float>   *TrackPFIsolationR01_sumChargedHadronPt;
   vector<float>   *TrackPFIsolationR01_sumNeutralHadronPt;
   vector<float>   *TrackPFIsolationR01_sumPhotonPt;
   vector<float>   *TrackPFIsolationR01_sumPUPt;
   vector<float>   *TrackPFIsolationR03_sumChargedHadronPt;
   vector<float>   *TrackPFIsolationR03_sumNeutralHadronPt;
   vector<float>   *TrackPFIsolationR03_sumPhotonPt;
   vector<float>   *TrackPFIsolationR03_sumPUPt;
   vector<float>   *TrackPFIsolationR05_sumChargedHadronPt;
   vector<float>   *TrackPFIsolationR05_sumNeutralHadronPt;
   vector<float>   *TrackPFIsolationR05_sumPhotonPt;
   vector<float>   *TrackPFIsolationR05_sumPUPt;
   vector<float>   *MuonPFIsolationR03_sumChargedHadronPt;
   vector<float>   *MuonPFIsolationR03_sumNeutralHadronPt;
   vector<float>   *MuonPFIsolationR03_sumPhotonPt;
   vector<float>   *MuonPFIsolationR03_sumPUPt;
   vector<float>   *Ih_noL1;
   vector<float>   *Ih_15drop;
   vector<float>   *Ih_StripOnly;
   vector<float>   *Ih_StripOnly_15drop;
   vector<float>   *Ih_PixelOnly_noL1;
   vector<float>   *Ih_SaturationCorrectionFromFits;
   vector<vector<float> > *clust_charge;
   vector<vector<float> > *clust_pathlength;
   vector<vector<unsigned int> > *clust_nstrip;
   vector<vector<bool> > *clust_sat254;
   vector<vector<bool> > *clust_sat255;
   vector<vector<unsigned int> > *clust_detid;
   vector<vector<bool> > *clust_isStrip;
   vector<vector<bool> > *clust_isPixel;
   vector<float>   *GenId;
   vector<float>   *GenCharge;
   vector<float>   *GenMass;
   vector<float>   *GenPt;
   vector<float>   *GenEta;
   vector<float>   *GenPhi;

   // List of branches
   TBranch        *b_Trig;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Lumi;   //!
   TBranch        *b_PileUp;   //!
   TBranch        *b_nofVtx;   //!
   TBranch        *b_Hscp;   //!
   TBranch        *b_nmuons;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_GeneratorWeight;   //!
   TBranch        *b_GeneratorBinningValues;   //!
   TBranch        *b_HLT_Mu50;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight;   //!
   TBranch        *b_HLT_PFHT500_PFMET100_PFMHT100_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_MET105_IsoTrk50;   //!
   TBranch        *b_RecoCaloMET;   //!
   TBranch        *b_RecoCaloMET_phi;   //!
   TBranch        *b_RecoCaloMET_sigf;   //!
   TBranch        *b_RecoPFMET;   //!
   TBranch        *b_RecoPFMET_phi;   //!
   TBranch        *b_RecoPFMET_sigf;   //!
   TBranch        *b_RecoPFMHT;   //!
   TBranch        *b_HLTCaloMET;   //!
   TBranch        *b_HLTCaloMET_phi;   //!
   TBranch        *b_HLTCaloMET_sigf;   //!
   TBranch        *b_HLTCaloMETClean;   //!
   TBranch        *b_HLTCaloMETClean_phi;   //!
   TBranch        *b_HLTCaloMETClean_sigf;   //!
   TBranch        *b_HLTCaloMHT;   //!
   TBranch        *b_HLTCaloMHT_phi;   //!
   TBranch        *b_HLTCaloMHT_sigf;   //!
   TBranch        *b_HLTPFMET;   //!
   TBranch        *b_HLTPFMET_phi;   //!
   TBranch        *b_HLTPFMET_sigf;   //!
   TBranch        *b_HLTPFMHT;   //!
   TBranch        *b_HLTPFMHT_phi;   //!
   TBranch        *b_HLTPFMHT_sigf;   //!
   TBranch        *b_matchedMuonWasFound;   //!
   TBranch        *b_Muon1_Pt;   //!
   TBranch        *b_Muon1_eta;   //!
   TBranch        *b_Muon1_phi;   //!
   TBranch        *b_Muon2_Pt;   //!
   TBranch        *b_Muon2_eta;   //!
   TBranch        *b_Muon2_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_pdgId;   //!
   TBranch        *b_Jet_et;   //!
   TBranch        *b_Jet_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_neutralEmEnergyFraction;   //!
   TBranch        *b_mT;   //!
   TBranch        *b_passCutPt55;   //!
   TBranch        *b_passPreselection;   //!
   TBranch        *b_passSelection;   //!
   TBranch        *b_Charge;   //!
   TBranch        *b_Pt;   //!
   TBranch        *b_PtErr;   //!
   TBranch        *b_Ias;   //!
   TBranch        *b_Ias_noTIBnoTIDno3TEC;   //!
   TBranch        *b_Ias_PixelOnly;   //!
   TBranch        *b_Ias_StripOnly;   //!
   TBranch        *b_Ias_PixelOnly_noL1;   //!
   TBranch        *b_Ih;   //!
   TBranch        *b_Ick;   //!
   TBranch        *b_Fmip;   //!
   TBranch        *b_ProbXY;   //!
   TBranch        *b_ProbXY_noL1;   //!
   TBranch        *b_ProbQ;   //!
   TBranch        *b_ProbQ_noL1;   //!
   TBranch        *b_Ndof;   //!
   TBranch        *b_Chi2;   //!
   TBranch        *b_QualityMask;   //!
   TBranch        *b_isHighPurity;   //!
   TBranch        *b_EoverP;   //!
   TBranch        *b_isMuon;   //!
   TBranch        *b_isPhoton;   //!
   TBranch        *b_isElectron;   //!
   TBranch        *b_isChHadron;   //!
   TBranch        *b_isNeutHadron;   //!
   TBranch        *b_isPfTrack;   //!
   TBranch        *b_isUndefined;   //!
   TBranch        *b_ECAL_energy;   //!
   TBranch        *b_HCAL_energy;   //!
   TBranch        *b_TOF;   //!
   TBranch        *b_TOFErr;   //!
   TBranch        *b_TOF_ndof;   //!
   TBranch        *b_DTTOF;   //!
   TBranch        *b_DTTOFErr;   //!
   TBranch        *b_DTTOF_ndof;   //!
   TBranch        *b_CSCTOF;   //!
   TBranch        *b_CSCTOFErr;   //!
   TBranch        *b_CSCTOF_ndof;   //!
   TBranch        *b_Mass;   //!
   TBranch        *b_MassErr;   //!
   TBranch        *b_dZ;   //!
   TBranch        *b_dXY;   //!
   TBranch        *b_dR;   //!
   TBranch        *b_p;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_NOH;   //!
   TBranch        *b_NOPH;   //!
   TBranch        *b_FOVH;   //!
   TBranch        *b_NOMH;   //!
   TBranch        *b_FOVHD;   //!
   TBranch        *b_NOM;   //!
   TBranch        *b_matchTrigMuon_minDeltaR;   //!
   TBranch        *b_matchTrigMuon_pT;   //!
   TBranch        *b_iso_TK;   //!
   TBranch        *b_iso_ECAL;   //!
   TBranch        *b_iso_HCAL;   //!
   TBranch        *b_track_genTrackMiniIsoSumPt;   //!
   TBranch        *b_PFMiniIso_relative;   //!
   TBranch        *b_PFMiniIso_wMuon_relative;   //!
   TBranch        *b_TrackPFIsolationR005_sumChargedHadronPt;   //!
   TBranch        *b_TrackPFIsolationR005_sumNeutralHadronPt;   //!
   TBranch        *b_TrackPFIsolationR005_sumPhotonPt;   //!
   TBranch        *b_TrackPFIsolationR005_sumPUPt;   //!
   TBranch        *b_TrackPFIsolationR01_sumChargedHadronPt;   //!
   TBranch        *b_TrackPFIsolationR01_sumNeutralHadronPt;   //!
   TBranch        *b_TrackPFIsolationR01_sumPhotonPt;   //!
   TBranch        *b_TrackPFIsolationR01_sumPUPt;   //!
   TBranch        *b_TrackPFIsolationR03_sumChargedHadronPt;   //!
   TBranch        *b_TrackPFIsolationR03_sumNeutralHadronPt;   //!
   TBranch        *b_TrackPFIsolationR03_sumPhotonPt;   //!
   TBranch        *b_TrackPFIsolationR03_sumPUPt;   //!
   TBranch        *b_TrackPFIsolationR05_sumChargedHadronPt;   //!
   TBranch        *b_TrackPFIsolationR05_sumNeutralHadronPt;   //!
   TBranch        *b_TrackPFIsolationR05_sumPhotonPt;   //!
   TBranch        *b_TrackPFIsolationR05_sumPUPt;   //!
   TBranch        *b_MuonPFIsolationR03_sumChargedHadronPt;   //!
   TBranch        *b_MuonPFIsolationR03_sumNeutralHadronPt;   //!
   TBranch        *b_MuonPFIsolationR03_sumPhotonPt;   //!
   TBranch        *b_MuonPFIsolationR03_sumPUPt;   //!
   TBranch        *b_Ih_noL1;   //!
   TBranch        *b_Ih_15drop;   //!
   TBranch        *b_Ih_StripOnly;   //!
   TBranch        *b_Ih_StripOnly_15drop;   //!
   TBranch        *b_Ih_PixelOnly_noL1;   //!
   TBranch        *b_Ih_SaturationCorrectionFromFits;   //!
   TBranch        *b_clust_charge;   //!
   TBranch        *b_clust_pathlength;   //!
   TBranch        *b_clust_nstrip;   //!
   TBranch        *b_clust_sat254;   //!
   TBranch        *b_clust_sat255;   //!
   TBranch        *b_clust_detid;   //!
   TBranch        *b_clust_isStrip;   //!
   TBranch        *b_clust_isPixel;   //!
   TBranch        *b_GenId;   //!
   TBranch        *b_GenCharge;   //!
   TBranch        *b_GenMass;   //!
   TBranch        *b_GenPt;   //!
   TBranch        *b_GenEta;   //!
   TBranch        *b_GenPhi;   //!

   HscpCandidates(TTree *tree=0);
   virtual ~HscpCandidates();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HscpCandidates_cxx
HscpCandidates::HscpCandidates(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmsxrootd.fnal.gov//store/group/lpchscp/ntuples/V1p0/MC_UL18/HSCPstopOnlyNeutral_M-100_TuneCP5_13TeV-pythia8/HSCPNtupler_2018_HSCPstopOnlyNeutral_M-100_wProbQ_v1p0_v1/221015_125013/0000/Histos_5.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://cmsxrootd.fnal.gov//store/group/lpchscp/ntuples/V1p0/MC_UL18/HSCPstopOnlyNeutral_M-100_TuneCP5_13TeV-pythia8/HSCPNtupler_2018_HSCPstopOnlyNeutral_M-100_wProbQ_v1p0_v1/221015_125013/0000/Histos_5.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://cmsxrootd.fnal.gov//store/group/lpchscp/ntuples/V1p0/MC_UL18/HSCPstopOnlyNeutral_M-100_TuneCP5_13TeV-pythia8/HSCPNtupler_2018_HSCPstopOnlyNeutral_M-100_wProbQ_v1p0_v1/221015_125013/0000/Histos_5.root:/HSCParticleAnalyzer/BaseName");
      dir->GetObject("HscpCandidates",tree);

   }
   Init(tree);
}

HscpCandidates::~HscpCandidates()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HscpCandidates::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HscpCandidates::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HscpCandidates::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Jet_pt = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_mass = 0;
   Jet_energy = 0;
   Jet_pdgId = 0;
   Jet_et = 0;
   Jet_chargedEmEnergyFraction = 0;
   Jet_neutralEmEnergyFraction = 0;
   mT = 0;
   passCutPt55 = 0;
   passPreselection = 0;
   passSelection = 0;
   Charge = 0;
   Pt = 0;
   PtErr = 0;
   Ias = 0;
   Ias_noTIBnoTIDno3TEC = 0;
   Ias_PixelOnly = 0;
   Ias_StripOnly = 0;
   Ias_PixelOnly_noL1 = 0;
   Ih = 0;
   Ick = 0;
   Fmip = 0;
   ProbXY = 0;
   ProbXY_noL1 = 0;
   ProbQ = 0;
   ProbQ_noL1 = 0;
   Ndof = 0;
   Chi2 = 0;
   QualityMask = 0;
   isHighPurity = 0;
   EoverP = 0;
   isMuon = 0;
   isPhoton = 0;
   isElectron = 0;
   isChHadron = 0;
   isNeutHadron = 0;
   isPfTrack = 0;
   isUndefined = 0;
   ECAL_energy = 0;
   HCAL_energy = 0;
   TOF = 0;
   TOFErr = 0;
   TOF_ndof = 0;
   DTTOF = 0;
   DTTOFErr = 0;
   DTTOF_ndof = 0;
   CSCTOF = 0;
   CSCTOFErr = 0;
   CSCTOF_ndof = 0;
   Mass = 0;
   MassErr = 0;
   dZ = 0;
   dXY = 0;
   dR = 0;
   p = 0;
   eta = 0;
   phi = 0;
   NOH = 0;
   NOPH = 0;
   FOVH = 0;
   NOMH = 0;
   FOVHD = 0;
   NOM = 0;
   matchTrigMuon_minDeltaR = 0;
   matchTrigMuon_pT = 0;
   iso_TK = 0;
   iso_ECAL = 0;
   iso_HCAL = 0;
   track_genTrackMiniIsoSumPt = 0;
   PFMiniIso_relative = 0;
   PFMiniIso_wMuon_relative = 0;
   TrackPFIsolationR005_sumChargedHadronPt = 0;
   TrackPFIsolationR005_sumNeutralHadronPt = 0;
   TrackPFIsolationR005_sumPhotonPt = 0;
   TrackPFIsolationR005_sumPUPt = 0;
   TrackPFIsolationR01_sumChargedHadronPt = 0;
   TrackPFIsolationR01_sumNeutralHadronPt = 0;
   TrackPFIsolationR01_sumPhotonPt = 0;
   TrackPFIsolationR01_sumPUPt = 0;
   TrackPFIsolationR03_sumChargedHadronPt = 0;
   TrackPFIsolationR03_sumNeutralHadronPt = 0;
   TrackPFIsolationR03_sumPhotonPt = 0;
   TrackPFIsolationR03_sumPUPt = 0;
   TrackPFIsolationR05_sumChargedHadronPt = 0;
   TrackPFIsolationR05_sumNeutralHadronPt = 0;
   TrackPFIsolationR05_sumPhotonPt = 0;
   TrackPFIsolationR05_sumPUPt = 0;
   MuonPFIsolationR03_sumChargedHadronPt = 0;
   MuonPFIsolationR03_sumNeutralHadronPt = 0;
   MuonPFIsolationR03_sumPhotonPt = 0;
   MuonPFIsolationR03_sumPUPt = 0;
   Ih_noL1 = 0;
   Ih_15drop = 0;
   Ih_StripOnly = 0;
   Ih_StripOnly_15drop = 0;
   Ih_PixelOnly_noL1 = 0;
   Ih_SaturationCorrectionFromFits = 0;
   clust_charge = 0;
   clust_pathlength = 0;
   clust_nstrip = 0;
   clust_sat254 = 0;
   clust_sat255 = 0;
   clust_detid = 0;
   clust_isStrip = 0;
   clust_isPixel = 0;
   GenId = 0;
   GenCharge = 0;
   GenMass = 0;
   GenPt = 0;
   GenEta = 0;
   GenPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Trig", &Trig, &b_Trig);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
   fChain->SetBranchAddress("PileUp", &PileUp, &b_PileUp);
   fChain->SetBranchAddress("nofVtx", &nofVtx, &b_nofVtx);
   fChain->SetBranchAddress("Hscp", &Hscp, &b_Hscp);
   fChain->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("GeneratorWeight", &GeneratorWeight, &b_GeneratorWeight);
   fChain->SetBranchAddress("GeneratorBinningValues", &GeneratorBinningValues, &b_GeneratorBinningValues);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight", &HLT_PFHT500_PFMET100_PFMHT100_IDTight, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_MET105_IsoTrk50", &HLT_MET105_IsoTrk50, &b_HLT_MET105_IsoTrk50);
   fChain->SetBranchAddress("RecoCaloMET", &RecoCaloMET, &b_RecoCaloMET);
   fChain->SetBranchAddress("RecoCaloMET_phi", &RecoCaloMET_phi, &b_RecoCaloMET_phi);
   fChain->SetBranchAddress("RecoCaloMET_sigf", &RecoCaloMET_sigf, &b_RecoCaloMET_sigf);
   fChain->SetBranchAddress("RecoPFMET", &RecoPFMET, &b_RecoPFMET);
   fChain->SetBranchAddress("RecoPFMET_phi", &RecoPFMET_phi, &b_RecoPFMET_phi);
   fChain->SetBranchAddress("RecoPFMET_sigf", &RecoPFMET_sigf, &b_RecoPFMET_sigf);
   fChain->SetBranchAddress("RecoPFMHT", &RecoPFMHT, &b_RecoPFMHT);
   fChain->SetBranchAddress("HLTCaloMET", &HLTCaloMET, &b_HLTCaloMET);
   fChain->SetBranchAddress("HLTCaloMET_phi", &HLTCaloMET_phi, &b_HLTCaloMET_phi);
   fChain->SetBranchAddress("HLTCaloMET_sigf", &HLTCaloMET_sigf, &b_HLTCaloMET_sigf);
   fChain->SetBranchAddress("HLTCaloMETClean", &HLTCaloMETClean, &b_HLTCaloMETClean);
   fChain->SetBranchAddress("HLTCaloMETClean_phi", &HLTCaloMETClean_phi, &b_HLTCaloMETClean_phi);
   fChain->SetBranchAddress("HLTCaloMETClean_sigf", &HLTCaloMETClean_sigf, &b_HLTCaloMETClean_sigf);
   fChain->SetBranchAddress("HLTCaloMHT", &HLTCaloMHT, &b_HLTCaloMHT);
   fChain->SetBranchAddress("HLTCaloMHT_phi", &HLTCaloMHT_phi, &b_HLTCaloMHT_phi);
   fChain->SetBranchAddress("HLTCaloMHT_sigf", &HLTCaloMHT_sigf, &b_HLTCaloMHT_sigf);
   fChain->SetBranchAddress("HLTPFMET", &HLTPFMET, &b_HLTPFMET);
   fChain->SetBranchAddress("HLTPFMET_phi", &HLTPFMET_phi, &b_HLTPFMET_phi);
   fChain->SetBranchAddress("HLTPFMET_sigf", &HLTPFMET_sigf, &b_HLTPFMET_sigf);
   fChain->SetBranchAddress("HLTPFMHT", &HLTPFMHT, &b_HLTPFMHT);
   fChain->SetBranchAddress("HLTPFMHT_phi", &HLTPFMHT_phi, &b_HLTPFMHT_phi);
   fChain->SetBranchAddress("HLTPFMHT_sigf", &HLTPFMHT_sigf, &b_HLTPFMHT_sigf);
   fChain->SetBranchAddress("matchedMuonWasFound", &matchedMuonWasFound, &b_matchedMuonWasFound);
   fChain->SetBranchAddress("Muon1_Pt", &Muon1_Pt, &b_Muon1_Pt);
   fChain->SetBranchAddress("Muon1_eta", &Muon1_eta, &b_Muon1_eta);
   fChain->SetBranchAddress("Muon1_phi", &Muon1_phi, &b_Muon1_phi);
   fChain->SetBranchAddress("Muon2_Pt", &Muon2_Pt, &b_Muon2_Pt);
   fChain->SetBranchAddress("Muon2_eta", &Muon2_eta, &b_Muon2_eta);
   fChain->SetBranchAddress("Muon2_phi", &Muon2_phi, &b_Muon2_phi);
   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
   fChain->SetBranchAddress("Jet_pdgId", &Jet_pdgId, &b_Jet_pdgId);
   fChain->SetBranchAddress("Jet_et", &Jet_et, &b_Jet_et);
   fChain->SetBranchAddress("Jet_chargedEmEnergyFraction", &Jet_chargedEmEnergyFraction, &b_Jet_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jet_neutralEmEnergyFraction", &Jet_neutralEmEnergyFraction, &b_Jet_neutralEmEnergyFraction);
   fChain->SetBranchAddress("mT", &mT, &b_mT);
   fChain->SetBranchAddress("passCutPt55", &passCutPt55, &b_passCutPt55);
   fChain->SetBranchAddress("passPreselection", &passPreselection, &b_passPreselection);
   fChain->SetBranchAddress("passSelection", &passSelection, &b_passSelection);
   fChain->SetBranchAddress("Charge", &Charge, &b_Charge);
   fChain->SetBranchAddress("Pt", &Pt, &b_Pt);
   fChain->SetBranchAddress("PtErr", &PtErr, &b_PtErr);
   fChain->SetBranchAddress("Ias", &Ias, &b_Ias);
   fChain->SetBranchAddress("Ias_noTIBnoTIDno3TEC", &Ias_noTIBnoTIDno3TEC, &b_Ias_noTIBnoTIDno3TEC);
   fChain->SetBranchAddress("Ias_PixelOnly", &Ias_PixelOnly, &b_Ias_PixelOnly);
   fChain->SetBranchAddress("Ias_StripOnly", &Ias_StripOnly, &b_Ias_StripOnly);
   fChain->SetBranchAddress("Ias_PixelOnly_noL1", &Ias_PixelOnly_noL1, &b_Ias_PixelOnly_noL1);
   fChain->SetBranchAddress("Ih", &Ih, &b_Ih);
   fChain->SetBranchAddress("Ick", &Ick, &b_Ick);
   fChain->SetBranchAddress("Fmip", &Fmip, &b_Fmip);
   fChain->SetBranchAddress("ProbXY", &ProbXY, &b_ProbXY);
   fChain->SetBranchAddress("ProbXY_noL1", &ProbXY_noL1, &b_ProbXY_noL1);
   fChain->SetBranchAddress("ProbQ", &ProbQ, &b_ProbQ);
   fChain->SetBranchAddress("ProbQ_noL1", &ProbQ_noL1, &b_ProbQ_noL1);
   fChain->SetBranchAddress("Ndof", &Ndof, &b_Ndof);
   fChain->SetBranchAddress("Chi2", &Chi2, &b_Chi2);
   fChain->SetBranchAddress("QualityMask", &QualityMask, &b_QualityMask);
   fChain->SetBranchAddress("isHighPurity", &isHighPurity, &b_isHighPurity);
   fChain->SetBranchAddress("EoverP", &EoverP, &b_EoverP);
   fChain->SetBranchAddress("isMuon", &isMuon, &b_isMuon);
   fChain->SetBranchAddress("isPhoton", &isPhoton, &b_isPhoton);
   fChain->SetBranchAddress("isElectron", &isElectron, &b_isElectron);
   fChain->SetBranchAddress("isChHadron", &isChHadron, &b_isChHadron);
   fChain->SetBranchAddress("isNeutHadron", &isNeutHadron, &b_isNeutHadron);
   fChain->SetBranchAddress("isPfTrack", &isPfTrack, &b_isPfTrack);
   fChain->SetBranchAddress("isUndefined", &isUndefined, &b_isUndefined);
   fChain->SetBranchAddress("ECAL_energy", &ECAL_energy, &b_ECAL_energy);
   fChain->SetBranchAddress("HCAL_energy", &HCAL_energy, &b_HCAL_energy);
   fChain->SetBranchAddress("TOF", &TOF, &b_TOF);
   fChain->SetBranchAddress("TOFErr", &TOFErr, &b_TOFErr);
   fChain->SetBranchAddress("TOF_ndof", &TOF_ndof, &b_TOF_ndof);
   fChain->SetBranchAddress("DTTOF", &DTTOF, &b_DTTOF);
   fChain->SetBranchAddress("DTTOFErr", &DTTOFErr, &b_DTTOFErr);
   fChain->SetBranchAddress("DTTOF_ndof", &DTTOF_ndof, &b_DTTOF_ndof);
   fChain->SetBranchAddress("CSCTOF", &CSCTOF, &b_CSCTOF);
   fChain->SetBranchAddress("CSCTOFErr", &CSCTOFErr, &b_CSCTOFErr);
   fChain->SetBranchAddress("CSCTOF_ndof", &CSCTOF_ndof, &b_CSCTOF_ndof);
   fChain->SetBranchAddress("Mass", &Mass, &b_Mass);
   fChain->SetBranchAddress("MassErr", &MassErr, &b_MassErr);
   fChain->SetBranchAddress("dZ", &dZ, &b_dZ);
   fChain->SetBranchAddress("dXY", &dXY, &b_dXY);
   fChain->SetBranchAddress("dR", &dR, &b_dR);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("NOH", &NOH, &b_NOH);
   fChain->SetBranchAddress("NOPH", &NOPH, &b_NOPH);
   fChain->SetBranchAddress("FOVH", &FOVH, &b_FOVH);
   fChain->SetBranchAddress("NOMH", &NOMH, &b_NOMH);
   fChain->SetBranchAddress("FOVHD", &FOVHD, &b_FOVHD);
   fChain->SetBranchAddress("NOM", &NOM, &b_NOM);
   fChain->SetBranchAddress("matchTrigMuon_minDeltaR", &matchTrigMuon_minDeltaR, &b_matchTrigMuon_minDeltaR);
   fChain->SetBranchAddress("matchTrigMuon_pT", &matchTrigMuon_pT, &b_matchTrigMuon_pT);
   fChain->SetBranchAddress("iso_TK", &iso_TK, &b_iso_TK);
   fChain->SetBranchAddress("iso_ECAL", &iso_ECAL, &b_iso_ECAL);
   fChain->SetBranchAddress("iso_HCAL", &iso_HCAL, &b_iso_HCAL);
   fChain->SetBranchAddress("track_genTrackMiniIsoSumPt", &track_genTrackMiniIsoSumPt, &b_track_genTrackMiniIsoSumPt);
   fChain->SetBranchAddress("PFMiniIso_relative", &PFMiniIso_relative, &b_PFMiniIso_relative);
   fChain->SetBranchAddress("PFMiniIso_wMuon_relative", &PFMiniIso_wMuon_relative, &b_PFMiniIso_wMuon_relative);
   fChain->SetBranchAddress("TrackPFIsolationR005_sumChargedHadronPt", &TrackPFIsolationR005_sumChargedHadronPt, &b_TrackPFIsolationR005_sumChargedHadronPt);
   fChain->SetBranchAddress("TrackPFIsolationR005_sumNeutralHadronPt", &TrackPFIsolationR005_sumNeutralHadronPt, &b_TrackPFIsolationR005_sumNeutralHadronPt);
   fChain->SetBranchAddress("TrackPFIsolationR005_sumPhotonPt", &TrackPFIsolationR005_sumPhotonPt, &b_TrackPFIsolationR005_sumPhotonPt);
   fChain->SetBranchAddress("TrackPFIsolationR005_sumPUPt", &TrackPFIsolationR005_sumPUPt, &b_TrackPFIsolationR005_sumPUPt);
   fChain->SetBranchAddress("TrackPFIsolationR01_sumChargedHadronPt", &TrackPFIsolationR01_sumChargedHadronPt, &b_TrackPFIsolationR01_sumChargedHadronPt);
   fChain->SetBranchAddress("TrackPFIsolationR01_sumNeutralHadronPt", &TrackPFIsolationR01_sumNeutralHadronPt, &b_TrackPFIsolationR01_sumNeutralHadronPt);
   fChain->SetBranchAddress("TrackPFIsolationR01_sumPhotonPt", &TrackPFIsolationR01_sumPhotonPt, &b_TrackPFIsolationR01_sumPhotonPt);
   fChain->SetBranchAddress("TrackPFIsolationR01_sumPUPt", &TrackPFIsolationR01_sumPUPt, &b_TrackPFIsolationR01_sumPUPt);
   fChain->SetBranchAddress("TrackPFIsolationR03_sumChargedHadronPt", &TrackPFIsolationR03_sumChargedHadronPt, &b_TrackPFIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("TrackPFIsolationR03_sumNeutralHadronPt", &TrackPFIsolationR03_sumNeutralHadronPt, &b_TrackPFIsolationR03_sumNeutralHadronPt);
   fChain->SetBranchAddress("TrackPFIsolationR03_sumPhotonPt", &TrackPFIsolationR03_sumPhotonPt, &b_TrackPFIsolationR03_sumPhotonPt);
   fChain->SetBranchAddress("TrackPFIsolationR03_sumPUPt", &TrackPFIsolationR03_sumPUPt, &b_TrackPFIsolationR03_sumPUPt);
   fChain->SetBranchAddress("TrackPFIsolationR05_sumChargedHadronPt", &TrackPFIsolationR05_sumChargedHadronPt, &b_TrackPFIsolationR05_sumChargedHadronPt);
   fChain->SetBranchAddress("TrackPFIsolationR05_sumNeutralHadronPt", &TrackPFIsolationR05_sumNeutralHadronPt, &b_TrackPFIsolationR05_sumNeutralHadronPt);
   fChain->SetBranchAddress("TrackPFIsolationR05_sumPhotonPt", &TrackPFIsolationR05_sumPhotonPt, &b_TrackPFIsolationR05_sumPhotonPt);
   fChain->SetBranchAddress("TrackPFIsolationR05_sumPUPt", &TrackPFIsolationR05_sumPUPt, &b_TrackPFIsolationR05_sumPUPt);
   fChain->SetBranchAddress("MuonPFIsolationR03_sumChargedHadronPt", &MuonPFIsolationR03_sumChargedHadronPt, &b_MuonPFIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("MuonPFIsolationR03_sumNeutralHadronPt", &MuonPFIsolationR03_sumNeutralHadronPt, &b_MuonPFIsolationR03_sumNeutralHadronPt);
   fChain->SetBranchAddress("MuonPFIsolationR03_sumPhotonPt", &MuonPFIsolationR03_sumPhotonPt, &b_MuonPFIsolationR03_sumPhotonPt);
   fChain->SetBranchAddress("MuonPFIsolationR03_sumPUPt", &MuonPFIsolationR03_sumPUPt, &b_MuonPFIsolationR03_sumPUPt);
   fChain->SetBranchAddress("Ih_noL1", &Ih_noL1, &b_Ih_noL1);
   fChain->SetBranchAddress("Ih_15drop", &Ih_15drop, &b_Ih_15drop);
   fChain->SetBranchAddress("Ih_StripOnly", &Ih_StripOnly, &b_Ih_StripOnly);
   fChain->SetBranchAddress("Ih_StripOnly_15drop", &Ih_StripOnly_15drop, &b_Ih_StripOnly_15drop);
   fChain->SetBranchAddress("Ih_PixelOnly_noL1", &Ih_PixelOnly_noL1, &b_Ih_PixelOnly_noL1);
   fChain->SetBranchAddress("Ih_SaturationCorrectionFromFits", &Ih_SaturationCorrectionFromFits, &b_Ih_SaturationCorrectionFromFits);
   fChain->SetBranchAddress("clust_charge", &clust_charge, &b_clust_charge);
   fChain->SetBranchAddress("clust_pathlength", &clust_pathlength, &b_clust_pathlength);
   fChain->SetBranchAddress("clust_nstrip", &clust_nstrip, &b_clust_nstrip);
   fChain->SetBranchAddress("clust_sat254", &clust_sat254, &b_clust_sat254);
   fChain->SetBranchAddress("clust_sat255", &clust_sat255, &b_clust_sat255);
   fChain->SetBranchAddress("clust_detid", &clust_detid, &b_clust_detid);
   fChain->SetBranchAddress("clust_isStrip", &clust_isStrip, &b_clust_isStrip);
   fChain->SetBranchAddress("clust_isPixel", &clust_isPixel, &b_clust_isPixel);
   fChain->SetBranchAddress("GenId", &GenId, &b_GenId);
   fChain->SetBranchAddress("GenCharge", &GenCharge, &b_GenCharge);
   fChain->SetBranchAddress("GenMass", &GenMass, &b_GenMass);
   fChain->SetBranchAddress("GenPt", &GenPt, &b_GenPt);
   fChain->SetBranchAddress("GenEta", &GenEta, &b_GenEta);
   fChain->SetBranchAddress("GenPhi", &GenPhi, &b_GenPhi);
   Notify();
}

Bool_t HscpCandidates::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HscpCandidates::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HscpCandidates::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HscpCandidates_cxx
