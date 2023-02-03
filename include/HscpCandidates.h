//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 10 11:43:44 2023 by ROOT version 6.14/09
// from TTree HscpCandidates/HscpCandidates
// found on file: root://cmsxrootd.fnal.gov//store/group/lpchscp/ntuples/V1p2/MC_UL18/HSCPgluino_M-1800_TuneCP5_13TeV-pythia8/HSCPNtupler_2018_HSCPgluino_M-1800_wProbQ_v1p2_v1/221224_214648/0000/Histos_numEvent1000_7.root
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
   vector<int>     *BunchXing;
   vector<int>     *nPU;
   vector<float>   *nPUmean;
   UInt_t          nofVtx;
   Int_t           npv;
   vector<float>   *pvX;
   vector<float>   *pvY;
   vector<float>   *pvZ;
   vector<float>   *pvRho;
   vector<int>     *pvNdof;
   vector<float>   *pvChi2;
   vector<float>   *pvSumPt2;
   UInt_t          Hscp;
   UInt_t          nMuons;
   UInt_t          njets;
   Float_t         Weight;
   Float_t         GeneratorWeight;
   Float_t         GeneratorBinningValues;
   vector<bool>    *triggerDecision;
   vector<int>     *triggerHLTPrescale;
   vector<vector<float> > *triggerObjectE;
   vector<vector<float> > *triggerObjectPt;
   vector<vector<float> > *triggerObjectEta;
   vector<vector<float> > *triggerObjectPhi;
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
   vector<int>     *gParticleId;
   vector<int>     *gParticleStatus;
   vector<float>   *gParticleE;
   vector<float>   *gParticlePt;
   vector<float>   *gParticlePz;
   vector<float>   *gParticleEta;
   vector<float>   *gParticlePhi;
   vector<float>   *gParticleBeta;
   vector<int>     *gParticleCharge;
   vector<float>   *gParticleProdVertexX;
   vector<float>   *gParticleProdVertexY;
   vector<float>   *gParticleProdVertexZ;
   vector<int>     *gParticleMotherId;
   vector<int>     *gParticleMotherIndex;
   vector<float>   *eleE;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleCharge;
   vector<float>   *eleE_SC;
   vector<float>   *eleEta_SC;
   vector<float>   *elePhi_SC;
   vector<float>   *eleSigmaIetaIeta;
   vector<float>   *eleFull5x5SigmaIetaIeta;
   vector<float>   *eleR9;
   vector<float>   *ele_dEta;
   vector<float>   *ele_dPhi;
   vector<float>   *ele_HoverE;
   vector<float>   *ele_d0;
   vector<float>   *ele_dZ;
   vector<float>   *ele_pileupIso;
   vector<float>   *ele_chargedIso;
   vector<float>   *ele_photonIso;
   vector<float>   *ele_neutralHadIso;
   vector<int>     *ele_MissHits;
   vector<bool>    *ele_passCutBasedIDVeto;
   vector<bool>    *ele_passCutBasedIDLoose;
   vector<bool>    *ele_passCutBasedIDMedium;
   vector<bool>    *ele_passCutBasedIDTight;
   vector<bool>    *ele_passMVAIsoIDWP80;
   vector<bool>    *ele_passMVAIsoIDWP90;
   vector<bool>    *ele_passMVAIsoIDWPHZZ;
   vector<bool>    *ele_passMVAIsoIDWPLoose;
   vector<bool>    *ele_passMVANoIsoIDWP80;
   vector<bool>    *ele_passMVANoIsoIDWP90;
   vector<bool>    *ele_passMVANoIsoIDWPLoose;
   vector<bool>    *ele_PassConvVeto;
   vector<float>   *ele_OneOverEminusOneOverP;
   vector<float>   *muonE;
   vector<float>   *muonPt;
   vector<float>   *muonEta;
   vector<float>   *muonPhi;
   vector<int>     *muonCharge;
   vector<bool>    *muonIsLoose;
   vector<bool>    *muonIsMedium;
   vector<bool>    *muonIsTight;
   vector<float>   *muon_d0;
   vector<float>   *muon_d0Err;
   vector<float>   *muon_dZ;
   vector<float>   *muon_ip3d;
   vector<float>   *muon_ip3dSignificance;
   vector<unsigned int> *muonType;
   vector<unsigned int> *muonQuality;
   vector<float>   *muon_pileupIso;
   vector<float>   *muon_chargedIso;
   vector<float>   *muon_photonIso;
   vector<float>   *muon_neutralHadIso;
   vector<float>   *muon_validFractionTrackerHits;
   vector<float>   *muTree_muon_normChi2onE;
   vector<float>   *muon_chi2LocalPosition;
   vector<float>   *muon_kinkFinder;
   vector<float>   *muon_segmentCompatability;
   vector<float>   *muon_trkIso;
   vector<float>   *muon_tuneP_Pt;
   vector<float>   *muon_tuneP_PtErr;
   vector<float>   *muon_tuneP_Eta;
   vector<float>   *muon_tuneP_Phi;
   vector<int>     *muon_tuneP_MuonBestTrackType;
   vector<bool>    *muon_isHighPtMuon;
   vector<bool>    *muon_isTrackerHighPtMuon;
   vector<float>   *Jet_pt;
   vector<float>   *Jet_eta;
   vector<float>   *Jet_phi;
   vector<float>   *Jet_mass;
   vector<float>   *Jet_energy;
   vector<float>   *Jet_pdgId;
   vector<float>   *Jet_et;
   vector<float>   *Jet_chargedEmEnergyFraction;
   vector<float>   *Jet_neutralEmEnergyFraction;
   vector<float>   *Jet_chargedHadronEnergyFraction;
   vector<float>   *Jet_neutralHadronEnergyFraction;
   vector<float>   *Jet_muonEnergyFraction;
   vector<int>     *Jet_chargedMultiplicity;
   vector<int>     *Jet_neutralMultiplicity;
   vector<float>   *Jet_jetArea;
   vector<float>   *Jet_pileupE;
   vector<float>   *mT;
   vector<bool>    *passCutPt55;
   vector<bool>    *passPreselection;
   vector<bool>    *passPreselectionSept8;
   vector<bool>    *passSelection;
   vector<bool>    *isPFMuon;
   vector<bool>    *PFMuonPt;
   vector<float>   *Charge;
   vector<float>   *Pt;
   vector<float>   *PtErr;
   vector<float>   *Is_StripOnly;
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
   vector<float>   *dZ_pv;
   vector<float>   *dXY_pv;
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
   vector<float>   *HSCP_tuneP_Pt;
   vector<float>   *HSCP_tuneP_PtErr;
   vector<float>   *HSCP_tuneP_Eta;
   vector<float>   *HSCP_tuneP_Phi;
   vector<int>     *HSCP_tuneP_MuonBestTrackType;
   vector<int>     *HSCP_ErrorHisto_bin;
   vector<int>     *HSCP_type;
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
   TBranch        *b_BunchXing;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_nPUmean;   //!
   TBranch        *b_nofVtx;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_pvX;   //!
   TBranch        *b_pvY;   //!
   TBranch        *b_pvZ;   //!
   TBranch        *b_pvRho;   //!
   TBranch        *b_pvNdof;   //!
   TBranch        *b_pvChi2;   //!
   TBranch        *b_pvSumPt2;   //!
   TBranch        *b_Hscp;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_GeneratorWeight;   //!
   TBranch        *b_GeneratorBinningValues;   //!
   TBranch        *b_triggerDecision;   //!
   TBranch        *b_triggerHLTPrescale;   //!
   TBranch        *b_triggerObjectE;   //!
   TBranch        *b_triggerObjectPt;   //!
   TBranch        *b_triggerObjectEta;   //!
   TBranch        *b_triggerObjectPhi;   //!
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
   TBranch        *b_gParticleId;   //!
   TBranch        *b_gParticleStatus;   //!
   TBranch        *b_gParticleE;   //!
   TBranch        *b_gParticlePt;   //!
   TBranch        *b_gParticlePz;   //!
   TBranch        *b_gParticleEta;   //!
   TBranch        *b_gParticlePhi;   //!
   TBranch        *b_gParticleBeta;   //!
   TBranch        *b_gParticleCharge;   //!
   TBranch        *b_gParticleProdVertexX;   //!
   TBranch        *b_gParticleProdVertexY;   //!
   TBranch        *b_gParticleProdVertexZ;   //!
   TBranch        *b_gParticleMotherId;   //!
   TBranch        *b_gParticleMotherIndex;   //!
   TBranch        *b_eleE;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleE_SC;   //!
   TBranch        *b_eleEta_SC;   //!
   TBranch        *b_elePhi_SC;   //!
   TBranch        *b_eleSigmaIetaIeta;   //!
   TBranch        *b_eleFull5x5SigmaIetaIeta;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_ele_dEta;   //!
   TBranch        *b_ele_dPhi;   //!
   TBranch        *b_ele_HoverE;   //!
   TBranch        *b_ele_d0;   //!
   TBranch        *b_ele_dZ;   //!
   TBranch        *b_ele_pileupIso;   //!
   TBranch        *b_ele_chargedIso;   //!
   TBranch        *b_ele_photonIso;   //!
   TBranch        *b_ele_neutralHadIso;   //!
   TBranch        *b_ele_MissHits;   //!
   TBranch        *b_ele_passCutBasedIDVeto;   //!
   TBranch        *b_ele_passCutBasedIDLoose;   //!
   TBranch        *b_ele_passCutBasedIDMedium;   //!
   TBranch        *b_ele_passCutBasedIDTight;   //!
   TBranch        *b_ele_passMVAIsoIDWP80;   //!
   TBranch        *b_ele_passMVAIsoIDWP90;   //!
   TBranch        *b_ele_passMVAIsoIDWPHZZ;   //!
   TBranch        *b_ele_passMVAIsoIDWPLoose;   //!
   TBranch        *b_ele_passMVANoIsoIDWP80;   //!
   TBranch        *b_ele_passMVANoIsoIDWP90;   //!
   TBranch        *b_ele_passMVANoIsoIDWPLoose;   //!
   TBranch        *b_ele_PassConvVeto;   //!
   TBranch        *b_ele_OneOverEminusOneOverP;   //!
   TBranch        *b_muonE;   //!
   TBranch        *b_muonPt;   //!
   TBranch        *b_muonEta;   //!
   TBranch        *b_muonPhi;   //!
   TBranch        *b_muonCharge;   //!
   TBranch        *b_muonIsLoose;   //!
   TBranch        *b_muonIsMedium;   //!
   TBranch        *b_muonIsTight;   //!
   TBranch        *b_muon_d0;   //!
   TBranch        *b_muon_d0Err;   //!
   TBranch        *b_muon_dZ;   //!
   TBranch        *b_muon_ip3d;   //!
   TBranch        *b_muon_ip3dSignificance;   //!
   TBranch        *b_muonType;   //!
   TBranch        *b_muonQuality;   //!
   TBranch        *b_muon_pileupIso;   //!
   TBranch        *b_muon_chargedIso;   //!
   TBranch        *b_muon_photonIso;   //!
   TBranch        *b_muon_neutralHadIso;   //!
   TBranch        *b_muon_validFractionTrackerHits;   //!
   TBranch        *b_muTree_muon_normChi2onE;   //!
   TBranch        *b_muon_chi2LocalPosition;   //!
   TBranch        *b_muon_kinkFinder;   //!
   TBranch        *b_muon_segmentCompatability;   //!
   TBranch        *b_muon_trkIso;   //!
   TBranch        *b_muon_tuneP_Pt;   //!
   TBranch        *b_muon_tuneP_PtErr;   //!
   TBranch        *b_muon_tuneP_Eta;   //!
   TBranch        *b_muon_tuneP_Phi;   //!
   TBranch        *b_muon_tuneP_MuonBestTrackType;   //!
   TBranch        *b_muon_isHighPtMuon;   //!
   TBranch        *b_muon_isTrackerHighPtMuon;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_pdgId;   //!
   TBranch        *b_Jet_et;   //!
   TBranch        *b_Jet_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_neutralEmEnergyFraction;   //!
   TBranch        *b_Jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jet_neutralHadronEnergyFraction;   //!
   TBranch        *b_Jet_muonEnergyFraction;   //!
   TBranch        *b_Jet_chargedMultiplicity;   //!
   TBranch        *b_Jet_neutralMultiplicity;   //!
   TBranch        *b_Jet_jetArea;   //!
   TBranch        *b_Jet_pileupE;   //!
   TBranch        *b_mT;   //!
   TBranch        *b_passCutPt55;   //!
   TBranch        *b_passPreselection;   //!
   TBranch        *b_passPreselectionSept8;   //!
   TBranch        *b_passSelection;   //!
   TBranch        *b_isPFMuon;   //!
   TBranch        *b_PFMuonPt;   //!
   TBranch        *b_Charge;   //!
   TBranch        *b_Pt;   //!
   TBranch        *b_PtErr;   //!
   TBranch        *b_Is_StripOnly;   //!
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
   TBranch        *b_dZ_pv;   //!
   TBranch        *b_dXY_pv;   //!
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
   TBranch        *b_HSCP_tuneP_Pt;   //!
   TBranch        *b_HSCP_tuneP_PtErr;   //!
   TBranch        *b_HSCP_tuneP_Eta;   //!
   TBranch        *b_HSCP_tuneP_Phi;   //!
   TBranch        *b_HSCP_tuneP_MuonBestTrackType;   //!
   TBranch        *b_HSCP_ErrorHisto_bin;   //!
   TBranch        *b_HSCP_type;   //!
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmsxrootd.fnal.gov//store/group/lpchscp/ntuples/V1p2/MC_UL18/HSCPgluino_M-1800_TuneCP5_13TeV-pythia8/HSCPNtupler_2018_HSCPgluino_M-1800_wProbQ_v1p2_v1/221224_214648/0000/Histos_numEvent1000_7.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://cmsxrootd.fnal.gov//store/group/lpchscp/ntuples/V1p2/MC_UL18/HSCPgluino_M-1800_TuneCP5_13TeV-pythia8/HSCPNtupler_2018_HSCPgluino_M-1800_wProbQ_v1p2_v1/221224_214648/0000/Histos_numEvent1000_7.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://cmsxrootd.fnal.gov//store/group/lpchscp/ntuples/V1p2/MC_UL18/HSCPgluino_M-1800_TuneCP5_13TeV-pythia8/HSCPNtupler_2018_HSCPgluino_M-1800_wProbQ_v1p2_v1/221224_214648/0000/Histos_numEvent1000_7.root:/HSCParticleAnalyzer/BaseName");
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
   BunchXing = 0;
   nPU = 0;
   nPUmean = 0;
   pvX = 0;
   pvY = 0;
   pvZ = 0;
   pvRho = 0;
   pvNdof = 0;
   pvChi2 = 0;
   pvSumPt2 = 0;
   triggerDecision = 0;
   triggerHLTPrescale = 0;
   triggerObjectE = 0;
   triggerObjectPt = 0;
   triggerObjectEta = 0;
   triggerObjectPhi = 0;
   gParticleId = 0;
   gParticleStatus = 0;
   gParticleE = 0;
   gParticlePt = 0;
   gParticlePz = 0;
   gParticleEta = 0;
   gParticlePhi = 0;
   gParticleBeta = 0;
   gParticleCharge = 0;
   gParticleProdVertexX = 0;
   gParticleProdVertexY = 0;
   gParticleProdVertexZ = 0;
   gParticleMotherId = 0;
   gParticleMotherIndex = 0;
   eleE = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleCharge = 0;
   eleE_SC = 0;
   eleEta_SC = 0;
   elePhi_SC = 0;
   eleSigmaIetaIeta = 0;
   eleFull5x5SigmaIetaIeta = 0;
   eleR9 = 0;
   ele_dEta = 0;
   ele_dPhi = 0;
   ele_HoverE = 0;
   ele_d0 = 0;
   ele_dZ = 0;
   ele_pileupIso = 0;
   ele_chargedIso = 0;
   ele_photonIso = 0;
   ele_neutralHadIso = 0;
   ele_MissHits = 0;
   ele_passCutBasedIDVeto = 0;
   ele_passCutBasedIDLoose = 0;
   ele_passCutBasedIDMedium = 0;
   ele_passCutBasedIDTight = 0;
   ele_passMVAIsoIDWP80 = 0;
   ele_passMVAIsoIDWP90 = 0;
   ele_passMVAIsoIDWPHZZ = 0;
   ele_passMVAIsoIDWPLoose = 0;
   ele_passMVANoIsoIDWP80 = 0;
   ele_passMVANoIsoIDWP90 = 0;
   ele_passMVANoIsoIDWPLoose = 0;
   ele_PassConvVeto = 0;
   ele_OneOverEminusOneOverP = 0;
   muonE = 0;
   muonPt = 0;
   muonEta = 0;
   muonPhi = 0;
   muonCharge = 0;
   muonIsLoose = 0;
   muonIsMedium = 0;
   muonIsTight = 0;
   muon_d0 = 0;
   muon_d0Err = 0;
   muon_dZ = 0;
   muon_ip3d = 0;
   muon_ip3dSignificance = 0;
   muonType = 0;
   muonQuality = 0;
   muon_pileupIso = 0;
   muon_chargedIso = 0;
   muon_photonIso = 0;
   muon_neutralHadIso = 0;
   muon_validFractionTrackerHits = 0;
   muTree_muon_normChi2onE = 0;
   muon_chi2LocalPosition = 0;
   muon_kinkFinder = 0;
   muon_segmentCompatability = 0;
   muon_trkIso = 0;
   muon_tuneP_Pt = 0;
   muon_tuneP_PtErr = 0;
   muon_tuneP_Eta = 0;
   muon_tuneP_Phi = 0;
   muon_tuneP_MuonBestTrackType = 0;
   muon_isHighPtMuon = 0;
   muon_isTrackerHighPtMuon = 0;
   Jet_pt = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_mass = 0;
   Jet_energy = 0;
   Jet_pdgId = 0;
   Jet_et = 0;
   Jet_chargedEmEnergyFraction = 0;
   Jet_neutralEmEnergyFraction = 0;
   Jet_chargedHadronEnergyFraction = 0;
   Jet_neutralHadronEnergyFraction = 0;
   Jet_muonEnergyFraction = 0;
   Jet_chargedMultiplicity = 0;
   Jet_neutralMultiplicity = 0;
   Jet_jetArea = 0;
   Jet_pileupE = 0;
   mT = 0;
   passCutPt55 = 0;
   passPreselection = 0;
   passPreselectionSept8 = 0;
   passSelection = 0;
   isPFMuon = 0;
   PFMuonPt = 0;
   Charge = 0;
   Pt = 0;
   PtErr = 0;
   Is_StripOnly = 0;
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
   dZ_pv = 0;
   dXY_pv = 0;
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
   HSCP_tuneP_Pt = 0;
   HSCP_tuneP_PtErr = 0;
   HSCP_tuneP_Eta = 0;
   HSCP_tuneP_Phi = 0;
   HSCP_tuneP_MuonBestTrackType = 0;
   HSCP_ErrorHisto_bin = 0;
   HSCP_type = 0;
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
   fChain->SetBranchAddress("BunchXing", &BunchXing, &b_BunchXing);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("nPUmean", &nPUmean, &b_nPUmean);
   fChain->SetBranchAddress("nofVtx", &nofVtx, &b_nofVtx);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("pvX", &pvX, &b_pvX);
   fChain->SetBranchAddress("pvY", &pvY, &b_pvY);
   fChain->SetBranchAddress("pvZ", &pvZ, &b_pvZ);
   fChain->SetBranchAddress("pvRho", &pvRho, &b_pvRho);
   fChain->SetBranchAddress("pvNdof", &pvNdof, &b_pvNdof);
   fChain->SetBranchAddress("pvChi2", &pvChi2, &b_pvChi2);
   fChain->SetBranchAddress("pvSumPt2", &pvSumPt2, &b_pvSumPt2);
   fChain->SetBranchAddress("Hscp", &Hscp, &b_Hscp);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("GeneratorWeight", &GeneratorWeight, &b_GeneratorWeight);
   fChain->SetBranchAddress("GeneratorBinningValues", &GeneratorBinningValues, &b_GeneratorBinningValues);
   fChain->SetBranchAddress("triggerDecision", &triggerDecision, &b_triggerDecision);
   fChain->SetBranchAddress("triggerHLTPrescale", &triggerHLTPrescale, &b_triggerHLTPrescale);
   fChain->SetBranchAddress("triggerObjectE", &triggerObjectE, &b_triggerObjectE);
   fChain->SetBranchAddress("triggerObjectPt", &triggerObjectPt, &b_triggerObjectPt);
   fChain->SetBranchAddress("triggerObjectEta", &triggerObjectEta, &b_triggerObjectEta);
   fChain->SetBranchAddress("triggerObjectPhi", &triggerObjectPhi, &b_triggerObjectPhi);
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
   fChain->SetBranchAddress("gParticleId", &gParticleId, &b_gParticleId);
   fChain->SetBranchAddress("gParticleStatus", &gParticleStatus, &b_gParticleStatus);
   fChain->SetBranchAddress("gParticleE", &gParticleE, &b_gParticleE);
   fChain->SetBranchAddress("gParticlePt", &gParticlePt, &b_gParticlePt);
   fChain->SetBranchAddress("gParticlePz", &gParticlePz, &b_gParticlePz);
   fChain->SetBranchAddress("gParticleEta", &gParticleEta, &b_gParticleEta);
   fChain->SetBranchAddress("gParticlePhi", &gParticlePhi, &b_gParticlePhi);
   fChain->SetBranchAddress("gParticleBeta", &gParticleBeta, &b_gParticleBeta);
   fChain->SetBranchAddress("gParticleCharge", &gParticleCharge, &b_gParticleCharge);
   fChain->SetBranchAddress("gParticleProdVertexX", &gParticleProdVertexX, &b_gParticleProdVertexX);
   fChain->SetBranchAddress("gParticleProdVertexY", &gParticleProdVertexY, &b_gParticleProdVertexY);
   fChain->SetBranchAddress("gParticleProdVertexZ", &gParticleProdVertexZ, &b_gParticleProdVertexZ);
   fChain->SetBranchAddress("gParticleMotherId", &gParticleMotherId, &b_gParticleMotherId);
   fChain->SetBranchAddress("gParticleMotherIndex", &gParticleMotherIndex, &b_gParticleMotherIndex);
   fChain->SetBranchAddress("eleE", &eleE, &b_eleE);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleE_SC", &eleE_SC, &b_eleE_SC);
   fChain->SetBranchAddress("eleEta_SC", &eleEta_SC, &b_eleEta_SC);
   fChain->SetBranchAddress("elePhi_SC", &elePhi_SC, &b_elePhi_SC);
   fChain->SetBranchAddress("eleSigmaIetaIeta", &eleSigmaIetaIeta, &b_eleSigmaIetaIeta);
   fChain->SetBranchAddress("eleFull5x5SigmaIetaIeta", &eleFull5x5SigmaIetaIeta, &b_eleFull5x5SigmaIetaIeta);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("ele_dEta", &ele_dEta, &b_ele_dEta);
   fChain->SetBranchAddress("ele_dPhi", &ele_dPhi, &b_ele_dPhi);
   fChain->SetBranchAddress("ele_HoverE", &ele_HoverE, &b_ele_HoverE);
   fChain->SetBranchAddress("ele_d0", &ele_d0, &b_ele_d0);
   fChain->SetBranchAddress("ele_dZ", &ele_dZ, &b_ele_dZ);
   fChain->SetBranchAddress("ele_pileupIso", &ele_pileupIso, &b_ele_pileupIso);
   fChain->SetBranchAddress("ele_chargedIso", &ele_chargedIso, &b_ele_chargedIso);
   fChain->SetBranchAddress("ele_photonIso", &ele_photonIso, &b_ele_photonIso);
   fChain->SetBranchAddress("ele_neutralHadIso", &ele_neutralHadIso, &b_ele_neutralHadIso);
   fChain->SetBranchAddress("ele_MissHits", &ele_MissHits, &b_ele_MissHits);
   fChain->SetBranchAddress("ele_passCutBasedIDVeto", &ele_passCutBasedIDVeto, &b_ele_passCutBasedIDVeto);
   fChain->SetBranchAddress("ele_passCutBasedIDLoose", &ele_passCutBasedIDLoose, &b_ele_passCutBasedIDLoose);
   fChain->SetBranchAddress("ele_passCutBasedIDMedium", &ele_passCutBasedIDMedium, &b_ele_passCutBasedIDMedium);
   fChain->SetBranchAddress("ele_passCutBasedIDTight", &ele_passCutBasedIDTight, &b_ele_passCutBasedIDTight);
   fChain->SetBranchAddress("ele_passMVAIsoIDWP80", &ele_passMVAIsoIDWP80, &b_ele_passMVAIsoIDWP80);
   fChain->SetBranchAddress("ele_passMVAIsoIDWP90", &ele_passMVAIsoIDWP90, &b_ele_passMVAIsoIDWP90);
   fChain->SetBranchAddress("ele_passMVAIsoIDWPHZZ", &ele_passMVAIsoIDWPHZZ, &b_ele_passMVAIsoIDWPHZZ);
   fChain->SetBranchAddress("ele_passMVAIsoIDWPLoose", &ele_passMVAIsoIDWPLoose, &b_ele_passMVAIsoIDWPLoose);
   fChain->SetBranchAddress("ele_passMVANoIsoIDWP80", &ele_passMVANoIsoIDWP80, &b_ele_passMVANoIsoIDWP80);
   fChain->SetBranchAddress("ele_passMVANoIsoIDWP90", &ele_passMVANoIsoIDWP90, &b_ele_passMVANoIsoIDWP90);
   fChain->SetBranchAddress("ele_passMVANoIsoIDWPLoose", &ele_passMVANoIsoIDWPLoose, &b_ele_passMVANoIsoIDWPLoose);
   fChain->SetBranchAddress("ele_PassConvVeto", &ele_PassConvVeto, &b_ele_PassConvVeto);
   fChain->SetBranchAddress("ele_OneOverEminusOneOverP", &ele_OneOverEminusOneOverP, &b_ele_OneOverEminusOneOverP);
   fChain->SetBranchAddress("muonE", &muonE, &b_muonE);
   fChain->SetBranchAddress("muonPt", &muonPt, &b_muonPt);
   fChain->SetBranchAddress("muonEta", &muonEta, &b_muonEta);
   fChain->SetBranchAddress("muonPhi", &muonPhi, &b_muonPhi);
   fChain->SetBranchAddress("muonCharge", &muonCharge, &b_muonCharge);
   fChain->SetBranchAddress("muonIsLoose", &muonIsLoose, &b_muonIsLoose);
   fChain->SetBranchAddress("muonIsMedium", &muonIsMedium, &b_muonIsMedium);
   fChain->SetBranchAddress("muonIsTight", &muonIsTight, &b_muonIsTight);
   fChain->SetBranchAddress("muon_d0", &muon_d0, &b_muon_d0);
   fChain->SetBranchAddress("muon_d0Err", &muon_d0Err, &b_muon_d0Err);
   fChain->SetBranchAddress("muon_dZ", &muon_dZ, &b_muon_dZ);
   fChain->SetBranchAddress("muon_ip3d", &muon_ip3d, &b_muon_ip3d);
   fChain->SetBranchAddress("muon_ip3dSignificance", &muon_ip3dSignificance, &b_muon_ip3dSignificance);
   fChain->SetBranchAddress("muonType", &muonType, &b_muonType);
   fChain->SetBranchAddress("muonQuality", &muonQuality, &b_muonQuality);
   fChain->SetBranchAddress("muon_pileupIso", &muon_pileupIso, &b_muon_pileupIso);
   fChain->SetBranchAddress("muon_chargedIso", &muon_chargedIso, &b_muon_chargedIso);
   fChain->SetBranchAddress("muon_photonIso", &muon_photonIso, &b_muon_photonIso);
   fChain->SetBranchAddress("muon_neutralHadIso", &muon_neutralHadIso, &b_muon_neutralHadIso);
   fChain->SetBranchAddress("muon_validFractionTrackerHits", &muon_validFractionTrackerHits, &b_muon_validFractionTrackerHits);
   fChain->SetBranchAddress("muTree_muon_normChi2onE", &muTree_muon_normChi2onE, &b_muTree_muon_normChi2onE);
   fChain->SetBranchAddress("muon_chi2LocalPosition", &muon_chi2LocalPosition, &b_muon_chi2LocalPosition);
   fChain->SetBranchAddress("muon_kinkFinder", &muon_kinkFinder, &b_muon_kinkFinder);
   fChain->SetBranchAddress("muon_segmentCompatability", &muon_segmentCompatability, &b_muon_segmentCompatability);
   fChain->SetBranchAddress("muon_trkIso", &muon_trkIso, &b_muon_trkIso);
   fChain->SetBranchAddress("muon_tuneP_Pt", &muon_tuneP_Pt, &b_muon_tuneP_Pt);
   fChain->SetBranchAddress("muon_tuneP_PtErr", &muon_tuneP_PtErr, &b_muon_tuneP_PtErr);
   fChain->SetBranchAddress("muon_tuneP_Eta", &muon_tuneP_Eta, &b_muon_tuneP_Eta);
   fChain->SetBranchAddress("muon_tuneP_Phi", &muon_tuneP_Phi, &b_muon_tuneP_Phi);
   fChain->SetBranchAddress("muon_tuneP_MuonBestTrackType", &muon_tuneP_MuonBestTrackType, &b_muon_tuneP_MuonBestTrackType);
   fChain->SetBranchAddress("muon_isHighPtMuon", &muon_isHighPtMuon, &b_muon_isHighPtMuon);
   fChain->SetBranchAddress("muon_isTrackerHighPtMuon", &muon_isTrackerHighPtMuon, &b_muon_isTrackerHighPtMuon);
   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
   fChain->SetBranchAddress("Jet_pdgId", &Jet_pdgId, &b_Jet_pdgId);
   fChain->SetBranchAddress("Jet_et", &Jet_et, &b_Jet_et);
   fChain->SetBranchAddress("Jet_chargedEmEnergyFraction", &Jet_chargedEmEnergyFraction, &b_Jet_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jet_neutralEmEnergyFraction", &Jet_neutralEmEnergyFraction, &b_Jet_neutralEmEnergyFraction);
   fChain->SetBranchAddress("Jet_chargedHadronEnergyFraction", &Jet_chargedHadronEnergyFraction, &b_Jet_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jet_neutralHadronEnergyFraction", &Jet_neutralHadronEnergyFraction, &b_Jet_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("Jet_muonEnergyFraction", &Jet_muonEnergyFraction, &b_Jet_muonEnergyFraction);
   fChain->SetBranchAddress("Jet_chargedMultiplicity", &Jet_chargedMultiplicity, &b_Jet_chargedMultiplicity);
   fChain->SetBranchAddress("Jet_neutralMultiplicity", &Jet_neutralMultiplicity, &b_Jet_neutralMultiplicity);
   fChain->SetBranchAddress("Jet_jetArea", &Jet_jetArea, &b_Jet_jetArea);
   fChain->SetBranchAddress("Jet_pileupE", &Jet_pileupE, &b_Jet_pileupE);
   fChain->SetBranchAddress("mT", &mT, &b_mT);
   fChain->SetBranchAddress("passCutPt55", &passCutPt55, &b_passCutPt55);
   fChain->SetBranchAddress("passPreselection", &passPreselection, &b_passPreselection);
   fChain->SetBranchAddress("passPreselectionSept8", &passPreselectionSept8, &b_passPreselectionSept8);
   fChain->SetBranchAddress("passSelection", &passSelection, &b_passSelection);
   fChain->SetBranchAddress("isPFMuon", &isPFMuon, &b_isPFMuon);
   fChain->SetBranchAddress("PFMuonPt", &PFMuonPt, &b_PFMuonPt);
   fChain->SetBranchAddress("Charge", &Charge, &b_Charge);
   fChain->SetBranchAddress("Pt", &Pt, &b_Pt);
   fChain->SetBranchAddress("PtErr", &PtErr, &b_PtErr);
   fChain->SetBranchAddress("Is_StripOnly", &Is_StripOnly, &b_Is_StripOnly);
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
   fChain->SetBranchAddress("dZ_pv", &dZ_pv, &b_dZ_pv);
   fChain->SetBranchAddress("dXY_pv", &dXY_pv, &b_dXY_pv);
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
   fChain->SetBranchAddress("HSCP_tuneP_Pt", &HSCP_tuneP_Pt, &b_HSCP_tuneP_Pt);
   fChain->SetBranchAddress("HSCP_tuneP_PtErr", &HSCP_tuneP_PtErr, &b_HSCP_tuneP_PtErr);
   fChain->SetBranchAddress("HSCP_tuneP_Eta", &HSCP_tuneP_Eta, &b_HSCP_tuneP_Eta);
   fChain->SetBranchAddress("HSCP_tuneP_Phi", &HSCP_tuneP_Phi, &b_HSCP_tuneP_Phi);
   fChain->SetBranchAddress("HSCP_tuneP_MuonBestTrackType", &HSCP_tuneP_MuonBestTrackType, &b_HSCP_tuneP_MuonBestTrackType);
   fChain->SetBranchAddress("HSCP_ErrorHisto_bin", &HSCP_ErrorHisto_bin, &b_HSCP_ErrorHisto_bin);
   fChain->SetBranchAddress("HSCP_type", &HSCP_type, &b_HSCP_type);
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
