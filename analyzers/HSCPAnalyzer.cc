#include "HSCPAnalyzer.h"
#include "RazorHelper.h"
#include "HSCPTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include <iostream>
#include <random>
//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"


using namespace std::chrono;
using namespace std;
using namespace ROOT::Math;
struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  float dZ;
  bool tightId;
  bool looseId;
  bool isGlobal;
  bool passVetoId;
  bool passLooseIso;
  bool passTightIso;
  bool passVTightIso;
  bool passVVTightIso;

};
struct largest_pt
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt;
void HSCPAnalyzer::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  //options format: MH/MX/ctau/condor: 1000/300/0/1


  int option = options%10;

  if( isData )
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }



  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;

  if (analysisTag == ""){
    analysisTag = "test";

  }


  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "HSCP_tree.root";
  TFile *outFile;
  outFile = new TFile(outfilename.c_str(), "RECREATE");


  HSCPTree *tree = new HSCPTree;
  tree->CreateTree();
  tree->tree_->SetAutoFlush(0);
  tree->InitTree();

  // for signals, need one output file for each signal point


  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);



  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1



  //--------------------------------
  //Initialize helper
  //--------------------------------
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData);


  //*************************************************************************
  //Look over Input File Events
  //*************************************************************************
  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();
  for (Long64_t jentry=0; jentry<fChain->GetEntries();jentry++) {

    //begin event
    if(jentry % 10000 == 0)
    {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      cout << "Processing entry " << jentry << endl;
      cout << "Time taken by program is : " << time_taken << endl;
      start = clock();
    }
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //GetEntry(ientry);
    nb = fChain->GetEntry(jentry); nbytes += nb;

    tree->InitVariables();



    //
    // //event info
    if (isData)
    {
      NEvents->Fill(1);
      tree->weight = 1;
    }
    else
    {
      //tree->weight = GeneratorWeight;
      //NEvents->Fill(1, GeneratorWeight);
      tree->weight = Weight;
      NEvents->Fill(1, Weight);
    }
    tree->runNum = Run;
    tree->lumiSec = Lumi;
    tree->evtNum = Event;
    tree->pileupWeight = Weight;
      //*************************************************************************
      // MET-related variables
      //*************************************************************************
      tree->met = RecoPFMET;
      tree->metPhi = RecoPFMET_phi;

      //Triggers

      for(unsigned int i = 0; i < triggerDecision->size(); i++){
        tree->HLTDecision[i] = triggerDecision->at(i);
      }
     if(isData && !tree->HLTDecision[196])continue;


     // gen particles

     if (!isData)
     {
       for (unsigned int i=0; i < gParticleStatus->size(); i++)
       {
         if (gParticleStatus->at(i)!=1)continue;
         //// stau

         int thePDGidForCandidate = abs(gParticleId->at(i));
         tree->genHSCP_type[tree->nGenHSCP] = 999;
         // R-hadrons
         if (   thePDGidForCandidate == 1000993 || thePDGidForCandidate == 1009113
             || thePDGidForCandidate == 1009223 || thePDGidForCandidate == 1009313
             || thePDGidForCandidate == 1009333 || thePDGidForCandidate == 1092114
             || thePDGidForCandidate == 1093214 || thePDGidForCandidate == 1093324
             || thePDGidForCandidate == 1000622 || thePDGidForCandidate == 1000642
             || thePDGidForCandidate == 1006113 || thePDGidForCandidate == 1006311
             || thePDGidForCandidate == 1006313 || thePDGidForCandidate == 1006333) {
           tree->genHSCP_type[tree->nGenHSCP] = 0;
         }
         // Single-charged HSCP
         else if (   thePDGidForCandidate == 1009213 || thePDGidForCandidate == 1009323
              || thePDGidForCandidate == 1091114 || thePDGidForCandidate == 1092214
              || thePDGidForCandidate == 1093114 || thePDGidForCandidate == 1093224
              || thePDGidForCandidate == 1093314 || thePDGidForCandidate == 1093334
              || thePDGidForCandidate == 1000612 || thePDGidForCandidate == 1000632
              || thePDGidForCandidate == 1000652 || thePDGidForCandidate == 1006211
              || thePDGidForCandidate == 1006213 || thePDGidForCandidate == 1006321
              || thePDGidForCandidate == 1006323 || thePDGidForCandidate == 1000015) {
           tree->genHSCP_type[tree->nGenHSCP] = 1;
         }
         // Double-charged R-hadrons
         else if (thePDGidForCandidate == 1092224 || thePDGidForCandidate == 1006223) {
           tree->genHSCP_type[tree->nGenHSCP] = 2;
         }
         // tau prime, could be single or multiple charged
         else if (thePDGidForCandidate == 17) {
           tree->genHSCP_type[tree->nGenHSCP] = 3;
         }
         else if (thePDGidForCandidate>=1000000)
         {
           tree->genHSCP_type[tree->nGenHSCP] = 4;
         }
         else {
           continue;
         }

         tree->genHSCP_pdgid[tree->nGenHSCP] = gParticleId->at(i);
         tree->genHSCP_pt[tree->nGenHSCP] = gParticlePt->at(i);
         tree->genHSCP_eta[tree->nGenHSCP] = gParticleEta->at(i);
         tree->genHSCP_phi[tree->nGenHSCP] = gParticlePhi->at(i);
         tree->genHSCP_e[tree->nGenHSCP] = gParticleE->at(i);

         tree->nGenHSCP++;
       }
     }// end if !isData


       // construct genlevel ZMass
      double gen_ZMass = -999;
      double gen_ZPt = -999;
      double gen_tmpDistToZPole = 9999;
      pair<uint,uint> gen_ZCandidateLeptonIndex;
      bool gen_foundZ = false;
      TLorentzVector gen_ZCandidate;
      for(int i = 0; i < tree->nGenHSCP; i++)
      {

        for(int j = i+1; j < tree->nGenHSCP; j++)
        {
          TLorentzVector hscp1, hscp2;
          hscp1.SetPtEtaPhiM(tree->genHSCP_pt[i], tree->genHSCP_eta[i], tree->genHSCP_phi[i], MU_MASS);
          hscp2.SetPtEtaPhiM(tree->genHSCP_pt[j], tree->genHSCP_eta[j], tree->genHSCP_phi[j], MU_MASS);
          double tmpMass = (hscp1 + hscp2).M();
          //select the pair closest to Z pole mass
          if ( fabs( tmpMass - Z_MASS) < gen_tmpDistToZPole)
          {
            gen_tmpDistToZPole = tmpMass;
            gen_ZCandidateLeptonIndex = pair<int,int>(i,j);
            gen_ZMass = tmpMass;
            gen_ZPt = (hscp1 + hscp2).Pt();
            gen_ZCandidate = hscp1 + hscp2;
            gen_foundZ = true;
          }
        }
      }
      if (gen_foundZ)
      {
        tree->genZMass = gen_ZMass;
        tree->genZPt   = gen_ZPt;
        tree->genZEta  = gen_ZCandidate.Eta();
        tree->genZPhi  = gen_ZCandidate.Phi();
        tree->genZleptonIndex1 = gen_ZCandidateLeptonIndex.first;
        tree->genZleptonIndex2 = gen_ZCandidateLeptonIndex.second;
      } // endif foundZ

      //*************************************************************************
    //Start Object Selection
    //*************************************************************************
    std::vector<leptons> Leptons;
   //-------------------------------
   //Muons
   //-------------------------------
   for( int i = 0; i < nMuons; i++ )
   {
     if (muonPt->at(i)<20)continue;
     if (fabs(muonEta->at(i))>2.4)continue;
     if (!muonIsLoose->at(i))continue;

     float muonIso = (muon_chargedIso->at(i) + fmax(0.0,  muon_photonIso->at(i) + muon_neutralHadIso->at(i) - 0.5*muon_pileupIso->at(i))) / muonPt->at(i);
     if (muonIso >= 0.25) continue;
     //remove overlaps
     bool overlap = false;
     for(auto& lep : Leptons)
     {
       if (RazorAnalyzer::deltaR(muonEta->at(i),muonPhi->at(i),lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
     }
     if(overlap) continue;

     leptons tmpMuon;
     tmpMuon.lepton.SetPtEtaPhiM(muonPt->at(i),muonEta->at(i), muonPhi->at(i), MU_MASS);
     tmpMuon.pdgId = 13 * -1 * muonCharge->at(i);
     tmpMuon.dZ = muon_dZ->at(i);
     tmpMuon.tightId = muonIsTight->at(i);
     tmpMuon.looseId = muonIsLoose->at(i);
     tmpMuon.passLooseIso = muonIso<0.25;
     tmpMuon.passTightIso = muonIso<0.15;
     tmpMuon.passVTightIso = muonIso<0.10;
     tmpMuon.passVVTightIso = muonIso<0.05;
     Leptons.push_back(tmpMuon);
   }

   sort(Leptons.begin(), Leptons.end(), my_largest_pt);
   for ( auto &tmp : Leptons )
   {
     tree->lepE[tree->nLeptons]      = tmp.lepton.E();
     tree->lepPt[tree->nLeptons]     = tmp.lepton.Pt();
     tree->lepEta[tree->nLeptons]    = tmp.lepton.Eta();
     tree->lepPhi[tree->nLeptons]    = tmp.lepton.Phi();
     tree->lepPdgId[tree->nLeptons]  = tmp.pdgId;
     tree->lepDZ[tree->nLeptons]     = tmp.dZ;
     tree->lepLooseId[tree->nLeptons] = tmp.looseId;
     tree->lepTightId[tree->nLeptons] = tmp.tightId;
     tree->lepPassLooseIso[tree->nLeptons] = tmp.passLooseIso;
     tree->lepPassTightIso[tree->nLeptons] = tmp.passTightIso;
     //tree->lepPassVTightIso[tree->nLeptons] = tmp.passVTightIso;
     //tree->lepPassVVTightIso[tree->nLeptons] = tmp.passVVTightIso;

     tree->nLeptons++;
   }


    /// HSCP stuff
    tree->nHSCP = 0;
    for(unsigned int i = 0; i < Pt->size(); i++)
    {
      // apply preselection
      if (isData && Pt->at(i) <= 55) continue;
      if (isData && fabs(eta->at(i)) > 1) continue;
      if (isData && NOPH->at(i) < 2) continue;
      if (isData && NOM->at(i) <= 9) continue;
      if (isData && !isHighPurity->at(i))continue;
      if (isData && FOVH->at(i) <= 0.8) continue;
      if (isData && EoverP->at(i) >= 0.3) continue;
      if (isData && Chi2->at(i)/Ndof->at(i) >= 5) continue;
      if (isData && track_genTrackMiniIsoSumPt->at(i) >= 15 )continue;
      if (isData && PFMiniIso_relative->at(i) >= 0.02) continue;
      if (isData && ProbQ_noL1->at(i) <= 0.0) continue;
      if (isData && ProbQ_noL1->at(i) > 0.7) continue;
      if (isData && PtErr->at(i)/(Pt->at(i)*Pt->at(i)) >= 0.0008)continue;
      //if (isData && fabs(dZ->at(i)) >= 0.1) continue;
      //if (isData && fabs(dXY->at(i)) >= 0.02) continue;
      tree->HSCP_pt[tree->nHSCP] = Pt->at(i);
      tree->HSCP_eta[tree->nHSCP] = eta->at(i);
      tree->HSCP_phi[tree->nHSCP] = phi->at(i);
      tree->HSCP_nPixelHit[tree->nHSCP] = NOPH->at(i);
      tree->HSCP_nHits[tree->nHSCP] = NOM->at(i);
      tree->HSCP_isHighPurity[tree->nHSCP] = isHighPurity->at(i);
      tree->HSCP_fracValidHits[tree->nHSCP] = FOVH->at(i);
      tree->HSCP_EOverP[tree->nHSCP] = EoverP->at(i);
      tree->HSCP_chi2[tree->nHSCP] = Chi2->at(i);
      tree->HSCP_nDof[tree->nHSCP] = Ndof->at(i);
      tree->HSCP_track_genTrackMiniIsoSumPt[tree->nHSCP] = track_genTrackMiniIsoSumPt->at(i);
      tree->HSCP_pfMiniIso_relative[tree->nHSCP] = PFMiniIso_relative->at(i);
      tree->HSCP_probQ[tree->nHSCP] = ProbQ_noL1->at(i);
      tree->HSCP_ptErr[tree->nHSCP] = PtErr->at(i);
      tree->HSCP_ias_StripOnly[tree->nHSCP] = Ias_StripOnly->at(i);
      tree->HSCP_mass[tree->nHSCP] = Mass->at(i);
      tree->HSCP_dZ[tree->nHSCP] = dZ->at(i);
      tree->HSCP_dXY[tree->nHSCP] = dXY->at(i);
      tree->HSCP_dZ_pv[tree->nHSCP] = dZ_pv->at(i);
      tree->HSCP_dXY_pv[tree->nHSCP] = dXY_pv->at(i);
      tree->HSCP_Ih[tree->nHSCP] = Ih_noL1->at(i);
      tree->HSCP_ErrorHisto_bin[tree->nHSCP] = HSCP_ErrorHisto_bin->at(i);
      tree->HSCP_type[tree->nHSCP] = HSCP_type->at(i);

      if (!isData)
      {

        float min_deltaR = 15.;
        int index = 999;
        for(int j = 0; j < tree->nGenHSCP;j++)
        {

          double current_delta_r = RazorAnalyzer::deltaR(tree->HSCP_eta[tree->nHSCP], tree->HSCP_phi[tree->nHSCP], tree->genHSCP_eta[j], tree->genHSCP_phi[j]);
          if (current_delta_r < min_deltaR)
          {
            min_deltaR = current_delta_r;
            index = j;
          }
        }
        if (min_deltaR < 0.4)tree->HSCP_match_genHSCP[tree->nHSCP] = true;
        else tree->HSCP_match_genHSCP[tree->nHSCP] = false;

        tree->HSCP_match_genHSCP_minDeltaR[tree->nHSCP] = min_deltaR;
        tree->HSCP_match_genHSCP_index[tree->nHSCP] = index;


      }
      tree->nHSCP++;
    }
     if (isData && tree->nHSCP < 1)continue;


     double ZMass = -999;
     double ZPt = -999;
     double tmpDistToZPole = 9999;
     pair<uint,uint> ZCandidateLeptonIndex;
     bool foundZ = false;
     TLorentzVector ZCandidate;
     for(unsigned int i = 0; i < Pt->size(); i++)
     {

       for(unsigned int j = i+1; j < Pt->size(); j++)
       {
         TLorentzVector hscp1, hscp2;
         hscp1.SetPtEtaPhiM(Pt->at(i), eta->at(i), phi->at(i), MU_MASS);
         hscp2.SetPtEtaPhiM(Pt->at(j), eta->at(j), phi->at(j), MU_MASS);
         double tmpMass = (hscp1 + hscp2).M();
         //select the pair closest to Z pole mass
         if ( fabs( tmpMass - Z_MASS) < tmpDistToZPole)
         {
           tmpDistToZPole = fabs( tmpMass - Z_MASS);
           ZCandidateLeptonIndex = pair<int,int>(i,j);
           ZMass = tmpMass;
           ZPt = (hscp1 + hscp2).Pt();
           ZCandidate = hscp1 + hscp2;
           foundZ = true;
         }
       }
     }
     if (foundZ)
   {
     tree->ZMass = ZMass;
     tree->ZPt   = ZPt;
     tree->ZEta  = ZCandidate.Eta();
     tree->ZPhi  = ZCandidate.Phi();
     tree->ZleptonIndex1 = ZCandidateLeptonIndex.first;
     tree->ZleptonIndex2 = ZCandidateLeptonIndex.second;
   } // endif foundZ


      tree->tree_->Fill();

    }// end of jentry loop
    if (!isData)
    {
       cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
       cout << "Writing output trees..." << endl;
       outFile->cd();
       tree->tree_->Write();
       NEvents->Write();
       outFile->Close();
    }


    else
    {
      cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
      cout << "Writing output trees..." << endl;
      outFile->cd();
      tree->tree_->Write();
      NEvents->Write();
      outFile->Close();
    }
}
