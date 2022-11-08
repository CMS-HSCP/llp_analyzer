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
      tree->weight = GeneratorWeight;
      NEvents->Fill(1, GeneratorWeight);
    }
    tree->runNum = Run;
    tree->lumiSec = Lumi;
    tree->evtNum = Event;
    tree->pileupWeight = Weight;
    //
    //
    //   for (int i=0; i < nBunchXing; i++)
    //   {
    //     if (BunchXing[i] == 0)
    //     {
    //       tree->npu = nPUmean[i];
    //     }
    //   }
    //   tree->pileupWeight = helper->getPileupWeight(tree->npu);
    //   tree->pileupWeightUp = helper->getPileupWeightUp(tree->npu) / tree->pileupWeight;
    //   tree->pileupWeightDown = helper->getPileupWeightDown(tree->npu) / tree->pileupWeight;
    //
    // }//end of isData
    //
    //   //get NPU
    //   tree->npv = nPV;
    //   tree->rho = fixedGridRhoFastjetAll;

      //*************************************************************************
      // MET-related variables
      //*************************************************************************
      tree->met = RecoPFMET;
      tree->metPhi = RecoPFMET_phi;


      //Triggers
      if (!HLT_Mu50) continue;
      tree->HLT_Mu50 = HLT_Mu50;
      tree->HLT_PFMET120_PFMHT120_IDTight = HLT_PFMET120_PFMHT120_IDTight;
      tree->HLT_PFHT500_PFMET100_PFMHT100_IDTight = HLT_PFHT500_PFMET100_PFMHT100_IDTight;
      tree->HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
      tree->HLT_MET105_IsoTrk50 = HLT_MET105_IsoTrk50;


    //*************************************************************************
    //Start Object Selection
    //*************************************************************************
    /// HSCP stuff
    tree->nHSCP = 0;
    for(unsigned int i = 0; i < Pt->size(); i++)
    {
      tree->HSCP_pt[tree->nHSCP] = Pt->at(i);
      tree->HSCP_eta[tree->nHSCP] = eta->at(i);
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

      tree->nHSCP++;
    }


    // muons


    // jets





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
