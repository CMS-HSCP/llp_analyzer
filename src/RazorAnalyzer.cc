#include "RazorAnalyzer.h"
#include "TLorentzVector.h"

using namespace std;

RazorAnalyzer::RazorAnalyzer(TTree *tree) : HscpCandidates(tree)
{
    //turn off all branches
    fChain->SetBranchStatus("*", 1);
}

RazorAnalyzer::~RazorAnalyzer()
{

}

void RazorAnalyzer::Analyze(bool isData, int option, string outputFileName, string label) {
    cout << "Analyze method called on base RazorAnalyzer instance.  Parameters were: " << isData << " " << option << " " << outputFileName << " " << label << endl;
}

//NOTE: the functions below need to be maintained by hand.  If variables are added or removed from the ntuple, these functions need to be updated to reflect the changes.

double RazorAnalyzer::deltaPhi(double phi1, double phi2) {
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}

double RazorAnalyzer::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
}

TLorentzVector RazorAnalyzer::makeTLorentzVector(double pt, double eta, double phi, double energy){
    TLorentzVector vec;
    vec.SetPtEtaPhiE(pt, eta, phi, energy);
    return vec;
}

TLorentzVector RazorAnalyzer::makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass){
    TLorentzVector vec;
    vec.SetPtEtaPhiM(pt, eta, phi, mass);
    return vec;
}
