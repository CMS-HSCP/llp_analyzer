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
