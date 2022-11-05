#ifndef DEF_HSCPAnalyzer
#define DEF_HSCPAnalyzer

#include "RazorAnalyzer.h"

class HSCPAnalyzer: public RazorAnalyzer {
    public: 
        HSCPAnalyzer(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
