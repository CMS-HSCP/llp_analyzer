RazorAnalyzer
=============

Class for analyzing the 2015 razor ntuples

Setup
-------------

    cmsrel CMSSW_10_6_30
    cd CMSSW_10_6_30/src
    cmsenv
    git clone git@github.com:CMS-HSCP/llp_analyzer.git
    cd llp_analyzer
    make
  
Defining a new analysis
-------------
1) Copy analyzers/DummyAnalyzer.cc and replace each instance of "DummyAnalyzer" with the name of your desired analyzer.
   Modify the body of the Analyze function to define your analyzer's behavior.
   DO NOT need to write a header file for the analyzer class; the Makefile will generate one for you automatically.  

2) Do `make`.  This will create an executable `bin/Run<name of your analyzer>`. You can execute your analysis using this program directly or by calling it via the `RazorRun` script. 

Running
------------
After compiling, 

    ./RazorRun <list of input files> <name of your analyzer> <options>
  
Example: to execute a HSCP analyzer on gluino 1.8 TeV:

    ./RazorRun lists/ntuples/V1p0/MC_UL18/HSCPgluino_M-1800_TuneCP5_13TeV-pythia8.txt HSCPAnalyzer

The "options" are the following:
    
    -d   --isData
    -f=  --outputFile=<output filename> (optional)
    -n=  --optionNumber=<option number> (optional)
    -l=  --optionLabel=<option Label> (optional)
    -h   --help


## Run the llp_analyzer
    ./RazorRun_T2 <list of input files> llp_vH -d=${isData} -n=${option} -f=${outputfile} -l=${tag}
* ```isData``` is ```yes``` or ```no```
* currently the options and tags don't do anything
* list of input files are stored in ```lists/ntuples``` 


submit_condor_caltech.py  submit_normalize_caltech.py

### Create input list
Run  ```scripts/make_input_list.py``` to create lists of input files for data and MC

### Submit condor jobs on Caltech tier2
Before submitting jobs, make sure proxy and CMSSW environment is setup.

* run the analyzer for signal, bkg or data:
	* ```scripts_condor/submit_condor_caltech.py```
	* make sure analyzer is correct
	* filesPerJob determines how the jobs are split, need to run test runs interactively to know how long the jobs would take
	* Check ```inputfilelist``` where the list you created in step 1 is stored, and ```output``` where you want the output do be stored


Normalizing the processed ntuples
------------
The NormalizeNtuple macro opens a specified set of files and adds a 'weight' branch to each TTree in each file.  The value of 'weight' is the same for all events in a tree and is equal to lumi * CrossSection/NEvents, where NEvents is the total number of events processed for the given dataset, and lumi is the luminosity normalized to.  The cross sections can be found in the file ```data/xSections.dat```.  To run NormalizeNtuple:

    ./NormalizeNtuple <input file list> [lumi (/pb)]

* input file list format, each line is:  <dataset name> <root file name>
* Make sure the dataset being processed have xSections in ```data/xSections.dat```
* Create input file list using ```scripts/create_normalize_txt.py```

* submit jobs to confor to hadd and normalize ntuples for all samples:
  * ```scripts_condor/submit_normalize_caltech.py```
  * Check ```outputDir``` and ```inputDir``` are correct
  
