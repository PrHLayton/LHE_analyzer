
#include "../include/SingleTopLHEAnalyzer.hpp"
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <sstream>

using namespace std;


int main()
{
  int nSamples = 1;
  int modes = 1;
  string* suffix = new string[nSamples];
  string* prefix = new string[modes];

prefix[0]="ST_t_channel_atLO_500kevents_cbw_m2_MG5_dim6top"; // Name of input File

suffix[0] ="";

  TFile** fInput = new TFile*[nSamples];
  TTree** tInput = new TTree*[nSamples];
  SingleTopLHEAnalyzer** singleTopLHEAnalyzer = new SingleTopLHEAnalyzer*[nSamples];

  for (int j=0; j<modes; j++)
  {
    for (int i=0; i<nSamples; i++)
    {
      string nb = to_string(i+1);
      suffix[i]= "";
      cout<<"Cooking "+prefix[j] + suffix[i]+".root"<<endl;

     ///////////////////////////Choose Input Path files//////////////////////////////////////
     string inputPath = "/gridgroup/cms/fillaudeau/eft_lhe_analyzer/data/madgraph/";
     string inputName = inputPath + prefix[j] + suffix[i] + ".root";
      
      
     string outputName = "output_" +prefix[j] + suffix[i] +  ".root";   
     fInput[i] = new TFile(inputName.c_str(),"READ");
     tInput[i] = (TTree*) fInput[i]->Get("LHEF");
     singleTopLHEAnalyzer[i] = new SingleTopLHEAnalyzer(tInput[i]);
     singleTopLHEAnalyzer[i]->Loop();
      
     ///////////////////////////Choose Output Path files//////////////////////////////////////
     string outputPath = "/gridgroup/cms/fillaudeau/eft_lhe_analyzer/data/madgraph/output/";
     string commandline = "mv output.root " + outputPath + outputName;

     system(commandline.c_str());
    }
  }
  
return 0;
}
