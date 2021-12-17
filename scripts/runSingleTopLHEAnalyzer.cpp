
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
  int nSamples = 4;
  int modes = 1;
  string* suffix = new string[nSamples];
  string* prefix = new string[modes];

prefix[0]="t_channel_";
//prefix[1]="madspin_";
//prefix[0]="EVENTS_t_channel_madspin_Rwt_ctwi_cbwi_";

suffix[0] = "cbw_m2";
suffix[1] = "cbw_m1";
suffix[2] = "cbw_p1";
suffix[3] = "cbw_p2";
// suffix[4] = "ctwi_p1";
// suffix[5] = "ctwi_p2";
// suffix[6] = "ctwi_n1";
// suffix[7] = "ctwi_n2";
// suffix[8] = "SM";
//suffix[9] = "cbwi_n1";
//suffix[10] = "cbwi_p1";
//suffix[11] = "cbwi_p2";
//suffix[12] = "cbwi_p5";
//suffix[13] = "ctwi_m1";
//suffix[14] = "ctwi_m2";
//suffix[15] = "ctwi_p1";



  TFile** fInput = new TFile*[nSamples];
  TTree** tInput = new TTree*[nSamples];
  SingleTopLHEAnalyzer** singleTopLHEAnalyzer = new SingleTopLHEAnalyzer*[nSamples];

  for (int j=0; j<modes; j++)
  {
    for (int i=0; i<nSamples; i++)
    {
      //string nb = to_string(i+1);
      //suffix[i]= nb;
      cout<<"Cooking "+prefix[j] + suffix[i]+".root"<<endl;

     ///////////////////////////Choose Input Path files//////////////////////////////////////
     string inputPath = "/gridgroup/cms/greenberg/Pheno/MG5_aMC_v2_6_7/top_EFT/LHEconvertedFiles/";
     string inputName = inputPath + prefix[j] + suffix[i] + ".root";
      
      
     string outputName = "output_" +prefix[j] + suffix[i] +  ".root";   
     fInput[i] = new TFile(inputName.c_str(),"READ");
     tInput[i] = (TTree*) fInput[i]->Get("LHEF");
     singleTopLHEAnalyzer[i] = new SingleTopLHEAnalyzer(tInput[i]);
     singleTopLHEAnalyzer[i]->Loop();
      
     ///////////////////////////Choose Output Path files//////////////////////////////////////
     string outputPath = "/gridgroup/cms/greenberg/Pheno/eft_lhe_analyzer/data/madgraph/CutData/Standalone/";
     string commandline = "mv output.root " + outputPath + outputName;

     system(commandline.c_str());
    }
  }



return 0;
}
