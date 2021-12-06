#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <sstream>
#include "TStyle.h"
#include <iostream>
#include <TStyle.h>
#include <string>
#include <TGraph.h>
#include <TF1.h>
#include <fstream>
#include <TSystem.h>
using namespace std;

bool reweight = false;

bool normalization = true; //Turn on normalization to 1
bool log_fity = false; //Turn on log scale in Y axis
bool diff_xsection = false; //Change this to have histograms normalized to Xsection
bool flows = false; //Change this to Active/Deactivate Overflow and Underflow for all Histograms

TH1D* GetHistoWeight(TTree* t, string variable, int nbins, double xmin, double xmax, string cut, string name)
{
        string sxmin, sxmax, snbins;
        stringstream ss[3];

        ss[0] << xmin;
        ss[0] >> sxmin;
        ss[1] << xmax;
        ss[1] >> sxmax;
        ss[2] << nbins;
        ss[2] >> snbins;

        string variablenew = variable + " >> h(" + snbins + "," + sxmin + "," + sxmax + ")";

        string cutnew = "1 * (" + cut + ")";
        //string cutnew = "(" + cut + ")";

//

        t->Draw(variablenew.c_str(), cutnew.c_str());
        TH1D *histo = (TH1D*)gDirectory->Get("h");
  if(flows)
  {
		if (histo->GetEntries()==0) return histo;
  
		double underflow = histo->GetBinContent(0);
		//cout << "underflow="<<underflow<<endl;
		double val = 0;
		if (underflow>0) {
			val = histo->GetBinContent(1);
      //cout<<"val= "<<val<<endl;
			histo->SetBinContent(1, val+underflow);
			 histo->SetBinContent(0, 0);
		}
		double overflow = histo->GetBinContent(nbins+1);
    //cout<<"overflow= "<<overflow<<endl;
		if (overflow>0) {
		  val = histo->GetBinContent(nbins);
		  histo->SetBinContent(nbins+1, 0);
		  histo->SetBinContent(nbins, val+overflow);
		}
  }
  //cout << "Area="<<histo->Integral()<<endl;
	//cout << "Nevents="<<histo->GetEntries()<<endl;
        histo->SetName(name.c_str());
        histo->SetTitle(name.c_str());

        return histo;
}

void Ratio_EFT_SM(TTree* t1, TTree* t2, TTree* t3, TTree* t4, TTree* t5, string variable, string EFT, int nbins, double xmin, double xmax, string selection1, string selection2, string selection3, string selection4, string selection5, string legendX, string legendY, string file_name, string lepton)
{
  TFile* file_output = new TFile(file_name.c_str(),"RECREATE");
  TTree* tree_file = new TTree("events","events");

  tree_file->Branch("");

  TH1D* Histo_SM = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection1, "Histo_SM");
  TH1D* Histo_EFT_p1 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection2, "Histo_EFT_p1");
  TH1D* Histo_EFT_p2 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection3, "Histo_EFT_p2");
  TH1D* Histo_EFT_m1 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection4, "Histo_EFT_m1");
  TH1D* Histo_EFT_m2 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection5, "Histo_EFT_m2");

  TH1D* Ratio_p2_SM = (TH1D*)Histo_EFT_p2->Clone("Ratio_p2_SM");
  TH1D* Ratio_p1_SM = (TH1D*)Histo_EFT_p1->Clone("Ratio_p1_SM");
  TH1D* Ratio_m1_SM = (TH1D*)Histo_EFT_m1->Clone("Ratio_m1_SM");
  TH1D* Ratio_m2_SM = (TH1D*)Histo_EFT_m2->Clone("Ratio_m2_SM");


  //X_sections values obtained with MG5
  double Xs_n5, Xs_n2, Xs_n1 ,Xs_SM = 35.2863;
  double Xs_p1, Xs_p2, Xs_p5;

  if(EFT == "cbwi")
  {
    Xs_n5 = 48.5315;
    Xs_n2 = 37.3699;
    Xs_n1 = 35.8251;
    Xs_p1 = 35.8091;
    Xs_p2 = 37.3955;
    Xs_p5 = 48.5184;
  }
  if(EFT == "ctwi")
  {
    Xs_n5 = 48.5683;
    Xs_n2 = 37.3838;
    Xs_n1 = 35.8029;
    Xs_p1 = 35.8126;
    Xs_p2 = 37.3893;
    Xs_p5 = 48.5945;
  }

  int N = 1000000; //sum des pts initiaux 
  Ratio_p2_SM->Scale(Xs_p2/N);
  Ratio_p1_SM->Scale(Xs_p1/N);
  Ratio_m2_SM->Scale(Xs_n2/N);
  Ratio_m1_SM->Scale(Xs_n1/N);
  Histo_SM->Scale(Xs_SM/N);

  //weights_sum[5] = {SM, p2, p1, m2, m1}

  //double weights_sum[5] = {77672, 82091.3, 79047.4, 82123.4, 79057.2}; //ctwi from ctwi2p5
  //double weights_sum[5] = {77672, 81820.2, 77911.8, 81826.8, 77912.2}; //cbwi from ctwi2p5
  //nb_weights = 230054;
  
  //double weights_sum[5] = {83125.7, 87952.1, 84560.8, 87933.9, 84570.7}; //ctwi from cbwi_p3
  //double weights_sum[5] = {83125.7, 84795.9, 83366.5, 84795.8, 83368.8}; //cbwi from cbwi_p3
  //nb_weights = 241499;

  //double weights_sum[5] = {73926.8, 77781.1, 74899.4, 77732.2, 74872.5}; //ctwi from cbwi2p5
  //double weights_sum[5] = {73926.8, 77618.7, 74852, 77619.4, 74852.4}; //cbwi from cbwi2p5
  //int nb_weights = 226646;
  
  //Histo_SM->Scale(weights_sum[0]/nb_weights);
  //Ratio_p2_SM->Scale(weights_sum[1]/nb_weights);
  //Ratio_p1_SM->Scale(weights_sum[2]/nb_weights);
  //Ratio_m2_SM->Scale(weights_sum[3]/nb_weights);
  //Ratio_m1_SM->Scale(weights_sum[4]/nb_weights);


  Ratio_p2_SM->Divide(Histo_SM);
  Ratio_p1_SM->Divide(Histo_SM);
  Ratio_m2_SM->Divide(Histo_SM);
  Ratio_m1_SM->Divide(Histo_SM);
  Histo_SM->Divide(Histo_SM);



  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  //double max = (Ratio_p2_SM->GetMaximum()>Ratio_p1_SM->GetMaximum()) ? Ratio_p2_SM->GetMaximum() : Ratio_p1_SM->GetMaximum();
  double max = Ratio_p2_SM->GetMaximum();
  //Ratio_p2_SM->SetAxisRange((2-max)*1.1, max*1.25, "Y");
  Ratio_p2_SM->SetAxisRange(2-max, max*1.1, "Y");
  Ratio_p2_SM->SetXTitle(legendX.c_str());
  Ratio_p2_SM->SetYTitle(legendY.c_str());
  Ratio_p2_SM->SetLineColor(kRed);
  Ratio_p2_SM->SetLineWidth(2);
  Ratio_p2_SM->Draw("");

  Ratio_p1_SM->SetLineColor(kBlue);
  Ratio_p1_SM->SetLineWidth(2);
  Ratio_p1_SM->Draw("SAME");

  Ratio_m2_SM->SetLineColor(kOrange);
  Ratio_m2_SM->SetLineWidth(2);
  Ratio_m2_SM->Draw("SAME");

  Ratio_m1_SM->SetLineColor(kGreen);
  Ratio_m1_SM->SetLineWidth(2);
  Ratio_m1_SM->Draw("SAME");

  Histo_SM->SetLineColor(kBlack);
  Histo_SM->SetLineWidth(2);
  Histo_SM->Draw("SAME");

  double lx0 = 0.6;
  double ly0 = 0.6;
  double lx1 = 0.99;
  double ly1 = 0.99;

  string legendtitle = "Value of the EFT";

  string eft_p2_legend = EFT + "/#\Lambda^{2} = 2 (TeV^{-2})";
  string eft_p1_legend = EFT + "/#\Lambda^{2} = 1 (TeV^{-2})";
  string eft_m1_legend = EFT + "/#\Lambda^{2} = -1 (TeV^{-2})";
  string eft_m2_legend = EFT + "/#\Lambda^{2} = -2 (TeV^{-2})";
  string SM_legend = "SM";

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.05);
  legend->AddEntry(Ratio_p2_SM->GetName(), eft_p2_legend.c_str(), "l");
  legend->AddEntry(Ratio_p1_SM->GetName(), eft_p1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m1_SM->GetName(), eft_m1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m2_SM->GetName(), eft_m2_legend.c_str(), "l");
  legend->AddEntry(Histo_SM->GetName(), SM_legend.c_str(), "l");

  legend->Draw("SAME");
  //Canvas->Print(file_name.c_str());


  //string TFile_name = "/results/ratio_madgraph/signal_proc_"+variable;
  //TFile* ratio_file = new TFile(("signal_proc_" + variable + "_" + EFT + "_" + lepton + "_TF1.root").c_str(), "RECREATE");
  TFile* ratio_file = new TFile(("results/EFT_vs_SM/Rwgt_cbwi2p5_XSect/signal_proc_" + variable + "_" + EFT + "_" + lepton + "_"+ to_string(nbins) + "bins" + "_Rwgt_cbwi2p5.root").c_str(), "RECREATE");
  ratio_file->cd();
  
  TGraph** ratio_Histo = new TGraph*[nbins];

  for (int i = 1 ; i<=nbins ; i++)
  {
    string number_plot = to_string(i);
    string file_name_eft = file_name;
    file_name_eft.insert(file_name.size()-4,"_"+number_plot);
    Canvas->Clear();
    //Canvas->SetLogy();
    ratio_Histo[i] = new TGraph(5);
    TF1* ratio_formula = new TF1(("bin_content_"+number_plot).c_str(),"[0]+[1]*x+[2]*x*x", -25 , 25);
    
    
     
    ratio_Histo[i]->SetPoint(0,-2,Ratio_m2_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(1,-1,Ratio_m1_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(2,0,Histo_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(3,1,Ratio_p1_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(4,2,Ratio_p2_SM->GetBinContent(i));
    
    ratio_Histo[i]->SetMarkerStyle(kStar);

    tree_file->Fill();
    ratio_Histo[i]->Fit(ratio_formula);
    ratio_Histo[i]->Draw();
    //Canvas->Print(("test.pdf"+number_plot).c_str());
    

    TLegend* legend2 = new TLegend(0.6, 0.7, 0.89, 0.89, "");
    legend2->SetTextSize(0.035);
    if(variable == "PhiStar" && nbins==5)
    {
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "0 < #phi^{*} < 1.25" ,"r");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "1.25 < #phi^{*} < 2.5" ,"r");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "2.5 < #phi^{*} < 4" ,"r");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "4 < #phi^{*} < 5.5" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "5.5 < #phi^{*} < 6.28" ,"r");
    }
    
    if(variable == "cosThetaStar" && nbins==5)
    {
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "-1 < cos(#theta^{*}) < -0.6" ,"r");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "-0.6 < cos(#theta^{*}) < -0.2" ,"r");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "-0.2 < cos(#theta^{*}) < 0.2" ,"r");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "0.2 < cos(#theta^{*}) < 0.6" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "0.6 < cos(#theta^{*}) < 1" ,"r");

    }
  
    legend2->Draw("SAME");   


    ratio_Histo[i]->GetYaxis()->SetRangeUser(0,Ratio_p2_SM->GetBinContent(i)*2);
    ratio_Histo[i]->GetYaxis()->SetTitle("EFT/SM");
    if(EFT == "ctwi") ratio_Histo[i]->GetXaxis()->SetTitle("C_{tW}^{I}");
    if(EFT == "cbwi") ratio_Histo[i]->GetXaxis()->SetTitle("C_{bW}^{I}");

    ratio_Histo[i]->SetTitle("");



    Canvas->Print(file_name_eft.c_str());

    //ratio_Histo[i]->Write(("bin_content_part_"+number_plot).c_str());
    
    //Canvas->Print(file_name_eft.c_str());

    ratio_formula->Write();
  }


  file_output->Write();
  ratio_file->Close();
}

void Compare_3Histos(TTree* t1, TTree* t2, TTree* t3, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string EFT, string Name){
  
  Name += variable;

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "Histo_3");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");

  if(log_fity) 
  {
    Canvas->SetLogy();
    Name += "_log"; 
    legendY += " log";
  }
  if(normalization) 
  {
    Double_t integral1 = 1/Histo_1->Integral();
    Histo_1->Scale(integral1);
    
    Double_t integral2 = 1/Histo_2->Integral();
    Histo_2->Scale(integral2);
    
    Double_t integral3 = 1/Histo_3->Integral();
    Histo_3->Scale(integral3);

    Name += "_normalized";
    legendY += " normalized";
  }
/*
  if(diff_xsection)
  {
    //int N1 = Histo_1->GetEntries() , N2 = Histo_2->GetEntries() , N3 = Histo_3->GetEntries();
    int N_initial = 1000000;
    double scale1 = 35.2863/N_initial , scale2 = 37.3838/N_initial , scale3 = 37.3893/N_initial ;
    Histo_1->Scale(scale1);
    Histo_2->Scale(scale2);
    Histo_3->Scale(scale3);

    Name += "_xsection" ;
    legendY += " #sigma";
  }

  if (legendX == "Top Mass (GeV)")
    {
      if (normalization == false) Histo_1->SetAxisRange(0.0001,max*1.5,"Y");
    }
  if (legendX == "W Mass (GeV)")
  {
    Histo_1->SetAxisRange(0.001,max*1.5,"Y");
  }

  if(log_fity==false)
  {
    if(variable == "PhiStar")
    {
      Histo_1->SetMinimum(0);
      double max_phi = Histo_1->GetMaximum();
      Histo_1->SetMaximum(max_phi*1.5);
    }
    if(variable == "cosThetaStar")
    {
      Histo_1->SetMinimum(0);
      double max_cos = Histo_1->GetMaximum();
      Histo_1->SetMaximum(max_cos*1.28);
    }
  }

  else
  {
    if(variable == "PhiStar")
    {
      Histo_1->SetMinimum(0.1);
      double max_1 = Histo_1->GetMaximum();
      Histo_1->SetMaximum(max_1*40);

    }
    if(variable =="cosThetaStar")
    {
      Histo_1->SetMinimum(0.1);     
      double max_1 = Histo_1->GetMaximum();
      Histo_1->SetMaximum(max_1*40);
    }
  }
  */
  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  Histo_1->SetMaximum(max*1.4 );

  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  
  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

  
  Histo_3->SetLineColor(kOrange);
  Histo_3->SetLineWidth(2);
  Histo_3->Draw("SAME");


 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.7;
	 lx1 = 0.6;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.6;
	 ly0 = 0.75; //Lower Y
	 lx1 = 0.99;
	 ly1 = 0.99; //Upper Y
  }
  
  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.04);
  if(EFT == "cbwi")
  {
    legend->AddEntry(Histo_1->GetName(), "C_{bW}^{I}/#\Lambda^{2} = 0 (TeV^{-2})","l");
    legend->AddEntry(Histo_2->GetName(), "C_{bW}^{I}/#\Lambda^{2} = -2 (TeV^{-2})","l");
    legend->AddEntry(Histo_3->GetName(), "C_{bW}^{I}/#\Lambda^{2} = 2 (TeV^{-2})","l");
  } 

  if(EFT == "ctwi")
  {
    legend->AddEntry(Histo_1->GetName(), "C_{tW}^{I}/#\Lambda^{2} = 0 (TeV^{-2})","l");
    legend->AddEntry(Histo_2->GetName(), "C_{tW}^{I}/#\Lambda^{2} = -2 (TeV^{-2})","l");
    legend->AddEntry(Histo_3->GetName(), "C_{tW}^{I}/#\Lambda^{2} = 2 (TeV^{-2})","l");
  }

  if(EFT == "")
  {
    legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(),"l");
    legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(),"l");
    legend->AddEntry(Histo_3->GetName(), legendEntry3.c_str(),"l");
  } 
  //legend->SetLegendSize(0.5);

  legend->Draw("SAME");

  Name += ".pdf";
  Canvas->Print(Name.c_str());
}

void Compare_1Histos(TTree* t1, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  Histo_1->SetStats(kFALSE);

  double max = Histo_1->GetMaximum();
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  double lx0, ly0, lx1, ly1;
   if (legendPlace=="legendUpLeft"){
 	 lx0 = 0.2;
 	 ly0 = 0.75;
 	 lx1 = 0.5;
 	 ly1 = 0.95;
   }
    if (legendPlace=="legendUpRight"){
 	 lx0 = 0.75;
 	 ly0 = 0.75;
 	 lx1 = 0.99;
 	 ly1 = 0.99;
   }

   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);
   legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");

   legend->Draw("SAME");

   Canvas->Print(Name.c_str());
}

void Compare_2Histos(TTree* t1, TTree* t2, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string Name){

  Name += "_" + variable;
  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);

  if(normalization)
  {
    double a = Histo_1->Integral();
    double b = Histo_2->Integral();
    Histo_1->Scale(1/a);
    Histo_2->Scale(1/b);
    Name += "_normalized";
  }
  //cout<<legendX<<endl;
  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.3,"Y");
  if (legendX == "Top Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.0001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  if (legendX == "W Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");
  legend->Draw("SAME");
  
  legend->SetTextSize(0.035);


  Name += ".pdf";
  Canvas->Print(Name.c_str());

  cout << "Histo1 mean: "<<Histo_1->GetMean()<<endl;
  cout << "Histo2 mean: "<<Histo_2->GetMean()<<endl;

}

void Compare_4Histos(TTree* t1, TTree* t2, TTree* t3, TTree* t4, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string legendEntry4, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "Histo_3");
  TH1D* Histo_4 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "Histo_4");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);
  Histo_4->SetStats(kFALSE);

  double a = Histo_1->Integral();
  double b = Histo_2->Integral();
  double c = Histo_3->Integral();
  double d = Histo_4->Integral();

  Histo_1->Scale(1/a);
  Histo_2->Scale(1/b);
  if (c>0) Histo_3->Scale(1/c);
  Histo_4->Scale(1/d);
  cout << "a="<<a<<" b="<<b<<" c="<<c<<endl;

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  if (c>0) max = (max>Histo_3->GetMaximum()) ? max : Histo_3->GetMaximum();
  if (d>0) max = (max>Histo_4->GetMaximum()) ? max : Histo_4->GetMaximum();
  cout << "max="<<max<<endl;

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  //Canvas->SetLogx();
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  if (legendX == "Top Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.0001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  if (legendX == "W Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

  if (c>0){
    Histo_3->SetLineColor(kOrange);
    Histo_3->SetLineWidth(2);
    Histo_3->Draw("SAME");
  }

  Histo_4->SetLineColor(kGreen);
  Histo_4->SetLineWidth(2);
  Histo_4->Draw("SAME");

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), "MadGraph EFT = 0", "l");
  legend->AddEntry(Histo_2->GetName(), "MadSpin EFT = 0", "l");
  legend->AddEntry(Histo_3->GetName(), "MadGraph EFT = 2", "l");
  legend->AddEntry(Histo_4->GetName(), "MadSpin EFT = 2", "l");
  legend->Draw("SAME");

  Canvas->Print(Name.c_str());

  cout << "Histo1 mean: "<<Histo_1->GetMean()<<endl;
  cout << "Histo2 mean: "<<Histo_2->GetMean()<<endl;
  if (c>0) cout << "Histo3 mean: "<<Histo_3->GetMean()<<endl;

}

//New Fcts

void Ratio_EFT_SM_7pts(TTree* t1, TTree* t2, TTree* t3, TTree* t4, TTree* t5, TTree* t6, TTree* t7, string variable, string EFT, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string file_name, string lepton)
{
  //Fcts to Plot EFT/SM for cbwi/ctwi = [-5;5]

  TFile* file_output = new TFile(file_name.c_str(),"RECREATE");
  TTree* tree_file = new TTree("events","events");

  tree_file->Branch("");

  TH1D* Histo_SM = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p1 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p2 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m1= GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m2 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p5 = GetHistoWeight(t6, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m5 = GetHistoWeight(t7, variable, nbins, xmin, xmax, selection, "");

  TH1D* Ratio_p2_SM = (TH1D*)Histo_EFT_p2->Clone("Ratio_p2_SM");
  TH1D* Ratio_p1_SM = (TH1D*)Histo_EFT_p1->Clone("Ratio_p1_SM");
  TH1D* Ratio_m1_SM = (TH1D*)Histo_EFT_m1->Clone("Ratio_m1_SM");
  TH1D* Ratio_m2_SM = (TH1D*)Histo_EFT_m2->Clone("Ratio_m2_SM");
  TH1D* Ratio_p5_SM = (TH1D*)Histo_EFT_p5->Clone("Ratio_p5_SM");
  TH1D* Ratio_m5_SM = (TH1D*)Histo_EFT_m5->Clone("Ratio_m5_SM");


  //Ratio_p2_SM->Divide(Histo_SM);
  //Ratio_p1_SM->Divide(Histo_SM);
  //Ratio_m2_SM->Divide(Histo_SM);
  //Ratio_m1_SM->Divide(Histo_SM);
  //Ratio_p5_SM->Divide(Histo_SM);
  //Ratio_m5_SM->Divide(Histo_SM);

  bool div_Xsection = true;
  if(div_Xsection)
  {
    double Xs_n5, Xs_n2, Xs_n1 ,Xs_SM = 35.2863;
    double Xs_p1, Xs_p2, Xs_p5;

    if(EFT == "cbwi")
    {
      Xs_n5 = 48.5315;
      Xs_n2 = 37.3699;
      Xs_n1 = 35.8251;
      Xs_p1 = 35.8091;
      Xs_p2 = 37.3955;
      Xs_p5 = 48.5184;
    }
    if(EFT == "ctwi")
    {
      Xs_n5 = 48.5683;
      Xs_n2 = 37.3838;
      Xs_n1 = 35.8029;
      Xs_p1 = 35.8126;
      Xs_p2 = 37.3893;
      Xs_p5 = 48.5945;
    }
    int N = 1000000;
    Ratio_p2_SM->Scale(Xs_p2/N);
    Ratio_p1_SM->Scale(Xs_p1/N);
    Ratio_m2_SM->Scale(Xs_n2/N);
    Ratio_m1_SM->Scale(Xs_n1/N);
    Ratio_p5_SM->Scale(Xs_p5/N);
    Ratio_m5_SM->Scale(Xs_n5/N);
    Histo_SM->Scale(Xs_SM/N);
  }

  bool ratio = true;
  if(ratio)
  {
    Ratio_p2_SM->Divide(Histo_SM);
    Ratio_p1_SM->Divide(Histo_SM);
    Ratio_m2_SM->Divide(Histo_SM);
    Ratio_m1_SM->Divide(Histo_SM);
    Ratio_p5_SM->Divide(Histo_SM);
    Ratio_m5_SM->Divide(Histo_SM);
    Histo_SM->Divide(Histo_SM);
  }
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  Ratio_p2_SM->SetLineColor(kRed);
  Ratio_p2_SM->SetLineWidth(2);

  Ratio_p1_SM->SetLineColor(kBlue);
  Ratio_p1_SM->SetLineWidth(2);

  Ratio_m2_SM->SetLineColor(kOrange);
  Ratio_m2_SM->SetLineWidth(2);

  Ratio_m1_SM->SetLineColor(kBlack);
  Ratio_m1_SM->SetLineWidth(2);

  Ratio_p5_SM->SetLineColor(kGreen);
  Ratio_p5_SM->SetLineWidth(2);

  Ratio_m5_SM->SetLineColor(kViolet);
  Ratio_m5_SM->SetLineWidth(2);


  Ratio_p2_SM->SetXTitle(legendX.c_str());
  Ratio_p2_SM->SetYTitle(legendY.c_str());

  double max = (Ratio_p5_SM->GetMaximum()>Ratio_m5_SM->GetMaximum()) ? Ratio_p5_SM->GetMaximum() : Ratio_m5_SM->GetMaximum();

  Ratio_p2_SM->SetAxisRange((2-max)*1.1, max*1.25, "Y");
  Ratio_p2_SM->Draw("");
  Ratio_p1_SM->Draw("SAME");
  Ratio_m1_SM->Draw("SAME");
  Ratio_m2_SM->Draw("SAME");
  Ratio_p5_SM->Draw("SAME");
  Ratio_m5_SM->Draw("SAME");

  double lx0 = 0.6;
  double ly0 = 0.6;
  double lx1 = 0.99;
  double ly1 = 0.99;
  string legendtitle = "Value of the EFT";

  string eft_p2_legend = EFT + "/#\Lambda^{2} = 2 (TeV^{-2})";
  string eft_p1_legend = EFT + "/#\Lambda^{2} = 1 (TeV^{-2})";
  string eft_m1_legend = EFT + "/#\Lambda^{2} = -1 (TeV^{-2})";
  string eft_m2_legend = EFT + "/#\Lambda^{2} = -2 (TeV^{-2})";
  string eft_p5_legend = EFT + "/#\Lambda^{2} = 5 (TeV^{-2})";
  string eft_m5_legend = EFT + "/#\Lambda^{2} = -5 (TeV^{-2})";



  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.05);
  legend->AddEntry(Ratio_p2_SM->GetName(), eft_p2_legend.c_str(), "l");
  legend->AddEntry(Ratio_p1_SM->GetName(), eft_p1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m1_SM->GetName(), eft_m1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m2_SM->GetName(), eft_m2_legend.c_str(), "l");
  legend->AddEntry(Ratio_p5_SM->GetName(), eft_p5_legend.c_str(), "l");
  legend->AddEntry(Ratio_m5_SM->GetName(), eft_m5_legend.c_str(), "l");

  legend->Draw("SAME");

  TGraph** ratio_Histo = new TGraph*[nbins];

    //string TFile_name = "/results/ratio_madgraph/signal_proc_"+variable;
  TFile* ratio_file = new TFile(("signal_proc_"+variable+"_"+EFT+"_"+lepton+"_TF1.root").c_str(),"RECREATE");
  ratio_file->cd();
  Canvas->Print(file_name.c_str());


  for (int i = 1 ; i<=nbins ; i++)
  {
    string number_plot = to_string(i); //Conversion of an int into a stream
    string file_name_eft = file_name;
    file_name_eft.insert(file_name.size()-4,"_"+number_plot);
    Canvas->Clear();
    //Canvas->SetLogy();
    ratio_Histo[i] = new TGraph(5);
    TF1* ratio_formula = new TF1(("bin_content_par1_"+number_plot).c_str(),"[0]+[1]*x+[2]*x*x", -3 , 3);
    
    bool Get_bins = false;
    if(Get_bins)
    {
      ratio_Histo[i]->SetPoint(0,-2,Histo_EFT_m2->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(1,-1,Histo_EFT_m1->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(2,0,Histo_SM->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(3,1,Histo_EFT_p1->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(4,2,Histo_EFT_p2->GetBinContent(i)/Histo_SM->GetBinContent(i));
    }
    else
    { 
      ratio_Histo[i]->SetPoint(0,-5,Ratio_m5_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(1,-2,Ratio_m2_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(2,-1,Ratio_m1_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(3,0,Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(4,1,Ratio_p1_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(5,2,Ratio_p2_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(6,5,Ratio_p5_SM->GetBinContent(i));
    }
    ratio_Histo[i]->SetMarkerStyle(kStar);

    tree_file->Fill();
    ratio_Histo[i]->Fit(ratio_formula);
    ratio_Histo[i]->Draw();
    //Canvas->Print(("test.pdf"+number_plot).c_str());

      TLegend* legend2 = new TLegend(0.49, 0.73, 0.75, 0.86, "");
      legend2->SetTextSize(0.035);
    if(variable == "PhiStar")
    {
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "0 < #phi^{*} < 1.25" ,"r");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "1.25 < #phi^{*} < 2.5" ,"r");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "2.5 < #phi^{*} < 4" ,"r");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "4 < #phi^{*} < 5.5" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "5.5 < #phi^{*} < 6.28" ,"r");
    }
    
    if(variable == "cosThetaStar")
    {
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "-1 < cos(#theta^{*}) < -0.6" ,"r");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "-0.6 < cos(#theta^{*}) < -0.2" ,"r");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "-0.2 < cos(#theta^{*}) < 0.2" ,"r");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "0.2 < cos(#theta^{*}) < 0.6" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "0.6 < cos(#theta^{*}) < 1" ,"r");

    }
    
    legend2->Draw("SAME");    

    ratio_Histo[i]->GetYaxis()->SetRangeUser(0,Ratio_p5_SM->GetBinContent(i)*2);
    ratio_Histo[i]->GetYaxis()->SetTitle("EFT/SM");
    if(EFT == "ctwi") ratio_Histo[i]->GetXaxis()->SetTitle("ctwi");
    if(EFT == "cbwi") ratio_Histo[i]->GetXaxis()->SetTitle("cbwi");
    ratio_Histo[i]->SetTitle("");
    Canvas->Print(file_name_eft.c_str());

    ratio_Histo[i]->Write(("bin_content_par1_"+number_plot).c_str());
    //ratio_formula->Write();
  }


  //file_output->Write();
  ratio_file->Close();




}

void Rwgt_vs(string selection, string variable, TTree* t1, TTree* t2,  int nbins, double xmin, double xmax, string EFT, string W_value, string legendtitle, string legendX, string legendY, string legendPlace, string legendEntry1, string legendEntry2, string Name){

  //Fcts to plot base MC Simulations KV data vs Weighted MC simulation KV Data
  Name += "_" + variable; 

  TH1D *Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D *Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, "1", "Histo_2");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);

  

 if(normalization)
 {
   double a = Histo_1->Integral();
   double b = Histo_2->Integral();
   Histo_1->Scale(1/a);
   Histo_2->Scale(1/b);

   Name = Name + "_normalized";
   legendY = legendY + "normalized";
 }

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");


  double lx0, ly0, lx1, ly1;
   if (legendPlace=="legendUpLeft"){
 	 lx0 = 0.2;
 	 ly0 = 0.71;
 	 lx1 = 0.5;
 	 ly1 = 0.95;
   }
    if (legendPlace=="legendUpRight"){
 	 lx0 = 0.75;
 	 ly0 = 0.75;
 	 lx1 = 0.99;
 	 ly1 = 0.99;
   }
   
   if(EFT=="ctwi")
   {
     legendEntry1 = "Rwgt C_{tW}^{I} = " + W_value + " TeV^{-2}";
     legendEntry2 = "C_{tW}^{I} = " + W_value + " TeV^{-2}";
   }
   
   if(EFT=="cbwi")
   {
     legendEntry1 = "Rwgt C_{bW}^{I} = " + W_value + " TeV^{-2}";
     legendEntry2 = "C_{bW}^{I} = " + W_value + " TeV^{-2}";
   }

    if(EFT == "SM")
   {
     legendEntry1 = "Rwgt SM";
     legendEntry2 = "SM";
   }

    if(EFT == "other")
   {
     legendEntry1 = "Dim6 = 2";
     legendEntry2 = "Dim6 = 1";
   }



   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);
   legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
   legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");
   legend->SetTextSize(0.035);

   legend->Draw("SAME");

   cout << "Histo1 Entries: "<<Histo_1->Integral()<<endl;
   cout << "Histo2 Entries: "<<Histo_2->Integral()<<endl;

   Name = Name + "_" +selection + ".pdf"; 

   Canvas->Print(Name.c_str());
}

void Compare_9Histos(int nbins, double xmin, double xmax, string selection1, string selection2, string selection3,string selection4, string selection5, string selection6, string selection7, string selection8, string selection9, TTree* t1, TTree* t2, TTree* t3, TTree* t4, TTree* t5, TTree* t6, TTree* t7, TTree* t8, TTree* t9, string variable1, string variable2, string variable3, string variable4, string variable5, string variable6, string variable7, string variable8, string variable9, string legendX, string legendY, string legendPlace, string legendtitle, string Name){
  
  //Fcts to Plot all weights in the same Canvas

  TH1D *Histo_1 = GetHistoWeight(t1, variable1, nbins, xmin, xmax, selection1, "Histo_1");
  TH1D *Histo_2 = GetHistoWeight(t2, variable2, nbins, xmin, xmax, selection2, "Histo_2");
  TH1D *Histo_3 = GetHistoWeight(t3, variable3, nbins, xmin, xmax, selection3, "Histo_3");
  TH1D *Histo_4 = GetHistoWeight(t4, variable4, nbins, xmin, xmax, selection4, "Histo_4");
  TH1D *Histo_5 = GetHistoWeight(t5, variable5, nbins, xmin, xmax, selection5, "Histo_5");
  TH1D *Histo_6 = GetHistoWeight(t6, variable6, nbins, xmin, xmax, selection6, "Histo_6");
  TH1D *Histo_7 = GetHistoWeight(t7, variable7, nbins, xmin, xmax, selection7, "Histo_7");
  TH1D *Histo_8 = GetHistoWeight(t8, variable8, nbins, xmin, xmax, selection8, "Histo_8");
  TH1D *Histo_9 = GetHistoWeight(t9, variable9, nbins, xmin, xmax, selection9, "Histo_9");


  double max = Histo_1->GetMaximum();
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*2.2,"Y");
  //Histo_SM->SetAxisRange(-0.3,0.5,"X");
  //Histo_SM->SetAxisRange(-10,10,"X");
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

  Histo_3->SetLineColor(kGreen);
  Histo_3->SetLineWidth(2);
  Histo_3->Draw("SAME");

  Histo_4->SetLineColor(kViolet);
  Histo_4->SetLineWidth(2);
  Histo_4->Draw("SAME");

  Histo_5->SetLineColor(kOrange);
  Histo_5->SetLineWidth(2);
  Histo_5->Draw("SAME");

  Histo_6->SetLineColor(kBlue);
  Histo_6->SetLineStyle(3);
  Histo_6->SetLineWidth(2);
  Histo_6->Draw("SAME");

  Histo_7->SetLineColor(kGreen);
  Histo_7->SetLineStyle(3);
  Histo_7->SetLineWidth(2);
  Histo_7->Draw("SAME");

  Histo_8->SetLineColor(kViolet);
  Histo_8->SetLineStyle(3);
  Histo_8->SetLineWidth(2);
  Histo_8->Draw("SAME");

  Histo_9->SetLineColor(kOrange);
  Histo_9->SetLineStyle(3);
  Histo_9->SetLineWidth(2);
  Histo_9->Draw("SAME");


  double lx0, ly0, lx1, ly1;
   if (legendPlace=="legendUpLeft"){
 	 lx0 = 0.2;
 	 ly0 = 0.75;
 	 lx1 = 0.5;
 	 ly1 = 0.95;
   }
    if (legendPlace=="legendUpRight"){
 	 lx0 = 0.75;
 	 ly0 = 0.60;
 	 lx1 = 0.99;
 	 ly1 = 0.99;
   }


   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);
   legend->AddEntry(Histo_1->GetName(), "SM", "l");
   legend->AddEntry(Histo_2->GetName(), "C_{bw}^{I}=-2", "l");
   legend->AddEntry(Histo_3->GetName(), "C_{bw}^{I}=-1", "l");
   legend->AddEntry(Histo_4->GetName(), "C_{bw}^{I}=1", "l");
   legend->AddEntry(Histo_5->GetName(), "C_{bw}^{I}=2", "l");
   legend->AddEntry(Histo_6->GetName(), "C_{tw}^{I}=-2", "l"); 
   legend->AddEntry(Histo_7->GetName(), "C_{tw}^{I}=-1", "l"); 
   legend->AddEntry(Histo_8->GetName(), "C_{tw}^{I}=1", "l"); 
   legend->AddEntry(Histo_9->GetName(), "C_{tw}^{I}=2", "l"); 

   legend->Draw("SAME");


   Name = Name + ".pdf";
   Canvas->Print(Name.c_str());

}


int main (){
int nbfiles = 21;
string suffix[nbfiles];

suffix[0] = "output_top_SM.root"; //SM

suffix[1] = "output_top_cbwi_n2.root"; // cbWi = -2
suffix[2] = "output_top_cbwi_n1.root"; //cbWi = -1
suffix[3] = "output_top_cbwi_p1.root"; //cbWi = 1
suffix[4] = "output_top_cbwi_p2.root"; //cbWi = 2

suffix[5] = "output_top_ctwi_n2.root"; //ctWi = -2
suffix[6] = "output_top_ctwi_n1.root"; //ctWi = -1
suffix[7] = "output_top_ctwi_p1.root"; //ctWi = 1
suffix[8] = "output_top_ctwi_p2.root"; //ctWi = 2

suffix[9] = "output_top_cbwi_n5.root"; //cbWi = -5
suffix[10] = "output_top_cbwi_p5.root"; //cbWi = 5
suffix[11] = "output_top_ctwi_n5.root"; //ctWi = -5
suffix[12] = "output_top_ctwi_p5.root"; //ctWi = 5

suffix[13] = "rwgt_ctwi_p3.root"; //Rwgt at ctwi = 3
suffix[14] = "rwgt_ctwi2p5.root"; //Rwgt at ctwi = 2.5
suffix[15] = "rwgt_cbwi_p3.root"; //Rwgt at cbwi = 3
suffix[16] = "rwgt_cbwi2p5.root"; //Rwgt at cbwi = 2.5

suffix[17] = "output_top_SM_1dim6.root"; //Dim6 = 1 ; SM
suffix[18] = "output_top_ctwi_m2_1dim6.root"; //Dim6 = 1 ; ctwi = -2
suffix[19] = "output_top_ctwi_p2_1dim6.root"; //Dim6 = 1 ; ctwi = 2
suffix[20] = "output_top_cbwi_m2_1dim6.root"; //Dim6 = 1 ; cbwi = -2
//suffix[21] = "output_top_cbwi_p2_1dim6.root"; //Dim6 = 1 ; cbwi = 2

TFile* fInput[nbfiles];
TTree* tInput[nbfiles];
string inputName;

for (int i=0; i<nbfiles; i++)
{
  inputName = "data/madgraph/output/" + suffix[i];
  //inputName = "data/madgraph/output/Streco_selection/" + suffix[i];
  fInput[i] = new TFile(inputName.c_str(),"READ");
  tInput[i] = (TTree*) fInput[i]->Get("Tree");
}


  //////////cbWi Plots//////////
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* [rad]", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "cbwi", "results/STreco_selection/cbWi_");
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "top_mass",20, 164, 180, "1", "Top Mass [GeV]", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "cbwi", "results/cbWi_");
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "cbwi", "results/STreco_selection/cbWi_");
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "top_pt", 20, 0, 400, "1", "Top Pt [GeV]", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "cbwi", "results/cbwi_");
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "W_pt", 20, 0, 120, "1", "W Pt [GeV]", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "cbwi", "results/cbwi_");
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "lepton_pt", 20, 0, 100, "1", "Lepton Pt [GeV]", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "cbwi", "results/cbwi_");
	//Compare_3Histos(tInput[0], tInput[1], tInput[4], "cosTheta", 20, -1, 1, "1", "cos(#theta)", "", "legendUpLeft", "dim6top", suffix[0], suffix[1], suffix[4], "cbwi", "results/STreco_selection/cbWi_");


  //////////ctWi Plots//////////
  //Compare_3Histos(tInput[0], tInput[5], tInput[8], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* [rad]", "", "legendUpRight", "dim6top", suffix[0], suffix[5], suffix[8], "ctwi", "results/STreco_selection/ctWi_");
  //Compare_3Histos(tInput[0], tInput[5], tInput[8], "top_mass",25, 150, 190, "1", "Top Mass [GeV]", "", "legendUpRight", "dim6top", suffix[0], suffix[5], suffix[8], "ctwi", "results/ctWi_");
  //Compare_3Histos(tInput[0], tInput[5], tInput[8], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "", "legendUpRight", "dim6top", suffix[0], suffix[5], suffix[8], "ctwi", "results/STreco_selection/ctWi_");
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "top_pt", 20, 0, 400, "1", "Top Pt [GeV]", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "ctwi", "results/ctwi_");
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "W_pt", 20, 0, 120, "1", "W Pt [GeV]", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "ctwi", "results/ctwi_");
  //Compare_3Histos(tInput[0], tInput[1], tInput[4], "lepton_pt", 20, 0, 100, "1", "Lepton Pt [GeV]", "", "legendUpRight", "dim6top", suffix[0], suffix[1], suffix[4], "ctwi", "results/ctwi_");
	//Compare_3Histos(tInput[0], tInput[5], tInput[8], "cosTheta", 20, -1, 1, "1", "cos(#theta)", "", "legendUpLeft", "dim6top", suffix[0], suffix[5], suffix[8], "ctwi", "results/STreco_selection/ctWi_");

  
  
  ////////////EFT VS SM////////////
  //costheta:
  //Ratio_EFT_SM(tInput[0],tInput[3],tInput[4],tInput[1],tInput[2],"cosTheta","cbwi",5,-1,1,"1","cos#theta","ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosTheta.pdf","lepton");
   
  //costhetaStar:
  //Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[2], tInput[1], "cosThetaStar","cbwi", 5, -1, 1 ,"nature_lepton==1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosThetaStar.pdf", "elec");
  //Ratio_EFT_SM(tInput[0], tInput[7], tInput[8], tInput[6], tInput[5], "cosThetaStar","ctwi", 5, -1, 1 ,"nature_lepton==1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_cosThetaStar.pdf", "elec");
  //Ratio_EFT_SM_7pts(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], tInput[10], tInput[9], "cosThetaStar","cbwi", 5, -1, 1 ,"nature_lepton==1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosThetaStar_pt7.pdf", "elec");

  //PhiStar:
  //Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[2], tInput[1], "PhiStar","cbwi", 5, 0, 2*TMath::Pi(),"1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_PhiStar.pdf","elec");
  //Ratio_EFT_SM(tInput[0], tInput[7], tInput[8], tInput[6], tInput[5], "PhiStar","ctwi", 5, 0, 2*TMath::Pi(),"1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_PhiStar.pdf","elec");
  //Ratio_EFT_SM_7pts(tInput[0], tInput[7], tInput[8], tInput[6], tInput[5], tInput[12], tInput[11], "PhiStar","ctwi", 5, 0, 2*TMath::Pi(),"1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_PhiStar.pdf","elec");
  
  //Weighted files
  // Ratio_EFT_SM(tInput[16], tInput[16], tInput[16], tInput[16], tInput[16], "PhiStar", "cbwi", 20, 0, 2*TMath::Pi(), "weight_SM", "weight_cbwi_p1","weight_cbwi_p2", "weight_cbwi_m1", "weight_cbwi_m2",  "C_{bW}^{I}", "EFT/SM", "results/EFT_vs_SM/Rwgt_cbwi2p5_XSect/PhiStar/Rwgt_cbwi2p5_cbwi_PhiStar_20bins.pdf", "elmu");
  // Ratio_EFT_SM(tInput[16], tInput[16], tInput[16], tInput[16], tInput[16], "cosThetaStar", "cbwi", 20, -1, 1, "weight_SM", "weight_cbwi_p1","weight_cbwi_p2", "weight_cbwi_m1", "weight_cbwi_m2", "C_{bW}^{I}", "EFT/SM", "results/EFT_vs_SM/Rwgt_cbwi2p5_XSect/cosThetaStar/Rwgt_cbwi2p5_cbwi_cosThetaStar_20bins.pdf", "elmu");
  // Ratio_EFT_SM(tInput[16], tInput[16], tInput[16], tInput[16], tInput[16], "PhiStar", "ctwi", 20, 0, 2*TMath::Pi(), "weight_SM", "weight_ctwi_p1","weight_ctwi_p2", "weight_ctwi_m1", "weight_ctwi_m2", "C_{tW}^{I}", "EFT/SM", "results/EFT_vs_SM/Rwgt_cbwi2p5_XSect/PhiStar/Rwgt_cbwi2p5_ctwi_PhiStar_20bins.pdf", "elmu");
  // Ratio_EFT_SM(tInput[16], tInput[16], tInput[16], tInput[16], tInput[16], "cosThetaStar", "ctwi", 20, -1, 1, "weight_SM", "weight_ctwi_p1","weight_ctwi_p2", "weight_ctwi_m1", "weight_ctwi_m2", "C_{tW}^{I}", "EFT/SM", "results/EFT_vs_SM/Rwgt_cbwi2p5_XSect/cosThetaStar/Rwgt_cbwi2p5_ctwi_cosThetaStar_20bins.pdf", "elmu");

if(reweight)
  {
    /////////////////Change the variables to get the Plot you want/////////////////

    string EFT = "";                 //EFT variable {ctwi, cbwi, SM, other} for legend name
    int W_value = 2;                   //Wilson coefficient Value for legend

    string weight = "weight_cbwi_p2";   //Cut value {weight_SM ; weight_ctwi_m2(m1,p1,p1) ; weight_cbwi_m2(m1,p1,p2)}
    

    //Rwgt_vs(weight, "PhiStar", tInput[16], tInput[19], 20, 0, 2*TMath::Pi(), EFT, to_string(W_value), "C_{bw}^{I} = 2", "#phi^{*} [rad]", "", "legendUpRight", "", "", "results/weighted/dim6");
    //Rwgt_vs(weight, "cosThetaStar", tInput[16], tInput[19], 20, -1, 1, EFT, to_string(W_value), "C_{bw}^{I} = 2", "cos(#theta^{*})", "", "legendUpRight", "", "", "results/weighted/dim6");
    //Rwgt_vs(weight, "top_mass", tInput[16], tInput[17], 40, 166, 178, EFT, to_string(W_value), "C_{tw}^{I} = -2", "M_{Top} [GeV]", "", "legendUpRight", "", "", "results/weighted/dim6");
    //Rwgt_vs(weight, "lepton_pt", tInput[16], tInput[17], 40, 0, 60, EFT, to_string(W_value), "C_{tW}^{I} = -2", "Pt_{lepton} [GeV]", "", "legendUpRight", "", "", "results/weighted/dim6");
    //Rwgt_vs(weight, "top_pt", tInput[16], tInput[17], 0, 20, 300, EFT, to_string(W_value), "C_{tW}^{I} = -2", "Pt_{Top} [GeV]", "", "legendUpRight", "", "", "results/weighted/dim6_zoom");
    //Rwgt_vs(weight, "W_pt", tInput[16], tInput[17], 40, 0, 100, EFT, to_string(W_value), "C_{tW}^{I} = -2", "Pt_{W} [GeV]", "", "legendUpRight", "", "", "results/weighted/dim6");


    //Compare_3Histos(tInput[14], tInput[14], tInput[14], "top_mass", 20, 165, 180, "", "M_{Top} [Gev]", "", "legendUpRight", "dim6top Base C_{tW}^{I} = 2.5 [TeV^{-2}]", "C_{tW}^{I} = -2 TeV^{-2}", "C_{tW}^{I} = 2 TeV^{-2}", "SM", "", "results/weighted/Rwgt_ctwi2p5_");
    //Compare_3Histos(tInput[15], tInput[15], tInput[15], "top_mass", 20, 165, 180, "", "M_{Top} [Gev]", "", "legendUpRight", "dim6top Base C_{bW}^{I} = 2.5 [TeV^{-2}]", "C_{tW}^{I} = -2 TeV^{-2}", "C_{tW}^{I} = 2 TeV^{-2}", "SM", "", "results/weighted/Rwgt_cbwi2p5_");
  }

  //////////Compare Weights//////////
  //base ctwi = 2.5
  //Compare_9Histos(100,0,0.02,"1","1","1","1","1","1","1","1","1", tInput[14],tInput[14],tInput[14],tInput[14],tInput[14],tInput[14],tInput[14],tInput[14],tInput[14],"weight_SM","weight_cbwi_m2","weight_cbwi_m1","weight_cbwi_p1","weight_cbwi_p2","weight_ctwi_m2","weight_ctwi_m1","weight_ctwi_p1","weight_ctwi_p2","C_{tW}^{I} = 2.5","","legendUpRight","dim6top EFT value [TeV^{-2}]","results/weighted/weights_EFT_ctwi2p5");
  //base ctwi = 3
  //Compare_9Histos(100,0,0.04,"1","1","1","1","1","1","1","1","1", tInput[13],tInput[13],tInput[13],tInput[13],tInput[13],tInput[13],tInput[13],tInput[13],tInput[13],"weight_SM","weight_cbwi_m2","weight_cbwi_m1","weight_cbwi_p1","weight_cbwi_p2","weight_ctwi_m2","weight_ctwi_m1","weight_ctwi_p1","weight_ctwi_p2","C_{tW}^{I} = 3","","legendUpRight","dim6top EFT value [TeV^{-2}]","results/weighted/weights_EFT_ctwi_p3");
  //base cbwi = 2.5
  //Compare_9Histos(100,0,0.02,"1","1","1","1","1","1","1","1","1", tInput[16],tInput[16],tInput[16],tInput[16],tInput[16],tInput[16],tInput[16],tInput[16],tInput[16],"weight_SM","weight_cbwi_m2","weight_cbwi_m1","weight_cbwi_p1","weight_cbwi_p2","weight_ctwi_m2","weight_ctwi_m1","weight_ctwi_p1","weight_ctwi_p2","C_{bW}^{I} = 2.5","","legendUpRight","dim6top EFT value [TeV^{-2}]","results/weighted/weights_EFT_cbwi2p5");
  //base cbwi = 3
  //Compare_9Histos(100,0,0.02,"1","1","1","1","1","1","1","1","1", tInput[15],tInput[15],tInput[15],tInput[15],tInput[15],tInput[15],tInput[15],tInput[15],tInput[15],"weight_SM","weight_cbwi_m2","weight_cbwi_m1","weight_cbwi_p1","weight_cbwi_p2","weight_ctwi_m2","weight_ctwi_m1","weight_ctwi_p1","weight_ctwi_p2","C_{bW}^{I} = 3","","legendUpRight","dim6top EFT value [TeV^{-2}]","results/weighted/weights_EFT_cbwi_p3");

  /////Made to compare Dim6=1 vs Dim6=2/////C_{tW}^{I} = -2
  //Compare_2Histos(tInput[1], tInput[20], "cosTheta", 20, -1, 1, "1", "cos(#theta)", "normalized", "legendUpRight", "C_{bW}^{I} = -2", "Dim6 = 2", "Dim6 = 1", "results/dim6_cbwi_m2");
  //Compare_2Histos(tInput[1], tInput[20], "cosThetaStar", 20, -1, 1, "1", "cos(#theta^{*})", "normalized", "legendUpRight", "C_{bW}^{I} = -2", "Dim6 = 2", "Dim6 = 1", "results/dim6_cbwi_m2");
  //Compare_2Histos(tInput[1], tInput[20], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi^{*}", "normalized", "legendUpRight", "C_{bW}^{I} = -2", "Dim6 = 2", "Dim6 = 1", "results/dim6_cbwi_m2");
  //Compare_2Histos(tInput[1], tInput[20], "top_pt", 20, 0, 300, "1", "Pt_{top} [GeV]", "normalized", "legendUpRight", "C_{bW}^{I} = -2", "Dim6 = 2", "Dim6 = 1", "results/dim6_cbwi_m2");
  //Compare_2Histos(tInput[1], tInput[20], "top_mass", 20, 166, 178, "1", "M_{top} [GeV]", "normalized", "legendUpRight", "C_{bW}^{I} = -2", "Dim6 = 2", "Dim6 = 1", "results/dim6_cbwi_m2");
  //Compare_2Histos(tInput[1], tInput[20], "lepton_pt", 20, 0, 60, "1", "Pt_{lepton} [GeV]", "normalized", "legendUpRight", "C_{bW}^{I} = -2", "Dim6 = 2", "Dim6 = 1", "results/dim6_cbwi_m2"); 
  //Compare_2Histos(tInput[1], tInput[20], "W_pt", 20, 0, 100, "1", "Pt_{W} [GeV]", "normalized", "legendUpRight", "C_{bW}^{I} = -2", "Dim6 = -2", "Dim6 = 1", "results/dim6_cbwi_m2");






/*
  //Plot_SM
  //Compare_2Histos(tInput[0], tInput[17], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_CosTheta.pdf"  );
  Compare_2Histos(tInput[0], tInput[17], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_CosThetaStar.pdf");
  //Compare_2Histos(tInput[0], tInput[17], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_PhiStar.pdf");
  //Compare_2Histos(tInput[0], tInput[17], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_top_pt.pdf");
  Compare_2Histos(tInput[0], tInput[17], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_W_pt.pdf");
  Compare_2Histos(tInput[0], tInput[17], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[17], suffix[17], "results/dim6top_compareSM_lepton_pt.pdf");
	Compare_2Histos(tInput[0], tInput[17], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_lepton_E_Wframe.pdf");
  Compare_2Histos(tInput[0], tInput[17], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_top_mass.pdf");
  Compare_2Histos(tInput[0], tInput[17], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_W_mass.pdf");

  //MadGraph

  //Plots ctW
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "cosTheta", 20, -1, 1, "1", "cos#theta", "a.u.", "legendUpLeft", "Operateur EFT", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_CosTheta.pdf");/*
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_CosThetaStar.pdf");
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_PhiStar.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_top_pt.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "Number of events", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_W_pt.pdf");
  /*Compare_3Histos(tInput[0], tInput[10], tInput[12], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_lepton_pt.pdf");
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "a.u.", "legendUpRight", "Operateur EFT", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[10], tInput[12], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_W_mass.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "W_transverse_mass",20, 50, 140, "1", "M_{T,W} (GeV)", "number of events", "legendUpRight", "Operateur de dimension 6", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_W_transverse_mass.pdf");

	//Plots ctWI*/
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_CosTheta.pdf");
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_CosThetaStar.pdf");
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_PhiStar.pdf");
  /*(tInput[0], tInput[14], tInput[16], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_top_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[14], tInput[16], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_W_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[14], tInput[16], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0],suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_lepton_E_Wframe.pdf");
  *///Compare_3Histos(tInput[0], tInput[14], tInput[16], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctwI_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[14], tInput[16], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_W_mass.pdf");

	//Plots cbWI
	//Compare_3Histos(tInput[0], tInput[2], tInput[4], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_CosTheta.pdf");
	Compare_3Histos(tInput[0], tInput[2], tInput[4], "cosThetaStar", 20, -1, 1, "1", "cos(#theta*)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_CosThetaStar.pdf");
	Compare_3Histos(tInput[0], tInput[2], tInput[4], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_PhiStar.pdf");
  Compare_3Histos(tInput[0], tInput[2], tInput[4], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_top_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[2], tInput[4], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_W_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[2], tInput[4], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[0], tInput[2], tInput[4], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_lepton_E_Wframe.pdf");
  *///Compare_3Histos(tInput[0], tInput[2], tInput[4], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbwI_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[2], tInput[4], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbwI_W_mass.pdf");

  //Plots cptbI
	//Compare_3Histos(tInput[0], tInput[6], tInput[8], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_CosTheta.pdf");
	*///Compare_3Histos(tInput[0], tInput[6], tInput[8], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_CosThetaStar.pdf");
	//Compare_3Histos(tInput[0], tInput[6], tInput[8], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_PhiStar.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_top_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_W_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_lepton_pt.pdf");
	//Compare_3Histos(tInput[0], tInput[6], tInput[8], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_lepton_E_Wframe.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[6], tInput[8], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_W_mass.pdf");

  //MadSpin

  //Plots ctW
	Compare_3Histos(tInput[17], tInput[27], tInput[29], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_CosTheta.pdf");
	Compare_3Histos(tInput[17], tInput[27], tInput[29], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[27], tInput[29], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_PhiStar.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27],suffix[29], "results/madspin_dim6top_ctW_top_pt.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[27], tInput[29], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[27], tInput[29], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_W_mass.pdf");


  //Plots ctWI
	Compare_3Histos(tInput[17], tInput[31], tInput[33], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[31], suffix[33], "results/madspin_dim6top_ctWI_CosTheta.pdf");
	Compare_3Histos(tInput[17], tInput[31], tInput[33], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0],  suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[31], tInput[33], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_PhiStar.pdf");
  Compare_3Histos(tInput[17], tInput[31], tInput[33], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0],  suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_top_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[31], tInput[33], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[31], tInput[33], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[31], tInput[33], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0],  suffix[31], suffix[33], "results/madspin_dim6top_ctWI_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[31], tInput[33], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31], suffix[33], "results/madspin_dim6top_ctWI_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[31], tInput[33], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31], suffix[33], "results/madspin_dim6top_ctWI_W_mass.pdf");


	//Plots cbWI
	//Compare_3Histos(tInput[17], tInput[19], tInput[21], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_CosTheta.pdf");
	Compare_3Histos(tInput[17], tInput[19], tInput[21], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19],suffix[21], "results/madspin_dim6top_cbWI_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[19], tInput[21], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_PhiStar.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_top_pt.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[19], tInput[21], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[19], tInput[21], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbwI_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbwI_W_mass.pdf");

  //Plots cptbI
	//Compare_3Histos(tInput[17], tInput[23], tInput[25], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_CosTheta.pdf");
	//Compare_3Histos(tInput[17], tInput[23], tInput[25], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top",  suffix[0],  suffix[23], suffix[25], "results/madspin_dim6top_cptbI_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[23], tInput[25], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_PhiStar.pdf");
  //Compare_3Histos(tInput[17], tInput[23], tInput[25], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_top_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[23], tInput[25], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[23],  suffix[25], "results/madspin_dim6top_cptbI_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[23], tInput[25], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[23], tInput[25], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[23], tInput[25], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[23], tInput[25], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_W_mass.pdf");

  //Madgraph + MadSpin

  //cbwi
  //Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_W_mass.pdf");

  //cptbI
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_W_mass.pdf");

  //ctw
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "Operateur de dimension 6", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_W_mass.pdf");


  //ctwI
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_W_mass.pdf");
*/

  //-----------------Ratio EFT/SM-------------//

  //ctWI

  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "cosTheta","ctwi", 5, -1, 1,"1", "cos#theta", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_cosTheta.pdf");
  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "PhiStar","ctwi",20, 0, 6.2831 ,"nature_lepton == 1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_PhiStar.pdf","elec");
  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "PhiStar","ctwi",20, 0, 6.2831 ,"nature_lepton == 2", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_PhiStar.pdf","muon");
  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "cosThetaStar","ctwi", 3, -1, 1 ,"1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_cosThetaStar.pdf");

/*  Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "cosTheta","cbwi", 5, -1, 1,"1", "cos#theta", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosTheta.pdf");
  Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "PhiStar","cbwi", 5, 0, 2*TMath::Pi(),"1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_PhiStar.pdf");*/
  //Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "cosThetaStar","cbwi", 20, -1, 1 ,"nature_lepton==1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosThetaStar.pdf", "elec");
  //Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "cosThetaStar","cbwi", 20, -1, 1 ,"nature_lepton==2", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosThetaStar.pdf", "muon");
 
 
  /*Compare_3Histos(tInput[0], tInput[13], tInput[14], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[13], suffix[14], "results/madgraph_dim6top_ctWI_CosTheta_test.pdf");*/

  return 0;
}