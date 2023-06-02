#define SingleTopLHEAnalyzer_cxx
#include "../include/SingleTopLHEAnalyzer.hpp"

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

// double weight_SM_size;
// double weight_ctwi_m5_size, weight_ctwi_m2_size, weight_ctwi_m1_size, weight_ctwi_p1_size, weight_ctwi_p2_size, weight_ctwi_p5_size;
// double weight_cbwi_m5_size, weight_cbwi_m2_size, weight_cbwi_m1_size, weight_cbwi_p1_size, weight_cbwi_p2_size, weight_cbwi_p5_size;

void SingleTopLHEAnalyzer::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L SingleTopLHEAnalyzer.C
//      Root > SingleTopLHEAnalyzer t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//.vscode/
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   float weight;
   float sinTheta, cosTheta;
   float sinThetaStar, cosThetaXStar, cosThetaYStar, cosThetaZStar;
   float sinPhiStar, cosPhiStar, PhiStar;
   float lepton_E_Wframe;
   float top_pt, W_pt, lepton_pt;
   float top_mass, W_mass, W_transverse_mass;
   float nature_lepton;
   double sumOf_LHEweights = 0.0;


   TFile* fOutput = new TFile("output.root","RECREATE");
   TTree* tOutput = new TTree("Tree","Tree");

   tOutput->Branch("nature_lepton",&nature_lepton,"nature_lepton/F");
   tOutput->Branch("weight",&weight,"weight/F");
   tOutput->Branch("sumOf_LHEweights",&sumOf_LHEweights,"sumOf_LHEweights/F");
   tOutput->Branch("cosTheta",&cosTheta,"cosTheta/F");
   tOutput->Branch("sinTheta",&sinTheta,"sinTheta/F");

   tOutput->Branch("cosThetaZStar",&cosThetaZStar,"cosThetaZStar/F");
   tOutput->Branch("cosThetaXStar",&cosThetaXStar,"cosThetaXStar/F");
   tOutput->Branch("cosThetaYStar",&cosThetaYStar,"cosThetaYStar/F");

   tOutput->Branch("sinThetaStar",&sinThetaStar,"sinThetaStar/F");
   tOutput->Branch("cosPhiStar",&cosPhiStar,"cosPhiStar/F");
   tOutput->Branch("sinPhiStar",&sinPhiStar,"sinPhiStar/F");
   tOutput->Branch("PhiStar",&PhiStar,"PhiStar/F");
   tOutput->Branch("lepton_E_Wframe",&lepton_E_Wframe,"lepton_E_Wframe/F");
   tOutput->Branch("top_pt",&top_pt,"top_pt/F");
   tOutput->Branch("W_pt",&W_pt,"W_pt/F");
   tOutput->Branch("lepton_pt",&lepton_pt,"lepton_pt/F");
   tOutput->Branch("W_mass", &W_mass, "W_mass/F");
   tOutput->Branch("top_mass", &top_mass, "top_mass/F");
   tOutput->Branch("W_transverse_mass", &W_transverse_mass, "W_transverse_mass/F");

  // Add New Branches in respect to your input file here


   if (fChain == 0) return;

	//Long64_t nentries = 100;
   Long64_t nentries = fChain->GetEntriesFast();
   //cout<<"nentries= "<<nentries<<endl;

   TLorentzVector Ptop;
   TLorentzVector Pb;
   TLorentzVector Pw;
   TLorentzVector Pl;
   TLorentzVector Pnu;
   TLorentzVector Pqspec;
   double Pl_ID;

     int value = 0;

	cout << "Type 1 for Z in the W direction or 2 for Z in the spectator quark direction" << endl;		//ask the user for the reference frame
	cin >> value;


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

	  //cout << "nParticles="<<Particle_<<endl;
	  for (int i=0; i<Particle_; i++){

      //If on wants to know the sum of the weights use a cout of this variables (sum of weights is linked to the Xsection)
       if (Rwgt_>0)
      {
        // Add reweight branches here

      }

		  if (TMath::Abs(Particle_PID[i])==24 && Particle_Status[i]==2)
			   Pw.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

		  if (TMath::Abs(Particle_PID[i])==5 && Particle_Status[i]==1)
			  Pb.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

		  if (TMath::Abs(Particle_PID[i])<=5 && Particle_Status[i]==1)
			  Pqspec.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

		  if ((TMath::Abs(Particle_PID[i])==11 || TMath::Abs(Particle_PID[i])==13) && Particle_Status[i]==1)
      {
        Pl.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
        Pl_ID = abs(Particle_PID[i]);
        if (Pl_ID == 11)
        {
          nature_lepton = 1;
        }
        if (Pl_ID == 13)
        {
          nature_lepton = 2;
        }
      }


      if ((TMath::Abs(Particle_PID[i])==12 || TMath::Abs(Particle_PID[i])==14) && Particle_Status[i]==1)
  		  Pnu.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

	  }

    Ptop = Pb + Pl + Pnu;

	  weight = Event_Weight[jentry];
	  sumOf_LHEweights += weight;

	  /* SELECTIONS */
    // M2 selection
  	if ((Pl.Pt()<=32 || TMath::Abs(Pl.Eta())>=2.1) && Pl_ID==11) continue; //Electron
    if ((Pl.Pt()<=30 || TMath::Abs(Pl.Eta())>=2.4) && Pl_ID==13) continue; //Muon
	  if (Pqspec.Pt()<=40 || TMath::Abs(Pqspec.Eta())>=4.7) continue; //Jet 1st selection
    if (abs(Pqspec.Eta()) >= 2.4 && Pqspec.Pt() <= 60) continue; //Jet 2nd selection
	  if (Pb.Pt()<=40 || TMath::Abs(Pb.Eta())>=2.5) continue;
	if (Pb.Pt()<=60 || TMath::Abs(Pb.Eta())>=2.4) continue;


    // STreco selection
  	// if ((Pl.Pt()<32 || (TMath::Abs(Pl.Eta()) > 1.4442 && TMath::Abs(Pl.Eta()) < 1.5660)) && Pl_ID==11) continue; //Electron Streco Selection
    // if ((Pl.Pt()<30 || TMath::Abs(Pl.Eta())>2.4) && Pl_ID==13) continue; //Muon Streco Selection [Pt(2016 and 2018: 26Gev ; 2017: 30Gev)]
    // if (Pqspec.Pt()<40 || TMath::Abs(Pqspec.Eta())>4.7) continue; //Jet 1st selection
    // if (TMath::Abs(Pqspec.Eta()) >= 2.4 && Pqspec.Pt() < 60) continue; //Jet Streco 2nd selection for |eta|>=2.4
	  // if (Pb.Pt()<40 || TMath::Abs(Pb.Eta())>2.5) continue; //Streco Selection [Eta for 2016 >2.4 and for 2017/2018 >2.5]
    

	  /* ANGLE RECONSTRUCTION */


	  TVector3 InvTopBoost;  
	  InvTopBoost.SetXYZ(-Ptop.Px()/Ptop.E(),-Ptop.Py()/Ptop.E(),-Ptop.Pz()/Ptop.E());


if(value == 1){				//calculation of cosines in the old reference frame

          Pw.Boost(InvTopBoost);
          Pl.Boost(InvTopBoost);
          Pqspec.Boost(InvTopBoost);

	  TVector3 Zdir = Pw.Vect().Unit();
   	  TVector3 qSpecUnit = Pqspec.Vect().Unit();
   	  TVector3 Ydir = qSpecUnit.Cross(Zdir).Unit();
   	  TVector3 Xdir = Ydir.Cross(Zdir);
   	  TVector3 leptonUnitary = Pl.Vect().Unit();

	  cosThetaZStar = leptonUnitary.Dot(Zdir);
	  cosThetaYStar = leptonUnitary.Dot(Ydir);
          cosThetaXStar = leptonUnitary.Dot(Xdir);
}

else{							//calculation of cosines in the new reference frame

          TVector3 lightQ = Pqspec.Vect().Unit() * Pqspec.Z();
	  
	  Pl.Boost(InvTopBoost);
	  Pqspec.Boost(InvTopBoost);
	  lightQ = lightQ + InvTopBoost;

	  TVector3 lightQunit = lightQ.Unit();

	  TVector3 Zdir = Pqspec.Vect().Unit();
	  TVector3 Ydir = -lightQunit.Cross(Zdir).Unit();
	  TVector3 Xdir = Ydir.Cross(Zdir);
	  TVector3 leptonUnitary = Pl.Vect().Unit();

	  cosThetaZStar = leptonUnitary.Dot(Zdir);
	  cosThetaYStar = leptonUnitary.Dot(Ydir);
          cosThetaXStar = leptonUnitary.Dot(Xdir);
}
   
    top_pt = Ptop.Pt();
    W_pt = Pw.Pt();
    lepton_pt = Pl.Pt();
    top_mass = Ptop.M();
    W_mass = Pw.M();
    W_transverse_mass = Pw.Mt();
    
    tOutput->Fill();
   }

   tOutput->Write();
   fOutput->Close();
}
