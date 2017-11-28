#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TVector2.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "/root/atlasstyle-00-03-05/AtlasStyle.h"
#include "/root/atlasstyle-00-03-05/AtlasStyle.C"
#include "/root/atlasstyle-00-03-05/AtlasUtils.h"
#include "/root/atlasstyle-00-03-05/AtlasUtils.C"
#include "/root/atlasstyle-00-03-05/AtlasLabels.h"
#include "/root/atlasstyle-00-03-05/AtlasLabels.C"
#include "TDatabasePDG.h"

void ReadNtuple()

{
  //  SetAtlasStyle();
  //TFile* f = TFile::Open("AnalysisResults.root","READ");
  //auto f = TFile::Open("GammaJet.root","READ");
  auto f = TFile::Open("16c2.root","READ");



  if(!f){
      printf("Error: cannot open ntuple.root");
      return;
  }
  f->Print();
 
  TTree* tree = (TTree*)f->Get("AliAnalysisTaskNTGJ/_tree_event");
  if(!tree){
    printf("Error: cannot find tree");
  }
  // tree->Print();
  TTreeReader myReader(tree);
  TTreeReaderArray<float> mc_truth_pt(myReader , "mc_truth_pt");
  TTreeReaderArray<float> track_pt(myReader , "track_pt");
  TTreeReaderArray<float> track_eta(myReader , "track_eta");
  TTreeReaderArray<float> track_phi(myReader , "track_phi");
  //TTreeReaderArray<char> track_charge(myReader , "track_charge");
  TTreeReaderArray<unsigned short> track_mc_truth_index(myReader, "track_mc_truth_index");
  TTreeReaderArray<unsigned char> track_quality(myReader,"track_quality");
  TTreeReaderArray<short> mc_truth_pdg_code(myReader,"mc_truth_pdg_code"); 
  TTreeReaderArray<float> mc_truth_eta(myReader, "mc_truth_eta");
  TTreeReaderArray<float> mc_truth_phi(myReader, "mc_truth_phi");
  TTreeReaderArray<char>  mc_truth_charge(myReader, "mc_truth_charge");

  TTreeReaderArray<float> track_its_ncluster(myReader,"track_its_ncluster");
  
  const Double_t bins[10] = {  0.1,          0.16681005,   0.27825594,   0.46415888 ,  0.77426368,
			      1.29154967,   2.15443469,   3.59381366,   5.9948425,   10.        };

  const int nbinseta = 30;
  Double_t etabins[nbinseta+1] = {};
  double etamin = -0.9;
  double etamax = 0.9;
  double etastep = (etamax-etamin)/nbinseta;
  for(int i=0; i<nbinseta+1; i++){
    etabins[i] = etamin + i*etastep;
  }

  const int nbinsphi = 100;
  Double_t phibins[nbinseta+1] = {};
  double phimin = -1.0*TMath::Pi();
  double phimax = 1.0*TMath::Pi();
  double phistep = (phimax-phimin)/nbinsphi;
  for(int i=0; i<nbinsphi+1; i++){
    phibins[i] = phimin + i*phistep;
  }

  auto hCorrelation   = new TH2F("hCorrelation", "", 80, 0, 10.0, 80, 0, 10.0);
  auto hRes_Pt   = new TH2F("hRes_Pt", "", 200, 0, 10.0, 80, -50, 50);
  auto hDen  = new TH1F("hDen", "", 9, bins);
  auto hNum  = new TH1F("hNum","", 9, bins);
  auto hFake = new TH1F("hFake", "", 9, bins);
  auto hReco = new TH1F("hReco","", 9,bins);
  auto hNum2Deta = new TH2F("hNum2Deta","", nbinseta, etabins, 9, bins);
  auto hDen2Deta = new TH2F("hDen2Deta","", nbinseta, etabins, 9, bins);
  auto hNum2Dphi = new TH2F("hNum2Dphi","", nbinsphi, phibins, 9, bins);
  auto hDen2Dphi = new TH2F("hDen2Dphi","", nbinsphi, phibins, 9, bins);
  auto hIsoCorrelation = new TH2F("hIsoCorrelation", "", 80, 0, 10.0, 80, 0, 10.0); 
  auto hIsoReco        = new TH1F("hIsoReco", "", 200, 0, 100);
  auto hIsoTrue        = new TH1F("hIsoTrue", "", 200, 0, 100);
  auto hCorrelation_dEtapT = new TH2F("hCorrelation_dEtapT", "" , 100, -0.1, 0.1, 100, 0.1, 10.0);
  auto hCorrelation_dPhipT = new TH2F("hCorrelation_dPhipT", "" , 100, -0.1, 0.1, 100, 0.1, 10.0);
  
  int nevent = 0; 
  const int TrackBit = 16; //ITSONLY==16; ITS--TPC with full-jet cuts

  std::cout << "nentries : " << tree->GetEntries() << std::endl;
  //myReader.Restart();
  while (myReader.Next()) {
    std::cout << "nentries : " << tree->GetEntries() << std::endl;
    nevent += 1;
    if(nevent%1==0) std::cout << nevent << std::endl;
    // if(nevent>3000) break;
    for (int n = 0;  n<track_pt.GetSize(); n++){
      auto index = track_mc_truth_index[n]; 
      if(index>65534) continue; //removing particles not associated with MC particle (i.e, secondaries or fakes)
      if((track_quality[n]&TrackBit)!=0){
          double resolution = 100*(track_pt[n]-mc_truth_pt[index])/(mc_truth_pt[index]);    
          hCorrelation->Fill(mc_truth_pt[index], track_pt[n]);
          hRes_Pt->Fill(mc_truth_pt[index], resolution);
          hNum->Fill(mc_truth_pt[index]);
          hNum2Deta->Fill(mc_truth_eta[index], mc_truth_pt[index]);
          hNum2Dphi->Fill(mc_truth_phi[index], mc_truth_pt[index]);
          double etares =  track_eta[n]-mc_truth_eta[index];
          hCorrelation_dEtapT->Fill(track_eta[n]-mc_truth_eta[index], mc_truth_pt[index]);
          hCorrelation_dPhipT->Fill( track_phi[n]-mc_truth_phi[index], mc_truth_pt[index]);
          //printf("eta %2.2f; mceta = %2.2f; res=%2.2f\n", track_eta[n], mc_truth_eta[index], etares);
      }
    }

    //Loop over MC particles (all are primaries), pick charged ones with |eta|<0.8
    for (int n = 0;  n< mc_truth_pt.GetSize(); n++){
      int pdgcode = mc_truth_pdg_code[n];
      // if(TDatabasePDG::Instance()->GetParticle(pdgcode)->Charge()==0) continue; // skip neutral particles
      if(int(mc_truth_charge[n])==0) continue;
      if(TMath::Abs(mc_truth_eta[n])>0.9) continue; //skip particles with |eta|<0.8
      hDen->Fill(mc_truth_pt[n]);
      hDen2Deta->Fill(mc_truth_eta[n], mc_truth_pt[n]);
      hDen2Dphi->Fill(mc_truth_phi[n], mc_truth_pt[n]);
    }

    //Loop to ver fraction of reconstructed ITS-only tracks that are not associated to any MC primary
    for (int n=0; n<track_pt.GetSize(); n++){ //loop over reconstructed tracks
      if( (track_quality[n]&TrackBit)!=0) {
      hReco->Fill(track_pt[n]);
      auto index = track_mc_truth_index[n]; 
      if(index==65535) hFake->Fill(track_pt[n]);     
     }
    }
   } //loop over events





   auto gausfit = new TF1("gaus","gaus", -20,20);
   gausfit->SetLineColor(kRed);
   auto g_mean = new TGraphErrors();
   auto g_sigma = new TGraphErrors();
 
   // auto c1 = new TCanvas();
   
   //Study of the ITS-only track pT resolution
   int nbins = 9;

   for(int i=0; i<nbins; i++){
     //double binwidth = 10.0/nbins;
     //double minpt = i*binwidth ;
     //double maxpt = (i+1)*binwidth;
     double minpt = bins[i];
     double maxpt = bins[i+1];
     double binwidth = maxpt-minpt;
     int minbin =  hRes_Pt->GetXaxis()->FindBin(minpt);
     int maxbin =  hRes_Pt->GetXaxis()->FindBin(maxpt);
     auto h1 = hRes_Pt->ProjectionY("h", minbin, maxbin);
     
     //h1->Draw();
     // h1->GetYaxis()->SetTitle("counts");
     h1->Fit(gausfit,"R");
     //gausfit->Draw("same");
     //myText(0.18, 0.8, kBlack, Form("%2.1f < p_{T}^{truth} < %2.1f GeV", minpt, maxpt));
     //myText(0.18, 0.74, kRed, Form("#mu = %2.1f [%]", gausfit->GetParameter(1)));
     //myText(0.18, 0.68, kRed, Form("#sigma = %2.1f [%]", gausfit->GetParameter(2))); 
     g_sigma->SetPoint(g_sigma->GetN(), (maxpt+minpt)/2.0, gausfit->GetParameter(2));     
     g_sigma->SetPointError(g_sigma->GetN()-1, binwidth/2.0, gausfit->GetParError(2));
     g_mean->SetPoint(g_mean->GetN(), (maxpt+minpt)/2.0, gausfit->GetParameter(1));
     g_mean->SetPointError(g_mean->GetN()-1, binwidth/2.0, gausfit->GetParError(1));
     //c1->SaveAs(Form("projecting%i.pdf", i));
   }
   
   auto p1 = new TF1("p1","[0] + [1]*x", 0.50, 10.0);
   p1->SetLineColor(kRed);   
   g_sigma->SetTitle("Relative resolution vs p_{T} ; p_{T}^{true} [GeV]; #sigma(p_{T})/p_{T} [%]"); 
   g_mean->SetTitle("; p_{T}^{true} [GeV]; Relative bias [%]");
   hCorrelation->SetTitle("; True p_{T} [GeV]; Reconstructed p_{T} [GeV]");
   hCorrelation_dEtapT->SetTitle("; #eta^{truth}-#eta^{reco}; p_{T}^{true} [GeV]");
   hCorrelation_dPhipT->SetTitle("; #phi^{truth}-#phi^{reco}; p_{T}^{true} [GeV]");
   
   TEfficiency* pEff = 0;
   TEfficiency* pFake = 0;
   auto fout = new TFile(Form("fout_%i.root",TrackBit), "RECREATE");
   //
   /*   if(TEfficiency::CheckConsistency(*hFake,*hDen))
   {
        std::cout << "calculating eff " << std::endl;
        pFake = new TEfficiency(*hFake,*hDen);
        pFake->Write("FakeRate");
    }


      if(TEfficiency::CheckConsistency(*hNum,*hDen))
   {
     std::cout << "calculating eff " << std::endl;
     pEff = new TEfficiency(*hNum,*hDen);
     pEff->Write("Efficiency");
   }
   */

   hFake->Write("hFake");
   hDen->Write("hDen");
   hReco->Write("hReco");
   g_sigma->Write("g_sigma");
   g_mean->Write("g_mean");

   //Writing efficiency
   hNum2Deta->Write("hNum2Deta");
   hDen2Deta->Write("hDen2Deta");
   hNum2Deta->Divide(hDen2Deta);
   hNum2Deta->SetMaximum(1.0);
   hNum2Deta->Write("hEff2Deta");

   hNum2Dphi->Write("hNum2Dphi");
   hDen2Dphi->Write("hDen2Dphi");
   hNum2Dphi->Divide(hDen2Dphi);
   hNum2Dphi->SetMaximum(1.0);
   hNum2Dphi->Write("hEff2Dphi");


   
     
   hCorrelation_dEtapT->Write("dEta_pt");
   hCorrelation_dPhipT->Write("dPhi_pt");
   hCorrelation->Write();
   fout->Close();

   auto c = new TCanvas();
   hCorrelation_dEtapT->Draw("colz");
   gPad->SetLogz();
   c->SaveAs(Form("Correlation_deta_pt%i.pdf",TrackBit));
   c->Clear();
   hCorrelation_dPhipT->Draw("colz");
   gPad->SetLogz();
   c->SaveAs(Form("Correlation_dphi_pt%i.pdf",TrackBit));


   
}
//  TTreeReaderArray<int> ncluster(myReader, "ncluster");
//Cluster variables
//  TTreeReaderArray<float>          cluster_pt(myReader, "cluster_pt");
//TTreeReaderArray<unsigned int>   cluster_nmc_truth(myReader, "cluster_nmc_truth");
//TTreeReaderArray<std::vector<unsigned short>> cluster_mc_truth_index(myReader, "cluster_mc_truth_index");
//TTreeReaderArray<float>          cluster_iso_tpc_04(myReader, "cluster_iso_tpc_04");
//TTreeReaderArray<float>          cluster_iso_04_truth(myReader, "cluster_iso_04_truth");
