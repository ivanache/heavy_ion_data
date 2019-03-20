#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>




int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");

  for (int iarg = 1; iarg < argc; iarg++) {
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
    TFile *file = TFile::Open((TString)argv[iarg]);

    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();

    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

    //TTree *_tree_event = dynamic_cast<TTree *>
    //    (dynamic_cast<TDirectoryFile *>
    //     (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
    if (_tree_event == NULL) {
      std::cout << " fail " << std::endl;
      exit(EXIT_FAILURE);
    }  
    //_tree_event->Print();

    TApplication application("", &dummyc, dummyv);
    TCanvas canvas("canvas", "");//, 960 + 4, 720 + 28);
    //canvas.SetLogy();

    const float ptmax = 30;
    const int nbins = 30; 
    TH1D histogram0("histogram0", "All Clusters", nbins, 0, ptmax);
    TH1D histogram1("histogram1", "DeepPion Selection", nbins, 0, ptmax);
    TH1D histogram2("histogram2", "lambda Selection", nbins, 0, ptmax);

    TH1D TrueSpectra("TrueSpectra", "True spectra", nbins, 0, ptmax);        
    TrueSpectra.Sumw2();
    TH1D eff_NN("eff_NN", "Efficiency with NN>0.85 using reco pT",nbins,0,ptmax);
    TH1D eff_Lambda("eff_Lambda", "Efficiency with Lambda>0.27 using reco pT",nbins,0,ptmax);
        

    TH1D num_purity_NN("num_purity_NN", "Numerator of purity with NN>0.85", nbins, 0,ptmax);
    TH1D num_purity_Lambda("num_purity_Lambda", "Numerator of purity with Lambda>0.27", nbins, 0,ptmax);

    const int nptbins = 11; 
    const float ptbins[nptbins+1] = {10.0, 12.0, 14.0, 16.0,18.0,20.0,22.0,24.0,26.0, 28.0,30.0,50.0};

    TH1F *hist[nptbins][5];
    THStack* hist_stack[nptbins];
    THStack* hist_stack_Lambda[nptbins];
    THStack* hist_stack_b5x5[nptbins];

    TH1F *hist_Lambda[nptbins][5];
    TH1F *hist_b5x5[nptbins][5];

    TH2F *hist_Lambda_NN[nptbins][5];

    char *histname = new char[30];

    const Color_t colors[6] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan};
       

    for (int ipt = 0; ipt<nptbins; ipt++){
      sprintf(histname, "hstack_NN_ptbin%d",ipt);
      hist_stack[ipt] = new THStack(histname,"");
      hist_stack[ipt]->SetTitle("; NN output (2 photons); entries");

      sprintf(histname, "hstack_Lambda_ptbin%d",ipt);
      hist_stack_Lambda[ipt] = new THStack(histname,"");
      hist_stack_Lambda[ipt]->SetTitle("; lambda_{0}^{2} ; entries");

      sprintf(histname, "hstack_b5x5_ptbin%d",ipt);
      hist_stack_b5x5[ipt] = new THStack(histname,"");
      hist_stack_b5x5[ipt]->SetTitle("; b5x5 ; entries");


      for(int iclass=0; iclass<5; iclass++){
	sprintf(histname, "h_NN_ptbin%d_class%d",ipt, iclass);
	hist[ipt][iclass]=new TH1F(histname,"",40, 0.0 , 1.0);
	hist[ipt][iclass]->SetLineColor(colors[iclass]);
	hist[ipt][iclass]->SetFillColor(colors[iclass]);
	hist[ipt][iclass]->SetLineWidth(2);
	       
	sprintf(histname, "h_Lambda_ptbin%d_class%d",ipt,iclass);
	hist_Lambda[ipt][iclass]=new TH1F(histname,"",40, 0.0 , 2.0);
	hist_Lambda[ipt][iclass]->SetLineColor(colors[iclass]);
	hist_Lambda[ipt][iclass]->SetFillColor(colors[iclass]);
	hist_Lambda[ipt][iclass]->SetLineWidth(2);
                        

	sprintf(histname, "h_b5x5_ptbin%d_class%d",ipt,iclass);
	hist_b5x5[ipt][iclass]=new TH1F(histname,"",40, 0.0 , 1.0);
	hist_b5x5[ipt][iclass]->SetLineColor(colors[iclass]);
	hist_b5x5[ipt][iclass]->SetFillColor(colors[iclass]);
	hist_b5x5[ipt][iclass]->SetLineWidth(2);

	if(iclass>0){
	  hist_stack[ipt]->Add(hist[ipt][iclass]);
	  hist_stack_Lambda[ipt]->Add(hist_Lambda[ipt][iclass]);
	  hist_stack_b5x5[ipt]->Add(hist_b5x5[ipt][iclass]);
	}
	sprintf(histname, "h_LambdaNN_ptbin%d_class%d",ipt,iclass);
	hist_Lambda_NN[ipt][iclass]=new TH2F(histname,"",40, 0.0 , 2.0,40,0,1.0);
	hist_Lambda_NN[ipt][iclass]->SetLineColor(iclass+1);
              
      }
    }

    histogram1.SetLineColor(kRed);

    Double_t primary_vertex[3];
       
    UInt_t ntrack;
    Float_t track_e[NTRACK_MAX];
    Float_t track_pt[NTRACK_MAX];
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    UChar_t track_quality[NTRACK_MAX];

    UInt_t ncluster;
    Float_t cluster_e[NTRACK_MAX];
    Float_t cluster_e_cross[NTRACK_MAX];
    Float_t cluster_pt[NTRACK_MAX];
    Float_t cluster_eta[NTRACK_MAX];
    Float_t cluster_phi[NTRACK_MAX]; 

    Float_t cluster_s_nphoton[NTRACK_MAX][4];
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    Int_t cluster_ncell[NTRACK_MAX];

    Float_t cluster_lambda_square[NTRACK_MAX][2];  
    Float_t  cluster_b5x5_lin[NTRACK_MAX]; 
    //MC
    unsigned int nmc_truth;
    Float_t mc_truth_pt[NTRACK_MAX];
    Float_t mc_truth_eta[NTRACK_MAX];  
    Float_t mc_truth_phi[NTRACK_MAX];
    short mc_truth_pdg_code[NTRACK_MAX];
    short mc_truth_first_parent_pdg_code[NTRACK_MAX];
    char mc_truth_charge[NTRACK_MAX];

    Float_t mc_truth_first_parent_e[NTRACK_MAX];
    Float_t mc_truth_first_parent_pt[NTRACK_MAX];
    Float_t mc_truth_first_parent_eta[NTRACK_MAX];
    Float_t mc_truth_first_parent_phi[NTRACK_MAX];
    UChar_t mc_truth_status[NTRACK_MAX];
 
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_e", track_e);
    _tree_event->SetBranchAddress("track_pt", track_pt);
    _tree_event->SetBranchAddress("track_eta", track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    _tree_event->SetBranchAddress("track_quality", track_quality);

    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_b5x5_lin", cluster_b5x5_lin);

    _tree_event->SetBranchAddress("nmc_truth",&nmc_truth);
    _tree_event->SetBranchAddress("mc_truth_pt",mc_truth_pt);
    _tree_event->SetBranchAddress("mc_truth_eta",mc_truth_eta);
    _tree_event->SetBranchAddress("mc_truth_phi",mc_truth_phi);
    _tree_event->SetBranchAddress("mc_truth_charge",mc_truth_charge);
    _tree_event->SetBranchAddress("mc_truth_pdg_code",mc_truth_pdg_code);
    _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code", mc_truth_first_parent_pdg_code);
    _tree_event->SetBranchAddress("mc_truth_first_parent_e", mc_truth_first_parent_e);
    _tree_event->SetBranchAddress("mc_truth_first_parent_pt", mc_truth_first_parent_pt);
    _tree_event->SetBranchAddress("mc_truth_first_parent_eta", mc_truth_first_parent_eta);
    _tree_event->SetBranchAddress("mc_truth_first_parent_phi", mc_truth_first_parent_phi);
    _tree_event->SetBranchAddress("mc_truth_status",  mc_truth_status);
    
    std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;


    for(Long64_t ievent = 0; ievent < _tree_event->GetEntries() ; ievent++){     
      // for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
      _tree_event->GetEntry(ievent);
      for (ULong64_t n = 0; n < ncluster; n++) {
	//event selection
	if(TMath::Abs(primary_vertex[2])>10) continue;
	//cluster selection
	if(cluster_pt[n]<8) continue;
	//if(cluster_pt[n]>20) continue;
	if(TMath::Abs(cluster_eta[n])>0.6) continue;
	if(cluster_ncell[n]<3) continue; 
	if(cluster_e_cross[n]/cluster_e[n]<0.05) continue;
              
	//std::cout << " cluster pt " << cluster_pt[n] << " NN output " << cluster_s_nphoton[n][2] << std::endl;
	int nphotons_pi0 = 0;
	int nphotons_eta = 0;
	int nelectrons_convertion = 0;
	unsigned short index_temp = 65535;
	unsigned short index_temp_e = 65535;

	for(int counter = 0 ; counter<32; counter++)
	  {
	    unsigned short index = cluster_mc_truth_index[n][counter];
	    if(index!=65535){
	      //if (mc_truth_first_parent_pdg_code[index]==22) continue;
	      // std::cout << "reco pt" << cluster_pt[n] << " true pt " << mc_truth_pt[index] << " pdgcode:" << mc_truth_pdg_code[index] << "  parent pdg:" << mc_truth_first_parent_pdg_code[index] << " first parent:" << mc_truth_first_parent_pt[index] << std::endl;
                   


	      if(mc_truth_pdg_code[index]==22){
		if(mc_truth_first_parent_pdg_code[index]==111 && index!=index_temp){
		  nphotons_pi0 += 1;
		  index_temp = index; 
		}
		else if(mc_truth_first_parent_pdg_code[index]==221){
		  nphotons_eta +=1;
		}               
	      }
	      else if(mc_truth_pdg_code[index]==11){
		if(mc_truth_first_parent_pdg_code[index]==22 && index!=index_temp_e){
		  nelectrons_convertion +=1;
		  index_temp_e = index;
		}
	      }
               
	    }
	  }//end loop on indices

	//std::cout << " truth photons from pi0 " << nphotons_pi0 << "; from eta " << nphotons_eta << " electrons from conv" << nelectrons_convertion << std::endl;
	//std::cout << std::endl;


	//Fill numerator of purity
	if(nphotons_pi0==2 || nphotons_pi0==1){
	  if(cluster_s_nphoton[n][2]>0.85){
	    num_purity_NN.Fill(cluster_pt[n]);
	    eff_NN.Fill(cluster_pt[n]);
	  }
	  if(cluster_lambda_square[n][0]>0.27){
	    num_purity_Lambda.Fill(cluster_pt[n]);
	    eff_Lambda.Fill(cluster_pt[n]);
	  }  
	}

	for(int ipt = 0; ipt<nptbins ; ipt++){
	  if(cluster_pt[n]>ptbins[ipt] && cluster_pt[n] < ptbins[ipt+1]){
	    hist[ipt][0]->Fill(cluster_s_nphoton[n][2]);
	    hist_Lambda[ipt][0]->Fill(cluster_lambda_square[n][0]);
	    hist_b5x5[ipt][0]->Fill(cluster_b5x5_lin[n]);
	    hist_Lambda_NN[ipt][0]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);
                
	    //start filling different classes 
	    if(nphotons_pi0==1){
	      hist[ipt][2]->Fill(cluster_s_nphoton[n][2]);
	      hist_Lambda[ipt][2]->Fill(cluster_lambda_square[n][0]);
	      hist_b5x5[ipt][2]->Fill(cluster_b5x5_lin[n]);
	      hist_Lambda_NN[ipt][2]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);

	    }
	    else if(nphotons_pi0==2){
	      hist[ipt][1]->Fill(cluster_s_nphoton[n][2]);
	      hist_Lambda[ipt][1]->Fill(cluster_lambda_square[n][0]);
	      hist_b5x5[ipt][1]->Fill(cluster_b5x5_lin[n]);

	      hist_Lambda_NN[ipt][1]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);
	    }
	    else if(nphotons_eta>0){
	      hist[ipt][3]->Fill(cluster_s_nphoton[n][2]);
	      hist_Lambda[ipt][3]->Fill(cluster_lambda_square[n][0]);
	      hist_b5x5[ipt][3]->Fill(cluster_b5x5_lin[n]);

	      hist_Lambda_NN[ipt][3]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);
	    }
	    else{
	      hist[ipt][4]->Fill(cluster_s_nphoton[n][2]);
	      hist_Lambda[ipt][4]->Fill(cluster_lambda_square[n][0]);
	      hist_b5x5[ipt][4]->Fill(cluster_b5x5_lin[n]);

	      hist_Lambda_NN[ipt][4]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);
	    }
	  }
	}//end for pt-differential results 

              
	histogram0.Fill(cluster_pt[n]);
	if(cluster_s_nphoton[n][2]>0.85) histogram1.Fill(cluster_pt[n]);
	if(cluster_lambda_square[n][0]>0.27) histogram2.Fill(cluster_pt[n]);   
      }//end loop on clusters. 
         
      //loop over truth particles
      for (ULong64_t nmc = 0; nmc < nmc_truth; nmc++) {
	if(mc_truth_pdg_code[nmc]==22 &&  mc_truth_first_parent_pdg_code[nmc]==111){
	  //std::cout << " pdg " << mc_truth_pdg_code[nmc] <<  " eta " << mc_truth_eta[nmc] << " " << mc_truth_first_parent_pdg_code[nmc] << 
	  //  " " << mc_truth_pt[nmc] << " parent " << mc_truth_first_parent_pt[nmc] << std::endl; 
	  if(TMath::Abs(mc_truth_eta[nmc])>0.6) continue;
	  if(TMath::Abs(mc_truth_first_parent_eta[nmc])<0.6 && mc_truth_first_parent_phi[nmc]>1.4){  
	    TrueSpectra.Fill(mc_truth_first_parent_pt[nmc],0.5);
	  }
	}
      }//end loop over truth particles


      if (ievent % 100000 == 0) {
	//hist[0][0]->Draw("e1x0");
	for(int iclass=1; iclass<5; iclass++){
	  hist_stack[5]->RecursiveRemove(hist[5][iclass]);
	  hist_stack_Lambda[5]->RecursiveRemove(hist_Lambda[3][iclass]);
	}
	for(int iclass=1; iclass<5; iclass++){
	  hist_stack[5]->Add(hist[5][iclass]);
	  hist_stack_Lambda[5]->Add(hist_Lambda[5][iclass]);
	}
              
	//hist_stack[5]->Draw();
	//canvas.Update();
	std::cout << "Event # " << ievent << " / " << _tree_event->GetEntries() << std::endl;
	// canvas.SaveAs("sub.eps");
      }
    } //end loop over events

    //TEfficiency* eff_Binomial_NN = 0;
    //TEfficiency* eff_Binomial_Lambda = 0;
    //TGraphAsymmErrors* eff_Binomial_NN = 0;
    //TGraphAsymmErrors* eff_Binomial_Lambda = 0;
    TFile* fout = new TFile("fout.root","RECREATE");
    
    histogram0.SetTitle("; cluster p_{T} [GeV]}; entries");
    histogram1.SetTitle("; cluster p_{T} [GeV]}; entries");
    histogram2.SetTitle("; cluster p_{T} [GeV]}; entries");
    TrueSpectra.SetTitle("; true p_{T} pi0; entries");
    TrueSpectra.Write("TrueSpectra");
    eff_NN.Write("Num_eff_NN");
    eff_NN.Divide(&TrueSpectra);
    eff_NN.Write("eff_NN");  
    eff_Lambda.Divide(&TrueSpectra);
    eff_Lambda.Write("eff_Lambda");
    TGraphAsymmErrors* eff_Binomial_NN = new TGraphAsymmErrors(&num_purity_NN, &TrueSpectra);
    eff_Binomial_NN->Write("eff_NN_bin");
    eff_Binomial_NN->Print();
    TGraphAsymmErrors*eff_Binomial_Lambda = new TGraphAsymmErrors(&num_purity_Lambda, &TrueSpectra);

    //eff_Binomial_NN = new TEfficiency(num_purity_NN,TrueSpectra);
    //eff_Binomial_Lambda = new TEfficiency(num_purity_Lambda,TrueSpectra);
    //eff_Binomial_NN->Write("eff_NN_bin");
    eff_Binomial_Lambda->Write("eff_Lambda_bin");

    eff_Binomial_Lambda->Print();

    num_purity_NN.SetTitle("; cluster p_{T} [GeV] with NN>0.85; entries");
    num_purity_Lambda.SetTitle("; cluster p_{T} [GeV] with lambda>0.27; entries");

    histogram0.Write("AllClusterSpectra");
    histogram1.Write("DeepPionSpectra");
    histograrity_Lambda.Sumw2();
    histogram2.Sumw2();
    num_purity_Lambda.Divide(&histogram2);
    num_purity_Lambda.SetTitle("; cluster p_{T} [GeV] with Lambda >0.27; Purity");
    num_purity_Lambda.Write("Purity_Lambda");

    histogram1.Divide(&num_purity_NN);
    histogram1.Divide(&eff_NN);
    histogram1.Divide(&TrueSpectra);
    histogram1.Write("ClosureTest_NN");
    histogram2.Divide(&num_purity_Lambda);
    histogram2.Divide(&eff_Lambda);
    histogram2.Divide(&TrueSpectra);
    histogram2.Write("ClosureTest_Lambda"_NN[ipt][iclass]->SetTitle("; #lambda^{2}_{0} ; NN output (2 photons)");
		     hist_Lambda_NN[ipt][iclass]->Write();
		     }
      hist_stack[ipt]->Write();
    hist_stack_Lambda[ipt]->Write();
    hist_stack_b5x5[ipt]->Write();
  }

  fout->Close();     
    
}//end loop over samples

std::cout << " ending " << std::endl;
return EXIT_SUCCESS;
}
