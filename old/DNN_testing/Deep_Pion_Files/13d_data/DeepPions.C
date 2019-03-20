#define DeepPions_cxx
#include "DeepPions.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void DeepPions::Loop()
{
//   In a ROOT session, you can do:
//      root> .L DeepPions.C
//      root> DeepPions t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

void DeepPions::Efficiency(){
    
    TCanvas canvas("canv", "");//, 960 + 4, 720 + 28);
    //canvas.SetLogy();
    
    TH1D histogram0("histogram0", "All Clusters", 12, 8.0, 20);
    TH1D histogram1("histogram1", "DeepPion Selection", 12, 8.0, 20);
    TH1D histogram2("histogram2", "lambda Selection", 12, 8.0, 20);
    
    const int nptbins = 6;
    const float ptbins[nptbins+1] = {8.0, 10.0, 12.0, 14.0, 16.0,18.0,20.0};
    const int nclass = 8;
    
    TH1F *hist[nptbins][nclass];
    THStack* hist_stack[nptbins];
    THStack* hist_stack_Lambda[nptbins];
    
    TH1F *hist_Lambda[nptbins][nclass];
    TH2F *hist_Lambda_NN[nptbins][nclass];
    
    char *histname = new char[30];
    
    const Color_t colors[8] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan, kOrange, kGray};
    
    
    for (int ipt = 0; ipt<nptbins; ipt++){
        sprintf(histname, "hstack_NN_ptbin%d",ipt);
        hist_stack[ipt] = new THStack(histname,"");
        hist_stack[ipt]->SetTitle("; NN output (2 photons); entries");
        
        sprintf(histname, "hstack_Lambda_ptbin%d",ipt);
        hist_stack_Lambda[ipt] = new THStack(histname,"");
        hist_stack_Lambda[ipt]->SetTitle("; lambda_{0}^{2} ; entries");
        
        for(int iclass=0; iclass<nclass; iclass++){
            sprintf(histname, "h_NN_ptbin%d_class%d",ipt, iclass);
            hist[ipt][iclass]=new TH1F(histname,"",10, 0.0 , 1.0);
            hist[ipt][iclass]->SetLineColor(colors[iclass]);
            hist[ipt][iclass]->SetFillColor(colors[iclass]);
            hist[ipt][iclass]->SetLineWidth(2);
            
            sprintf(histname, "h_Lambda_ptbin%d_class%d",ipt,iclass);
            hist_Lambda[ipt][iclass]=new TH1F(histname,"",20, 0.0 , 2.0);
            hist_Lambda[ipt][iclass]->SetLineColor(colors[iclass]);
            hist_Lambda[ipt][iclass]->SetFillColor(colors[iclass]);
            hist_Lambda[ipt][iclass]->SetLineWidth(2);
            
            
            if(iclass>0){
                hist_stack[ipt]->Add(hist[ipt][iclass]);
                hist_stack_Lambda[ipt]->Add(hist_Lambda[ipt][iclass]);
            }
            sprintf(histname, "h_LambdaNN_ptbin%d_class%d",ipt,iclass);
            hist_Lambda_NN[ipt][iclass]=new TH2F(histname,"",20, 0.0 , 2.0,10,0,1.0);
            hist_Lambda_NN[ipt][iclass]->SetLineColor(iclass+1);
            
        }
    }
    
    histogram1.SetLineColor(kRed);
    
    
    //-------------------------------------Main Loop---------------------------------------------------
    
    
    for(Long64_t ievent = 0; ievent < fChain->GetEntriesFast() ; ievent++){
        //for(Long64_t ievent = 0; ievent < 100000 ; ievent++){
        fChain->GetEntry(ievent);
        for (ULong64_t n = 0; n < ncluster; n++) {
            //event selection
            if(TMath::Abs(primary_vertex[2])>10) continue;
            //cluster selection
            if(cluster_pt[n]<8) continue;
            if(cluster_ncell[n]<3) continue;
            
            //std::cout << " cluster pt " << cluster_pt[n] << " NN output " << cluster_s_nphoton[n][2] << std::endl;
            int nphotons_pi0 = 0;
            int nphotons_eta = 0;
            int nelectrons_conversion = 0;
            int nchargedhadrons = 0; //aren't from pi, eta, or el, and have charge
            unsigned short index_temp = 65535;
            unsigned short index_temp_e = 65535;
            
            for(int counter = 0 ; counter<32; counter++)
            {
                unsigned short index = cluster_mc_truth_index[n][counter];
                if(index!=65535){
                    
                    //std::cout << "true pt " << mc_truth_pt[index] << " " << mc_truth_pdg_code[index] << " " << mc_truth_first_parent_pdg_code[index] << " " << mc_truth_first_parent_pt[index] << std::endl;
                    if(mc_truth_pdg_code[index]==22) //photon
                    {
                        if(mc_truth_first_parent_pdg_code[index]==111 && index!=index_temp){ //pi0
                            nphotons_pi0 += 1;
                            index_temp = index;
                            continue; }
                        
                        else if(mc_truth_first_parent_pdg_code[index]==221){ //eta
                            nphotons_eta +=1;
                            continue;}
                    }
                    
                    else if(mc_truth_pdg_code[index]==11)
                    {
                        if(mc_truth_first_parent_pdg_code[index]==22 && index!=index_temp_e){
                            nelectrons_conversion +=1; //does this include positrons?
                            index_temp_e = index;
                            continue;}
                    }
                    
                    else if(mc_truth_charge[index] != 0){//the particle is neither an electron nor photon
                        nchargedhadrons +=1;
                    }
                }
            }//end loop on indices
            
            // std::cout << " truth photons from pi0 " << nphotons_pi0 << "; from eta " << nphotons_eta << " electrons from conv" << nelectrons_conversion << std::endl;
            //std::cout << std::endl;
            
            for(int ipt = 0; ipt<nptbins ; ipt++){
                if(cluster_pt[n]>ptbins[ipt] && cluster_pt[n] < ptbins[ipt+1]){
                    hist[ipt][0]->Fill(cluster_s_nphoton[n][2]);
                    hist_Lambda[ipt][0]->Fill(cluster_lambda_square[n][0]);
                    hist_Lambda_NN[ipt][0]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);
                    
                    if(nphotons_pi0==1){
                        hist[ipt][1]->Fill(cluster_s_nphoton[n][2]);
                        hist_Lambda[ipt][1]->Fill(cluster_lambda_square[n][0]);
                        hist_Lambda_NN[ipt][1]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);
                        
                        if(nchargedhadrons > 0){
                            hist[ipt][1]->Fill(cluster_s_nphoton[n][2]);
                            hist_Lambda[ipt][6]->Fill(cluster_lambda_square[n][0]);
                            hist_Lambda_NN[ipt][6]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]); }
                    }
                    
                    else if(nphotons_pi0==2){
                        hist[ipt][2]->Fill(cluster_s_nphoton[n][2]);
                        hist_Lambda[ipt][2]->Fill(cluster_lambda_square[n][0]);
                        hist_Lambda_NN[ipt][2]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);
                        
                        if(nchargedhadrons > 0){
                            hist[ipt][1]->Fill(cluster_s_nphoton[n][2]);
                            hist_Lambda[ipt][7]->Fill(cluster_lambda_square[n][0]);
                            hist_Lambda_NN[ipt][7]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]); }
                    }
                    
                    else if(nphotons_eta>0){
                        hist[ipt][3]->Fill(cluster_s_nphoton[n][2]);
                        hist_Lambda[ipt][3]->Fill(cluster_lambda_square[n][0]);
                        hist_Lambda_NN[ipt][3]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);
                    }
                    else if(nchargedhadrons > 0){
                        hist[ipt][5]->Fill(cluster_s_nphoton[n][2]);
                        hist_Lambda[ipt][5]->Fill(cluster_lambda_square[n][0]);
                        hist_Lambda_NN[ipt][5]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);
                    }
                    else {
                        hist[ipt][4]->Fill(cluster_s_nphoton[n][2]);
                        hist_Lambda[ipt][4]->Fill(cluster_lambda_square[n][0]);
                        hist_Lambda_NN[ipt][4]->Fill(cluster_lambda_square[n][0], cluster_s_nphoton[n][2]);
                    }
                }
            }//end for pt-differential results
            
            
            histogram0.Fill(cluster_pt[n]);
            if(cluster_s_nphoton[n][2]>0.8) histogram1.Fill(cluster_pt[n]);
            if(cluster_lambda_square[n][0]>0.3) histogram2.Fill(cluster_pt[n]);
        }//end loop on clusters.
        
        for (ULong64_t nmc = 0; nmc < nmc_truth; nmc++) {
            if(mc_truth_pdg_code[nmc]==111 && TMath::Abs(mc_truth_eta[nmc])<0.8){
                histogram0.Fill(mc_truth_pt[nmc]);
            }
        }
        
        
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
            //hist_stack_Lambda[3]->Draw();
            hist_stack[5]->Draw();
            
            
            canvas.Update();
            std::cout << "Event # " << ievent << " / " << fChain->GetEntries() << std::endl;
            canvas.SaveAs("sub.eps");
        }
    } //end loop over events
    
    TFile* fout = new TFile("fout_16c3b_small.root","RECREATE");
    
    histogram0.SetTitle("; cluster p_{T} [GeV}; entries");
    histogram1.SetTitle("; cluster p_{T} [GeV}; entries");
    histogram2.SetTitle("; cluster p_{T} [GeV}; entries");
    
    histogram0.Write("AllClusterSpectra");
    histogram1.Write("DeepPionSpectra");
    histogram2.Write("MergedPionLambdaSpectra");
    
    for (int ipt = 0; ipt<nptbins; ipt++){
        for(int iclass=0; iclass<nclass; iclass++){
            hist[ipt][iclass]->SetTitle("; NN output (2 photons); Entries");
            hist[ipt][iclass]->Write();
            hist_Lambda[ipt][iclass]->SetTitle("; #lambda^{2}_{0}; Entries");
            hist_Lambda[ipt][iclass]->Write();
            hist_Lambda_NN[ipt][iclass]->SetTitle("; #lambda^{2}_{0} ; NN output (2 photons)");
            hist_Lambda_NN[ipt][iclass]->Write();
        }
        hist_stack[ipt]->Write();
        hist_stack_Lambda[ipt]->Write();
    }
    
    fout->Close();
    
    //end loop over samples
    
    std::cout << " ending " << std::endl;
    //    return EXIT_SUCCESS;
}
