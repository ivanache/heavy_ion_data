/**
   This program produces energy response plots from Monte-Carlo simulations.
*/
// Original author: Miguel Arratia
// Taken and adapted by Ivan Chernyshev, June 12, 2018
// Syntax: takes in all arguments as ROOT files to fill histograms from, and fills a single set of histograms, and each histograms has input from all arguments

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <THStack.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <TGraphAsymmErrors.h>

#define NTRACK_MAX (1U << 15)

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
    TApplication application("", &dummyc, dummyv);
    std::cout <<" Number of arguments " << argc << std::endl;
    
    const int xjbins = 20;
    const int phibins = 20;
    const int etabins = 20;
    
    TH1D hjet_Xj("hjet_Xj", "Xj distribution, Jets", xjbins, 0.0,2.0);
    TH1D hjet_Xj_truth("hjet_Xj_truth", "True Xj distribution, Jets", xjbins, 0.0,2.0);
    
    
    TH1D hjet_dPhi("hjet_dPhi", "delta phi jet-jet signal region", phibins, 0, TMath::Pi());
    TH1D hjet_dPhi_truth("h_dPhi_truth", "delta phi jet-jet truth MC", phibins, 0, TMath::Pi());
    
    TH1D hjet_pTD("hjet_pTD", "pTD distribution, Jets", 40, 0.0,1.0);
    TH1D hjet_Multiplicity("hjet_Multiplicity", "Secondary Jet Multiplicity distribution", 20, -0.5 , 19.5);
    TH1D hjet_jetwidth("hjet_jetwidth", "jet width distribution", 20, 0.0,0.01);
    
    
    TH1D hjet_dEta("hjet_dEta", "delta eta jet-jet", etabins, -1.2, 1.2);
    TH1D hjet_dEta_truth("hjet_dEta_truth", "delta eta jet-jet, truth", etabins, -1.2, 1.2);
    
    TH1D hjet_AvgEta("hjet_AvgEta", "Average eta jet-jet", 2*etabins, -1.2, 1.2);
    TH1D hjet_AvgEta_truth("hjet_AvgEta_truth", "Average eta jet-jet, truth", 2*etabins, -1.2, 1.2);
    
    TH1D hjet_leadingpT("hjet_leadingpT", "Leading jet pT distribution", 25, 5, 30);
    TH1D hjet_leadingpT_truth("hjet_leadingpT_truth", "True Leading jet pT distribution, Jets", 25, 5, 30);
    
    TH1D hjet_leadingEta("hjet_leadingEta", "Leading jet eta distribution", 2*etabins, -1.2, 1.2);
    TH1D hjet_leadingEta_truth("hjet_leadingEta_truth", "True Leading jet eta distribution", 2*etabins, -1.2, 1.2);
    
    TH1D hjet_leadingPhi("hjet_leadingPhi", "Leading jet phi distribution", 2*phibins, -TMath::Pi(), TMath::Pi());
    TH1D hjet_leadingPhi_truth("hjet_leadingPhi_truth", "True Leading jet phi distribution", 2*phibins, -TMath::Pi(), TMath::Pi());
    
    hjet_Xj.Sumw2();
    hjet_Xj_truth.Sumw2();
    hjet_dPhi.Sumw2();
    hjet_dPhi_truth.Sumw2();
    hjet_pTD.Sumw2();
    hjet_Multiplicity.Sumw2();
    hjet_dEta.Sumw2();
    hjet_dEta_truth.Sumw2();
    hjet_AvgEta.Sumw2();
    hjet_AvgEta_truth.Sumw2();
    
    hjet_Xj.SetTitle("; X_{j} ; 1/N_{#gamma} dN_{J#gamma}/dX_{j}");
    hjet_Xj_truth.SetTitle("; X_{j}^{true} ; counts");
    
    double sum_of_weights = 0;
    double weightsum_individual_jets = 0;
    int N_individual_jets = 0;
    int N_eventpassed = 0;
    
    for (int iarg = 1; iarg < argc; iarg++) {
        std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
        TFile *file = TFile::Open((TString)argv[iarg]);
        
        if (file == NULL) {
            std::cout << " fail; could not open file" << std::endl;
            exit(EXIT_FAILURE);
        }
        file->Print();
        
        // Get all the TTree variables from the file to open
        TTree *_tree_event = NULL;
        std::cout << " About to try getting the ttree" << std::endl;
        _tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
        if (_tree_event == NULL) {
            std::cout << "First try did not got trying again" << std::endl;
            _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
            if (_tree_event == NULL) {
                std::cout << " fail; could not find _tree_event " << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
        if (_tree_event == NULL) {
            std::cout << " fail; the _tree_event is NULL " << std::endl;
            exit(EXIT_FAILURE);
        }
        
        TCanvas* canvas = new TCanvas();
        
        //you define variables
        Double_t primary_vertex[3];
        Bool_t is_pileup_from_spd_5_08;
        Bool_t is_pileup_from_spd_3_08;
        
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
        Float_t cluster_iso_tpc_04[NTRACK_MAX];
        Float_t cluster_iso_its_04[NTRACK_MAX];
        Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
        Float_t cluster_frixione_its_04_02[NTRACK_MAX];
        Float_t cluster_s_nphoton[NTRACK_MAX][4];
        UChar_t cluster_nlocal_maxima[NTRACK_MAX];
        Float_t cluster_distance_to_bad_channel[NTRACK_MAX];
        
        unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
        Int_t cluster_ncell[NTRACK_MAX];
        UShort_t  cluster_cell_id_max[NTRACK_MAX];
        Float_t cluster_lambda_square[NTRACK_MAX][2];
        Float_t cell_e[17664];
        
        //Jets
        UInt_t njet_ak04its;
        Float_t jet_ak04its_pt_raw[NTRACK_MAX];
        Float_t jet_ak04its_eta_raw[NTRACK_MAX];
        Float_t jet_ak04its_phi[NTRACK_MAX];
        
        Float_t jet_ak04its_width_sigma_raw[NTRACK_MAX][2];
        Float_t jet_ak04its_ptd_raw[NTRACK_MAX];
        UShort_t jet_ak04its_multiplicity[NTRACK_MAX];
        
        Float_t jet_ak04its_pt_truth[NTRACK_MAX];
        Float_t jet_ak04its_eta_truth[NTRACK_MAX];
        Float_t jet_ak04its_phi_truth[NTRACK_MAX];
        
        
        //Truth Jets
        UInt_t njet_truth_ak04;
        Float_t jet_truth_ak04_pt[NTRACK_MAX];
        Float_t jet_truth_ak04_eta[NTRACK_MAX];
        Float_t jet_truth_ak04_phi[NTRACK_MAX];
        
        //Int_t eg_ntrial;
        
        Float_t eg_cross_section;
        Int_t   eg_ntrial;
        
        //MC
        unsigned int nmc_truth;
        Float_t mc_truth_pt[NTRACK_MAX];
        Float_t mc_truth_eta[NTRACK_MAX];
        Float_t mc_truth_phi[NTRACK_MAX];
        short mc_truth_pdg_code[NTRACK_MAX];
        short mc_truth_first_parent_pdg_code[NTRACK_MAX];
        char mc_truth_charge[NTRACK_MAX];
        UChar_t mc_truth_status[NTRACK_MAX];
        
        Float_t mc_truth_first_parent_e[NTRACK_MAX];
        Float_t mc_truth_first_parent_pt[NTRACK_MAX];
        Float_t mc_truth_first_parent_eta[NTRACK_MAX];
        Float_t mc_truth_first_parent_phi[NTRACK_MAX];
        
        
        
        // Set the branch addresses of the branches in the TTrees
        //    _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
        _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
        _tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
        _tree_event->SetBranchAddress("is_pileup_from_spd_3_08", &is_pileup_from_spd_3_08);
        
        _tree_event->SetBranchAddress("ntrack", &ntrack);
        _tree_event->SetBranchAddress("track_e", track_e);
        _tree_event->SetBranchAddress("track_pt", track_pt);
        _tree_event->SetBranchAddress("track_eta", track_eta);
        _tree_event->SetBranchAddress("track_phi", track_phi);
        _tree_event->SetBranchAddress("track_quality", track_quality);
        
        _tree_event->SetBranchAddress("ncluster", &ncluster);
        _tree_event->SetBranchAddress("cluster_e", cluster_e);
        _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
        _tree_event->SetBranchAddress("cluster_pt", cluster_pt); // here
        _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
        _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
        _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton); // here
        _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
        _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
        _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
        _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
        _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
        _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
        _tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);
        _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);
        
        _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
        _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
        _tree_event->SetBranchAddress("cell_e", cell_e);
        
        _tree_event->SetBranchAddress("nmc_truth", &nmc_truth);
        _tree_event->SetBranchAddress("mc_truth_pdg_code", mc_truth_pdg_code);
        _tree_event->SetBranchAddress("mc_truth_pt", mc_truth_pt);
        _tree_event->SetBranchAddress("mc_truth_phi", mc_truth_phi);
        _tree_event->SetBranchAddress("mc_truth_eta", mc_truth_eta);
        _tree_event->SetBranchAddress("mc_truth_status", mc_truth_status);
        _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code",mc_truth_first_parent_pdg_code);
        
        _tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
        _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
        
        
        //jets
        _tree_event->SetBranchAddress("njet_ak04its", &njet_ak04its);
        _tree_event->SetBranchAddress("jet_ak04its_pt_raw", jet_ak04its_pt_raw);
        _tree_event->SetBranchAddress("jet_ak04its_eta_raw", jet_ak04its_eta_raw);
        _tree_event->SetBranchAddress("jet_ak04its_phi", jet_ak04its_phi);
        _tree_event->SetBranchAddress("jet_ak04its_pt_truth", jet_ak04its_pt_truth);
        _tree_event->SetBranchAddress("jet_ak04its_eta_truth", jet_ak04its_eta_truth);
        _tree_event->SetBranchAddress("jet_ak04its_phi_truth", jet_ak04its_phi_truth);
        
        _tree_event->SetBranchAddress("jet_ak04its_ptd_raw", jet_ak04its_ptd_raw);
        _tree_event->SetBranchAddress("jet_ak04its_multiplicity_raw", jet_ak04its_multiplicity);
        _tree_event->SetBranchAddress("jet_ak04its_width_sigma_raw", jet_ak04its_width_sigma_raw);
        
        //truth jets
        _tree_event->SetBranchAddress("njet_truth_ak04", &njet_truth_ak04);
        _tree_event->SetBranchAddress("jet_truth_ak04_pt", jet_truth_ak04_pt);
        _tree_event->SetBranchAddress("jet_truth_ak04_phi", jet_truth_ak04_phi);
        _tree_event->SetBranchAddress("jet_truth_ak04_eta", jet_truth_ak04_eta);
        // Loop over events
        
        Bool_t isRealData = true;
        _tree_event->GetEntry(1);
        if(nmc_truth>0) isRealData= false;
        
        for(Long64_t ievent = 0; ievent < _tree_event->GetEntries() && ievent < 3000000 ; ievent++){
            //for(Long64_t ievent = 0; ievent < 500000 ; ievent++){
            _tree_event->GetEntry(ievent);
            //Eevent Selection:
            if(not( TMath::Abs(primary_vertex[2])<10.0)) continue; //vertex z position
            if(is_pileup_from_spd_5_08) continue; //removes pileup
            N_eventpassed +=1;
            
            std::string filestring = argv[iarg];
            double weight = 1.0;
            if(eg_ntrial>0) weight = eg_cross_section/(double)eg_ntrial;
            //17g6a3 weights
            if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a3/17g6a3_pthat1.root"){
                if(ievent == 0)
                    std::cout << "PtHat 1 of 17g6a3 series opened" << std::endl;
                weight = 4.47e-11;
            }
            if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a3/17g6a3_pthat2.root"){
                if(ievent == 0)
                    std::cout << "PtHat 2 of 17g6a3 series opened" << std::endl;
                weight = 9.83e-11;
            }
            if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a3/17g6a3_pthat3.root"){
                if(ievent == 0)
                    std::cout << "PtHat 3 of 17g6a3 series opened" << std::endl;
                weight = 1.04e-10;
            }
            if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a3/17g6a3_pthat4.root"){
                if(ievent == 0)
                    std::cout << "PtHat 4 of 17g6a3 series opened" << std::endl;
                weight = 1.01e-11;
            }
            if(filestring == "/project/projectdirs/alice/NTuples/MC/16c3b/16c3b_pthat1.root"){
                if(ievent == 0)
                    std::cout << "PtHat 1 of 16c3b series opened" << std::endl;
                weight = 6.071458e-03;
            }
            if(filestring == "/project/projectdirs/alice/NTuples/MC/16c3b/16c3b_pthat2.root"){
                if(ievent == 0)
                    std::cout << "PtHat 2 of 16c3b series opened" << std::endl;
                weight = 3.941701e-03;
            }
            if(filestring == "/project/projectdirs/alice/NTuples/MC/16c3b/16c3b_pthat3.root"){
                if(ievent == 0)
                    std::cout << "PtHat 3 of 16c3b series opened" << std::endl;
                weight = 2.001984e-03;
            }
            if(filestring == "/project/projectdirs/alice/NTuples/MC/16c3b/16c3b_pthat4.root"){
                if(ievent == 0)
                    std::cout << "PtHat 4 of 16c3b series opened" << std::endl;
                weight = 9.862765e-04;
            }
            sum_of_weights += weight;
            
            //loop over jets
            for (ULong64_t kjet = 0; kjet < njet_ak04its; kjet++) {
                if(not (jet_ak04its_pt_raw[kjet]>5)) continue;
                if(not (jet_ak04its_pt_raw[kjet]<30)) continue;
                if(not (TMath::Abs(jet_ak04its_eta_raw[kjet])<0.5)) continue;
                
                N_individual_jets++;
                weightsum_individual_jets += weight;
                
                //start second jet loop
                
                for (ULong64_t ijet = kjet + 1; ijet < njet_ak04its; ijet++) { //start loop over jets
                    if(not (jet_ak04its_pt_raw[ijet]>5)) continue;
                    if(not (jet_ak04its_pt_raw[ijet]<30)) continue;
                    if(not (TMath::Abs(jet_ak04its_eta_raw[ijet])<0.5)) continue;
                    //N_jets += 1;
                    
                    // If it's a Monte-Carlo, then test phi to see if there is an "NaN" for either jet in the pair. If it is, cut out this pair
                    if (!isRealData)
                        if(isnan(jet_ak04its_phi_truth[ijet]) || isnan(jet_ak04its_phi_truth[kjet]))
                            continue;
                    
                    Float_t dphi = TMath::Abs(TVector2::Phi_mpi_pi(jet_ak04its_phi[ijet] - jet_ak04its_phi[kjet]));
                    Float_t dEta = TMath::Abs(jet_ak04its_eta_raw[ijet] - jet_ak04its_eta_raw[kjet]);
                    Float_t AvgEta = (jet_ak04its_eta_raw[ijet] + jet_ak04its_eta_raw[kjet])/2;
                    
                    hjet_dPhi.Fill(dphi,weight);
                    hjet_dEta.Fill(dEta, weight);
                    hjet_AvgEta.Fill(AvgEta, weight);
                    
                    if (jet_ak04its_pt_raw[ijet] < jet_ak04its_pt_raw[kjet]) {
                        hjet_leadingEta.Fill(jet_ak04its_eta_raw[kjet], weight);
                        hjet_leadingPhi.Fill(jet_ak04its_phi[kjet], weight);
                    }
                    else {
                        hjet_leadingEta.Fill(jet_ak04its_eta_raw[ijet], weight);
                        hjet_leadingPhi.Fill(jet_ak04its_phi[ijet], weight);
                    }
                    hjet_leadingpT.Fill(jet_ak04its_pt_raw[kjet], weight);
                    
                    
                    if (!isRealData){
                        Float_t dphi_truth = TMath::Abs(TVector2::Phi_mpi_pi(jet_ak04its_phi_truth[ijet] - jet_ak04its_phi_truth[kjet]));
                        Float_t dEta_truth = TMath::Abs(jet_ak04its_eta_truth[ijet] - jet_ak04its_eta_truth[kjet]);
                        Float_t AvgEta_truth = (jet_ak04its_eta_truth[ijet] + jet_ak04its_eta_truth[kjet])/2;
                        
                        hjet_dPhi_truth.Fill(dphi_truth,weight);
                        hjet_dEta_truth.Fill(dEta_truth, weight);
                        hjet_AvgEta_truth.Fill(AvgEta_truth, weight);
                        
                        if (jet_ak04its_pt_raw[ijet] < jet_ak04its_pt_raw[kjet]) {
                            hjet_leadingEta_truth.Fill(jet_ak04its_eta_truth[kjet], weight);
                            hjet_leadingPhi_truth.Fill(jet_ak04its_phi_truth[kjet], weight);
                        }
                        else {
                            hjet_leadingEta_truth.Fill(jet_ak04its_eta_truth[ijet], weight);
                            hjet_leadingPhi_truth.Fill(jet_ak04its_phi_truth[ijet], weight);
                        }
                    }
                    
                    if(not (dphi>(TMath::Pi()/2))) continue;
                    //std::cout <<"truthptjet: " << jet_ak04its_pt_truth[ijet] << "reco pt jet" << jet_ak04its_pt_raw[ijet] << std::endl;
                    
                    Float_t xj = 0;
                    Float_t ptD = 0;
                    Float_t mult = 0;
                    Float_t width = 0;
                    
                    if (jet_ak04its_pt_raw[ijet] < jet_ak04its_pt_raw[kjet]) {
                        xj = jet_ak04its_pt_raw[ijet]/jet_ak04its_pt_raw[kjet];
                        ptD = jet_ak04its_ptd_raw[ijet];
                        mult = jet_ak04its_multiplicity[ijet];
                        width =  jet_ak04its_width_sigma_raw[ijet][0];
                    }
                    else {
                        xj = jet_ak04its_pt_raw[kjet]/jet_ak04its_pt_raw[ijet];
                        ptD = jet_ak04its_ptd_raw[kjet];
                        mult = jet_ak04its_multiplicity[kjet];
                        width =  jet_ak04its_width_sigma_raw[ijet][0];
                    }
                    hjet_Xj.Fill(xj, weight);
                    hjet_pTD.Fill(ptD, weight);
                    hjet_Multiplicity.Fill(mult, weight);
                    hjet_jetwidth.Fill(width, weight);
                    //std::cout<<" xj " << xj << " " << " xj_truth "<< xj_truth << std::endl;
                    
                    // Fill truth events
                    if (!isRealData){
                        Float_t xj_truth = 0;
                        
                        if (jet_ak04its_pt_raw[ijet] < jet_ak04its_pt_raw[kjet]) {
                            xj_truth = jet_ak04its_pt_truth[ijet]/jet_ak04its_pt_truth[kjet];
                            hjet_leadingpT_truth.Fill(jet_ak04its_pt_truth[kjet], weight);
                        }
                        else {
                            xj_truth = jet_ak04its_pt_truth[kjet]/jet_ak04its_pt_truth[ijet];
                            hjet_leadingpT_truth.Fill(jet_ak04its_pt_truth[ijet], weight);
                        }
                        hjet_Xj_truth.Fill(xj_truth,weight);
                    }
                }//end second loop over jets
            
                
            }//end loop on jets
            //loop over truth particles
            
            //std::cout << " about to loop over mc truth particles " << std::endl;
            /**
            if (!isRealData){
            for (ULong64_t nmc = 0; nmc < njet_truth_ak04; nmc++) {
                if(not(jet_truth_ak04_pt[nmc]>5)) continue;
                if(not(jet_truth_ak04_pt[nmc]<30)) continue;
                if(not(TMath::Abs(jet_truth_ak04_eta[nmc])<0.5)) continue;
                
                hjet_leadingEta_truth.Fill(jet_truth_ak04_eta[nmc], weight);
                hjet_leadingPhi_truth.Fill(jet_truth_ak04_phi[nmc], weight);
                hjet_leadingpT_truth.Fill(jet_truth_ak04_pt[nmc], weight);
                
                N_individual_truth++;
                weightsum_individual_truth += weight;
                
                    //std::cout << " number of truth jets " << njet_truth_ak04 << std::endl;
                    for (ULong64_t ijet = nmc+1; ijet < njet_truth_ak04; ijet++) {
                        if(not(jet_truth_ak04_pt[ijet]>5.0)) continue;
                        if(not(jet_truth_ak04_pt[ijet]<30.0)) continue;
                        if(not(TMath::Abs(jet_truth_ak04_eta[ijet])<0.5)) continue;
                        
                        N_truth +=1;
                        weightsum_truth += weight;
                        
                        Float_t dphi_truth = TMath::Abs(TVector2::Phi_mpi_pi(jet_truth_ak04_phi[ijet] - jet_truth_ak04_phi[nmc]));
                        Float_t dEta_truth = TMath::Abs(jet_truth_ak04_eta[ijet] - jet_truth_ak04_eta[nmc]);
                        Float_t AvgEta_truth = (jet_truth_ak04_eta[ijet] + jet_truth_ak04_eta[nmc])/2;
                        //std::cout<< dphi_truth << std::endl;
                        
                        hjet_dPhi_truth.Fill(dphi_truth,weight);
                        hjet_dEta_truth.Fill(dEta_truth, weight);
                        hjet_AvgEta_truth.Fill(AvgEta_truth, weight);
                        
                        if( not(dphi_truth>(TMath::Pi()/2))) continue;
                        
                        Float_t xj_truth = 0;
                        Float_t leading_pT_truth = 0;
                        if (jet_ak04its_pt_raw[ijet] < jet_ak04its_pt_raw[nmc]) {
                            xj_truth = jet_truth_ak04_pt[ijet]/jet_truth_ak04_pt[nmc];
                        }
                        else {
                            xj_truth = jet_truth_ak04_pt[nmc]/jet_truth_ak04_pt[ijet];
                        }
                        hjet_Xj_truth.Fill(xj_truth,weight);
                        
                    }//end loop over truth jets
            }//end loop over mc particles
            }
            */
             
            // Create the file label, to be used within the filenames, to represent the source file
        
            
            if (ievent % 10000 == 0) {
                //SR_Xj.Draw("e1x0nostack");
                hjet_dPhi.Draw("e1x0nostack");
                std::cout << ievent << " " << _tree_event->GetEntries() << std::endl;
                std::cout << "Weight: " << weight << std::endl;
                canvas->Update();
            }
            
        }//end over events
        
    }//end of arguments
    
    std::cout << " Numbers of events passing selection " << N_eventpassed << std::endl;
    std::cout << " Number of individual jets " << N_individual_jets << std::endl;
    //std::cout << " Number of truth jets " << N_truth << std::endl;
    
    std::string opened_files = "";
    for (int iarg = 1; iarg < argc; iarg++) {
        std::string filepath = argv[iarg];
        opened_files += "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
    }
    
    TFile* fout = new TFile(Form("JetJet_config%s.root", opened_files.c_str()),"RECREATE");
    fout->Print();
    hjet_Xj.Scale(1.0/weightsum_individual_jets);
    hjet_dPhi.Scale(1.0/weightsum_individual_jets);
    hjet_pTD.Scale(1.0/weightsum_individual_jets);
    hjet_Multiplicity.Scale(1.0/weightsum_individual_jets);
    hjet_jetwidth.Scale(1.0/weightsum_individual_jets);
    hjet_dEta.Scale(1.0/weightsum_individual_jets);
    hjet_AvgEta.Scale(1.0/weightsum_individual_jets);
    hjet_leadingpT.Scale(1.0/weightsum_individual_jets);
    hjet_leadingEta.Scale(1.0/weightsum_individual_jets);
    hjet_leadingPhi.Scale(1.0/weightsum_individual_jets);
    
    hjet_dPhi_truth.Scale(1.0/weightsum_individual_jets);
    hjet_Xj_truth.Scale(1.0/weightsum_individual_jets);
    hjet_dEta_truth.Scale(1.0/weightsum_individual_jets);
    hjet_AvgEta_truth.Scale(1.0/weightsum_individual_jets);
    hjet_leadingpT_truth.Scale(1.0/weightsum_individual_jets);
    hjet_leadingEta_truth.Scale(1.0/weightsum_individual_jets);
    hjet_leadingPhi_truth.Scale(1.0/weightsum_individual_jets);
    
    hjet_Xj.Write("hjet_Xj");
    hjet_Xj_truth.Write("hjet_Xj_truth");
    
    hjet_dPhi.Write("hjet_dPhi");
    hjet_dPhi_truth.Write("hjet_dPhi_truth");
    
    hjet_pTD.Write("hjet_ptD");
    hjet_Multiplicity.Write("hjet_Multiplicity");
    hjet_jetwidth.Write("hjet_jetwidth");
    
    hjet_dEta.Write("hjet_dEta");
    hjet_dEta_truth.Write("hjet_dEta_truth");
    
    hjet_AvgEta.Write("hjet_AvgEta");
    hjet_AvgEta_truth.Write("hjet_AvgEta_truth");
    
    hjet_leadingpT.Write("hjet_leading_pT");
    hjet_leadingpT_truth.Write("hjet_leading_pT_truth");
    
    hjet_leadingEta.Write("hjet_leading_Eta");
    hjet_leadingEta_truth.Write("hjet_leading_Eta_truth");
    
    hjet_leadingPhi.Write("hjet_leading_Phi");
    hjet_leadingPhi_truth.Write("hjet_leadingPhi_truth");
    
    std::cout << "Sum of weights over events: " << sum_of_weights << "\nover Nindividualjets: " << weightsum_individual_jets << std::endl;
    std::cout << " ending " << std::endl;
    
    return EXIT_SUCCESS;
    
}
