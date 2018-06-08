// This file selects detected and generated photons from MC files to determine efficiency, which it graphs against momentum as a 1D plot and against phi and eta as a 2D plot

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

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;

//enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

// 1D histogram dividing function
// Precondition: the two histograms must have the same x-dimensions
TH1D* divide_histograms1D(TH1D* graph1, TH1D* graph2){
    // Make the 1D histogram to contain the quotient and find minimum and maximum bins along both axes
    TH1D* quotient = new TH1D(*graph2);
    double x_bin_min = quotient->GetXaxis()->FindBin(quotient->GetXaxis()->GetXmin());
    double x_bin_max = quotient->GetXaxis()->FindBin(quotient->GetXaxis()->GetXmax());
    
    // Loop over all bins, divide the element from graph 1 by its counterpart on graph 2. Manually propogate the error
    for(int i = x_bin_min; i <= x_bin_max; i++){
            if (graph2->GetBinContent(i) != 0)
                quotient->SetBinContent(i, (graph1->GetBinContent(i))/(graph2->GetBinContent(i)));
            // Failsafe
            else
                quotient->SetBinContent(i, 0);
            
            // Error propogation
            if (graph2->GetBinContent(i) == 0) {
                if (graph1->GetBinContent(i) == 0)
                    quotient->SetBinError(i, 0);
                else
                    quotient->SetBinError(i, (quotient->GetBinContent(i)*TMath::Sqrt(( (graph1->GetBinError(i)/graph1->GetBinContent(i)) * (graph1->GetBinError(i)/graph1->GetBinContent(i)) ) )));
            }
            else if (graph1->GetBinContent(i) == 0)
                quotient->SetBinError(i, (quotient->GetBinContent(i)*TMath::Sqrt(( (graph2->GetBinError(i)/graph2->GetBinContent(i)) * (graph2->GetBinError(i)/graph2->GetBinContent(i)) ) )));
            else
                quotient->SetBinError(i, (quotient->GetBinContent(i)*TMath::Sqrt(( (graph1->GetBinError(i)/graph1->GetBinContent(i)) * (graph1->GetBinError(i)/graph1->GetBinContent(i)) ) + ( (graph2->GetBinError(i)/graph2->GetBinContent(i)) * (graph2->GetBinError(i)/graph2->GetBinContent(i)) ) )));
        }
    
    
    return quotient;
}

// 2D histogram dividing function
// Precondition: the two histograms must have the same y- and x-dimensions
TH2D* divide_histograms2D(TH2D* graph1, TH2D* graph2){
    // Make the 2D histogram to contain the quotient and find minimum and maximum bins along both axes
    TH2D* quotient = new TH2D(*graph1);
    double x_bin_min = quotient->GetXaxis()->FindBin(quotient->GetXaxis()->GetXmin());
    double x_bin_max = quotient->GetXaxis()->FindBin(quotient->GetXaxis()->GetXmax());
    double y_bin_min = quotient->GetYaxis()->FindBin(quotient->GetYaxis()->GetXmin());
    double y_bin_max = quotient->GetYaxis()->FindBin(quotient->GetYaxis()->GetXmax());
    
    // Loop over all bins, divide the element from graph 1 by its counterpart on graph 2. Manually propogate the error
    for(int i = x_bin_min; i <= x_bin_max; i++)
        for(int j = y_bin_min; j <= y_bin_max; j++) {
            if (graph2->GetBinContent(i, j) != 0)
                quotient->SetBinContent(i, j, (graph1->GetBinContent(i, j))/(graph2->GetBinContent(i, j)));
            // Failsafe
            else
                quotient->SetBinContent(i, j, 0);
            
            // Error propogation
            if (graph2->GetBinContent(i, j) == 0) {
                if (graph1->GetBinContent(i, j) == 0)
                    quotient->SetBinError(i, j, 0);
                else
                    quotient->SetBinError(i, j, (quotient->GetBinContent(i, j)*TMath::Sqrt(( (graph1->GetBinError(i, j)/graph1->GetBinContent(i, j)) * (graph1->GetBinError(i, j)/graph1->GetBinContent(i, j)) ) )));
            }
            else if (graph1->GetBinContent(i, j) == 0)
                quotient->SetBinError(i, j, (quotient->GetBinContent(i, j)*TMath::Sqrt(( (graph2->GetBinError(i, j)/graph2->GetBinContent(i, j)) * (graph2->GetBinError(i, j)/graph2->GetBinContent(i, j)) ) )));
            else
                quotient->SetBinError(i, j, (quotient->GetBinContent(i, j)*TMath::Sqrt(( (graph1->GetBinError(i, j)/graph1->GetBinContent(i, j)) * (graph1->GetBinError(i, j)/graph1->GetBinContent(i, j)) ) + ( (graph2->GetBinError(i, j)/graph2->GetBinContent(i, j)) * (graph2->GetBinError(i, j)/graph2->GetBinContent(i, j)) ) )));
        }
    
    
    return quotient;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }
    int dummyc = 1;
    char **dummyv = new char *[1];
    
    dummyv[0] = strdup("main");
    
    // Loop over files to open, I guess. I may be wrong though.
    for (int iarg = 1; iarg < argc; iarg++) {
        std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
        TFile *file = TFile::Open((TString)argv[iarg]);
        
        if (file == NULL) {
            std::cout << " fail" << std::endl;
            exit(EXIT_FAILURE);
        }
        file->Print();
        
        // Get all the TTree variables from the file to open, I guess
        TTree *_tree_event = NULL;
        _tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
        
        if (_tree_event == NULL) {
            std::cout << "First try did not got (AliAnalysisTaskNTGJ does not exist), trying again" << std::endl;
            _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
            if (_tree_event == NULL) {
                std::cout << " fail " << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        //_tree_event->Print();
        
        TApplication application("", &dummyc, dummyv);
        
        // 1D Histograms
        TH1D* hist_measured = new TH1D("hist_measured", "", 50, 0, 50);
        TH1D* hist_generated = new TH1D("hist_generated", "", 50, 0, 50);
        TH1D* hist_ratio = new TH1D("hist_ratio", "", 50, 0, 50);
        
        // 2D Histograms
        TH2D* phietamap_measured = new TH2D("phi_eta_map_measured", "", 40, 1.2, 3.2, 32, -0.8, 0.8);
        TH2D* phietamap_generated = new TH2D("phi_eta_map_generated", "", 40, 1.2, 3.2, 32, -0.8, 0.8);
        TH2D* phietamap_ratio = new TH2D("phi_eta_map_ratio", "", 40, 1.2, 3.2, 32, -0.8, 0.8);
        
        TCanvas* canvas = new TCanvas();
        
        //you define variables
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
        
        //Int_t eg_ntrial;
        
        Float_t eg_cross_section;
        Int_t   eg_ntrial;
        
        // Set the branch addresses of the branches in the TTrees
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
        
        std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;
        
        // Loop over events
        for(Long64_t ievent = 0; ievent < _tree_event->GetEntries() ; ievent++){
            //for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
            _tree_event->GetEntry(ievent);
            
            // Weight the event by simulation pT bin
            double weight = 1.0;
            if(eg_ntrial>0) weight = eg_cross_section/(double)eg_ntrial;
            
            //loop over clusters
            for (ULong64_t n = 0; n < ncluster; n++) {
                // Apply cuts
                
                // Photon Selection
                //if( not(cluster_pt[n]>10)) continue; //select pt of photons
                if( not(cluster_ncell[n]>2)) continue;   //removes clusters with 1 or 2 cells
                if( not(cluster_e_cross[n]/cluster_e[n]>0.05)) continue; //removes "spiky" clusters
                if( not(cluster_nlocal_maxima[n]<= 2)) continue; //require to have at most 2 local maxima.
                if( not(cluster_distance_to_bad_channel[n]>=2.0)) continue;
                
                //Isolation and shower shape selection
                //if( not(cluster_s_nphoton[n][1]<=0.85)) continue; // NN max: deep photons
                //if( not(cluster_s_nphoton[n][1]>=0.55)) continue; // NN min: deep photons
                if( not(cluster_iso_its_04[n] < 1.0)) continue;
                if( not(cluster_lambda_square[n][0]<0.27)) continue;
                
                // Truth-matching cut
                Bool_t isTruePhoton = false;
                //Float_t truth_pt = -999.0;
                for(int counter = 0 ; counter<32; counter++){
                    unsigned short index = cluster_mc_truth_index[n][counter];
                    
                    if(isTruePhoton) break;
                    if(index==65535) continue;
                    //std::cout<<"truth, pt: " << mc_truth_pt[index] << "phi " << mc_truth_phi[index] << " eta " << mc_truth_eta[index] << " code: " << mc_truth_pdg_code[index] << " status " << int(mc_truth_status[index]) << " parentpdg " << mc_truth_first_parent_pdg_code[index] << std::endl;
                    if(mc_truth_pdg_code[index]!=22) continue;
                    if(mc_truth_first_parent_pdg_code[index]!=22) continue;
                    if( not (mc_truth_status[index] >0)) continue;
                    isTruePhoton = true;
                    //truth_pt = mc_truth_pt[index];
                }//end loop over indices
                
                // Fill the measured histogram bins
                if(isTruePhoton){
                    double phi = cluster_phi[n]/(TMath::Pi());
                    while(phi < 1.2)
                        phi += 2;
                    while (phi > 3.2)
                        phi -= 2;
                    double eta = cluster_eta[n];
                    if(TMath::Abs(eta) > 0.7) continue;
                    hist_measured->Fill(cluster_pt[n], weight);
                    phietamap_measured->Fill(phi, cluster_eta[n], weight);
                }
                
            }//end loop on clusters.
            
            int mctruths_rejected = 0;
            int mctruths_accepted = 0;
            
            //loop over mc truth
            for (unsigned int m = 0; m < nmc_truth; m++) {
                // Apply cuts
                
                //std::cout << mc_truth_pt[m] << "phi " << mc_truth_phi[m] << " eta " << mc_truth_eta[m] << " code: " << mc_truth_pdg_code[m] << " status " << int(mc_truth_status[m]) << " parentpdg " << mc_truth_first_parent_pdg_code[m] << std::endl;
                if( not(mc_truth_pdg_code[m]==22 && int(mc_truth_status[m])>0 &&  mc_truth_first_parent_pdg_code[m]==22))
                    {mctruths_rejected++;  continue; }
                mctruths_accepted++;
                
                // Fill the measured histogram bin
                double phi = mc_truth_phi[m]/(TMath::Pi());
                while(phi < 1.2)
                    phi += 2;
                while (phi > 3.2)
                    phi -= 2;
                double eta = mc_truth_eta[m];
                // Cut out non-DCal and non-EMCal data
                if(TMath::Abs(eta) > 0.7) continue;
                if(phi < 1.45) continue;
                if(phi >= 1.45 && phi < 1.75 && TMath::Abs(eta) < 0.2) continue;
                if(phi > 1.85 && phi < 2.45) continue;
                if(phi > 3.05) continue;
                hist_generated->Fill(mc_truth_pt[m], weight);
                phietamap_generated->Fill(phi, mc_truth_eta[m], weight);
                
            }// end loop on mc truth
            
            if (ievent % 10000 == 0)
                std::cout<< "Event no. " << ievent << std::endl;
            
        } //end loop over events
        
        // Calculate ratio
        hist_ratio = divide_histograms1D(hist_measured, hist_generated);
        phietamap_ratio = divide_histograms2D(phietamap_measured, phietamap_generated);
        
        // Create the file label, to be used within the filenames, to represent the source file
        std::string opened_files = "";
        for (int iarg = 1; iarg < argc; iarg++) {
            std::string filepath = argv[iarg];
            
            opened_files += "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
        }
        
        // Draw all graphs
        TFile* efficiencyOut = new TFile(Form("fout_efficiency%s.root", opened_files.c_str()), "RECREATE");
        hist_measured->Write("measured_photons");
        hist_generated->Write("generated_photons");
        hist_ratio->Write("ratio_photons");
        phietamap_measured->Write("phi_eta_map_measured");
        phietamap_generated->Write("phi_eta_map_generated");
        phietamap_ratio->Write("phi_eta_map_ratio");
        efficiencyOut->Close();
        
        // Save the ratio graphs
        hist_measured->SetTitle("Measured Photons; P_{T} (GeV); # of photons");
        hist_measured->GetYaxis()->SetTitleOffset(1.5);
        hist_measured->Draw();
        canvas->SaveAs(Form("measured_photons%s.png", opened_files.c_str()));
        canvas->Clear();
        
        hist_generated->SetTitle("Generated Photons; P_{T} (GeV); # of photons");
        hist_generated->GetYaxis()->SetTitleOffset(1.5);
        hist_generated->Draw();
        canvas->SaveAs(Form("generated_photons%s.png", opened_files.c_str()));
        canvas->Clear();
        
        hist_ratio->SetTitle("Efficiency; P_{T} (GeV); efficiency");
        hist_ratio->GetYaxis()->SetTitleOffset(1.5);
        hist_ratio->Draw();
        canvas->SaveAs(Form("efficiency%s.png", opened_files.c_str()));
        canvas->Clear();
        
        phietamap_measured->SetTitle("Measured Photons; #phi (#frac{rad}{#pi}); #eta");
        phietamap_measured->GetXaxis()->SetTitleOffset(1.5);
        phietamap_measured->GetYaxis()->SetTitleOffset(1.5);
        phietamap_measured->Draw("COLZ");
        canvas->SaveAs(Form("measured_photons_phietamap_%s.png", opened_files.c_str()));
        canvas->Clear();
        
        phietamap_generated->SetTitle("Generated Photons; #phi (#frac{rad}{#pi}); #eta");
        phietamap_generated->GetXaxis()->SetTitleOffset(1.5);
        phietamap_generated->GetYaxis()->SetTitleOffset(1.5);
        phietamap_generated->Draw("COLZ");
        canvas->SaveAs(Form("generated_photons_phietamap_%s.png", opened_files.c_str()));
        canvas->Clear();
        
        phietamap_ratio->SetTitle("Efficiency; #phi (#frac{rad}{#pi}); #eta");
        phietamap_ratio->GetXaxis()->SetTitleOffset(1.5);
        phietamap_ratio->GetYaxis()->SetTitleOffset(1.5);
        phietamap_ratio->Draw("COLZ");
        canvas->SaveAs(Form("efficiency_phietamap_%s.png", opened_files.c_str()));
        canvas->Clear();
        
    }//end loop over samples
    
    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}

