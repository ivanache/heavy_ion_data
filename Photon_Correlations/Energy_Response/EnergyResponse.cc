/**
   This program produces energy response plots from Monte-Carlo simulations
*/

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
        TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
        
        if (_tree_event == NULL) {
            std::cout << " fail " << std::endl;
            exit(EXIT_FAILURE);
        }
        
        // The histograms
        TH2D* energyresponse_scatter = new TH2D("measured_pt_chart", "", 22, 10, 21, 26, 10, 23);
        TH1D* energyresolution = new TH1D("energy_resolution", "", 80, -.4, .4);
        
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
        
        Float_t mc_truth_first_parent_e[NTRACK_MAX];
        Float_t mc_truth_first_parent_pt[NTRACK_MAX];
        Float_t mc_truth_first_parent_eta[NTRACK_MAX];
        Float_t mc_truth_first_parent_phi[NTRACK_MAX];
        UChar_t mc_truth_status[NTRACK_MAX];
        
        
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
        
        _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
        _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
        _tree_event->SetBranchAddress("cell_e", cell_e);
        
        _tree_event->SetBranchAddress("nmc_truth", &nmc_truth);
        _tree_event->SetBranchAddress("mc_truth_pdg_code", mc_truth_pdg_code);
        _tree_event->SetBranchAddress("mc_truth_pt", mc_truth_pt);
        _tree_event->SetBranchAddress("mc_truth_phi", mc_truth_phi);
        _tree_event->SetBranchAddress("mc_truth_eta", mc_truth_eta);
        
        
        // Loop over events
        for(Long64_t ievent = 0; ievent < _tree_event->GetEntries() ; ievent++){
            //for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
            _tree_event->GetEntry(ievent);
            
            double detectedpT = -1;
            double generatedpT = -1;
            
            //loop over clusters
            for (ULong64_t n = 0; n < ncluster; n++) {
                // Apply cuts
                if( not(cluster_pt[n]>10)) {continue;} //select pt of photons
                if( not(cluster_s_nphoton[n][1]<=0.85)) {continue;} // NN max: deep photons
                if( not(cluster_s_nphoton[n][1]>=0.55)) {continue;} // NN min: deep photons
                
                // Keep the cluster pt of the cluster within the event that passed the cuts (there should be only one)
                detectedpT = cluster_pt[n];
                break;
                
            }//end loop on clusters
            
            //loop over mc truth
            for (unsigned int m = 0; m < nmc_truth; m++) {
                // Apply cuts
                if( not(mc_truth_pt[m]>10)) {continue;} //select pt of photons
                if( not(mc_truth_pdg_code[m]==22)) {continue;} // MC truth pT
                
                // Keep the truth pt of the cluster within the event that passed the cuts (there should be only one, any more than one is probably a bug)
                generatedpT = mc_truth_pt[m];
                
                
            }// end loop on mc truth
            
            // -1 (see above) is an indicator value for no results
            if (detectedpT != -1 && generatedpT != -1) {
                energyresponse_scatter->Fill(generatedpT, detectedpT);
                energyresolution->Fill((generatedpT - detectedpT)/(generatedpT));
            }
            
            if(ievent % 10000 == 0)
                std::cout << "Event " << ievent << " has been processed" << std::endl;
            
        }//end loop on events
        
        // Create the file label, to be used within the filenames, to represent the source file
        std::string opened_files = "";
        for (int iarg = 1; iarg < argc; iarg++) {
            std::string filepath = argv[iarg];
            
            opened_files += "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
        }
        
        // Draw all graphs
        TFile* energyresponseOut = new TFile(Form("fout_energyresponse%s.root", opened_files.c_str()), "RECREATE");
        energyresponse_scatter->SetTitle("p_{T} response; True p_{T} (GeV); Measured p_{T} (GeV)");
        energyresolution->SetTitle("p_{T} resolution; #frac{true p_{T} - meas p_{T}}{true p_{T}}; Number of Photons");
        energyresponse_scatter->Write("energy_response");
        energyresolution->Write("energy_resolution");
        energyresponseOut->Close();
        
        energyresponse_scatter->Draw("COLZ");
        canvas->SaveAs(Form("energy_response_scatter_%s.png", opened_files.c_str()));
        canvas->Clear();
        
        energyresolution->Draw();
        canvas->SaveAs(Form("energy_response_resolution_%s.png", opened_files.c_str()));
        canvas->Clear();
        
        canvas->Close();
    } // End loop over samples
    
    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
