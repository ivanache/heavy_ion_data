/*
  This code is a variation of Fernando's Skeleton_Mix_Correlations.cc code, this time for events
*/
// Author: Ivan Chernyshev; Creator of template code: Fernando Torales-Acosta

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
#include "H5Cpp.h"

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

using namespace H5;

int main(int argc, char *argv[])
{
    if (argc < 2) {
      fprintf(stderr,"Please Indicate ROOT file");
        exit(EXIT_FAILURE);
    }
    
    // Number of events
    int nevents = 0;
    
    //HISTOGRAMS
    TCanvas canvas("canvas", "");
    
    TH1D* z_Vertices_individual = new TH1D("Primary_Vertex_root", "Z-vertex (ROOT)", 240, -12, 12);
    
    TH1D* Multiplicity_individual = new TH1D("Multiplicity_root", "Multiplicity (ROOT)", 1000, 0, 1000);
    
    TH1D* jet_pT_ROOT = new TH1D("jet_pt_distribution", "Jet p_{T} distribution (ROOT)", 60, -15, 15);
    TH1D* jet_eta_ROOT = new TH1D("jet_eta_distribution", "Jet #eta distribution (ROOT)", 60, -0.8, 0.8);
    TH1D* jet_phi_ROOT = new TH1D("jet_phi_distribution", "Jet #phi distribution (ROOT)", 60, -3.1415926, 3.1415926);
    TH1D* jet_pTD_ROOT = new TH1D("jet_ptD_distribution", "Jet p_{T}D distribution (ROOT)", 60, -1, 1);
    TH1D* jet_Multiplicity_ROOT = new TH1D("jet_ptD_distribution", "Jet p_{T}D distribution (ROOT)", 60, -15, 15);
    
    for (int iarg = 1; iarg < argc; iarg++) {
        std::string filestring = (std::string)argv[iarg];
        std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
        TFile *file = TFile::Open((TString)argv[iarg]);
        
        if (file == NULL) {
            std::cout << " fail; could not open file" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        // Get all the TTree variables from the file to open, I guess
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
        
        //
        
        std::cout <<"_tree_event->GetEntries() " << _tree_event->GetEntries() << std::endl;
        
        //you define variables
        Double_t primary_vertex[3];
        Bool_t is_pileup_from_spd_5_08;
        Bool_t is_pileup_from_spd_3_08;
        Float_t ue_estimate_its_const;
        Float_t ue_estimate_tpc_const;
        
        UInt_t ntrack;
        Float_t track_e[NTRACK_MAX];
        Float_t track_pt[NTRACK_MAX];
        Float_t track_eta[NTRACK_MAX];
        Float_t track_phi[NTRACK_MAX];
        UChar_t track_quality[NTRACK_MAX];
        
        UInt_t ncluster;
        Float_t cluster_e[NTRACK_MAX];
        Float_t cluster_e_cross[NTRACK_MAX];
        Float_t cluster_e_max[NTRACK_MAX];
        Float_t cluster_pt[NTRACK_MAX];
        Float_t cluster_eta[NTRACK_MAX];
        Float_t cluster_phi[NTRACK_MAX];
        
        Float_t cluster_iso_tpc_04[NTRACK_MAX];
        Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
        
        Float_t cluster_iso_its_04[NTRACK_MAX];
        
        Float_t cluster_iso_its_04_ue[NTRACK_MAX];
        
        
        Float_t cluster_frixione_its_04_02[NTRACK_MAX];
        Float_t cluster_s_nphoton[NTRACK_MAX][4];
        UChar_t cluster_nlocal_maxima[NTRACK_MAX];
        Float_t cluster_distance_to_bad_channel[NTRACK_MAX];
        
        unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
        Int_t cluster_ncell[NTRACK_MAX];
        UShort_t  cluster_cell_id_max[NTRACK_MAX];
        Float_t cluster_lambda_square[NTRACK_MAX][2];
        Float_t cell_e[17664];
        
        //Jets reco
        UInt_t njet_ak04its;
        Float_t jet_ak04its_pt_raw[NTRACK_MAX];
        Float_t jet_ak04its_eta_raw[NTRACK_MAX];
        Float_t jet_ak04its_phi[NTRACK_MAX];
        
        Float_t jet_ak04its_pt_truth[NTRACK_MAX];
        Float_t jet_ak04its_eta_truth[NTRACK_MAX];
        Float_t jet_ak04its_phi_truth[NTRACK_MAX];
        
        //The z_reco is defined as the fraction of the true jet that ended up in this reco jet
        //There are two entries and indices, the first is the best.
        Int_t   jet_ak04its_truth_index_z_reco[NTRACK_MAX][2];
        Float_t jet_ak04its_truth_z_reco[NTRACK_MAX][2];
        Float_t jet_ak04its_ptd_raw[NTRACK_MAX];
        Float_t jet_ak04its_width_sigma[NTRACK_MAX][2];
        UShort_t jet_ak04its_multiplicity[NTRACK_MAX];
        
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
        
        ULong64_t trigger_mask[2];
        
        // Set the branch addresses of the branches in the TTrees
        //    _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
        _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
        _tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
        _tree_event->SetBranchAddress("is_pileup_from_spd_3_08", &is_pileup_from_spd_3_08);
        _tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);
        _tree_event->SetBranchAddress("ue_estimate_tpc_const", &ue_estimate_tpc_const);
        
        _tree_event->SetBranchAddress("trigger_mask", &trigger_mask);
        
        
        _tree_event->SetBranchAddress("ntrack", &ntrack);
        _tree_event->SetBranchAddress("track_e", track_e);
        _tree_event->SetBranchAddress("track_pt", track_pt);
        _tree_event->SetBranchAddress("track_eta", track_eta);
        _tree_event->SetBranchAddress("track_phi", track_phi);
        _tree_event->SetBranchAddress("track_quality", track_quality);
        
        _tree_event->SetBranchAddress("ncluster", &ncluster);
        _tree_event->SetBranchAddress("cluster_e", cluster_e);
        _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
        _tree_event->SetBranchAddress("cluster_e_max", cluster_e_max);
        _tree_event->SetBranchAddress("cluster_pt", cluster_pt); // here
        _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
        _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
        _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton); // here
        _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
        _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
        
        _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
        _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
        
        _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
        _tree_event->SetBranchAddress("cluster_iso_its_04_ue",cluster_iso_its_04_ue);
        
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
        
        //quark-gluon discriminator variables
        _tree_event->SetBranchAddress("jet_ak04its_ptd_raw", jet_ak04its_ptd_raw);
        _tree_event->SetBranchAddress("jet_ak04its_width_sigma", jet_ak04its_width_sigma);
        _tree_event->SetBranchAddress("jet_ak04its_multiplicity_raw", jet_ak04its_multiplicity);
        
        
        
        _tree_event->SetBranchAddress("jet_ak04its_truth_index_z_reco",     jet_ak04its_truth_index_z_reco);
        _tree_event->SetBranchAddress("jet_ak04its_truth_z_reco", jet_ak04its_truth_z_reco);
        
        //truth jets
        _tree_event->SetBranchAddress("njet_truth_ak04", &njet_truth_ak04);
        _tree_event->SetBranchAddress("jet_truth_ak04_pt", jet_truth_ak04_pt);
        _tree_event->SetBranchAddress("jet_truth_ak04_phi", jet_truth_ak04_phi);
        _tree_event->SetBranchAddress("jet_truth_ak04_eta", jet_truth_ak04_eta);
        // Loop over events
        
        _tree_event->GetEntry(1);
        
        if( not(nevents>0)){
            nevents = _tree_event->GetEntries();
        }
        
        for(Long64_t ievent = 0; ievent < nevents ; ievent++){
            if (ievent % 100000 == 0) std::cout << " event " << ievent << std::endl;
            _tree_event->GetEntry(ievent);
            
            // Loop over jets
            for (Long64_t ijet = 0; ijet < njet_ak04its; ijet++) {
                double jet_pT = -9000;
                double jet_phi = -9000;
                double jet_eta = -9000;
                double jet_pTD = -9000;
                double jet_Multiplicity = -9000;
                
                jet_pT = jet_ak04its_pt_raw[ijet];
                jet_phi = jet_ak04its_phi[ijet];
                jet_eta = jet_ak04its_eta_raw[ijet];
                jet_pTD = jet_ak04its_ptd_raw[ijet];
                jet_Multiplicity = jet_ak04its_multiplicity[ijet];
                
                // Cut on pT jet = 0 and pT jet NaN
                if(TMath::IsNaN(jet_pT)) continue;
                if(TMath::IsNaN(jet_eta)) continue;
                if(TMath::IsNaN(jet_phi)) continue;
                if(TMath::IsNaN(jet_pTD)) continue;
                if(TMath::IsNaN(jet_Multiplicity)) continue;
                //if(not(jet_data_out[0][ijet][0] == 0)) {continue;}
                
                // Fill the histogram
                jet_pT_ROOT->Fill(jet_pT);
                jet_eta_ROOT->Fill(jet_eta);
                jet_phi_ROOT->Fill(jet_phi);
                jet_pTD_ROOT->Fill(jet_pTD);
                jet_Multiplicity_ROOT->Fill(jet_Multiplicity);
            }
        }//end loop over events
    }
        // if(ievent % 10000 == 0)
        //     std::cout << "Event " << ievent << std::endl;

    
    //very particular about file names to ease scripting
    // Write to fout
    std::string filepath = argv[1];
    std::string opened_files = "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
    //std::string rawname = std::string(argv[1]);
    TFile* fout = new TFile("ROOT_jetpT_Distribution.root","RECREATE");
    std::cout<< "Created ROOT file " <<std::endl;
    
    
    //Write histograms here
    z_Vertices_individual->Write();
    Multiplicity_individual->Write();
    jet_pT_ROOT->Write();
    jet_eta_ROOT->Write();
    jet_phi_ROOT->Write();
    jet_pTD_ROOT->Write();
    jet_Multiplicity_ROOT->Write();
    
    TCanvas* c = new TCanvas();
    jet_pT_ROOT->Draw();
    c->SaveAs("ROOT_jetpT_Distribution.pdf");
    c->Clear();
    
    jet_eta_ROOT->Draw();
    c->SaveAs("ROOT_jeteta_Distribution.pdf");
    c->Clear();
    
    jet_phi_ROOT->Draw();
    c->SaveAs("ROOT_jetphi_Distribution.pdf");
    c->Clear();
    
    jet_pTD_ROOT->Draw();
    c->SaveAs("ROOT_jetpTD_Distribution.pdf");
    c->Clear();
    
    jet_Multiplicity_ROOT->Draw();
    c->SaveAs("ROOT_jetMultiplicity_Distribution.pdf");
    c->Clear();
    
    fout->Close();
    
    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
