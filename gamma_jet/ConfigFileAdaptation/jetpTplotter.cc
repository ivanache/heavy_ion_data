/**
   This program plots the number of gamma-jet pairs in a data sample with respect to the jet pT, both in a 2D plot where cluster pT is also included and in a 1D plot where the data is cut to a certain relevant cluster pT interval
*/
// Author: Ivan Chernyshev; Date: 9/27/2018

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
#include <set>

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};
enum photon_IDVARS {LAMBDA_0, DNN, EMAX_OVER_ECLUSTER};

int main(int argc, char *argv[]) {
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }
    
    // Read configuration file for the cut variables
    FILE* config = fopen("GammaJet_config.yaml", "r");
    if (config == NULL)  std::cout<<"no config"<<std::endl;
    // Default values of various variables used in the file (actual values are to be determined by the configuration file)
    // Cut variables
    double primary_vertex_max = 0.0;
    double SIG_DNN_min = 0.0;
    double SIG_DNN_max = 0.0;
    double BKG_DNN_min = 0.0;
    double BKG_DNN_max = 0.0;
    double SIG_lambda_min = 0.0;
    double SIG_lambda_max = 0.0;
    double BKG_lambda_min = 0.0;
    double BKG_lambda_max = 0.0;
    double SIG_Emax_over_Ecluster_min = 0.0;
    double SIG_Emax_over_Ecluster_max = 0.0;
    double BKG_Emax_over_Ecluster_min = 0.0;
    double BKG_Emax_over_Ecluster_max = 0.0;
    double clus_pT_min = 0;
    double clus_pT_max = 0;
    double track_pT_max = 0;
    double jet_pT_min = 0.0;
    double Eta_max = 0.0;
    double Cluster_ncell_min = 0;
    double Cluster_locmaxima_max = 0.0;
    double Cluster_distobadchannel = 0.0;
    double EcrossoverE_min = 0.00;
    
    // The bounds for the events to fal into the isolation and nonisolation areas
    double iso_max = 0.0;
    double noniso_min = 0.0;
    double noniso_max = 0.0;
    
    // Delta eta
    double deta_max = 0.5;
    
    // Number of bins in correlation functions
    int xjbins = 10;
    int phibins = 5;
    int etabins = 20;
    
    // Which branch should be used to determine whether a cluster should fall into iso, noniso, or neither
    isolationDet determiner = CLUSTER_ISO_ITS_04;
    photon_IDVARS photon_identifier = LAMBDA_0;
    
    // Truth cuts
    int rightpdgcode = 22;
    int rightparentpdgcode = 22;
    
    // Number of events
    int nevents = 0;
    // Loop through config file
    char line[MAX_INPUT_LENGTH];
    while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
        if (line[0] == '#') {
            continue;
        }
        
        // Declare char arrays needed to read the line
        char key[MAX_INPUT_LENGTH];
        char dummy[MAX_INPUT_LENGTH];
        char value[MAX_INPUT_LENGTH];
        
        // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays
        key[0] = '\0';
        value[0] = '\0';
        sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);
        
        // Use if statements to detect, based on key, which variable the line's content should be used to fill and fill that variable
        if (strcmp(key, "primary_vertex_max") == 0) {
            // Assign primary_vertex_max to the double-converted version of value
            primary_vertex_max = atof(value);
            std::cout << "primary_vertex_max is " << primary_vertex_max << std::endl;
        }
        else if (strcmp(key, "SIG_DNN_min") == 0) {
            // Assign SIG_DNN_min to the double-converted version of value
            SIG_DNN_min = atof(value);
            std::cout << "SIG_DNN_min is " << SIG_DNN_min << std::endl;
        }
        else if (strcmp(key, "SIG_DNN_max") == 0) {
            // Assign SIG_DNN_max to the double-converted version of value
            SIG_DNN_max = atof(value);
            std::cout << "SIG_DNN_max is " << SIG_DNN_max << std::endl;
        }
        else if (strcmp(key, "BKG_DNN_min") == 0) {
            // Assign BKG_DNN_min to the double-converted version of value
            BKG_DNN_min = atof(value);
            std::cout << "BKG_DNN_min is " << BKG_DNN_min << std::endl;
        }
        else if (strcmp(key, "BKG_DNN_max") == 0) {
            // Assign BKG_DNN_max to the double-converted version of value
            BKG_DNN_max = atof(value);
            std::cout << "BKG_DNN_max is " << BKG_DNN_max << std::endl;
        }
        else if (strcmp(key, "SIG_lambda_min") == 0) {
            // Assign SIG_lambda_min to the double-converted version of value
            SIG_lambda_min = atof(value);
            std::cout << "SIG_lambda_min is " << SIG_lambda_min << std::endl;
        }
        else if (strcmp(key, "SIG_lambda_max") == 0) {
            // Assign SIG_lambda_max to the double-converted version of value
            SIG_lambda_max = atof(value);
            std::cout << "SIG_lambda_max is " << SIG_lambda_max << std::endl;
        }
        else if (strcmp(key, "BKG_lambda_min") == 0) {
            // Assign BKG_lambda_min to the double-converted version of value
            BKG_lambda_min = atof(value);
            std::cout << "BKG_lambda_min is " << BKG_lambda_min << std::endl;
        }
        else if (strcmp(key, "BKG_lambda_max") == 0) {
            // Assign BKG_lambda_max to the double-converted version of value
            BKG_lambda_max = atof(value);
            std::cout << "BKG_lambda_max is " << BKG_lambda_max << std::endl;
        }
        else if (strcmp(key, "SIG_Emax_over_Ecluster_min") == 0) {
            // Assign SIG_lambda_min to the double-converted version of value
            SIG_Emax_over_Ecluster_min = atof(value);
            std::cout << "SIG_Emax_over_Ecluster_min is " << SIG_Emax_over_Ecluster_min << std::endl;
        }
        else if (strcmp(key, "SIG_Emax_over_Ecluster_max") == 0) {
            // Assign SIG_lambda_max to the double-converted version of value
            SIG_Emax_over_Ecluster_max = atof(value);
            std::cout << "SIG_Emax_over_Ecluster_max is " << SIG_Emax_over_Ecluster_max << std::endl;
        }
        else if (strcmp(key, "BKG_Emax_over_Ecluster_min") == 0) {
            // Assign BKG_lambda_min to the double-converted version of value
            BKG_Emax_over_Ecluster_min = atof(value);
            std::cout << "BKG_Emax_over_Ecluster_min is " << BKG_Emax_over_Ecluster_min << std::endl;
        }
        else if (strcmp(key, "BKG_Emax_over_Ecluster_max") == 0) {
            // Assign BKG_lambda_max to the double-converted version of value
            BKG_Emax_over_Ecluster_max = atof(value);
            std::cout << "BKG_Emax_over_Ecluster_max is " << BKG_Emax_over_Ecluster_max << std::endl;
        }
        else if (strcmp(key, "clus_pT_min") == 0) {
            clus_pT_min = atof(value);
            std::cout << "clus_pT_min is " << clus_pT_min << std::endl;
        }
        else if (strcmp(key, "clus_pT_max") == 0) {
            clus_pT_max = atof(value);
            std::cout << "clus_pT_max is " << clus_pT_max << std::endl;
        }
        else if (strcmp(key, "track_pT_max") == 0) {
            track_pT_max = atof(value);
            std::cout << "track_pT_max is " << track_pT_max << std::endl;
        }
        else if (strcmp(key, "jet_pT_min") == 0) {
            jet_pT_min = atof(value);
            std::cout << "jet_pT_min is " << jet_pT_min << std::endl;
        }
        else if (strcmp(key, "Eta_max") == 0) {
            Eta_max = atof(value);
            std::cout << "Eta_max is " << Eta_max << std::endl;
        }
        else if (strcmp(key, "Cluster_ncell_min") == 0) {
            Cluster_ncell_min = atof(value);
            std::cout << "Cluster_ncell_min is " << Cluster_ncell_min << std::endl;
        }
        else if (strcmp(key, "Cluster_locmaxima_max") == 0) {
            Cluster_locmaxima_max = atof(value);
            std::cout << "Cluster_locmaxima_max is " << Cluster_locmaxima_max << std::endl;
        }
        else if (strcmp(key, "Cluster_distobadchannel") == 0) {
            Cluster_distobadchannel = atof(value);
            std::cout << "Cluster_distobadchannel is " << Cluster_distobadchannel << std::endl;
        }
        else if (strcmp(key, "EcrossoverE_min") == 0) {
            EcrossoverE_min = atof(value);
            std::cout << "EcrossoverE_min is " << EcrossoverE_min << std::endl;
        }
        else if (strcmp(key, "iso_max") == 0) {
            iso_max = atof(value);
            std::cout << "iso_max is " << iso_max << std::endl;
        }
        else if (strcmp(key, "noniso_min") == 0) {
            noniso_min = atof(value);
            std::cout << "noniso_min is " << noniso_min << std::endl;
        }
        else if (strcmp(key, "noniso_max") == 0) {
            noniso_max = atof(value);
            std::cout << "noniso_max is " << noniso_max << std::endl;
        }
        else if (strcmp(key, "deta_max") == 0) {
            deta_max = atof(value);
            std::cout << "deta_max is " << deta_max << std::endl;
        }
        else if (strcmp(key, "phi_func_bins") == 0) {
            phibins = atoi(value);
            std::cout << "Bins in a phi function: " << phibins << std::endl;
        }
        else if (strcmp(key, "eta_func_bins") == 0) {
            etabins = atoi(value);
            std::cout << "Bins in an eta function: " << etabins << std::endl;
        }
        else if (strcmp(key, "xj_func_bins") == 0) {
            xjbins = atoi(value);
            std::cout << "Bins in an xj function: " << xjbins << std::endl;
        }
        else if (strcmp(key, "photon_idvar") == 0) {
            if (strcmp(value, "lambda_0") == 0){
                photon_identifier = LAMBDA_0;
                std::cout << "lambda_0 will determine photon selection" << std::endl;
            }
            else if (strcmp(value, "DNN") == 0){
                photon_identifier = DNN;
                std::cout << "Deep Neural Net will determine photon selection" << std::endl;
            }
            else if (strcmp(value, "Emax_over_Ecluster") == 0){
                photon_identifier = EMAX_OVER_ECLUSTER;
                std::cout << "#frac{E_{max}}{E_{cluster}} will determine photon selection" << std::endl;
            }
            else {
                std::cout << "ERROR: Photon selection determinant in configuration file must be \"lambda_0\", \"DNN\", or \"Emax_over_Ecluster\"" << std::endl << "Aborting the program" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(key, "Cluster_isolation_determinant") == 0) {
            if (strcmp(value, "cluster_iso_tpc_04") == 0){
                determiner = CLUSTER_ISO_TPC_04;
                std::cout << "cluster_iso_tpc_04 will determine the isolation and non-isolation placement" << std::endl;
            }
            else if (strcmp(value, "cluster_iso_its_04") == 0){
                determiner = CLUSTER_ISO_ITS_04;
                std::cout << "cluster_iso_its_04 will determine the isolation and non-isolation placement" << std::endl;
            }
            else if (strcmp(value, "cluster_frixione_tpc_04_02") == 0){
                determiner = CLUSTER_FRIXIONE_TPC_04_02;
                std::cout << "cluster_frixione_tpc_04_02 will determine the isolation and non-isolation placement" << std::endl;
            }
            else if (strcmp(value, "cluster_frixione_its_04_02") == 0){
                determiner = CLUSTER_FRIXIONE_ITS_04_02;
                std::cout << "cluster_frixione_its_04_02 will determine the isolation and non-isolation placement" << std::endl;
            }
            else {
                std::cout << "ERROR: Cluster_isolation_determinant in configuration file must be \"cluster_iso_tpc_04\", \"cluster_iso_its_04\", \"cluster_frixione_tpc_04_02\", or \"cluster_frixione_its_04_02\"" << std::endl << "Aborting the program" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(key, "pdg_code") == 0) {
            rightpdgcode = atoi(value);
            std::cout << "Right pdg_code: " << rightpdgcode << std::endl;
        }
        else if (strcmp(key, "parent_pdg_code") == 0) {
            rightparentpdgcode = atoi(value);
            std::cout << "Right parent_pdg_code: " << rightparentpdgcode << std::endl;
        }
        else if (strcmp(key, "Num_events") == 0) {
            nevents = atoi(value);
            std::cout << "Num_events: " << nevents << std::endl;
        }
        else {
            std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
        }
    }
    fclose(config);
    
    int dummyc = 1;
    char **dummyv = new char *[1];
    
    
    std::cout << " Number of events requested " << nevents << std::endl;
    std::cout << " ptmin " << clus_pT_min << " ptmax = " << clus_pT_max << std::endl;
    
    dummyv[0] = strdup("main");
    TApplication application("", &dummyc, dummyv);
    std::cout <<" Number of arguments " << argc << std::endl;
    
    // Declare histograms
    TH1D hSR_ak04_jetpT("sig_ak04_jetpT", Form("Reco Jet pT from ak04 (%2.2f GeV < p_{T}^{clus} < %2.2f GeV), signal region; p_{T}^{jet} (GeV); #frac{dN}{dp_{T}^{jet}}", clus_pT_min, clus_pT_max), 20, 0, 20);
    TH1D hBR_ak04_jetpT("bkg_ak04_jetpT", Form("Reco Jet pT from ak04 (%2.2f GeV < p_{T}^{clus} < %2.2f GeV), background region; p_{T}^{jet} (GeV); #frac{dN}{dp_{T}^{jet}}", clus_pT_min, clus_pT_max), 20, 0, 20);
    TH1D hSR_ak03_jetpT("sig_ak03_jetpT", Form("Reco Jet pT from ak03 (%2.2f GeV < p_{T}^{clus} < %2.2f GeV), signal region; p_{T}^{jet} (GeV); #frac{dN}{dp_{T}^{jet}}", clus_pT_min, clus_pT_max), 20, 0, 20);
    TH1D hBR_ak03_jetpT("bkg_ak03_jetpT", Form("Reco Jet pT from ak03 (%2.2f GeV < p_{T}^{clus} < %2.2f GeV), background region; p_{T}^{jet} (GeV); #frac{dN}{dp_{T}^{jet}}", clus_pT_min, clus_pT_max), 20, 0, 20);
    
    TH2D hSR_cluspT_ak04_jetpT("sig_cluspT_ak04_jetpT", "Cluster pT and reco Jet pT from ak04, signal region; p_{T}^{#gamma} (GeV); p_{T}^{jet} (GeV)", 30, 0, 30, 20, 0, 20);
    TH2D hBR_cluspT_ak04_jetpT("bkg_cluspT_ak04_jetpT", "Cluster pT and reco Jet pT from ak04, background region; p_{T}^{#gamma} (GeV); p_{T}^{jet} (GeV)", 30, 0, 30, 20, 0, 20);
    TH2D hSR_cluspT_ak03_jetpT("sig_cluspT_ak03_jetpT", "Cluster pT and reco Jet pT from ak03, signal region; p_{T}^{#gamma} (GeV); p_{T}^{jet} (GeV)", 30, 0, 30, 20, 0, 20);
    TH2D hBR_cluspT_ak03_jetpT("bkg_cluspT_ak03_jetpT", "Cluster pT and reco Jet pT from ak03, background region; p_{T}^{#gamma} (GeV); p_{T}^{jet} (GeV)", 30, 0, 30, 20, 0, 20);
    
    hSR_ak04_jetpT.Sumw2();
    hBR_ak04_jetpT.Sumw2();
    hSR_ak03_jetpT.Sumw2();
    hBR_ak03_jetpT.Sumw2();
    
    hSR_cluspT_ak04_jetpT.Sumw2();
    hBR_cluspT_ak04_jetpT.Sumw2();
    hSR_cluspT_ak03_jetpT.Sumw2();
    hBR_cluspT_ak03_jetpT.Sumw2();
    
    // Loop over all files on which the variable is called
    Bool_t isRealData = true;
    for (int iarg = 1; iarg < argc; iarg++) {
        std::string filestring = (std::string)argv[iarg];
        std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
        TFile *file = TFile::Open((TString)argv[iarg]);
        
        if (file == NULL) {
            std::cout << " fail; could not open file" << std::endl;
            exit(EXIT_FAILURE);
        }
        file->Print();
        
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
        
        std::cout <<"_tree_event->GetEntries() " << _tree_event->GetEntries() << std::endl;
        
        // Define the jet variables
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
        
        //Jets reco (ak04)
        UInt_t njet_ak04its;
        Float_t jet_ak04its_pt_raw[NTRACK_MAX];
        Float_t jet_ak04its_eta_raw[NTRACK_MAX];
        Float_t jet_ak04its_phi[NTRACK_MAX];
        
        Float_t jet_ak04its_pt_truth[NTRACK_MAX];
        Float_t jet_ak04its_eta_truth[NTRACK_MAX];
        Float_t jet_ak04its_phi_truth[NTRACK_MAX];
        
        // Jets reco (ak03)
        UInt_t njet_ak03its;
        Float_t jet_ak03its_pt_raw[NTRACK_MAX];
        Float_t jet_ak03its_eta_raw[NTRACK_MAX];
        Float_t jet_ak03its_phi[NTRACK_MAX];
        
        Float_t jet_ak03its_pt_truth[NTRACK_MAX];
        Float_t jet_ak03its_eta_truth[NTRACK_MAX];
        Float_t jet_ak03its_phi_truth[NTRACK_MAX];
        
        //The z_reco is defined as the fraction of the true jet that ended up in this reco jet
        //There are two entries and indices, the first is the best.
        Int_t   jet_ak04its_truth_index_z_reco[NTRACK_MAX][2];
        Float_t jet_ak04its_truth_z_reco[NTRACK_MAX][2];
        Float_t jet_ak04its_ptd_raw[NTRACK_MAX];
        Float_t jet_ak04its_width_sigma[NTRACK_MAX][2];
        UShort_t jet_ak04its_multiplicity[NTRACK_MAX];
        
        //Truth Jets (ak04)
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
        //ak03
        _tree_event->SetBranchAddress("njet_ak03its", &njet_ak03its);
        _tree_event->SetBranchAddress("jet_ak03its_pt_raw", jet_ak03its_pt_raw);
        _tree_event->SetBranchAddress("jet_ak03its_eta_raw", jet_ak03its_eta_raw);
        _tree_event->SetBranchAddress("jet_ak03its_phi", jet_ak03its_phi);
        _tree_event->SetBranchAddress("jet_ak03its_pt_truth", jet_ak03its_pt_truth);
        _tree_event->SetBranchAddress("jet_ak03its_eta_truth", jet_ak03its_eta_truth);
        _tree_event->SetBranchAddress("jet_ak03its_phi_truth", jet_ak03its_phi_truth);
        
        //ak04
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
        
        _tree_event->GetEntry(1);
        if(nmc_truth>0) isRealData= false;
        else isRealData = true;
        
        // Loop over events
        std::cout<<" About to start looping over events to get weights" << std::endl;
        
        if( not(nevents>0)){
            nevents = _tree_event->GetEntries();
        }
        for(Long64_t ievent = 0; ievent < nevents ; ievent++){
            if (ievent % 100000 == 0) std::cout << " event " << ievent << std::endl;
            
            _tree_event->GetEntry(ievent);
            if(not( TMath::Abs(primary_vertex[2])<primary_vertex_max)) continue; //vertex z position
            if(not (primary_vertex[2]!=0.00 )) continue; //removes default of vertex z = 0
            if(is_pileup_from_spd_5_08) continue; //removes pileup
            
            ULong64_t one1 = 1;
            ULong64_t triggerMask_13data = (one1 << 17) | (one1 << 18) | (one1 << 19) | (one1 << 20); //EG1 or EG2 or EJ1 or EJ2
            //if(triggerMask_13data & trigger_mask[0] == 0) continue; //trigger selection
            
            // Determine weights
            double weight = 1.0;
            if(not isRealData){
                if(eg_cross_section>0 and eg_ntrial>0){
                    weight = eg_cross_section/(double)eg_ntrial;
                }
                //17g6a1 weights
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat1_ptmin12.0_Nevent_500000.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 1 of 17g6a1 series opened" << std::endl;
                    weight = 1.60e-11;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat2_ptmin12.0_Nevent_500000.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 2 of 17g6a1 series opened" << std::endl;
                    weight = 2.72e-12;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat3_ptmin12.0_Nevent_500000.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 3 of 17g6a1 series opened" << std::endl;
                    weight = 3.69e-13;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat4_ptmin12.0_Nevent_500000.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 4 of 17g6a1 series opened" << std::endl;
                    weight = 6.14e-14;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat5_ptmin12.0_Nevent_500000.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 5 of 17g6a1 series opened" << std::endl;
                    weight = 1.27e-14;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat1_ptmin12.0_Nevent_300000.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 1 of 17g6a1 series opened" << std::endl;
                    weight = 1.60e-11;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a2_pthat1_ptmin12.0_Nevent_300000.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 2 of 17g6a1 series opened" << std::endl;
                    weight = 2.72e-12;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a3_pthat1_ptmin12.0_Nevent_300000.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 3 of 17g6a1 series opened" << std::endl;
                    weight = 3.69e-13;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a4_pthat1_ptmin12.0_Nevent_300000.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 4 of 17g6a1 series opened" << std::endl;
                    weight = 6.14e-14;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a5_pthat1_ptmin12.0_Nevent_300000.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 5 of 17g6a1 series opened" << std::endl;
                    weight = 1.27e-14;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/17g6a1_pthat1.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 1 of 17g6a1 series opened" << std::endl;
                    weight = 1.60e-11;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/17g6a1_pthat2.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 2 of 17g6a1 series opened" << std::endl;
                    weight = 2.72e-12;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/17g6a1_pthat3.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 3 of 17g6a1 series opened" << std::endl;
                    weight = 3.69e-13;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/17g6a1_pthat4.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 4 of 17g6a1 series opened" << std::endl;
                    weight = 6.14e-14;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/17g6a1_pthat5.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 5 of 17g6a1 series opened" << std::endl;
                    weight = 1.27e-14;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat2_4L_allruns_ptmin15.0.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 2 of 17g6a1 series opened" << std::endl;
                    weight = 2.72e-12;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat3_4L_allruns_ptmin15.0.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 3 of 17g6a1 series opened" << std::endl;
                    weight = 3.69e-13;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat4_4L_allruns_ptmin15.0.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 4 of 17g6a1 series opened" << std::endl;
                    weight = 6.14e-14;
                }
                if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat5_4L_allruns_ptmin15.0.root"){
                    if(ievent == 0)
                        std::cout << "PtHat 5 of 17g6a1 series opened" << std::endl;
                    weight = 1.27e-14;
                }
            }
            
            // Loop over clusters
            for (ULong64_t n = 0; n < ncluster; n++) {
                
                // Cuts
                double isolation;
                //remove UE subtraction
                if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n] + cluster_iso_its_04_ue[n];
                else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n];
                else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n] + cluster_iso_its_04_ue[n];
                else isolation = cluster_frixione_its_04_02[n] + cluster_iso_its_04_ue[n];
                
                isolation = isolation - ue_estimate_its_const*0.4*0.4*TMath::Pi(); //Use rhoxA subtraction
                
                if( not(cluster_ncell[n]>Cluster_ncell_min)) continue;   //removes clusters with 1 or 2 cells
                if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters
                if( not(cluster_nlocal_maxima[n]<= Cluster_locmaxima_max)) continue; //require to have at most 2 local maxima.
                if( not(cluster_distance_to_bad_channel[n]>=Cluster_distobadchannel)) continue;
                if( not(isolation < iso_max)) continue;
                
                // Shower shape cuts
                Bool_t inSignalRegion;
                Bool_t inBkgRegion;
                
                if (photon_identifier == DNN) {
                    inSignalRegion = ((cluster_s_nphoton[n][1] > SIG_DNN_min) and (cluster_s_nphoton[n][1]<SIG_DNN_max));
                    inBkgRegion    = ((cluster_s_nphoton[n][1]>BKG_DNN_min) and (cluster_s_nphoton[n][1]<BKG_DNN_max));
                }
                else if (photon_identifier == LAMBDA_0) {
                    inSignalRegion = ((cluster_lambda_square[n][0]>SIG_lambda_min) and (cluster_lambda_square[n][0]<SIG_lambda_max));
                    inBkgRegion    = ((cluster_lambda_square[n][0]>BKG_lambda_min) and (cluster_lambda_square[n][0]<BKG_lambda_max));
                }
                else {
                    float eratio = cluster_e_max[n]/cluster_e[n];
                    inSignalRegion = (( eratio > SIG_Emax_over_Ecluster_min) and (eratio < SIG_Emax_over_Ecluster_max));
                    inBkgRegion    = ((eratio > BKG_Emax_over_Ecluster_min) and (eratio < BKG_Emax_over_Ecluster_max));
                }
                
                // Loop over jets
                for (ULong64_t ijet = 0; ijet < njet_ak04its; ijet++) {
                    if(not (TMath::Abs(jet_ak04its_eta_raw[ijet]) <Eta_max)) continue; // Jet cut
                    
                    // Plot 2D histograms
                    if (inSignalRegion){
                        hSR_cluspT_ak04_jetpT.Fill(cluster_pt[n], jet_ak04its_pt_raw[ijet], weight);
                        hSR_cluspT_ak03_jetpT.Fill(cluster_pt[n], jet_ak03its_pt_raw[ijet], weight);
                    }
                    else if(inBkgRegion) {
                        hBR_cluspT_ak04_jetpT.Fill(cluster_pt[n], jet_ak04its_pt_raw[ijet], weight);
                        hBR_cluspT_ak03_jetpT.Fill(cluster_pt[n], jet_ak03its_pt_raw[ijet], weight);
                    }
                    
                    if( not(cluster_pt[n]>clus_pT_min)) continue; //select pt of photons
                    if( not(cluster_pt[n]<clus_pT_max)) continue;
                    
                    // Plot 1D histograms
                    if (inSignalRegion){
                        hSR_ak04_jetpT.Fill(jet_ak04its_pt_raw[ijet], weight);
                        hSR_ak03_jetpT.Fill(jet_ak03its_pt_raw[ijet], weight);
                    }
                    else if(inBkgRegion) {
                        hBR_ak04_jetpT.Fill(jet_ak04its_pt_raw[ijet], weight);
                        hBR_ak03_jetpT.Fill(jet_ak03its_pt_raw[ijet], weight);
                    }
                } // End loop over jets
            } // End loop over clusters
        } // End loop over events
    } // End loop over files
    
    std::string opened_files = "";
    for (int iarg = 1; iarg < argc; iarg++) {
        std::string filepath = argv[iarg];
        opened_files += "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
    }
    
    std::string photonselectionvar;
    if (photon_identifier == DNN) {
        photonselectionvar = "DNN";
    }
    else if (photon_identifier == LAMBDA_0) {
        photonselectionvar = "Lambda0";
    }
    else {
        photonselectionvar = "EmaxOverEcluster";
    }
    TFile* fout = new TFile(Form("GammaJet_jetpT_charts_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_%s_PHOTONSELECT_%s.root", clus_pT_min, clus_pT_max, jet_pT_min, opened_files.c_str(), photonselectionvar.c_str()),"RECREATE");
    //TFile* fout = new TFile(Form("GammaJet_config_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_MC17g6a1_PHOTONSELECT_%s.root", clus_pT_min, clus_pT_max, jet_pT_min, photonselectionvar.c_str()),"RECREATE");
    //TFile* fout = new TFile(Form("GammaJet_config_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_MCdijet_PHOTONSELECT_%s.root", clus_pT_min, clus_pT_max, jet_pT_min, photonselectionvar.c_str()),"RECREATE");
    //TFile* fout = new TFile(Form("GammaJet_config_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_MCgammajet_PHOTONSELECT_%s.root", clus_pT_min, clus_pT_max, jet_pT_min, photonselectionvar.c_str()),"RECREATE");
    fout->Print();
    
    // Load the plot to ROOT file
    hSR_ak04_jetpT.Write("sig_ak04_jetpT");
    hBR_ak04_jetpT.Write("bkg_ak04_jetpT");
    hSR_ak03_jetpT.Write("sig_ak03_jetpT");
    hBR_ak03_jetpT.Write("bkg_ak03_jetpT");
    
    hSR_cluspT_ak04_jetpT.Write("sig_cluspT_ak04_jetpT");
    hBR_cluspT_ak04_jetpT.Write("bkg_cluspT_ak04_jetpT");
    hSR_cluspT_ak03_jetpT.Write("sig_cluspT_ak03_jetpT");
    hBR_cluspT_ak03_jetpT.Write("bkg_cluspT_ak03_jetpT");
    
    // Graph the plots for ak04
    TCanvas canvas1("canvas1", "");
    canvas1.SetLogy();
    hSR_ak04_jetpT.Draw();
    canvas1.SaveAs(Form("hSR_ak04_jetpT_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_%s_PHOTONSELECT_%s.pdf", clus_pT_min, clus_pT_max, jet_pT_min, opened_files.c_str(), photonselectionvar.c_str()));
    canvas1.Clear();
    hBR_ak04_jetpT.Draw();
    canvas1.SaveAs(Form("hBR_ak04_jetpT_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_%s_PHOTONSELECT_%s.pdf", clus_pT_min, clus_pT_max, jet_pT_min, opened_files.c_str(), photonselectionvar.c_str()));
    canvas1.Clear();
    
    TCanvas canvas2("canvas2", "");
    canvas2.SetLogz();
    hSR_cluspT_ak04_jetpT.SetStats(0);
    hSR_cluspT_ak04_jetpT.Draw("COLZ");
    canvas2.SaveAs(Form("hSR_cluspT_ak04_jetpT_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_%s_PHOTONSELECT_%s.pdf", clus_pT_min, clus_pT_max, jet_pT_min, opened_files.c_str(), photonselectionvar.c_str()));
    canvas2.Clear();
    hBR_cluspT_ak04_jetpT.SetStats(0);
    hBR_cluspT_ak04_jetpT.Draw("COLZ");
    canvas2.SaveAs(Form("hBR_cluspT_ak04_jetpT_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_%s_PHOTONSELECT_%s.pdf", clus_pT_min, clus_pT_max, jet_pT_min, opened_files.c_str(), photonselectionvar.c_str()));
    canvas2.Clear();
    
    
    canvas1.Close();
    canvas2.Close();
    std::cout << " ending " << std::endl;
    fout->Close();
    //end of arguments
    return EXIT_SUCCESS;
}
