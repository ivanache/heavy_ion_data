/*
  This code is a variation of Fernando's Skeleton_Mix_Correlations.cc code, this time for events
*/
// Author: Ivan Chernyshev; Creator of template code: Fernando Torales-Acosta

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1D.h>

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

// Energy of lead, in GeV
const double EPb = 1560;

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

using namespace H5;

double calculatebinwidth(int numofbins, double binmin, double binmax){
    return (binmax - binmin)/numofbins;
}

int main(int argc, char *argv[])
{
    if (argc < 9) {
        fprintf(stderr,"Batch Syntax is [Gamma-Triggered Paired Root], [Min-Bias HDF5] [Mix Start] [Mix End] [Track Skim GeV] [Min cluster pT] [Max cluster pT] [Min jet pT]");
        exit(EXIT_FAILURE);
    }
    
    int dummyc = 1;
    char **dummyv = new char *[1];
    
    dummyv[0] = strdup("main");
    
    TString root_file = (TString)argv[1];
    std::cout << "Opening: " << (TString)argv[1] << std::endl;
    
    const H5std_string hdf5_file_name(argv[2]);
    TString hdf5_file = (TString)argv[2];
    fprintf(stderr,hdf5_file);
    
    size_t mix_start = atoi(argv[3]);
    size_t mix_end = atoi(argv[4]);
    
    int GeV_Track_Skim = atoi(argv[5]);
    std::cout<<"mix start is "<<mix_start<<std::endl;
    std::cout<<"mix end is "<<mix_end<<std::endl;
    fprintf(stderr,"Using %iGeV Track Skimmed from batch Script \n",GeV_Track_Skim);
    
    double cluspTmin = atof(argv[6]);
    double cluspTmax = atof(argv[7]);
    std::cout<< "Cluster pT min: " << cluspTmin << "; max: " << cluspTmax << std::endl;
    
    double jetpTmin = atof(argv[8]);
    std::cout << "Minimum jet pT: " << jetpTmin << std::endl;
    
    size_t nmix = 300;
    fprintf(stderr,"Number of Mixed Events: %i \n",nmix);
    
    // Declare histograms
    TH1D* SIGcluster_pt_dist = new TH1D("sig_Cluster_pT", "Signal Cluster p_{T} distribution; cluster p_{T} (GeV); #frac{dN}{N_{#gamma}*N_{minbias}}", 5, cluspTmin, cluspTmax);
    TH1D* SIGjet_pt_dist = new TH1D("sig_Jet_pT", "Signal Jet p_{T} distribution; jet p_{T} (GeV); #frac{dN}{N_{#gamma}*N_{minbias}}", 25, 5, 30);
    TH1D* SIGpt_diff_dist = new TH1D("sig_clusjet_pT_diff", "Signal p_{T}^{cluster}-p_{T}^{jet} distribution; #Delta p_{T} (GeV); #frac{dN}{N_{#gamma}*N_{minbias}}", 20, 0, 20);
    double SIGcluster_pt_dist_binwidth = calculatebinwidth(5, cluspTmin, cluspTmax);
    double SIGjet_pt_dist_binwidth = calculatebinwidth(25, 5, 30);
    double SIGpt_diff_dist_binwidth = calculatebinwidth(25, 0, 20);
    SIGcluster_pt_dist->Sumw2();
    SIGjet_pt_dist->Sumw2();
    SIGpt_diff_dist->Sumw2();
    
    TH1D* SIGdPhi = new TH1D("sig_dPhi", "Signal #Delta #phi distribution; #Delta #phi (rads); #frac{dN}{N_{#gamma}*N_{minbias}}", 7, 0, TMath::Pi());
    TH1D* SIGclusterPhi = new TH1D("sig_clusterPhi", "Signal #phi_{cluster} distribution; #phi (rads); #frac{dN}{N_{#gamma}*N_{minbias}}", 14, -TMath::Pi(), TMath::Pi());
    TH1D* SIGjetPhi = new TH1D("sig_jetPhi", "Signal #phi_{jet} distribution; #phi (rads); #frac{dN}{N_{#gamma}*N_{minbias}}", 14, -TMath::Pi(), TMath::Pi());
    double SIGdPhi_binwidth = calculatebinwidth(7, 0, TMath::Pi());
    double SIGclusterPhi_binwidth = calculatebinwidth(14, -TMath::Pi(), TMath::Pi());
    double SIGjetPhi_binwidth = calculatebinwidth(14, -TMath::Pi(), TMath::Pi());
    SIGdPhi->Sumw2();
    SIGclusterPhi->Sumw2();
    SIGjetPhi->Sumw2();
    
    TH1D* SIGdEta = new TH1D("sig_dEta", "Signal #Delta #eta distribution; #Delta #eta; #frac{dN}{N_{#gamma}*N_{minbias}}", 40, -2.4, 2.4);
    TH1D* SIGclusterEta = new TH1D("sig_clusterEta", "#Signal eta_{cluster} distribution; #eta; #frac{dN}{N_{#gamma}*N_{minbias}}", 20, -1.2, 1.2);
    TH1D* SIGjetEta = new TH1D("sig_jetEta", "Signal #eta_{jet} distribution; #eta; #frac{dN}{N_{#gamma}*N_{minbias}}", 20, -1.2, 1.2);
    double SIGdEta_binwidth = calculatebinwidth(40, -2.4, 2.4);
    double SIGclusterEta_binwidth = calculatebinwidth(20, -1.2, 1.2);
    double SIGjetEta_binwidth = calculatebinwidth(20, -1.2, 1.2);
    SIGdEta->Sumw2();
    SIGclusterEta->Sumw2();
    SIGjetEta->Sumw2();
    
    TH1D* SIGXj = new TH1D("sig_Xj", "Signal Xj distribution; Xj; #frac{dN}{N_{#gamma}*N_{minbias}}", 10, 0.0,2.0);
    TH1D* SIGpTD = new TH1D("sig_pTD", "Signal Jet pTD distribution; p_{T}D (GeV); #frac{dN}{N_{#gamma}*N_{minbias}}", 5, 0.0,1.0);
    TH1D* SIGMultiplicity = new TH1D("sig_Multiplicity", "Signal Jet Multiplicity; Multiplicity; #frac{dN}{N_{#gamma}*N_{minbias}}", 10, 0.0 , 20.0);
    TH1D* SIGXobsPb = new TH1D("sig_XobsPb", "x_{pPb}^{obs} distribution: signal region; x_{pPb}^{obs}; #frac{d #sigma}{dx^{obs}_{pPb}}", 5, 0.004, 0.024);
    double SIGXj_binwidth = calculatebinwidth(10, 0.0,2.0);
    double SIGpTD_binwidth = calculatebinwidth(5, 0.0,1.0);
    double SIGMultiplicity_binwidth = calculatebinwidth(10, 0.0 , 20.0);
    double SIGXobsPb_binwidth = calculatebinwidth(5, 0.004, 0.024);
    SIGXj->Sumw2();
    SIGpTD->Sumw2();
    SIGMultiplicity->Sumw2();
    SIGXobsPb->Sumw2();
    
    TH1D* BKGcluster_pt_dist = new TH1D("bkg_Cluster_pT", "Background Cluster p_{T} distribution; cluster p_{T} (GeV); #frac{dN}{N_{#gamma}*N_{minbias}}", 7, cluspTmin, cluspTmax);
    TH1D* BKGjet_pt_dist = new TH1D("bkg_Jet_pT", "Background Jet p_{T} distribution; jet p_{T} (GeV); #frac{dN}{N_{#gamma}*N_{minbias}}", 25, 5, 30);
    TH1D* BKGpt_diff_dist = new TH1D("bkg_clusjet_pT_diff", "Background p_{T}^{cluster}-p_{T}^{jet} distribution; #Delta p_{T} (GeV); #frac{dN}{N_{#gamma}*N_{minbias}}", 20, 0, 20);
    double BKGcluster_pt_dist_binwidth = calculatebinwidth(7, cluspTmin, cluspTmax);
    double BKGjet_pt_dist_binwidth = calculatebinwidth(25, 5, 30);
    double BKGpt_diff_dist_binwidth = calculatebinwidth(20, 0, 20);
    BKGcluster_pt_dist->Sumw2();
    BKGjet_pt_dist->Sumw2();
    BKGpt_diff_dist->Sumw2();
    
    TH1D* BKGdPhi = new TH1D("bkg_dPhi", "Background #Delta #phi distribution; #Delta #phi (rads); #frac{dN}{N_{#gamma}*N_{minbias}}", 7, 0, TMath::Pi());
    TH1D* BKGclusterPhi = new TH1D("bkg_clusterPhi", "Background #phi_{cluster} distribution; #phi (rads); #frac{dN}{N_{#gamma}*N_{minbias}}", 14, -TMath::Pi(), TMath::Pi());
    TH1D* BKGjetPhi = new TH1D("bkg_jetPhi", "Background #phi_{jet} distribution; #phi (rads); #frac{dN}{N_{#gamma}*N_{minbias}}", 14, -TMath::Pi(), TMath::Pi());
    double BKGdPhi_binwidth = calculatebinwidth(7, 0, TMath::Pi());
    double BKGclusterPhi_binwidth = calculatebinwidth(14, -TMath::Pi(), TMath::Pi());
    double BKGjetPhi_binwidth = calculatebinwidth(14, -TMath::Pi(), TMath::Pi());
    BKGdPhi->Sumw2();
    BKGclusterPhi->Sumw2();
    BKGjetPhi->Sumw2();
    
    TH1D* BKGdEta = new TH1D("bkg_dEta", "Background #Delta #eta distribution; #Delta #eta; #frac{dN}{N_{#gamma}*N_{minbias}}", 40, -2.4, 2.4);
    TH1D* BKGclusterEta = new TH1D("bkg_clusterEta", "Background #eta_{cluster} distribution; #eta; #frac{dN}{N_{#gamma}*N_{minbias}}", 20, -1.2, 1.2);
    TH1D* BKGjetEta = new TH1D("bkg_jetEta", "Background #eta_{jet} distribution; #eta; #frac{dN}{N_{#gamma}*N_{minbias}}", 20, -1.2, 1.2);
    double BKGdEta_binwidth = calculatebinwidth(40, -2.4, 2.4);
    double BKGclusterEta_binwidth = calculatebinwidth(20, -1.2, 1.2);
    double BKGjetEta_binwidth = calculatebinwidth(20, -1.2, 1.2);
    BKGdEta->Sumw2();
    BKGclusterEta->Sumw2();
    BKGjetEta->Sumw2();
    
    TH1D* BKGXj = new TH1D("bkg_Xj", "Background Xj distribution; Xj; #frac{dN}{N_{#gamma}*N_{minbias}}", 10, 0.0,2.0);
    TH1D* BKGpTD = new TH1D("bkg_pTD", "Background Jet pTD distribution; p_{T}D (GeV); #frac{dN}{N_{#gamma}*N_{minbias}}", 5, 0.0,1.0);
    TH1D* BKGMultiplicity = new TH1D("bkg_Multiplicity", "Background Jet Multiplicity distribution; Multiplicity; #frac{dN}{N_{#gamma}*N_{minbias}}", 10, 0.0 , 20.0);
    TH1D* BKGXobsPb = new TH1D("bkg_XobsPb", "x_{pPb}^{obs} distribution: background region; x_{pPb}^{obs}; #frac{d #sigma}{dx^{obs}_{pPb}}", 5, 0.004, 0.024);
    double BKGXj_binwidth = calculatebinwidth(10, 0.0,2.0);
    double BKGpTD_binwidth = calculatebinwidth(5, 0.0,1.0);
    double BKGMultiplicity_binwidth = calculatebinwidth(10, 0.0 , 20.0);
    double BKGXobsPb_binwidth = calculatebinwidth(5, 0.004, 0.024);
    BKGXj->Sumw2();
    BKGpTD->Sumw2();
    BKGMultiplicity->Sumw2();
    BKGXobsPb->Sumw2();
    
    TH1D* z_Vertices_individual = new TH1D("z_Vertices_individual", "Z-vertex (ROOT)", 50, 0, 25);
    TH1D* z_Vertices_hdf5 = new TH1D("z_Vertices_hdf5", "Z-vertex (hdf5)", 50, 0, 25);
    TH1D* z_Vertices = new TH1D("z_Vertices", "Z-vertex difference distribution", 50, 0, 25);
    double z_Vertices_individual_binwidth = calculatebinwidth(50, 0, 25);
    double z_Vertices_hdf5_binwidth = calculatebinwidth(50, 0, 25);
    double z_Vertices_binwidth = calculatebinwidth(50, 0, 25);
    z_Vertices_individual->Sumw2();
    z_Vertices_hdf5->Sumw2();
    z_Vertices->Sumw2();
    
    TH1D* Multiplicity_individual = new TH1D("mult_Vertices_individual", "Multiplicity (ROOT)", 427, 0, 1281);
    TH1D* Multiplicity_hdf5 = new TH1D("mult_Vertices_hdf5", "Multiplicity (hdf5)", 427, 0, 1281);
    TH1D* Multiplicity = new TH1D("mult_Vertices", "Multiplicity differnce distribution", 1281, 0, 1281);
    double Multiplicity_individual_binwidth = calculatebinwidth(50, 0, 25);
    double Multiplicity_hdf5_binwidth = calculatebinwidth(50, 0, 25);
    double Multiplicity_binwidth = calculatebinwidth(50, 0, 25);
    Multiplicity_individual->Sumw2();
    Multiplicity_hdf5->Sumw2();
    Multiplicity->Sumw2();
    
    //Config File ---------------------------------------------------------------------------
    
    //Declaration and Initialize Variables TO BE SET BY [Corr_config.yaml]
    FILE* config = fopen("Corr_config.yaml", "r");
    double DNN_min = 0;
    double DNN_max = 0;
    double pT_min = 0;
    double pT_max = 0;
    double Eta_max = 0;
    double Cluster_min = 0;
    double EcrossoverE_min = 0;
    int Track_Cut_Bit = 0;
    double iso_max = 0;
    double noniso_min = 0;
    double noniso_max = 0;
    double deta_max = 0;
    isolationDet determiner = CLUSTER_ISO_ITS_04;
    int n_eta_bins = 0;
    int n_phi_bins = 0;
    
    // zT & pT bins
    int nztbins = 7;
    float* ztbins;
    ztbins = new float[nztbins+1];
    ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
    
    int nptbins = 3;
    float* ptbins;
    ptbins = new float[nptbins+1];
    ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;
    
    
    //READ CONFIG
    char line[MAX_INPUT_LENGTH];
    while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
        if (line[0] == '#') continue;
        
        char key[MAX_INPUT_LENGTH];
        char dummy[MAX_INPUT_LENGTH];
        char value[MAX_INPUT_LENGTH];
        
        // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays
        key[0] = '\0';
        value[0] = '\0';
        sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);
        
        //Read Config File: Detect Keys
        if (strcmp(key, "DNN_min") == 0) {
            DNN_min = atof(value);
            std::cout << "DNN_min: " << DNN_min << std::endl; }
        
        else if (strcmp(key, "DNN_max") == 0) {
            DNN_max = atof(value);
            std::cout << "DNN_max: " << DNN_max << std::endl; }
        
        else if (strcmp(key, "pT_min") == 0) {
            pT_min = atof(value);
            std::cout << "pT_min: " << pT_min << std::endl; }
        
        else if (strcmp(key, "pT_max") == 0) {
            pT_max = atof(value);
            std::cout << "pT_max: " << pT_max << std::endl; }
        
        else if (strcmp(key, "Eta_max") == 0) {
            Eta_max = atof(value);
            std::cout << "Eta_max: " << Eta_max << std::endl;
        }
        else if (strcmp(key, "Cluster_min") == 0) {
            Cluster_min = atof(value);
            std::cout << "Cluster_min: " << Cluster_min << std::endl; }
        
        else if (strcmp(key, "EcrossoverE_min") == 0) {
            EcrossoverE_min = atof(value);
            std::cout << "EcrossoverE_min; " << EcrossoverE_min << std::endl; }
        
        else if (strcmp(key, "iso_max") == 0) {
            iso_max = atof(value);
            std::cout << "iso_max: " << iso_max << std::endl; }
        
        else if (strcmp(key, "noniso_min") == 0) {
            noniso_min = atof(value);
            std::cout << "noniso_min: " << noniso_min << std::endl; }
        
        else if (strcmp(key, "noniso_max") == 0) {
            noniso_max = atof(value);
            std::cout << "noniso_max: " << noniso_max << std::endl; }
        
        else if (strcmp(key, "deta_max") == 0) {
            deta_max = atof(value);
            std::cout << "deta_max: " << deta_max << std::endl; }
        
        else if (strcmp(key, "N_Phi_Bins") == 0) {
            n_phi_bins = atoi(value);
            std::cout << "Number of Phi Bins: " << n_phi_bins << std::endl; }
        
        else if (strcmp(key, "N_Eta_Bins") == 0) {
            n_eta_bins = atoi(value);
            std::cout << "Number of Eta Bins: " << n_eta_bins << std::endl; }
        
        else if (strcmp(key, "Track_Cut_Bit") == 0) {
            Track_Cut_Bit = atoi(value);
            std::cout << "Track Cut Bit: " << Track_Cut_Bit << std::endl; }
        
        else if (strcmp(key, "Zt_bins") == 0) {
            nztbins = -1;
            for (const char *v = value; *v != ']';) {
                while (*v != ']' && !isdigit(*v)) v++;
                nztbins++;
                while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
            
            ztbins = new float[nztbins + 1];
            int i = 0;
            for (const char *v = value; *v != ']' ;) {
                while (*v != ']' && !isdigit(*v)) v++;
                ztbins[i] = atof(v);
                i++;
                while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
            
            std::cout << "Number of Zt bins: " << nztbins << std::endl << "Zt bins: {";
            for (int i = 0; i <= nztbins; i++)
                std::cout << ztbins[i] << ", ";
            std::cout << "}\n";
        }
        
        else if (strcmp(key, "Pt_bins") == 0) {
            nptbins = -1;
            for (const char *v = value; *v != ']';) {
                while (*v != ']' && !isdigit(*v)) v++;
                nptbins++;
                while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
            
            ptbins = new float[nptbins + 1];
            int i = 0;
            for (const char *v = value; *v != ']' ;) {
                while (*v != ']' && !isdigit(*v))  v++;
                ptbins[i] = atof(v);
                i++;
                while (*v != ']' && (isdigit(*v) || *v == '.')) v++; }
            
            std::cout << "Number of Pt bins: " << nptbins << std::endl << "Pt bins: {";
            for (int i = 0; i <= nptbins; i++)
                std::cout << ptbins[i] << ", ";
            std::cout << "}\n";
        }
        
        else if (strcmp(key, "Cluster_isolation_determinant") == 0) {
            if (strcmp(value, "cluster_iso_tpc_04") == 0){
                determiner = CLUSTER_ISO_TPC_04;
                std::cout << "Isolation Variable: cluster_iso_tpc_04" << std::endl; }
            
            else if (strcmp(value, "cluster_iso_its_04") == 0){
                determiner = CLUSTER_ISO_ITS_04;
                std::cout << "Isolation Variable: cluster_iso_its_04" << std::endl; }
            
            else if (strcmp(value, "cluster_frixione_tpc_04_02") == 0){
                determiner = CLUSTER_FRIXIONE_TPC_04_02;
                std::cout << "Isolation Variable: cluster_frixione_tpc_04_02" << std::endl; }
            
            else if (strcmp(value, "cluster_frixione_its_04_02") == 0){
                determiner = CLUSTER_FRIXIONE_ITS_04_02;
                std::cout << "Isolation Variable: cluster_frixione_its_04_02" << std::endl; }
            
            else {
                std::cout << "ERROR: Cluster_isolation_determinant in configuration file must be \"cluster_iso_tpc_04\", \"cluster_iso_its_04\", \"cluster_frixione_tpc_04_02\", or \"cluster_frixione_its_04_02\"" << std::endl << "Aborting the program" << std::endl;
                exit(EXIT_FAILURE); }
        }
        
        else std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
        
    }
    //end Config Loop
    
    fclose(config);
    
    for (int i = 0; i <= nztbins; i++)
        std::cout << "zt bound: " << ztbins[i] << std::endl;
    for (int i = 0; i <= nptbins; i++)
        std::cout << "pt bound: " << ptbins[i] << std::endl;
    
    
    //HISTOGRAMS
    TCanvas canvas("canvas", "");
    
    //TH2D* Signal_pT_Dist = new TH2D("Signal_pT_Dist","Cluster Pt Spectrum For Isolation (its_04) bins 0.55 < DNN < 0.85",59,0.5,30,59,0.5,30);
    //For this example, fill with pt and energy
    
    //ROOT --------------------------------------------------------------------------------------
    
    TFile *file = TFile::Open(root_file);
    
    if (file == NULL) {
        std::cout << " fail" << std::endl;
        exit(EXIT_FAILURE);
    }
    file->Print();
    
    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
    if (_tree_event == NULL) {
        _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
        if (_tree_event == NULL) {
            std::cout << " tree fail " << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    
    //variables
    UInt_t nevent;
    std::vector<Double_t> primary_vertex(3, NAN);
    std::vector<Float_t> multiplicity_v0(64, NAN);//64 channels for v0 detector, to be summed
    
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
    
    Long64_t mix_events[300];
    
    _tree_event->SetBranchAddress("primary_vertex", &primary_vertex[0]);
    _tree_event->SetBranchAddress("multiplicity_v0", &multiplicity_v0[0]);
    
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
    
    //jets
    _tree_event->SetBranchAddress("njet_ak04its", &njet_ak04its);
    _tree_event->SetBranchAddress("jet_ak04its_pt_raw", jet_ak04its_pt_raw);
    _tree_event->SetBranchAddress("jet_ak04its_eta_raw", jet_ak04its_eta_raw);
    _tree_event->SetBranchAddress("jet_ak04its_phi", jet_ak04its_phi);
    
    _tree_event->SetBranchAddress("mixed_events", mix_events);
    
    std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;
    
    
    //Using low level hdf5 API -------------------------------------------------------------------------------
    
    //open hdf5: Define size of data from file, explicitly allocate memory in hdf5 space and array size
    const H5std_string event_ds_name( "event" );
    H5File h5_file( hdf5_file_name, H5F_ACC_RDONLY ); //hdf5_file_name from argv[2]
    DataSet event_dataset = h5_file.openDataSet( event_ds_name );
    DataSpace event_dataspace = event_dataset.getSpace();
    //Load the dimensions of dataset from file, to be used in array/hyperslab
    const int event_ndims = event_dataspace.getSimpleExtentNdims();
    hsize_t event_maxdims[event_ndims];
    hsize_t eventdims[event_ndims];
    event_dataspace.getSimpleExtentDims(eventdims, event_maxdims);
    //UInt_t nevent_max = eventdims[1];
    UInt_t NEvent_Vars = eventdims[1];
    fprintf(stderr, "\n%s:%d: n track variables\n", __FILE__, __LINE__, NEvent_Vars);
    
    //Define array hyperslab will be fed into
    float event_data_out[1][NEvent_Vars];
    
    //Define hyperslab size and offset in  FILE;
    hsize_t event_offset[2] = {0, 0};
    hsize_t event_count[2] = {1, NEvent_Vars};
    
    /*
     The Offset is how we iterate over the entire hdf5 file.
     For example, To obtain data for event 68, set the
     offset's to {68, ntrack_max, NTrack_Vars}.
     */
    
    
    event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");
    
    //Define the memory dataspace in which to place hyperslab
    const int RANK_OUT = 2; //# of Dimensions
    DataSpace event_memspace( RANK_OUT, eventdims );
    
    //Define memory offset for hypreslab starting at begining:
    hsize_t event_offset_out[2] = {0};
    
    //define Dimensions of array, for writing slab to array
    hsize_t event_count_out[2] = {1, NEvent_Vars};
    
    //define space in memory for hyperslab, then write from file to memory
    event_memspace.selectHyperslab( H5S_SELECT_SET, event_count_out, event_offset_out );
    std::cout << "Made it to line 394" <<  std::endl;
    event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "event dataset read into array: OK");
    
    // Get jet
    const H5std_string jet_ds_name( "jet" );
    H5File h5_file_jet( hdf5_file_name, H5F_ACC_RDONLY ); //hdf5_file_name from argv[2]
    DataSet jet_dataset = h5_file_jet.openDataSet( jet_ds_name );
    DataSpace jet_dataspace = jet_dataset.getSpace();
    
    //Load the dimensions of dataset from file, to be used in array/hyperslab
    const int jet_ndims = jet_dataspace.getSimpleExtentNdims();
    hsize_t jet_maxdims[jet_ndims];
    hsize_t jetdims[jet_ndims];
    jet_dataspace.getSimpleExtentDims(jetdims, jet_maxdims);
    
    UInt_t nEvents = jetdims[0];
    UInt_t njet_max = jetdims[1];
    UInt_t Njet_Vars = jetdims[2];
    fprintf(stderr, "\n%s:%d: n jet variables\n", __FILE__, __LINE__, Njet_Vars);
    
    //Define array hyperslab will be fed into
    float jet_data_out[1][njet_max][Njet_Vars];
    
    //Define hyperslab size and offset in  FILE;
    hsize_t jet_offset[3] = {0, 0, 0};
    hsize_t jet_count[3] = {1, njet_max, Njet_Vars};
    
    /*
     The Offset is how we iterate over the entire hdf5 file.
     For example, To obtain data for jet 68, set the
     offset's to {68, njet_max, Njet_Vars}.
     */
    
    
    jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");
    
    //Define the memory dataspace in which to place hyperslab
    const int RANK_OUT_jet = 3; //# of Dimensions
    DataSpace jet_memspace( RANK_OUT_jet, jetdims );
    
    //Define memory offset for hypreslab starting at begining:
    hsize_t jet_offset_out[3] = {0};
    
    //define Dimensions of array, for writing slab to array
    hsize_t jet_count_out[3] = {1, njet_max, Njet_Vars};
    std::cout<< "Made it to line 398" << std::endl;
    
    //define space in memory for hyperslab, then write from file to memory
    jet_memspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
    jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "jet dataset read into array: OK");
    
    //MONEY MAKING LOOP
    Long64_t nentries = _tree_event->GetEntries();
    
    int N_SR = 0;
    int N_BR = 0;
    
    for(Long64_t ievent = 0; ievent < nentries ; ievent++){
        _tree_event->GetEntry(ievent);

        //Cuts/Variables from the ROOT file go here
        for (Long64_t imix = mix_start; imix < mix_end+1; imix++){
            Long64_t mix_event = mix_events[imix];
            //fprintf(stderr,"\n %s:%d: Mixed event = %lu",__FILE__,__LINE__,mix_event);
            
            //if (mix_event == ievent) continue; //not needed for gamma-MB pairing: Different Triggers
            if(mix_event >= 9999999) continue;
            event_dataspace.selectHyperslab( H5S_SELECT_SET, event_count, event_offset );
            event_dataset.read( event_data_out, PredType::NATIVE_FLOAT, event_memspace, event_dataspace );
            
            //adjust offset for next mixed event
            jet_offset[0]=mix_event;
            jet_dataspace.selectHyperslab( H5S_SELECT_SET, jet_count, jet_offset );
            jet_dataset.read( jet_data_out, PredType::NATIVE_FLOAT, jet_memspace, jet_dataspace );
            
            double cluspT = -9000;
            double jet_pT = -9000;
            double clusphi = -9000;
            double jet_phi = -9000;
            double cluseta = -9000;
            double jet_eta = -9000;
            double jet_pTD = -9000;
            double jet_multiplicity = -9000;
            
            for(Long64_t icluster = 0; icluster < ncluster; icluster++) {
                if(not(cluster_pt[icluster] > cluspTmin)) {continue;}
                if(not(cluster_pt[icluster] < cluspTmax)) {continue;}
                if( not(cluster_ncell[icluster]>2)) {continue;}   //removes clusters with 1 or 2 cells
                if( not(cluster_e_cross[icluster]/cluster_e[icluster]>0.03)) {continue;} //removes "spiky" clusters
                if( not(cluster_nlocal_maxima[icluster]<= 2)) {continue;} //require to have at most 2 local maxima.
                if( not(cluster_distance_to_bad_channel[icluster]>=2.0)) {continue;}
                
                cluspT = cluster_pt[icluster];
                clusphi = cluster_phi[icluster];
                cluseta = cluster_eta[icluster];
                
                while(clusphi >= TMath::Pi()) clusphi -= (2*TMath::Pi());
                while(clusphi <= -TMath::Pi()) clusphi += (2*TMath::Pi());
                
               // After cluster cuts, increment number of triggers and loop over jets
                if((cluster_lambda_square[icluster][0] > 0.05) && (cluster_lambda_square[icluster][0] < 0.3)) {
                    N_SR++;
                }
                if((cluster_lambda_square[icluster][0] > 0.4) && (cluster_lambda_square[icluster][0] < 1.0)) {
                    N_BR++;
                }
                for(Long64_t ijet = 0; ijet < njet_ak04its; ijet++){
                    if(TMath::IsNaN(jet_data_out[0][ijet][0])) continue;
                    if(not(jet_data_out[0][ijet][0] > jetpTmin)) {continue;}
                    std::cout << "icluster " << icluster << " ijet " << ijet << " has jet pT " << jet_data_out[0][ijet][0] << std::endl;
                    // After the jet cuts, fill histograms
                    jet_pT = jet_data_out[0][ijet][0];
                    jet_phi = jet_data_out[0][ijet][2];
                    jet_eta = jet_data_out[0][ijet][1];
                    jet_pTD = jet_data_out[0][ijet][3];
                    jet_multiplicity = jet_data_out[0][ijet][4];
                    if (not(TMath::Abs(jet_eta) < 0.5)) {continue;}
                    
                    std::cout << "icluster " << icluster << " ijet " << ijet << " passed all jet cuts " << std::endl;
                    
                    while(jet_phi >= TMath::Pi()) jet_phi -= (2*TMath::Pi());
                    while(jet_phi <= -TMath::Pi()) jet_phi += (2*TMath::Pi());
                    
                    if((cluster_lambda_square[icluster][0] > 0.05) && (cluster_lambda_square[icluster][0] < 0.3)) {
                        SIGcluster_pt_dist->Fill(cluspT);
                        SIGjet_pt_dist->Fill(jet_pT);
                        SIGpt_diff_dist->Fill(TMath::Abs(cluspT-jet_pT));
                        
                        double dphinum = clusphi-jet_phi;
                        while(dphinum < -TMath::Pi()) dphinum += 2*TMath::Pi();
                        while(dphinum > TMath::Pi()) dphinum -= 2*TMath::Pi();
                        SIGdPhi->Fill(TMath::Abs(dphinum));
                        if(not(dphinum > TMath::Pi()/2)) continue;
                        SIGclusterPhi->Fill(clusphi);
                        SIGjetPhi->Fill(jet_phi);
                        
                        SIGdEta->Fill(jet_eta-cluseta);
                        SIGclusterEta->Fill(cluseta);
                        SIGjetEta->Fill(jet_eta);
                        
                        SIGXj->Fill(jet_pT/cluspT);
                        SIGpTD->Fill(jet_pTD);
                        SIGMultiplicity->Fill(jet_multiplicity);
                        SIGXobsPb->Fill(((cluspT*TMath::Exp(-cluseta))+(jet_pT*TMath::Exp(-jet_eta)))/(2*EPb));
                        
                        z_Vertices->Fill(TMath::Abs(event_data_out[0][0] - primary_vertex[2]));
                        z_Vertices_individual->Fill(primary_vertex[2]);
                        z_Vertices_hdf5->Fill(event_data_out[0][0]);
                        
                        float multiplicity_sum = 0;
                        for (int k = 0; k < 64; k++)  multiplicity_sum += multiplicity_v0[k];
                        //std::cout << "Multiplicity difference " << TMath::Abs(event_data_out[0][1] - multiplicity_sum) << std::endl;
                        Multiplicity->Fill(TMath::Abs(event_data_out[0][1] - multiplicity_sum));
                        Multiplicity_individual->Fill(multiplicity_sum);
                        Multiplicity_hdf5->Fill(event_data_out[0][1]);
                    }
                    
                    if((cluster_lambda_square[icluster][0] > 0.4) && (cluster_lambda_square[icluster][0] < 1.0)) {
                        BKGcluster_pt_dist->Fill(cluspT);
                        BKGjet_pt_dist->Fill(jet_pT);
                        BKGpt_diff_dist->Fill(TMath::Abs(cluspT-jet_pT));
                        
                        double dphinum = clusphi-jet_phi;
                        while(dphinum < -TMath::Pi()) dphinum += 2*TMath::Pi();
                        while(dphinum > TMath::Pi()) dphinum -= 2*TMath::Pi();
                        BKGdPhi->Fill(TMath::Abs(dphinum));
                        if(not(dphinum > TMath::Pi()/2)) continue;
                        BKGclusterPhi->Fill(clusphi);
                        BKGjetPhi->Fill(jet_phi);
                        
                        BKGdEta->Fill(jet_eta-cluseta);
                        BKGclusterEta->Fill(cluseta);
                        BKGjetEta->Fill(jet_eta);
                        
                        BKGXj->Fill(jet_pT/cluspT);
                        BKGpTD->Fill(jet_pTD);
                        BKGMultiplicity->Fill(jet_multiplicity);
                        BKGXobsPb->Fill(((cluspT*TMath::Exp(-cluseta))+(jet_pT*TMath::Exp(-jet_eta)))/(2*EPb));
                        
                    }
                    
                }
            }
            
        }//end loop over mixed events
        if(ievent % 10000 == 0)
            std::cout << "Event " << ievent << std::endl;
    } //end loop over events
    
    //very particular about file names to ease scripting
    // Write to fout
    
    std::string rawname = ((std::string)root_file).substr(((std::string)root_file).find_last_of("/")+1, ((std::string)root_file).find_last_of(".")-((std::string)root_file).find_last_of("/")-1);
    //std::string rawname = std::string(argv[1]);
    TFile* fout = new TFile(Form("New_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.root",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin),"RECREATE");
    std::cout<< "Created datafile: " << Form("New_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.root",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin) << std::endl;
    // Normalize
    SIGcluster_pt_dist->Scale(1.0/(N_SR*SIGcluster_pt_dist_binwidth));
    SIGjet_pt_dist->Scale(1.0/(N_SR*SIGjet_pt_dist_binwidth));
    SIGpt_diff_dist->Scale(1.0/(N_SR*SIGpt_diff_dist_binwidth));
    
    SIGdPhi->Scale(1.0/(N_SR*SIGdPhi_binwidth));
    SIGclusterPhi->Scale(1.0/(N_SR*SIGclusterPhi_binwidth));
    SIGjetPhi->Scale(1.0/(N_SR*SIGjetPhi_binwidth));
    
    SIGdEta->Scale(1.0/(N_SR*SIGdEta_binwidth));
    SIGclusterEta->Scale(1.0/(N_SR*SIGclusterEta_binwidth));
    SIGjetEta->Scale(1.0/(N_SR*SIGjetEta_binwidth));
    
    SIGXj->Scale(1.0/(N_SR*SIGXj_binwidth));
    SIGpTD->Scale(1.0/(N_SR*SIGpTD_binwidth));
    SIGMultiplicity->Scale(1.0/(N_SR*SIGMultiplicity_binwidth));
    SIGXobsPb->Scale(1.0/(N_SR*SIGXobsPb_binwidth));
    
    BKGcluster_pt_dist->Scale(1.0/(N_BR*BKGcluster_pt_dist_binwidth));
    BKGjet_pt_dist->Scale(1.0/(N_BR*BKGjet_pt_dist_binwidth));
    BKGpt_diff_dist->Scale(1.0/(N_BR*BKGpt_diff_dist_binwidth));
    
    BKGdPhi->Scale(1.0/(N_BR*BKGdPhi_binwidth));
    BKGclusterPhi->Scale(1.0/(N_BR*BKGclusterPhi_binwidth));
    BKGjetPhi->Scale(1.0/(N_BR*BKGjetPhi_binwidth));
    
    BKGdEta->Scale(1.0/(N_BR*BKGdEta_binwidth));
    BKGclusterEta->Scale(1.0/(N_BR*BKGclusterEta_binwidth));
    BKGjetEta->Scale(1.0/(N_BR*BKGjetEta_binwidth));
    
    BKGXj->Scale(1.0/(N_BR*BKGXj_binwidth));
    BKGpTD->Scale(1.0/(N_BR*BKGpTD_binwidth));
    BKGMultiplicity->Scale(1.0/(N_BR*BKGMultiplicity_binwidth));
    BKGXobsPb->Scale(1.0/(N_BR*BKGXobsPb_binwidth));
    
    // Set minima
    SIGcluster_pt_dist->SetMinimum(0);
    SIGjet_pt_dist->SetMinimum(0);
    SIGpt_diff_dist->SetMinimum(0);
    
    SIGdPhi->SetMinimum(0);
    SIGclusterPhi->SetMinimum(0);
    SIGjetPhi->SetMinimum(0);
    
    SIGdEta->SetMinimum(0);
    SIGclusterEta->SetMinimum(0);
    SIGjetEta->SetMinimum(0);
    
    SIGXj->SetMinimum(0);
    SIGpTD->SetMinimum(0);
    SIGMultiplicity->SetMinimum(0);
    SIGXobsPb->SetMinimum(0);
    
    BKGcluster_pt_dist->SetMinimum(0);
    BKGjet_pt_dist->SetMinimum(0);
    BKGpt_diff_dist->SetMinimum(0);
    
    BKGdPhi->SetMinimum(0);
    BKGclusterPhi->SetMinimum(0);
    BKGjetPhi->SetMinimum(0);
    
    BKGdEta->SetMinimum(0);
    BKGclusterEta->SetMinimum(0);
    BKGjetEta->SetMinimum(0);
    
    BKGXj->SetMinimum(0);
    BKGpTD->SetMinimum(0);
    BKGMultiplicity->SetMinimum(0);
    BKGXobsPb->SetMinimum(0);
    
    //Write histograms here
    
    SIGcluster_pt_dist->Write();
    SIGjet_pt_dist->Write();
    SIGpt_diff_dist->Write();
    
    SIGdPhi->Write();
    SIGclusterPhi->Write();
    SIGjetPhi->Write();
    
    SIGdEta->Write();
    SIGclusterEta->Write();
    SIGjetEta->Write();
    
    SIGXj->Write();
    SIGpTD->Write();
    SIGMultiplicity->Write();
    SIGXobsPb->Write();
    
    BKGcluster_pt_dist->Write();
    BKGjet_pt_dist->Write();
    BKGpt_diff_dist->Write();
    
    BKGdPhi->Write();
    BKGclusterPhi->Write();
    BKGjetPhi->Write();
    
    BKGdEta->Write();
    BKGclusterEta->Write();
    BKGjetEta->Write();
    
    BKGXj->Write();
    BKGpTD->Write();
    BKGMultiplicity->Write();
    BKGXobsPb->Write();
    
    z_Vertices->Write();
    Multiplicity->Write();
    z_Vertices_individual->Write();
    z_Vertices_hdf5->Write();
    Multiplicity_individual->Write();
    Multiplicity_hdf5->Write();
    
    // Commented out due to segfaults
    /*
    SIGcluster_pt_dist->Draw();
    canvas.SaveAs(Form("signal_cluster_pT_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    SIGjet_pt_dist->Draw();
    canvas.SaveAs(Form("signal_jet_pt_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f._minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    SIGpt_diff_dist->Draw();
    canvas.SaveAs(Form("signal_pt_diff_differences_individual_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    
    SIGclusterPhi->Draw();
    canvas.SaveAs(Form("signal_cluster_phi_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    SIGjetPhi->Draw();
    canvas.SaveAs(Form("signal_jet_phi_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    SIGdPhi->Draw();
    canvas.SaveAs(Form("signal_delta_phi_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    
    SIGdEta->Draw();
    canvas.SaveAs(Form("signal_delta_eta_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    SIGclusterEta->Draw();
    canvas.SaveAs(Form("signal_cluster_eta_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    SIGjetEta->Draw();
    canvas.SaveAs(Form("signal_jet_eta_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    
    SIGXj->Draw();
    canvas.SaveAs(Form("signal_Xj_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    SIGpTD->Draw();
    canvas.SaveAs(Form("signal_pTD_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    SIGMultiplicity->Draw();
    canvas.SaveAs(Form("signal_multiplicity_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    
    BKGcluster_pt_dist->Draw();
    canvas.SaveAs(Form("background_cluster_pT_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    BKGjet_pt_dist->Draw();
    canvas.SaveAs(Form("background_jet_pt_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    BKGpt_diff_dist->Draw();
    canvas.SaveAs(Form("background_pt_diff_differences_individual_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    
    BKGclusterPhi->Draw();
    canvas.SaveAs(Form("background_cluster_phi_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    BKGjetPhi->Draw();
    canvas.SaveAs(Form("background_jet_phi_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    BKGdPhi->Draw();
    canvas.SaveAs(Form("background_delta_phi_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    
    BKGdEta->Draw();
    canvas.SaveAs(Form("background_delta_eta_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    BKGclusterEta->Draw();
    canvas.SaveAs(Form("background_cluster_eta_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    BKGjetEta->Draw();
    canvas.SaveAs(Form("background_jet_eta_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    
    BKGXj->Draw();
    canvas.SaveAs(Form("background_Xj_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    BKGpTD->Draw();
    canvas.SaveAs(Form("background_pTD_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    BKGMultiplicity->Draw();
    canvas.SaveAs(Form("background_multiplicity_distribution_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu_cluspT_%1.1f_to_%1.1f_minjetpT_%2.1f.png",rawname.data(),GeV_Track_Skim,mix_start,mix_end, cluspTmin, cluspTmax, jetpTmin));
    canvas.Clear();
    
    z_Vertices->Draw();
   
    std::string filepath = argv[1];
    std::string opened_files = "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
    //std::string rawname = std::string(argv[1]);
    canvas.SaveAs(Form("z_Vertices_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu.png",opened_files.c_str(),GeV_Track_Skim,mix_start,mix_end));
    canvas.Clear();
    Multiplicity->Draw();
    canvas.SaveAs(Form("Multiplicity_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu.png",opened_files.c_str(),GeV_Track_Skim,mix_start,mix_end));
    //Signal_pT_Dist->Write();
    canvas.Clear();
    z_Vertices_individual->Draw();
    canvas.SaveAs(Form("z_Vertices_individual_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu.png",opened_files.c_str(),GeV_Track_Skim,mix_start,mix_end));
    canvas.Clear();
    z_Vertices_hdf5->Draw();
    canvas.SaveAs(Form("z_Vertices_hdf5l_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu.png",opened_files.c_str(),GeV_Track_Skim,mix_start,mix_end));
    canvas.Clear();
    Multiplicity_individual->Draw();
    canvas.SaveAs(Form("Multiplicity_individual_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu.png",opened_files.c_str(),GeV_Track_Skim,mix_start,mix_end));
    canvas.Clear();
    Multiplicity_hdf5->Draw();
    canvas.SaveAs(Form("Multiplicity_hdf5_%s_%luGeVTracks_Correlation_%1.1lu_to_%1.1lu.png",opened_files.c_str(),GeV_Track_Skim,mix_start,mix_end));
    canvas.Clear(); */
    
    canvas.Close();
    
    fout->Close();
    
    std::cout << " ending; num of signal triggers is " << N_SR << " background " << N_BR << std::endl;
    return EXIT_SUCCESS;
}
