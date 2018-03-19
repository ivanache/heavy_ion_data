// This file makes combined direct-photon correlations between same-event and mixed-event samples
// Author: Ivan Chernyshev; Date: 3/15/2018
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

// 2D histogram dividing function, fitted with parameters for constructing the quotient histogram
// Precondition: the two histograms must have the same y- and x-dimensions, and so must the quotient
TH2D* divide_histograms2D(TH2D* graph1, TH2D* graph2, const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup){
    // Make the 2D histogram to contain the quotient and find minimum and maximum bins along both axes
    TH2D* quotient = new TH2D(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
    quotient->Sumw2();
    quotient->SetMinimum(0.);
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

// 1D histogram dividing function, fitted with parameters for constructing the quotient histogram
// Precondition: the two histograms must have the same x-dimensions, and the quotient must have the same x-dimensions
TH1F* divide_histograms1D(TH1F* graph1, TH1F* graph2, const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup){
    // Make the 1D histogram to contain the quotient and find minimum and maximum bins along both axes
    TH1F* quotient = new TH1F(name, title, nbinsx, xlow, xup);
    quotient->Sumw2();
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

int main(int argc, char *argv[])
{
    
    if (argc < 3) {
        std::cout<<"Temporary Syntax for Mixing is [Command] [mixed_root_file] [hdf5_file] (repeat for mult. samples)"<<std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Read configuration file for the DNN_min and DNN_max variables
    FILE* config = fopen("Corr_config.yaml", "r");
    if (config == NULL)  std::cout<<"no config"<<std::endl;
    // Default values of various variables used in the file (actual values are to be determined by the configuration file)
    // Cut variables
    double DNN_min = 0.55;
    double DNN_max = 0.85;
    double pT_min = 10;
    double pT_max = 16;
    double Eta_max = 0.6;
    double Cluster_min = 2;
    double EcrossoverE_min = 0.03;
    int selection_number = 3;
    
    // The bounds for the events to fal into the isolation and nonisolation areas
    double iso_max = 2.0;
    double noniso_min = 5.0;
    double noniso_max = 10.0;
    
    // Delta eta
    double deta_max = 0.6;
    
    // Zt bins
    int nztbins = 7;
    float* ztbins;
    ztbins = new float[nztbins+1];
    ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
    
    // Number of bins in correlation functions
    int n_correlationbins = 18;
    
    // Which branch should be used to determine whether a cluster should fall into iso, noniso, or neither
    isolationDet determiner = CLUSTER_ISO_ITS_04;
    
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
        if (strcmp(key, "DNN_min") == 0) {
            // Assign DNN_min to the double-converted version of value
            DNN_min = atof(value);
            std::cout << "DNN_min is " << DNN_min << std::endl;
        }
        else if (strcmp(key, "DNN_max") == 0) {
            // Assign DNN_max to the double-converted version of value
            DNN_max = atof(value);
            std::cout << "DNN_max is " << DNN_max << std::endl;
        }
        else if (strcmp(key, "pT_min") == 0) {
            pT_min = atof(value);
            std::cout << "pT_min is " << pT_min << std::endl;
        }
        else if (strcmp(key, "pT_max") == 0) {
            pT_max = atof(value);
            std::cout << "pT_max is " << pT_max << std::endl;
        }
        else if (strcmp(key, "Eta_max") == 0) {
            Eta_max = atof(value);
            std::cout << "Eta_max is " << Eta_max << std::endl;
        }
        else if (strcmp(key, "Cluster_min") == 0) {
            Cluster_min = atof(value);
            std::cout << "Cluster_min is " << Cluster_min << std::endl;
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
        else if (strcmp(key, "Correlation_func_bins") == 0) {
            n_correlationbins = atoi(value);
            std::cout << "Bins in a correlation function: " << n_correlationbins << std::endl;
        }
        else if (strcmp(key, "track_selection_num") == 0) {
            selection_number = atoi(value);
            std::cout << "Number of Selection that tracks must pass: " << selection_number << std::endl;
        }
        else if (strcmp(key, "Zt_bins") == 0) {
            nztbins = -1;
            for (const char *v = value; *v != ']';) {
                while (*v != ']' && !isdigit(*v)) {
                    v++;
                }
                
                nztbins++;
                
                while (*v != ']' && (isdigit(*v) || *v == '.')) {
                    v++;
                }
            }
            ztbins = new float[nztbins + 1];
            int i = 0;
            for (const char *v = value; *v != ']' ;) {
                while (*v != ']' && !isdigit(*v)) {
                    v++;
                }
                ztbins[i] = atof(v);
                i++;
                while (*v != ']' && (isdigit(*v) || *v == '.')) {
                    v++;
                }
            }
            std::cout << "Number of Zt bins: " << nztbins << std::endl << "Zt bins: {";
            for (int i = 0; i <= nztbins; i++)
                std::cout << ztbins[i] << ", ";
            std::cout << "}\n";
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
        else {
            std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
        }
    }
    fclose(config);
    
    for (int i = 0; i <= nztbins; i++){
        std::cout << "zt bound: " << ztbins[i] << std::endl;
    }
    
    // Create the TCanvas and the histograms
    TCanvas canvas("canvas", "");
    TH1D histogram0("histogram0", "", 16, 8.0, 16.0);
    //TH2D histogram1("histogram1", "", 30, -1.5, 1.5, 18, -0.5, 1.5);
    //TH1D histogram2("histogram2", "", 18, -0.5,1.5);
    TH1D histogram3("histogram3", "", 18, -0.5,1.5);
    TH1D h_ntrig("h_ntrig", "", 2, -0.5,1.0);
    
    // Create the histogram for the 2D plots
    // TH2D histogram2D0("histogram2D0", "", );
    
    // Zt bins
    //const int nztbins = 7;
    //const float ztbins[nztbins+1] = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
    
    // Function declarations of all graphs that will be included here, for the mixed, same-event, and correlation function
    TH1F* h_dPhi_iso_same[nztbins];
    TH1F* h_dPhi_noniso_same[nztbins];
    TH2D* Map_same[nztbins];
    TH2D* IsoMap_same[nztbins];
    TH2D* AntiIsoMap_same[nztbins];
    
    TH1F* h_dPhi_iso_mixed[nztbins];
    TH1F* h_dPhi_noniso_mixed[nztbins];
    TH2D* Map_mixed[nztbins];
    TH2D* IsoMap_mixed[nztbins];
    TH2D* AntiIsoMap_mixed[nztbins];
    
    TH1F* h_dPhi_iso_corr[nztbins];
    TH1F* h_dPhi_noniso_corr[nztbins];
    TH2D* Map_corr[nztbins];
    TH2D* IsoMap_corr[nztbins];
    TH2D* AntiIsoMap_corr[nztbins];
    
    int mixed_clusters_passed_iso[nztbins];
    int mixed_clusters_passed_Antiiso[nztbins];
    int same_clusters_passed_iso[nztbins];
    int same_clusters_passed_Antiiso[nztbins];
    int corr_clusters_passed_iso[nztbins];
    int corr_clusters_passed_Antiiso[nztbins];
    
    int n_eta_bins = 14;
    int n_phi_bins = 24;
    
    for (int izt = 0; izt<nztbins; izt++){
        // Same events
        h_dPhi_iso_same[izt] = new TH1F(Form("dPhi_iso_same_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]) ,"", n_correlationbins,-0.5,1.5);
        h_dPhi_noniso_same[izt] = new TH1F(Form("dPhi_noniso_same_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt], 10*ztbins[izt+1]), "", n_correlationbins,-0.5,1.5);
        h_dPhi_iso_same[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
        h_dPhi_noniso_same[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
        h_dPhi_iso_same[izt]->Sumw2();
        h_dPhi_noniso_same[izt]->Sumw2();
        
        
        Map_same[izt] = new TH2D(Form("Map_same_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),
                                 "Same Event #gamma-H [all] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        Map_same[izt]->Sumw2();
        Map_same[izt]->SetMinimum(0.);
        
        IsoMap_same[izt] = new TH2D(Form("IsoMap_same_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),
                                    "Same Event #gamma-H [Iso] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        IsoMap_same[izt]->Sumw2();
        IsoMap_same[izt]->SetMinimum(0.);
        
        AntiIsoMap_same[izt] = new TH2D(Form("AntiIsoMap_same_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),
                                        "Same Event #gamma-H [AntiIso] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        AntiIsoMap_same[izt]->Sumw2();
        AntiIsoMap_same[izt]->SetMinimum(0.);
        
        same_clusters_passed_iso[izt] = 0;
        same_clusters_passed_Antiiso[izt] = 0;
        
        // Mixed events
        
        h_dPhi_iso_mixed[izt] = new TH1F(Form("dPhi_iso_mixed_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]) ,"", n_correlationbins,-0.5,1.5);
        h_dPhi_noniso_mixed[izt] = new TH1F(Form("dPhi_noniso_mixed_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt], 10*ztbins[izt+1]), "", n_correlationbins,-0.5,1.5);
        h_dPhi_iso_mixed[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
        h_dPhi_noniso_mixed[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
        h_dPhi_iso_mixed[izt]->Sumw2();
        h_dPhi_noniso_mixed[izt]->Sumw2();
        
        
        Map_mixed[izt] = new TH2D(Form("Map_mixed_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),
                                  "Mixed Event #gamma-H [all] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        Map_mixed[izt]->Sumw2();
        Map_mixed[izt]->SetMinimum(0.);
        
        IsoMap_mixed[izt] = new TH2D(Form("IsoMap_mixed_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),
                                     "Mixed Event #gamma-H [Iso] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        IsoMap_mixed[izt]->Sumw2();
        IsoMap_mixed[izt]->SetMinimum(0.);
        
        AntiIsoMap_mixed[izt] = new TH2D(Form("AntiIsoMap_mixed_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),
                                         "Mixed Event #gamma-H [AntiIso] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        AntiIsoMap_mixed[izt]->Sumw2();
        AntiIsoMap_mixed[izt]->SetMinimum(0.);
        
        mixed_clusters_passed_iso[izt] = 0;
        mixed_clusters_passed_Antiiso[izt] = 0;
        
        /*
        // Correlations
        h_dPhi_iso_corr[izt] = new TH1F(Form("dPhi_iso_corr_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]) ,"", n_correlationbins,-0.5,1.5);
        h_dPhi_noniso_corr[izt] = new TH1F(Form("dPhi_noniso_corr_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt], 10*ztbins[izt+1]), "", n_correlationbins,-0.5,1.5);
        h_dPhi_iso_corr[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
        h_dPhi_noniso_corr[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
        h_dPhi_iso_corr[izt]->Sumw2();
        h_dPhi_noniso_corr[izt]->Sumw2();
        
        
        Map_corr[izt] = new TH2D(Form("Map_corr_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),
                                 "Correlation #gamma-H [all] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        Map_corr[izt]->Sumw2();
        Map_corr[izt]->SetMinimum(0.);
        
        IsoMap_corr[izt] = new TH2D(Form("IsoMap_corr_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),
                                    "Corrrelation #gamma-H [Iso] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        IsoMap_corr[izt]->Sumw2();
        IsoMap_corr[izt]->SetMinimum(0.);
        
        AntiIsoMap_corr[izt] = new TH2D(Form("AntiIsoMap_corr_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]),
                                        "Correlation #gamma-H [AntiIso] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        AntiIsoMap_corr[izt]->Sumw2();
        AntiIsoMap_corr[izt]->SetMinimum(0.);
        
        corr_clusters_passed_iso[izt] = 0;
        corr_clusters_passed_Antiiso[izt] = 0;*/
    }
    
    //histogram2.Sumw2();
    //histogram3.Sumw2();
    
    //Loop over samples
    for (int iarg = 1; iarg < argc; iarg+=2) {
        std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
        TFile *file = TFile::Open((TString)argv[iarg]);
        
        if (file == NULL) {
            std::cout << " fail" << std::endl;
            exit(EXIT_FAILURE);
        }
        file->Print();
        
        TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));
        if (_tree_event == NULL) {
            std::cout << " fail " << std::endl;
            exit(EXIT_FAILURE);
        }
        
        //variables
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
        
        Long64_t Mix_Events[50];
        
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
        _tree_event->SetBranchStatus("*mc*", 0);
        
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
        _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
        _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
        _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
        _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
        
        _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
        _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
        _tree_event->SetBranchAddress("cell_e", cell_e);
        
        _tree_event->SetBranchAddress("Mix_Events", Mix_Events);
        //_tree_event->SetBranchAddress("LimitUse_Mixed_Events", Mix_Events);
        
        std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;
        
        //     UInt_t ntrack_max = 553;
        //     UInt_t ncluster_max = 150;
        
        UInt_t ntrack_max = 0;
        UInt_t ncluster_max = 0;
        
        fprintf(stderr, "\r%s:%d: %s\n", __FILE__, __LINE__, "Determining ntrack_max and ncluster_max needed for hdf5 hyperslab");
        for (Long64_t i = 0; i < _tree_event->GetEntries(); i++) {
            _tree_event->GetEntry(i);
            ntrack_max = std::max(ntrack_max, ntrack);
            ncluster_max = std::max(ncluster_max, ncluster);
            fprintf(stderr, "\r%s:%d: %llu", __FILE__, __LINE__, i);
        }
        fprintf(stderr, "\n%s:%d: maximum tracks:%i maximum clusters:%i\n", __FILE__, __LINE__, ntrack_max,ncluster_max);
        
        //open hdf5
        const H5std_string hdf5_file_name(argv[iarg+1]);
        //FIXME: Can parse the string s.t. remove .root and simply add .hdf5
        //FIXME: This will obviously require stringent naming conventions
        
        const H5std_string track_ds_name( "track" );
        H5File h5_file( hdf5_file_name, H5F_ACC_RDONLY );
        DataSet track_dataset = h5_file.openDataSet( track_ds_name );
        DataSpace track_dataspace = track_dataset.getSpace();
        
        const H5std_string cluster_ds_name( "cluster" );
        DataSet cluster_dataset = h5_file.openDataSet( cluster_ds_name );
        DataSpace cluster_dataspace = cluster_dataset.getSpace();
        
        //Define array hyperslab will be read into
        float track_data_out[1][ntrack_max][7];
        float cluster_data_out[1][ncluster_max][5];
        
        //Define hyperslab size and offset in  FILE;
        hsize_t track_offset[3] = {0, 0, 0};
        hsize_t track_count[3] = {1, ntrack_max, 7};
        hsize_t cluster_offset[3] = {0, 0, 0};
        hsize_t cluster_count[3] = {1, ncluster_max, 5};
        
        track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
        cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");
        
        //Define the memory dataspace to place hyperslab
        const int RANK_OUT = 3;
        hsize_t track_dimsm[3] = {1, ntrack_max, 7};      //memory space dimensions
        DataSpace track_memspace( RANK_OUT, track_dimsm );
        hsize_t cluster_dimsm[3] = {1, ncluster_max, 5};      //memory space dimensions
        DataSpace cluster_memspace( RANK_OUT, cluster_dimsm );
        
        //Define memory offset for hypreslab->array
        hsize_t track_offset_out[3] = {0};
        hsize_t cluster_offset_out[3] = {0};
        
        //define 2D memory hyperslab
        hsize_t track_count_out[3] = {1, ntrack_max, 7};    // size of the hyperslab in memory
        hsize_t cluster_count_out[3] = {1, ncluster_max, 5};
        
        //define space in memory for hyperslab, then write from file to memory
        track_memspace.selectHyperslab( H5S_SELECT_SET, track_count_out, track_offset_out );
        track_dataset.read( track_data_out, PredType::NATIVE_FLOAT, track_memspace, track_dataspace );
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "track dataset read into array: OK");
        
        cluster_memspace.selectHyperslab( H5S_SELECT_SET, cluster_count_out, cluster_offset_out );
        cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "cluster dataset read into array: OK");
        
        Long64_t nentries = _tree_event->GetEntries();
        
        for(Long64_t ievent = 0; ievent < nentries ; ievent++){
            //for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
            _tree_event->GetEntry(ievent);
            fprintf(stderr, "\r%s:%d: %llu / %llu", __FILE__, __LINE__, ievent, nentries);
            //fprintf(stderr, "%s:%d: %llu / %llu ENTRIES\n", __FILE__, __LINE__, ievent, nentries);
            
            for (ULong64_t n = 0; n < ncluster; n++) {
                if( not(cluster_pt[n]>pT_min and cluster_pt[n]<pT_max)) continue; //select pt of photons
                if( not(cluster_s_nphoton[n][1]>DNN_min and cluster_s_nphoton[n][1]<DNN_max)) continue; //select deep-photons
                if( not(TMath::Abs(cluster_eta[n])<Eta_max)) continue; //cut edges of detector
                if( not(cluster_ncell[n]>Cluster_min)) continue;   //removes clusters with 1 or 2 cells
                if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters
                
                
                //determiner: which frixione variable based on Corr_config.yaml
                float isolation;
                if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n];
                else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n];
                else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n];
                else isolation = cluster_frixione_its_04_02[n];
                
                //isolation
                if(isolation<iso_max){
                    histogram0.Fill(cluster_pt[n]);
                    h_ntrig.Fill(0);
                }
                if(isolation>noniso_min && isolation<noniso_max){
                    h_ntrig.Fill(0.5);
                }
                
                // Mixed events
                for (Long64_t imix = 0; imix < 50; imix++){
                    Long64_t mix_event = Mix_Events[imix];
                    if (mix_event == ievent) continue;
                    
                    //adjust hyper slab offset for next mixed event
                    track_offset[0]=mix_event;
                    track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
                    track_dataset.read( track_data_out, PredType::NATIVE_FLOAT, track_memspace, track_dataspace );
                    
                    cluster_offset[0]=mix_event;
                    cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
                    cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
                    
                    const int TrackCutBit =16;
                    for (ULong64_t itrack = 0; itrack < ntrack_max; itrack++) {
                        if (std::isnan(track_data_out[0][itrack][1])) break;
                        if((track_quality[itrack]&selection_number)==0) continue; //pass 3 cut
                        //if (track_data_out[0][itrack][1] < 0.15) continue;
                        if (track_data_out[0][itrack][1] < 2) continue; //less than 2GeV
                        
                        //veto charged particles from mixed event tracks
                        //the isolation takes care of the initial culster, but does nothing for mixed event track
                        bool MixTrack_HasMatch = false;
                        for (unsigned int l = 0; l < ncluster_max; l++){
                            if (std::isnan(cluster_data_out[0][l][0])) break;
                            if (TMath::Abs(cluster_data_out[0][l][2] - track_data_out[0][itrack][5]) < 0.015  &&
                                TMath::Abs(cluster_data_out[0][l][3] - track_data_out[0][itrack][6]) < 0.015) {
                                MixTrack_HasMatch = true;
                                break;
                            }
                        }
                        if (MixTrack_HasMatch) continue;
                        //fprintf(stderr, "%s:%d: Mixed Event: %llu Track: %llu\n", __FILE__, __LINE__, mix_event, itrack);
                        
                        //FIXME: Lazy implementation from past code. Will use this repositories ∆'s soon
                        Float_t DeltaPhi = cluster_phi[n] - track_data_out[0][itrack][3];
                        if (DeltaPhi < -M_PI/2){DeltaPhi += 2*M_PI;}  //if less then -pi/2 add 2pi
                        if (DeltaPhi > 3*M_PI/2){DeltaPhi =DeltaPhi -2*M_PI;}
                        Float_t DeltaEta = cluster_eta[n] - track_data_out[0][itrack][2];
                        if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue;
                        
                        Double_t zt = track_data_out[0][itrack][1]/cluster_pt[n];
                        Float_t deta =  cluster_eta[n]-track_data_out[0][itrack][2];
                        Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_data_out[0][itrack][3]);
                        dphi = dphi/TMath::Pi();
                        //if(!(TMath::Abs(deta)<0.6)) continue; //deta cut
                        if(dphi<-0.5) dphi +=2;
                        
                        // Loop over zt bins
                        for(int izt = 0; izt<nztbins ; izt++){
                            if(zt>ztbins[izt] and  zt<ztbins[izt+1])
                            {
                                // Where the  h_dPhi_iso and h_dPhi_noniso bins are filled
                                if(isolation< iso_max){
                                    h_dPhi_iso_mixed[izt]->Fill(DeltaPhi);
                                    IsoMap_mixed[izt]->Fill(DeltaPhi,DeltaEta);
                                }
                                if(isolation> noniso_min && isolation<noniso_max){
                                    h_dPhi_noniso_mixed[izt]->Fill(DeltaPhi);
                                    AntiIsoMap_mixed[izt]->Fill(DeltaPhi,DeltaEta);
                                }
                                Map_mixed[izt]->Fill(DeltaPhi,DeltaEta);
                            }
                        } // end loop over zt bins
                    }//end loop over tracks
                }//end loop over mixed events
                
                // Loop over tracks: for the same events
                for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {
                    if((track_quality[itrack]&selection_number)==0) continue; //select only tracks that pass selection 3
                    
                    // Fernando's new dphi and deta, the one I'll be using
                    Float_t DeltaPhi = cluster_phi[n] - track_phi[itrack];
                    if (DeltaPhi < -M_PI/2){DeltaPhi += 2*M_PI;}  //if less then -pi/2 add 2pi
                    if (DeltaPhi > 3*M_PI/2){DeltaPhi =DeltaPhi -2*M_PI;}
                    Float_t DeltaEta = cluster_eta[n] - track_eta[itrack];
                    if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue;
                    
                    // Find the track's zt, delta eta, and delta phi with the cluster
                    Double_t zt = track_pt[itrack]/cluster_pt[n];
                    Float_t deta =  cluster_eta[n]-track_eta[itrack];
                    Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_phi[itrack])/TMath::Pi();
                    //if(!(TMath::Abs(deta)<deta_max)) continue; // delta eta cut
                    if(dphi<-0.5) dphi +=2;
                    
                    // Loop over zt bins
                    for(int izt = 0; izt<nztbins ; izt++){
                        if(zt>ztbins[izt] and  zt<ztbins[izt+1])
                        {
                            // Where the  h_dPhi_iso and h_dPhi_noniso bins are filled, along with all of the maps
                            if(isolation< iso_max){
                                h_dPhi_iso_same[izt]->Fill(DeltaPhi);
                                IsoMap_same[izt]->Fill(DeltaPhi,DeltaEta);
                            }
                            if(isolation> noniso_min && isolation<noniso_max){
                                h_dPhi_noniso_same[izt]->Fill(DeltaPhi);
                                AntiIsoMap_same[izt]->Fill(DeltaPhi,DeltaEta);
                            }
                            Map_same[izt]->Fill(DeltaPhi,DeltaEta);
                        }
                    } // end loop over bins
                }//end loop over tracks (same events)
                
            }//end loop on clusters.
            if (ievent % 25000 == 0) {
                histogram0.Draw("e1x0");
                canvas.Update();
            }
        } //end loop over events
        
    }//end loop over samples
    
    // Solve for the correlations
    for (int izt = 0; izt<nztbins; izt++){
        h_dPhi_iso_corr[izt] = divide_histograms1D(h_dPhi_iso_same[izt], h_dPhi_iso_mixed[izt], Form("dPhi_iso_corr_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]) ,"", n_correlationbins,-0.5,1.5);
        h_dPhi_iso_corr[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
        h_dPhi_noniso_corr[izt] = divide_histograms1D(h_dPhi_noniso_same[izt], h_dPhi_noniso_mixed[izt], Form("dPhi_noniso_corr_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt], 10*ztbins[izt+1]), "", n_correlationbins,-0.5,1.5);
        h_dPhi_noniso_corr[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
        Map_corr[izt] = divide_histograms2D(Map_same[izt], Map_mixed[izt], Form("Map_corr_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]), "Correlation #gamma-H [all] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        IsoMap_corr[izt] = divide_histograms2D(IsoMap_same[izt], IsoMap_mixed[izt], Form("IsoMap_corr_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]), "Corrrelation #gamma-H [Iso] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
        AntiIsoMap_corr[izt] = divide_histograms2D(AntiIsoMap_same[izt], AntiIsoMap_mixed[izt], Form("AntiIsoMap_corr_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]), "Correlation #gamma-H [AntiIso] Map", n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
    }
    
    // Write to fout
    //TFile* fout = new TFile(Form("fout_Corr_config%s.root", opened_files.c_str()),"RECREATE");
    TFile* fout = new TFile("fout_all_frixione.root","RECREATE");
    histogram0.Write("DeepPhotonSpectra");
    h_ntrig.Write("ntriggers");
    
    
    for (int izt = 0; izt<nztbins; izt++){
        h_dPhi_iso_same[izt]->SetMinimum(0.0);
        h_dPhi_iso_same[izt]->Write();
        h_dPhi_iso_same[izt]->Draw();
        canvas.SaveAs(Form("dPhi_iso_same_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
        
        h_dPhi_noniso_same[izt]->SetMinimum(0.0);
        h_dPhi_noniso_same[izt]->Write();
        h_dPhi_noniso_same[izt]->Draw();
        canvas.SaveAs(Form("dPhi_noniso_same_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt], 10*ztbins[izt+1]));
        canvas.Clear();
    }
    for (int izt = 0; izt<nztbins; izt++){
        Map_same[izt]->Write();
        Map_same[izt]->Draw("COLZ");
        canvas.SaveAs(Form("Map_same_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
        IsoMap_same[izt]->Write();
        IsoMap_same[izt]->Draw("COLZ");
        canvas.SaveAs(Form("IsoMap_same_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
        AntiIsoMap_same[izt]->Write();
        AntiIsoMap_same[izt]->Draw("COLZ");
        canvas.SaveAs(Form("AntiIsoMap_same_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
    }
    
    for (int izt = 0; izt<nztbins; izt++){
        h_dPhi_iso_mixed[izt]->SetMinimum(0.0);
        h_dPhi_iso_mixed[izt]->Write();
        h_dPhi_iso_mixed[izt]->Draw();
        canvas.SaveAs(Form("dPhi_iso_mixed_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
        h_dPhi_noniso_mixed[izt]->SetMinimum(0.0);
        h_dPhi_noniso_mixed[izt]->Write();
        h_dPhi_noniso_mixed[izt]->Draw();
        canvas.SaveAs(Form("dPhi_noniso_mixed_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt], 10*ztbins[izt+1]));
        canvas.Clear();
    }
    
    for (int izt = 0; izt<nztbins; izt++){
        Map_mixed[izt]->Write();
        Map_mixed[izt]->Draw("COLZ");
        canvas.SaveAs(Form("Map_mixed_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
        IsoMap_mixed[izt]->Write();
        IsoMap_mixed[izt]->Draw("COLZ");
        canvas.SaveAs(Form("IsoMap_mixed_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
        AntiIsoMap_mixed[izt]->Write();
        AntiIsoMap_mixed[izt]->Draw("COLZ");
        canvas.SaveAs(Form("AntiIsoMap_mixed_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
    }
    
    for (int izt = 0; izt<nztbins; izt++){
        h_dPhi_iso_corr[izt]->SetMinimum(0.0);
        h_dPhi_iso_corr[izt]->Write();
        h_dPhi_iso_corr[izt]->Draw();
        canvas.SaveAs(Form("dPhi_iso_corr_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
        h_dPhi_noniso_corr[izt]->SetMinimum(0.0);
        h_dPhi_noniso_corr[izt]->Write();
        h_dPhi_noniso_corr[izt]->Draw();
        canvas.SaveAs(Form("dPhi_noniso_corr_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt], 10*ztbins[izt+1]));
        canvas.Clear();
    }
    
    for (int izt = 0; izt<nztbins; izt++){
        Map_corr[izt]->Write();
        Map_corr[izt]->Draw("COLZ");
        canvas.SaveAs(Form("Map_corr_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
        IsoMap_corr[izt]->Write();
        IsoMap_corr[izt]->Draw("COLZ");
        canvas.SaveAs(Form("IsoMap_corr_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
        AntiIsoMap_corr[izt]->Write();
        AntiIsoMap_corr[izt]->Draw("COLZ");
        canvas.SaveAs(Form("AntiIsoMap_corr_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        canvas.Clear();
    }
    fout->Close();
    
    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
