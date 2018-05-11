// This program is a combined same-event and mixed-event correlation program that has been given updates from Fernando Torales-Acosta's code
// Author: Ivan Chernyshev
// Date: 5/2/2018

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
TH1D* divide_histograms1D(TH1D* graph1, TH1D* graph2, const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup){
    // Make the 1D histogram to contain the quotient and find minimum and maximum bins along both axes
    TH1D* quotient = new TH1D(name, title, nbinsx, xlow, xup);
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

// A function for dividing a 2D histogram by a scalar
// Precondition: the scalar cannot be 0 and shouldn't have an error
void divide_histograms2D_byscalar(TH2D* graph, double scalar){
    // Find minimum and maximum bins along both axes
    double x_bin_min = graph->GetXaxis()->FindBin(graph->GetXaxis()->GetXmin());
    double x_bin_max = graph->GetXaxis()->FindBin(graph->GetXaxis()->GetXmax());
    double y_bin_min = graph->GetYaxis()->FindBin(graph->GetYaxis()->GetXmin());
    double y_bin_max = graph->GetYaxis()->FindBin(graph->GetYaxis()->GetXmax());
    
    // Loop over all bins, divide the element from graph 1 by its counterpart on graph 2. Manually propogate the error
    for(int i = x_bin_min; i <= x_bin_max; i++)
        for(int j = y_bin_min; j <= y_bin_max; j++) {
            if (scalar != 0)
                graph->SetBinContent(i, j, (graph->GetBinContent(i, j))/scalar);
            // Failsafe
            else
                graph->SetBinContent(i, j, 0);
            
            // Error propogation
            graph->SetBinError(i, j, graph->GetBinError(i, j)/scalar);
        }
    
}

// A function for dividing a 1D histogram by a scalar
// Precondition: the scalar cannot be 0 and shouldn't have an error
void divide_histograms1D_byscalar(TH1D* graph, double scalar){
    // Find minimum and maximum bins along both axes
    double x_bin_min = graph->GetXaxis()->FindBin(graph->GetXaxis()->GetXmin());
    double x_bin_max = graph->GetXaxis()->FindBin(graph->GetXaxis()->GetXmax());
    
    // Loop over all bins, divide the element from graph 1 by its counterpart on graph 2. Manually propogate the error
    for(int i = x_bin_min; i <= x_bin_max; i++){
        if (scalar != 0)
            graph->SetBinContent(i, (graph->GetBinContent(i))/(scalar));
        // Failsafe
        else
            graph->SetBinContent(i, 0);
        
        // Error propogation
        graph->SetBinError(i, graph->GetBinError(i)/scalar);
    }
    
}

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

using namespace H5;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }
    int dummyc = 1;
    char **dummyv = new char *[1];
    
    dummyv[0] = strdup("main");
    
    //Config File
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
    isolationDet determiner = CLUSTER_ISO_ITS_04; //replaced by config file. Check on Print
    int n_eta_bins = 0;
    int n_phi_bins = 0;
    
    // Zt bins
    //FIXME: Will have to likely set nztbins first, then initialize array
    int nztbins = 7;
    float* ztbins;
    ztbins = new float[nztbins+1];
    ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
    
    int nptbins = 3;
    float* ptbins;
    ptbins = new float[nptbins+1];
    ptbins[0] = 10.0; ptbins[1] = 11; ptbins[2] = 12.5; ptbins[3] = 16;
    
    
    // Loop through config file
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
    
    //TH1D Mixed_h_ntrig("h_ntrig_MixedEvent", "", 2, -0.5,1.0);
    
    //TH2D* Mixed_Signal_pT_Dist = new TH2D("Signal_pT_Dist_MixedEvent","Mixed Event Cluster Pt Spectrum For Isolation (its_04) bins 0.55 < DNN < 0.85",24,10,16,5,-0.5,2);
    //TH2D* Mixed_BKGD_pT_Dist = new TH2D("BKGD_pT_Dist_MixedEvent","Mixed Event Cluster Pt Spectrum For Isolation (its_04) bins 0.0 < DNN < 0.3",24,10,16,5,-0.5,2);
    
    // Same event
    TH2D* Same_Corr[nztbins*nptbins];
    TH2D* Same_IsoCorr[nztbins*nptbins];
    TH2D* Same_BKGD_IsoCorr[nztbins*nptbins];
    
    TH1D* Same_dPhi_Corr[nztbins*nptbins];
    TH1D* Same_dPhi_IsoCorr[nztbins*nptbins];
    TH1D* Same_dPhi_BKGD_IsoCorr[nztbins*nptbins];
    
    float Same_N_Total_Triggers = 0;
    float Same_N_Signal_Triggers = 0;
    float Same_N_BKGD_Triggers = 0;
    
    // Mixed Event
    TH2D* Mixed_Corr[nztbins*nptbins];
    TH2D* Mixed_IsoCorr[nztbins*nptbins];
    TH2D* Mixed_BKGD_IsoCorr[nztbins*nptbins];
    
    TH1D* Mixed_dPhi_Corr[nztbins*nptbins];
    TH1D* Mixed_dPhi_IsoCorr[nztbins*nptbins];
    TH1D* Mixed_dPhi_BKGD_IsoCorr[nztbins*nptbins];
    
    //Overall Correlations
    TH2D* Corr[nztbins*nptbins];
    TH2D* IsoCorr[nztbins*nptbins];
    TH2D* BKGD_IsoCorr[nztbins*nptbins];
    
    TH1D* dPhi_Corr[nztbins*nptbins];
    TH1D* dPhi_IsoCorr[nztbins*nptbins];
    TH1D* dPhi_BKGD_IsoCorr[nztbins*nptbins];
    
    //TH1D* Mixed_H_Signal_Triggers[nptbins];
    //TH1D* Mixed_H_BKGD_Triggers[nptbins];
    //float Mixed_N_Signal_Triggers = 0;
    //float Mixed_N_BKGD_Triggers = 0;
    
    //FIXME: Add to config file
    
    for (int ipt = 0; ipt <nptbins; ipt++) {
        //Mixed_H_Signal_Triggers[ipt] = new TH1D(Form("N_DNN%i_Triggers_MixedEvent_pT%1.0f_%1.0f",1,ptbins[ipt],ptbins[ipt+1]), "Number of Isolated Photon Triggers", 2, -0.5,1.0);
        
        //Mixed_H_BKGD_Triggers[ipt] = new TH1D(Form("N_DNN%i_Triggers_MixedEvent_pT%1.0f_%1.0f",2,ptbins[ipt],ptbins[ipt+1]),"Number of Isolated Low DNN Photon Triggers", 2, -0.5,1.0);
        
        for (int izt = 0; izt<nztbins; izt++){
            // 2D Histograms
            Same_Corr[izt+ipt*nztbins] = new TH2D(Form("Correlation_SameEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#gamma-H [all] [SameEvent] Correlation,  %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; #Delta #eta; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
            
            Same_Corr[izt+ipt*nztbins]->Sumw2();
            Same_Corr[izt+ipt*nztbins]->SetMinimum(0.);
            
            Same_IsoCorr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation_SameEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f", 1,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#gamma-H [Iso] [SameEvent] Correlation,  %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; #Delta #eta; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2,n_eta_bins, -1.4, 1.4);
            
            Same_IsoCorr[izt+ipt*nztbins]->Sumw2();
            Same_IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);
            
            Same_BKGD_IsoCorr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation_SameEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#gamma-H [AntiIso] [SameEvent] Correlation, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; #Delta #eta; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
            
            Same_BKGD_IsoCorr[izt+ipt*nztbins]->Sumw2();
            Same_BKGD_IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);
            
            Mixed_Corr[izt+ipt*nztbins] = new TH2D(Form("Correlation__pT%1.0f_%1.0f_MixedEvent_zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#gamma-H [all] [Mixed-Event] Correlation,  %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; #Delta #eta; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
            
            Mixed_Corr[izt+ipt*nztbins]->Sumw2();
            Mixed_Corr[izt+ipt*nztbins]->SetMinimum(0.);
            
            Mixed_IsoCorr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",1,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#gamma-H [Iso] [Mixed-Event] Correlation,  %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; #Delta #eta; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2,n_eta_bins, -1.4, 1.4);
            
            Mixed_IsoCorr[izt+ipt*nztbins]->Sumw2();
            Mixed_IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);
            
            Mixed_BKGD_IsoCorr[izt+ipt*nztbins] = new TH2D(Form("DNN%i_Correlation_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#gamma-H [AntiIso] [Mixed-Event] Correlation,  %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; #Delta #eta; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
            
            Mixed_BKGD_IsoCorr[izt+ipt*nztbins]->Sumw2();
            Mixed_BKGD_IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);
            
            // 1D Histograms
            Same_dPhi_Corr[izt+ipt*nztbins] = new TH1D(Form("DeltaPhiMap_SameEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#Delta #phi [all] [SameEvent] Map, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2);
            
            Same_dPhi_Corr[izt+ipt*nztbins]->Sumw2();
            Same_dPhi_Corr[izt+ipt*nztbins]->SetMinimum(0.);
            
            Same_dPhi_IsoCorr[izt+ipt*nztbins] = new TH1D(Form("DNN%i_DeltaPhiMap_SameEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",1,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#Delta #phi [Iso] [SameEvent] Map, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2);
            
            Same_dPhi_IsoCorr[izt+ipt*nztbins]->Sumw2();
            Same_dPhi_IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);
            
            Same_dPhi_BKGD_IsoCorr[izt+ipt*nztbins] = new TH1D(Form("DNN%i_DeltaPhiMap_SameEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#Delta #phi [AntiIso] [SameEvent] Map, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2);
            
            Same_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->Sumw2();
            Same_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);
            
            Mixed_dPhi_Corr[izt+ipt*nztbins] = new TH1D(Form("DeltaPhiMap_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#Delta #phi [all] [MixedEvent] Map, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2);
            
            Mixed_dPhi_Corr[izt+ipt*nztbins]->Sumw2();
            Mixed_dPhi_Corr[izt+ipt*nztbins]->SetMinimum(0.);
            
            Mixed_dPhi_IsoCorr[izt+ipt*nztbins] = new TH1D(Form("DNN%i_DeltaPhiMap_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",1,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#Delta #phi [Iso] [MixedEvent] Map, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2);
            
            Mixed_dPhi_IsoCorr[izt+ipt*nztbins]->Sumw2();
            Mixed_dPhi_IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);
            
            Mixed_dPhi_BKGD_IsoCorr[izt+ipt*nztbins] = new TH1D(Form("DNN%i_DeltaPhiMap_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#Delta #phi [AntiIso] [MixedEvent] Map, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2);
            
            Mixed_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->Sumw2();
            Mixed_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->SetMinimum(0.);
        }//zt bins
    }//pt bins
    
    
    //LOOP OVER SAMPLES
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
            _tree_event = dynamic_cast<TTree *>(file->Get("AliAnalysisTaskNTGJ/_tree_event"));
            if (_tree_event == NULL) {
                std::cout << " fail " << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
        //variables
        Double_t primary_vertex[3];
        UInt_t ntrack;
        Float_t track_e[NTRACK_MAX];
        Float_t track_pt[NTRACK_MAX];
        Float_t track_eta[NTRACK_MAX];
        Float_t track_phi[NTRACK_MAX];
        UChar_t track_quality[NTRACK_MAX];
        Float_t track_eta_emcal[NTRACK_MAX];
        Float_t track_phi_emcal[NTRACK_MAX];
        UChar_t track_its_ncluster[NTRACK_MAX];
        Float_t track_its_chi_square[NTRACK_MAX];
        Float_t track_dca_xy[NTRACK_MAX];
        
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
        
        UInt_t ntrack_max = 0;
        UInt_t ncluster_max = 0;
        
        fprintf(stderr, "\r%s:%d: %s\n", __FILE__, __LINE__, "Determining ntrack_max and ncluster_max needed for hdf5 hyperslab");
        for (Long64_t i = 0; i < _tree_event->GetEntries(); i++) {
            _tree_event->GetEntry(i);
            ncluster_max = std::max(ncluster_max, ncluster);
            ntrack_max = std::max(ntrack_max, ntrack);
            fprintf(stderr, "\r%s:%d: %llu", __FILE__, __LINE__, i);
        }
        ntrack_max = 414;
        ncluster_max = 23;
        
        //    UInt_t ntrack_max = 414;
        //    UInt_t ncluster_max = 23;
        
        fprintf(stderr, "\n%s:%d: %s", __FILE__, __LINE__, "USING HARDCODED HDF5 DIMENSIONS");
        fprintf(stderr, "\n%s:%d: maximum tracks:%i maximum clusters:%i\n", __FILE__, __LINE__, ntrack_max,ncluster_max);
        
        //open hdf5: Define size of data from file, explicitly allocate memory in hdf5 space and array size
        const H5std_string hdf5_file_name(argv[iarg+1]);
        
        const H5std_string track_ds_name( "track" );
        H5File h5_file( hdf5_file_name, H5F_ACC_RDONLY );
        DataSet track_dataset = h5_file.openDataSet( track_ds_name );
        DataSpace track_dataspace = track_dataset.getSpace();
        
        const H5std_string cluster_ds_name( "cluster" );
        DataSet cluster_dataset = h5_file.openDataSet( cluster_ds_name );
        DataSpace cluster_dataspace = cluster_dataset.getSpace();
        
        //Define array hyperslab will be read into
        float track_data_out[1][ntrack_max][10];
        float cluster_data_out[1][ncluster_max][5];
        
        //Define hyperslab size and offset in  FILE;
        hsize_t track_offset[3] = {0, 0, 0};
        hsize_t track_count[3] = {1, ntrack_max, 10};
        hsize_t cluster_offset[3] = {0, 0, 0};
        hsize_t cluster_count[3] = {1, ncluster_max, 5};
        
        track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
        cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "select Hyperslab OK");
        
        //Define the memory dataspace to place hyperslab
        const int RANK_OUT = 3; //# of Dimensions
        hsize_t track_dimsm[3] = {1, ntrack_max, 10};
        DataSpace track_memspace( RANK_OUT, track_dimsm );
        hsize_t cluster_dimsm[3] = {1, ncluster_max, 5};
        DataSpace cluster_memspace( RANK_OUT, cluster_dimsm );
        
        //Define memory offset for hypreslab starting at begining:
        hsize_t track_offset_out[3] = {0};
        hsize_t cluster_offset_out[3] = {0};
        
        //define Dimensions of array, for writing slab to array
        hsize_t track_count_out[3] = {1, ntrack_max, 10};
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
            
            for (ULong64_t n = 0; n < ncluster; n++) {
                if( not(cluster_pt[n]>pT_min and cluster_pt[n]<pT_max)) continue; //select pt of photons
                if( not(TMath::Abs(cluster_eta[n])<Eta_max)) continue; //cut edges of detector
                if( not(cluster_ncell[n]>Cluster_min)) continue;   //removes clusters with 1 or 2 cells
                if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters
                
                float isolation;
                if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n];
                else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n];
                else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n];
                else isolation = cluster_frixione_its_04_02[n];
                if (isolation>iso_max) continue;
                /*
                for (Long64_t imix = 0; imix < 50; imix++){
                    Long64_t mix_event = Mix_Events[imix];
                    if (mix_event == ievent) continue;
                    
                    //adjust offset for next mixed event
                    track_offset[0]=mix_event;
                    track_dataspace.selectHyperslab( H5S_SELECT_SET, track_count, track_offset );
                    track_dataset.read( track_data_out, PredType::NATIVE_FLOAT, track_memspace, track_dataspace );
                    
                    cluster_offset[0]=mix_event;
                    cluster_dataspace.selectHyperslab( H5S_SELECT_SET, cluster_count, cluster_offset );
                    cluster_dataset.read( cluster_data_out, PredType::NATIVE_FLOAT, cluster_memspace, cluster_dataspace );
                    
                    //MIXED associated
                    const int TrackCutBit =16;
                    for (ULong64_t itrack = 0; itrack < ntrack_max; itrack++) {
                        if (std::isnan(track_data_out[0][itrack][1])) break;
                        //if ((int(track_data_out[0][itrack][4]+0.5)&selection_number)==0) continue;
                        if ((int(track_data_out[0][itrack][4]+ 0.5)&TrackCutBit)==0) continue; //selection 16
                        if (track_data_out[0][itrack][1] < 0.5) continue; //less than 1GeV
                        if (track_data_out[0][itrack][1] > 30) continue; //less than 1GeV
                        if (abs(track_data_out[0][itrack][2]) > 0.8) continue;
                        if (track_data_out[0][itrack][7] < 4) continue;
                        if ((track_data_out[0][itrack][8]/track_data_out[0][itrack][7]) > 36) continue;
                        if( not(TMath::Abs(track_data_out[0][itrack][9])<0.0231+0.0315/TMath::Power(track_data_out[0][itrack][4],1.3 ))) continue;
                        
                        double dRmin = 0.02;
                        //veto charged particles from mixed event tracks
                        bool MixTrack_HasMatch = false;
                        for (unsigned int l = 0; l < ncluster_max; l++){
                            if (std::isnan(cluster_data_out[0][l][0])) break;
                            float dphi = TMath::Abs(cluster_data_out[0][l][2] - track_data_out[0][itrack][5]);
                            float deta = TMath::Abs(cluster_data_out[0][l][3] - track_data_out[0][itrack][6]);
                            float dR = sqrt(dphi*dphi + deta*deta);
                            if (dR < dRmin)    MixTrack_HasMatch = true;
                            break;
                        }
                        if (MixTrack_HasMatch) continue;
                        
                        //fprintf(stderr, "%s:%d: Mixed Event: %llu Track: %llu\n", __FILE__, __LINE__, mix_event, itrack);
                        
                        Float_t DeltaPhi = cluster_phi[n] - track_data_out[0][itrack][3];
                        if (DeltaPhi < -M_PI/2){DeltaPhi += 2*M_PI;}  //if less then -pi/2 add 2pi
                        if (DeltaPhi > 3*M_PI/2){DeltaPhi =DeltaPhi -2*M_PI;}
                        Float_t DeltaEta = cluster_eta[n] - track_data_out[0][itrack][2];
                        if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue;
                        
                        Double_t zt = track_data_out[0][itrack][1]/cluster_pt[n];
                        
                        for (int ipt = 0; ipt < nptbins; ipt++){
                            if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1]){
                                for(int izt = 0; izt<nztbins ; izt++){
                                    if(zt>ztbins[izt] and  zt<ztbins[izt+1]){
                                        
                                        if(isolation< iso_max){
                                            Mixed_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
                                            Mixed_dPhi_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi);
                                        }
                                        if (cluster_s_nphoton[n][1] > 0.0 && cluster_s_nphoton[n][1] < 0.3){ //sel deep photons
                                            Mixed_BKGD_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
                                            Mixed_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi);
                                        }
                                        
                                        Mixed_Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
                                        Mixed_dPhi_Corr[izt+ipt*nztbins]->Fill(DeltaPhi);
                                    }//if in zt bin
                                } // end loop over zt bins
                            }//end if in pt bin
                        }//end pt loop bin
                    }//end loop over tracks
                }//end loop over mixed events
                */
                
                // Now same events
                Same_N_Total_Triggers += 1;
                
                //High DNN Trigger Signal
                if ((cluster_s_nphoton[n][1] > DNN_min) && (cluster_s_nphoton[n][1]<DNN_max)){
                    Same_N_Signal_Triggers += 1;
                }
                //Low DNN Trigger BKGD
                if ((cluster_s_nphoton[n][1]>0.0) && (cluster_s_nphoton[n][1]<0.3)){
                    Same_N_BKGD_Triggers += 1;
                }
                
                //Track Loop
                for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {
                    if(track_pt[itrack] < 0.5) continue; //500 MeV Tracks
                    if(track_pt[itrack] > 30) continue;
                    if((track_quality[itrack]&Track_Cut_Bit)==0) continue; //select only tracks that pass selection
                    if(abs(track_eta[itrack]) > 0.8) continue;
                    if( not(track_its_ncluster[itrack]>4)) continue;
                    if( not(track_its_chi_square[itrack]/track_its_ncluster[itrack] <36)) continue;
                    if( not(TMath::Abs(track_dca_xy[itrack])<0.0231+0.0315/TMath::Power(track_pt[itrack],1.3 ))) continue;
                    
                    //Electron Veto for associated tracks outside of isolation cone
                    double dRmin = 0.02;
                    bool Track_HasMatch = false;
                    for (ULong64_t c = 0; c < ncluster; c++){
                        Float_t deta =  cluster_eta[n]-track_eta_emcal[itrack];
                        Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[c]-track_phi_emcal[itrack])/TMath::Pi();
                        float dR = sqrt(dphi*dphi + deta*deta);
                        if (dR < dRmin) Track_HasMatch = true;
                        break;
                    }
                    if (Track_HasMatch) continue;
                    
                    //Observables:
                    Double_t zt = track_pt[itrack]/cluster_pt[n];
                    Float_t DeltaPhi = cluster_phi[n] - track_phi[itrack];
                    if (DeltaPhi < -M_PI/2) DeltaPhi += 2*M_PI;
                    if (DeltaPhi > 3*M_PI/2) DeltaPhi =DeltaPhi -2*M_PI;
                    Float_t DeltaEta = cluster_eta[n] - track_eta[itrack];
                    if ((TMath::Abs(DeltaPhi) < 0.005) && (TMath::Abs(DeltaEta) < 0.005)) continue; //Match Mixing Cut
                    
                    //fprintf(stderr,"%f   ",track_pt[itrack]);
                    
                    for (int ipt = 0; ipt < nptbins; ipt++){
                        if (cluster_pt[n] >ptbins[ipt] && cluster_pt[n] <ptbins[ipt+1]){
                            for(int izt = 0; izt<nztbins ; izt++){
                                if(zt>ztbins[izt] and  zt<ztbins[izt+1]){
                                    
                                    //2 DNN Regions
                                    if (cluster_s_nphoton[n][1] > DNN_min && cluster_s_nphoton[n][1] < DNN_max) {
                                        Same_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
                                        Same_dPhi_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi);
                                    }
                                    if (cluster_s_nphoton[n][1] > 0.0 && cluster_s_nphoton[n][1] < 0.3) {//sel deep photons
                                        Same_BKGD_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
                                        Same_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->Fill(DeltaPhi);
                                    }
                                    Same_Corr[izt+ipt*nztbins]->Fill(DeltaPhi,DeltaEta);
                                    Same_dPhi_Corr[izt+ipt*nztbins]->Fill(DeltaPhi);
                                }//if in zt bin
                            } // end loop over zt bins
                        }//end if in pt bin
                    }//end pt loop bin
                }//end loop over tracks
                
            }//end loop on clusters.
            
        } //end loop over events
        
    }//end loop over samples
    
    // Normalize all same event correlation functions with number of triggers
    for (int ipt = 0; ipt<nptbins; ipt++){
        for (int izt = 0; izt<nztbins; izt++) {
            divide_histograms2D_byscalar(Same_IsoCorr[izt+ipt*nztbins], Same_N_Signal_Triggers);
            divide_histograms1D_byscalar(Same_dPhi_IsoCorr[izt+ipt*nztbins], Same_N_Signal_Triggers);
            
            divide_histograms2D_byscalar(Same_BKGD_IsoCorr[izt+ipt*nztbins], Same_N_BKGD_Triggers);
            divide_histograms1D_byscalar(Same_dPhi_BKGD_IsoCorr[izt+ipt*nztbins], Same_N_BKGD_Triggers);
            
            divide_histograms2D_byscalar(Same_Corr[izt+ipt*nztbins], Same_N_Total_Triggers);
            divide_histograms1D_byscalar(Same_dPhi_Corr[izt+ipt*nztbins], Same_N_Total_Triggers);
        }
    }
    
    // Normalize all correlation mixed event functions with the origin
    for (int ipt = 0; ipt<nptbins; ipt++){
        for (int izt = 0; izt<nztbins; izt++) {
            divide_histograms1D_byscalar(Mixed_dPhi_Corr[izt+ipt*nztbins], Mixed_dPhi_Corr[izt+ipt*nztbins]->GetBinContent(Mixed_dPhi_Corr[izt+ipt*nztbins]->GetXaxis()->FindBin(0.0)));
            divide_histograms1D_byscalar(Mixed_dPhi_IsoCorr[izt+ipt*nztbins], Mixed_dPhi_IsoCorr[izt+ipt*nztbins]->GetBinContent(Mixed_dPhi_IsoCorr[izt+ipt*nztbins]->GetXaxis()->FindBin(0.0)));
            divide_histograms1D_byscalar(Mixed_dPhi_BKGD_IsoCorr[izt+ipt*nztbins], Mixed_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->GetBinContent(Mixed_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->GetXaxis()->FindBin(0.0) ));
            
            divide_histograms2D_byscalar(Mixed_Corr[izt+ipt*nztbins], Mixed_Corr[izt+ipt*nztbins]->GetBinContent(Mixed_Corr[izt+ipt*nztbins]->GetXaxis()->FindBin(0.0), Mixed_Corr[izt+ipt*nztbins]->GetYaxis()->FindBin(0.0)));
            divide_histograms2D_byscalar(Mixed_IsoCorr[izt+ipt*nztbins], Mixed_IsoCorr[izt+ipt*nztbins]->GetBinContent(Mixed_IsoCorr[izt+ipt*nztbins]->GetXaxis()->FindBin(0.0), Mixed_IsoCorr[izt+ipt*nztbins]->GetYaxis()->FindBin(0.0)));
            divide_histograms2D_byscalar(Mixed_BKGD_IsoCorr[izt+ipt*nztbins], Mixed_BKGD_IsoCorr[izt+ipt*nztbins]->GetBinContent(Mixed_BKGD_IsoCorr[izt+ipt*nztbins]->GetXaxis()->FindBin(0.0), Mixed_BKGD_IsoCorr[izt+ipt*nztbins]->GetYaxis()->FindBin(0.0)));
        }
    }
    
    // Divided same event by mixed event in order to obtain correlation functions
    for (int ipt = 0; ipt<nptbins; ipt++){
        for (int izt = 0; izt<nztbins; izt++) {
            // Divide
            Corr[izt+ipt*nztbins] = divide_histograms2D(Same_Corr[izt+ipt*nztbins], Mixed_Corr[izt+ipt*nztbins], Form("Correlation_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#gamma-H [all] Correlation,  %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; #Delta #eta; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
            IsoCorr[izt+ipt*nztbins] = divide_histograms2D(Same_IsoCorr[izt+ipt*nztbins], Mixed_IsoCorr[izt+ipt*nztbins], Form("DNN%i_Correlation_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f", 1,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#gamma-H [Iso] Correlation,  %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; #Delta #eta; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
            BKGD_IsoCorr[izt+ipt*nztbins] = divide_histograms2D(Same_BKGD_IsoCorr[izt+ipt*nztbins], Mixed_BKGD_IsoCorr[izt+ipt*nztbins], Form("DNN%i_Correlation_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]), Form("#gamma-H [AntiIso] Correlation, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; #Delta #eta; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2, n_eta_bins, -1.4, 1.4);
            
            dPhi_Corr[izt+ipt*nztbins] = divide_histograms1D(Same_dPhi_Corr[izt+ipt*nztbins], Mixed_dPhi_Corr[izt+ipt*nztbins], Form("DeltaPhiMap_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#Delta #phi [all] Map, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2);
            dPhi_IsoCorr[izt+ipt*nztbins] = divide_histograms1D(Same_dPhi_IsoCorr[izt+ipt*nztbins], Mixed_dPhi_IsoCorr[izt+ipt*nztbins], Form("DNN%i_DeltaPhiMap_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",1,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#Delta #phi [Iso] Map, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2);
            dPhi_BKGD_IsoCorr[izt+ipt*nztbins] = divide_histograms1D(Same_dPhi_BKGD_IsoCorr[izt+ipt*nztbins], Mixed_dPhi_BKGD_IsoCorr[izt+ipt*nztbins], Form("DNN%i_DeltaPhiMap_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f",2,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]),Form("#Delta #phi [AntiIso] Map, %2.1f < zt < %2.1f, %2.0f < pt < %2.0f; #Delta #phi [rad]; entries", ztbins[izt], ztbins[izt+1], ptbins[ipt], ptbins[ipt+1]), n_phi_bins,-M_PI/2,3*M_PI/2);
            
        }
    }
    
    // Write to fout
    TFile* fout = new TFile("Combined_Correlation_13def.root","RECREATE");
    for (int ipt = 0; ipt<nptbins; ipt++){
        for (int izt = 0; izt<nztbins; izt++){
            Same_Corr[izt+ipt*nztbins]->Write();
        }
        for (int izt = 0; izt<nztbins; izt++){
            Same_IsoCorr[izt+ipt*nztbins]->Write();
        }
        for (int izt = 0; izt<nztbins; izt++){
            Same_BKGD_IsoCorr[izt+ipt*nztbins]->Write();
        }
    }
    for (int ipt = 0; ipt<nptbins; ipt++){
        for (int izt = 0; izt<nztbins; izt++)
            Same_dPhi_Corr[izt+ipt*nztbins]->Write();
        
        for (int izt = 0; izt<nztbins; izt++)
            Same_dPhi_IsoCorr[izt+ipt*nztbins]->Write();
        
        for (int izt = 0; izt<nztbins; izt++)
            Same_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->Write();
    }
    
    
    for (int ipt = 0; ipt<nptbins; ipt++){
        for (int izt = 0; izt<nztbins; izt++){
            Mixed_Corr[izt+ipt*nztbins]->Write();
        }
        for (int izt = 0; izt<nztbins; izt++){
            Mixed_IsoCorr[izt+ipt*nztbins]->Write();
        }
        for (int izt = 0; izt<nztbins; izt++){
            Mixed_BKGD_IsoCorr[izt+ipt*nztbins]->Write();
        }
    }
    for (int ipt = 0; ipt<nptbins; ipt++){
        for (int izt = 0; izt<nztbins; izt++)
            Mixed_dPhi_Corr[izt+ipt*nztbins]->Write();
        
        for (int izt = 0; izt<nztbins; izt++)
            Mixed_dPhi_IsoCorr[izt+ipt*nztbins]->Write();
        
        for (int izt = 0; izt<nztbins; izt++)
            Mixed_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->Write();
    }
    
    for (int ipt = 0; ipt<nptbins; ipt++){
        for (int izt = 0; izt<nztbins; izt++){
            Corr[izt+ipt*nztbins]->Write();
        }
        for (int izt = 0; izt<nztbins; izt++){
            IsoCorr[izt+ipt*nztbins]->Write();
        }
        for (int izt = 0; izt<nztbins; izt++){
            BKGD_IsoCorr[izt+ipt*nztbins]->Write();
        }
    }
    
    for (int ipt = 0; ipt<nptbins; ipt++){
        for (int izt = 0; izt<nztbins; izt++)
            dPhi_Corr[izt+ipt*nztbins]->Write();
        
        for (int izt = 0; izt<nztbins; izt++)
            dPhi_IsoCorr[izt+ipt*nztbins]->Write();
        
        for (int izt = 0; izt<nztbins; izt++)
            dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->Write();
    }
    
    fout->Close();
    
    //Put Correlation functions into PNG files
    for (int ipt = 0; ipt<nptbins; ipt++){
        /**for (int izt = 0; izt<nztbins; izt++) {
         Mixed_Corr[izt+ipt*nztbins]->Draw("SURF2");
         canvas.SaveAs(Form("Correlation_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f.png",ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));
         canvas.Clear();
         Mixed_dPhi_Corr[izt+ipt*nztbins]->Draw();
         canvas.SaveAs(Form("DeltaPhiMap_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f.png",ptbins[ipt],ptbins[ipt+1],10*ztbins[izt],10*ztbins[izt+1]));
         canvas.Clear();
         
         Mixed_IsoCorr[izt+ipt*nztbins]->Draw("SURF2");
         canvas.SaveAs(Form("DNN%i_Correlation_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f.png", 1,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]));
         canvas.Clear();
         Mixed_dPhi_IsoCorr[izt+ipt*nztbins]->Draw();
         canvas.SaveAs(Form("DNN%i_DeltaPhiMap_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f.png",1,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]));
         canvas.Clear();
         
         Mixed_BKGD_IsoCorr[izt+ipt*nztbins]->Draw("SURF2");
         canvas.SaveAs(Form("DNN%i_Correlation_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f.png",2,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]));
         canvas.Clear();
         Mixed_dPhi_BKGD_IsoCorr[izt+ipt*nztbins]->Draw();
         canvas.SaveAs(Form("DNN%i_DeltaPhiMap_MixedEvent_pT%1.0f_%1.0f__zT%1.0f_zT%1.0f.png",2,ptbins[ipt],ptbins[ipt+1], 10*ztbins[izt],10*ztbins[izt+1]));
         canvas.Clear();
         }*/
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Same_Corr[izt+ipt*nztbins]->Draw("SURF2");
        }
        canvas.SaveAs(Form("Correlation_SameEvent_pT%1.0f_%1.0f.png",ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Same_dPhi_Corr[izt+ipt*nztbins]->Draw();
        }
        canvas.SaveAs(Form("DeltaPhiMap_SameEvent_pT%1.0f_%1.0f.png",ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Same_IsoCorr[izt+ipt*nztbins]->Draw("SURF2");
        }
        canvas.SaveAs(Form("DNN%i_Correlation_SameEvent_pT%1.0f_%1.0f.png", 1,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Same_dPhi_Corr[izt+ipt*nztbins]->Draw();
        }
        canvas.SaveAs(Form("DNN%i_DeltaPhiMap_SameEvent_pT%1.0f_%1.0f.png",1,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Same_BKGD_IsoCorr[izt+ipt*nztbins]->Draw("SURF2");
        }
        canvas.SaveAs(Form("DNN%i_Correlation_SameEvent_pT%1.0f_%1.0f.png",2,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Same_dPhi_Corr[izt+ipt*nztbins]->Draw();
        }
        canvas.SaveAs(Form("DNN%i_DeltaPhiMap_SameEvent_pT%1.0f_%1.0f.png",2,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Mixed_Corr[izt+ipt*nztbins]->Draw("SURF2");
        }
        canvas.SaveAs(Form("Correlation_MixedEvent_pT%1.0f_%1.0f.png",ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Mixed_dPhi_Corr[izt+ipt*nztbins]->Draw();
        }
        canvas.SaveAs(Form("DeltaPhiMap_MixedEvent_pT%1.0f_%1.0f.png",ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Mixed_IsoCorr[izt+ipt*nztbins]->Draw("SURF2");
        }
        canvas.SaveAs(Form("DNN%i_Correlation_MixedEvent_pT%1.0f_%1.0f.png", 1,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Mixed_dPhi_Corr[izt+ipt*nztbins]->Draw();
        }
        canvas.SaveAs(Form("DNN%i_DeltaPhiMap_MixedEvent_pT%1.0f_%1.0f.png",1,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Mixed_BKGD_IsoCorr[izt+ipt*nztbins]->Draw("SURF2");
        }
        canvas.SaveAs(Form("DNN%i_Correlation_MixedEvent_pT%1.0f_%1.0f.png",2,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Mixed_dPhi_Corr[izt+ipt*nztbins]->Draw();
        }
        canvas.SaveAs(Form("DNN%i_DeltaPhiMap_MixedEvent_pT%1.0f_%1.0f.png",2,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            Corr[izt+ipt*nztbins]->Draw("SURF2");
        }
        canvas.SaveAs(Form("Correlation_pT%1.0f_%1.0f.png",ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            dPhi_Corr[izt+ipt*nztbins]->Draw();
        }
        canvas.SaveAs(Form("DeltaPhiMap_pT%1.0f_%1.0f.png",ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            IsoCorr[izt+ipt*nztbins]->Draw("SURF2");
        }
        canvas.SaveAs(Form("DNN%i_Correlation_pT%1.0f_%1.0f.png", 1,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            dPhi_Corr[izt+ipt*nztbins]->Draw();
        }
        canvas.SaveAs(Form("DNN%i_DeltaPhiMap_pT%1.0f_%1.0f.png",1,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            BKGD_IsoCorr[izt+ipt*nztbins]->Draw("SURF2");
        }
        canvas.SaveAs(Form("DNN%i_Correlation_pT%1.0f_%1.0f.png",2,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
        
        canvas.Divide(4, 2, 0.01, 0.01);
        for (int izt = 0; izt<nztbins; izt++) {
            canvas.cd(izt+1);
            dPhi_Corr[izt+ipt*nztbins]->Draw();
        }
        canvas.SaveAs(Form("DNN%i_DeltaPhiMap_pT%1.0f_%1.0f.png",2,ptbins[ipt],ptbins[ipt+1]));
        canvas.Clear();
    }
    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}

