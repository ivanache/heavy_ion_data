// This program fuses results from 3 different data samples
// Author: Ivan Chernyshev


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

// Adds the corresponding elements of the 3 input histograms
// Precondition: the histograms must have the same x-dimensions
TH1D* addition_histograms1D(TH1D* graph1, int size1, TH1D* graph2, int size2, TH1D* graph3, int size3){
    // Make the 1D histogram to contain the result and find minimum and maximum bins along both axes
    TH1D* result = new TH1D(*graph1);
    double x_bin_min = result->GetXaxis()->FindBin(result->GetXaxis()->GetXmin());
    double x_bin_max = result->GetXaxis()->FindBin(result->GetXaxis()->GetXmax());
    int totalsize = size1 + size2 + size3;
    
    // Loop over all bins, set the result's bins to the sum. Manually propogate the error
    for(int i = x_bin_min; i <= x_bin_max; i++){
        double val1 = graph1->GetBinContent(i);
        double val2 = graph2->GetBinContent(i);
        double val3 = graph3->GetBinContent(i);
        result->SetBinContent(i, ((val1*size1) + (val2*size2) + (val3*size3))/totalsize);
        
        // Error propogation
        double delta1_scaled = ((graph1->GetBinError(i))*size1)/totalsize;
        double delta2_scaled = ((graph2->GetBinError(i))*size2)/totalsize;
        double delta3_scaled = ((graph3->GetBinError(i))*size3)/totalsize;
        double deltares = TMath::Sqrt((delta1_scaled*delta1_scaled) + (delta2_scaled*delta2_scaled) + (delta3_scaled*delta3_scaled));
        result->SetBinError(i, deltares);
    }
    
    
    return result;
}

void rootresultcombiner(std::string sample1, int size1_sig, int size1_bkg, std::string sample2, int size2_sig, int size2_bkg, std::string sample3, int size3_sig, int size3_bkg, std::string resname) {
    // Get the 3 files
    TFile *file1 = new TFile(sample1.c_str(), "READ");
    TFile *file2 = new TFile(sample2.c_str(), "READ");
    TFile *file3 = new TFile(sample3.c_str(), "READ");
    
    // Get the histograms
    TH1D* dPhiSignal1 = 0;
    TH1D* dPhiBackground1 = 0;
    file1->GetObject("sig_dPhi", dPhiSignal1);
    file1->GetObject("bkg_dPhi", dPhiBackground1);
    TH1D* dPhiSignal2 = 0;
    TH1D* dPhiBackground2 = 0;
    file2->GetObject("sig_dPhi", dPhiSignal2);
    file2->GetObject("bkg_dPhi", dPhiBackground2);
    TH1D* dPhiSignal3 = 0;
    TH1D* dPhiBackground3 = 0;
    file3->GetObject("sig_dPhi", dPhiSignal3);
    file3->GetObject("bkg_dPhi", dPhiBackground3);
    
    // Add the dPhi plots
    TH1D* dPhiSignal_tot = addition_histograms1D(dPhiSignal1, size1_sig, dPhiSignal2, size2_sig, dPhiSignal3, size3_sig);
    TH1D* dPhiBackground_tot = addition_histograms1D(dPhiBackground1, size1_bkg, dPhiBackground2, size2_bkg, dPhiBackground3, size3_bkg);
    
    TFile* fout = new TFile(resname.c_str(), "RECREATE");
    dPhiSignal_tot->Write();
    dPhiBackground_tot->Write();
    
    file1->Close();
    file2->Close();
    file3->Close();
    fout->Close();
}
