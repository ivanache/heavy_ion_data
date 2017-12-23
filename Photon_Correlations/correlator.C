// This program takes the output from the 13e_clusv1_small.root NTuple and creates photon-hadron correlations from it
// Programmer: Ivan Chernyshev; Date created: 12/21/2017

#include <TFile.h>
#include <TCanvas.h>
#include <TROOT.h>

#include "TH1F.h"

#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

// 1D histogram dividing function
// Precondition: the two histograms must have the same x- and y- dimensions
TH1F* divide_correlations(TH1F* dividend. TH1F* divisor) {
    // Create the quotient-containing histogram and find the minimum and maximum bins along the x-axis
    TH1F* quotient = new TH1F(*dividend);
    double x_bin_min = quotient->GetXaxis()->FindBin(quotient->GetXaxis()->GetXmin());
    double x_bin_max = quotient->GetXaxis()->FindBin(quotient->GetXaxis()->GetXmax());
    
    // Loop over all bins, divide the element from the dividend histogram by the element from the divisor histogram
    for (int i = x_bin_min; i <= x_bin_max; i++) {
        if (divisor->GetBinContent(i) != 0)
            quotient->SetBinContent(i, (dividend->GetBinContent(i))/(divisor->GetBinContent(i)));
    // Failsafe
        else
            quotient->SetBinContent(i, 0);
    
        // Error propogation
        if (divisor->GetBinContent(i) == 0) {
            if (dividend->GetBinContent(i) == 0)
                quotient->SetBinError(i, 0);
            else
                quotient->SetBinError(i, (quotient->GetBinContent(i)*TMath::Sqrt(( (dividend->GetBinError(i)/dividend->GetBinContent(i)) * (dividend->GetBinError(i)/dividend->GetBinContent(i)) ) )));
        }
        else if (dividend->GetBinContent(i) == 0)
            quotient->SetBinError(i, (quotient->GetBinContent(i)*TMath::Sqrt(( (divisor->GetBinError(i)/divisor->GetBinContent(i)) * (divisor->GetBinError(i)/divisor->GetBinContent(i)) ) )));
        else
            quotient->SetBinError(i, (quotient->GetBinContent(i)*TMath::Sqrt(( (dividend->GetBinError(i)/dividend->GetBinContent(i)) * (dividend->GetBinError(i)/dividend->GetBinContent(i)) ) + ( (divisor->GetBinError(i)/divisor->GetBinContent(i)) * (divisor->GetBinError(i)/divisor->GetBinContent(i)) ) )));

    }
    
    return quotient;
}

void correlator() {
    // Open the NTuple output file
    TFile* correlation_source = new TFile("fout_13eclusv1.root", "READ");
    
    // Get Atlas Style and define the Canvas
    TCanvas* c = new TCanvas();
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Zt bins
    const int nztbins = 7;
    const float ztbins[nztbins][2] = {{0.0, 0.1}, {0.1, 0.2} , {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 1.0}, {1.0, 1.2}};
    
    // Declare correlation functions
    TH1F* dPhi_correlations[nztbins];
    
    // For each zt bin, create the correlations by dividing the within-isolation signal by non-isolation background
    for (int zti = 0; zti < nztbins; zti++) {
        // Get the isolation signal and non-isolation background
        TH1F* dPhi_iso = 0;
        correlation_source->GetObject(Form("dPhi_iso_ztmin%i_ztmax%i", (int) ztbins[zti][0]*10, (int) ztbins[zti][1]*10), dPhi_iso);
        TH1F* dPhi_noniso = 0;
        correlation_source->GetObject(Form("dPhi_noniso_ztmin%i_ztmax%i", (int) ztbins[zti][0]*10, (int) ztbins[zti][1]*10), dPhi_iso);
        
        
    }
}
