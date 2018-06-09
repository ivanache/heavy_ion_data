/**
 my_code is a data-processing macro meant for processing THnSparses root files, specifically the h_Pion one in THnSparses_LHC13d_090717.root, THnSparses_LHC13e_090717.root, THnSparses_LHC13f_090717.root
 Programmer: Ivan Chernyshev
 */
#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"
#include "TAxis.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

// The cutting function
void SetCut(THnSparse* h, const int axis, double min, double max){
    //make a selection on the chosen variable
    double width = h->GetAxis(axis)->GetBinWidth(1);
    int binmin = h->GetAxis(axis)->FindBin(min);
    int binmax = h->GetAxis(axis)->FindBin(max);
    h->GetAxis(axis)->SetRange(binmin, binmax - 1);
    return;
}

// The Gaussian fit-function, for the residual distributions
double gaussian_peak(Double_t *x, Double_t *par) {
    double A = par[0];
    double mean = par[1];
    double sigma = par[2];
    
    double arg;
    if (sigma != 0)
        arg = (x[0] - mean)/sigma;
    
    double fitval = A*TMath::Exp(-0.5*arg*arg)/TMath::Sqrt(2*TMath::Pi()*sigma*sigma);
    return fitval;
}

/**
 Main function
 */
void my_code() {
    
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
   
    TCanvas* graphcanvas = new TCanvas();
    
    // Grab the file and the energy response histogram
    TFile* fout = new TFile("fout_energyresponse_18b10b_calo_pthat1to4_v1.root", "READ");
    TH1D* fitgraph = 0;
    fout->GetObject("energy_resolution", fitgraph);
    
    // Make a copy of the energy response histogram and cut it. The copy will be fitted to.
    TH1D fittedcopy = TH1D(*fitgraph);
    int binmax = fittedcopy.GetXaxis()->FindBin(0.08);
    fittedcopy.GetXaxis()->SetRange(0, binmax);
    
    // Fit and plot the Energy Response distribution
    graphcanvas->Clear();
    TF1* gaussian_fit = new TF1("fit", gaussian_peak, -5, 5, 3);
    gaussian_fit->SetLineColor(2);
    gaussian_fit->SetParameters(50,  0, 1);
    gaussian_fit->SetParLimits(0, 0, 10000);
    gaussian_fit->SetParLimits(1, -10, 10);
    gaussian_fit->SetParLimits(2, 0, 10);
    fittedcopy.Fit(gaussian_fit);
    fitgraph->Draw();
    gaussian_fit->Draw("same");
    myText(.40, .94, kBlack, "#scale[1.5]{p_{T} Response, p_{T}>10 GeV}");
    myText(.2, .75, kBlack, Form("#scale[0.75]{Mean Deviation: %1.5f+/-%1.5f}", gaussian_fit->GetParameter(1), gaussian_fit->GetParError(1)));
    myText(.2, .7, kBlack, Form("#scale[0.75]{Sigma: %1.5f+/-%1.5f}", gaussian_fit->GetParameter(2), gaussian_fit->GetParError(2)));
    graphcanvas->SaveAs("EnergyResponseFit.png");
    
    graphcanvas->Close();
    return;
}
