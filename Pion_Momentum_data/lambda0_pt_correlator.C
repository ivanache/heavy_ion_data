/**
   This macro graphs lambda0 vs momentum correlation from the pion data THnSparses
   Author: Ivan Chernyshev; Date: 6/28/2017
*/
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TCanvas.h>
#include <iostream>
#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

//variables of hPion
const int axis_pion_Cen           = 0;
const int axis_pion_Zvtx          = 1;
const int axis_pionMass           = 2;
const int axis_pionPt             = 3;
const int axis_photon1E           = 7;
const int axis_photon2E           = 8;
const int axis_asymmetry          = 9;
const int axis_pionAngle          = 16;
const int axis_pionLambda1        = 17;
const int axis_pionLambda2        = 18;
const int axis_pionNcells1        = 19;
const int axis_pionNcells2        = 20;
const int axis_pionMatchedTracks1 = 21;
const int axis_pionMatchedTracks2 = 22;
const int axis_pionDisToBorder1   = 27;
const int axis_pionDisToBorder2   = 28;

// The cutting function
void SetCut(THnSparse* h, const int axis, double min, double max){
    //make a selection on the chosen variable
    double width = h->GetAxis(axis)->GetBinWidth(1);
    int binmin = h->GetAxis(axis)->FindBin(min);
    int binmax = h->GetAxis(axis)->FindBin(max);
    h->GetAxis(axis)->SetRange(binmin, binmax - 1);
    return;
}

/**
 Main Function
*/
void lambda0_pt_correlator() {
    TCanvas* canvas = new TCanvas();
    
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Load the data from the THnSparses and apply appropriate data cuts
    TFile* fIn = new TFile("THnSparses_062717.root", "READ");
    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion", h_Pion);
    SetCut(h_Pion, axis_pionMatchedTracks1, -1.5, -0.75);
    SetCut(h_Pion, axis_pionMatchedTracks2, -1.5, -0.75);
    SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
    SetCut(h_Pion, axis_pionAngle, 0.015, 0.5);
    SetCut(h_Pion, axis_pionNcells1, 1.0, 30.0);
    SetCut(h_Pion, axis_pionNcells2, 1.0, 30.0);
    SetCut(h_Pion, axis_pionDisToBorder1, 2.0, 6.0);
    SetCut(h_Pion, axis_pionDisToBorder1, 2.0, 6.0);
    TH2D* pt_lambda = h_Pion->Projection(axis_pionPt, axis_pionLambda1);
    
    // Plot the data
    pt_lambda->SetTitle("Pi0 transverse momentum vs leading photon lambda0; lambda0; Transverse Momentum (GeV/c)");
    pt_lambda->GetYaxis()->SetRangeUser(6.0, 50.0);
    pt_lambda->Draw("COLZ");
    myText(0.20, 0.95, kBlack, "Pi0 transverse momentum vs leading photon lambda0");
    canvas->SaveAs("Pt_vs_Lambda_correlation.png");
    
    // Now, import the fit from the mass data and repeat what you did above for various intervals away fro mthe mean
    TFile* fFit = new TFile("Pion5CutsSparsesOutput.root", "READ");
    TH1D* massdata = 0;
    fFit->GetObject("mass_pion", massdata);
    TF1* model = (TF1*) massdata->GetListOfFunctions()->FindObject("mass peak");
    double mean = model->GetParameter(1);
    double sigma = model->GetParameter(2);
    
    int num_of_intervals = 4;
    for (int i = 1; i <= num_of_intervals; i++) {
        cout << "\nLower bound: " << mean - i*sigma << "\nUpper bound: " << mean + i*sigma << std::endl;
        SetCut(h_Pion, axis_pionMass, mean - i*sigma, mean + i*sigma);
        pt_lambda = h_Pion->Projection(axis_pionPt, axis_pionLambda1);
        
        // Plot the data
        pt_lambda->SetTitle("Pi0 transverse momentum vs leading photon lambda0 (within 1 sigma of mean mass); lambda0; Transverse Momentum (GeV/c)");
        pt_lambda->GetYaxis()->SetRangeUser(6.0, 50.0);
        pt_lambda->Draw("COLZ");
        myText(0.03, 0.95, kBlack, Form("Pi0 transverse momentum vs leading photon lambda0 (within 1 sigma of mean mass)", i));
        canvas->SaveAs(Form("Pt_vs_Lambda_correlation_%isigma.png", i));

    }
    canvas->Close();
}
