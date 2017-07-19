// This macro graphs a 2-D histogram of lambda vs energy readings
// Author: Ivan Chernyshev; Date: July 18, 2017
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

//const int axis_photonEnergy       =  2;
const int axis_photonPt           =  3;
const int axis_photonPseudoRapidity= 4;
const int axis_photonLambda       =  6;
const int axis_photonNcells       =  7;
const int axis_photonDisToBorder  =  9;
const int axis_photonDisToBadCell = 10;
const int axis_photonDisToCharged = 11;
const int axis_photonExoticity    = 12;
const int axis_photonTime         = 13;
const int axis_photonIsolation    = 15;

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
 Main function
 */
void photon_Elambda() {
    // Load ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    TCanvas* canvas = new TCanvas();
    
    // Load data from THnSparses_071617.root
    TFile* fIn = new TFile("THnSparses_071617.root", "READ");
    THnSparse* hPhoton = 0;
    fIn->GetObject("h_Cluster", hPhoton);
    
    // Do the baseline cuts
    //SetCut(hPhoton, axis_photonPt, 10, 50);
    SetCut(hPhoton, axis_photonDisToCharged, 0.02, 0.15);
    SetCut(hPhoton, axis_photonDisToBorder, 0.9, 6);
    SetCut(hPhoton, axis_photonDisToBadCell, 1.9, 10);
    SetCut(hPhoton, axis_photonNcells, 1.9, 30.0);
    SetCut(hPhoton, axis_photonExoticity, 0, 0.97);
    SetCut(hPhoton, axis_photonTime, -30, 30);
    SetCut(hPhoton, axis_photonLambda, 0, 1.54);
    //SetCut(hPhoton, axis_photonPseudoRapidity, -0.27, 0.27);
    
    // Create the 2D histogram and graph
    TH2D* e_lambda_hist = hPhoton->Projection(axis_photonIsolation, axis_photonLambda);
    e_lambda_hist->GetYaxis()->SetTitle("#Sigma E^{cone}_{T} [GeV]");
    e_lambda_hist->Draw("COLZ");
    gPad->SetLogz(kTRUE);
    myText(.35, .92, kBlack, "#scale[1.5]{Cone energy vs. Lambda}");
    canvas->SaveAs("e_lambda_correlation.png");
    
    canvas->Close();
}
