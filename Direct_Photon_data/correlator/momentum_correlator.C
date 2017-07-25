// This macro graphs photon momentum/energy vs isolation and photon momentum/energy vs. UE
// Author: Ivan Chernyshev; Date: 7/24/2017

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

const int axis_photonPt                  =  3;
const int axis_photonPseudoRapidity      =  4;
const int axis_photonLambda              =  6;
const int axis_photonNcells              =  7;
const int axis_photonDisToBorder         =  9;
const int axis_photonDisToBadCell        = 10;
const int axis_photonDisToCharged        = 11;
const int axis_photonExoticity           = 14;
const int axis_photonTime                = 15;
const int axis_photonIsolation_R_track   = 18;
const int axis_photonIsolation_R_cluster = 22;
const int axis_photonUE                  = 20;


// The cutting function
void SetCut(THnSparse* h, const int axis, double min, double max){
    //make a selection on the chosen variable
    double width = h->GetAxis(axis)->GetBinWidth(1);
    int binmin = h->GetAxis(axis)->FindBin(min);
    int binmax = h->GetAxis(axis)->FindBin(max);
    h->GetAxis(axis)->SetRange(binmin, binmax - 1);
    return;
}

void momentum_correlator() {
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    TCanvas* canvas = new TCanvas();
    
    // Load isolation vs. momentum data from THnSparses_072417.root
    TFile* fIn = new TFile("THnSparses_072417.root", "READ");
    THnSparse* hPhoton = 0;
    fIn->GetObject("h_Cluster", hPhoton);
    
    // Cut the data
    SetCut(hPhoton, axis_photonPt, 6, 50);
    SetCut(hPhoton, axis_photonLambda, 0.1, 0.4);
    SetCut(hPhoton, axis_photonDisToCharged, 0.02, 0.15);
    SetCut(hPhoton, axis_photonDisToBorder, 0.9, 6);
    SetCut(hPhoton, axis_photonDisToBadCell, 1.9, 10);
    SetCut(hPhoton, axis_photonNcells, 1.9, 30.0);
    SetCut(hPhoton, axis_photonExoticity, 0, 0.97);
    SetCut(hPhoton, axis_photonTime, -30, 30);
    SetCut(hPhoton, axis_photonPseudoRapidity, -0.27, 0.27);
    
    // Graph the data
    TH2D* pt_isolation_hist = hPhoton->Projection(axis_photonIsolation_R_track, axis_photonPt);
    pt_isolation_hist->SetTitle("Isolation-energy correlation; Photon energy (GeV); #Sigma E^{cone}_{T} [GeV], track");
    gPad->SetLogz(kTRUE);
    pt_isolation_hist->Draw("COLZ");
    myText(.20, .92, kBlack, "Isolation-energy correlation");
    canvas->SaveAs("isolation_energy_correlation.png");
    
    // Repeat for UE vs. momentum, but scale each bin by 1/0.57.
    SetCut(hPhoton, axis_photonPt, 6, 25);
    SetCut(hPhoton, axis_photonUE, 0, 25);
    canvas->Clear();
    TH2D* pt_UE_hist = hPhoton->Projection(axis_photonUE, axis_photonPt);
    int xbinmin = pt_UE_hist->GetXaxis()->
    int xbinmax = pt_UE_hist->GetXaxis()->
    int ybinmin = pt_UE_hist->GetYaxis()->
    int ybinmax = pt_UE_hist->GetYaxis()->
    pt_UE_hist->SetTitle("UE-energy correlation; Photon energy (GeV); UE (GeV)");
    pt_UE_hist->Draw("COLZ");
    myText(.20, .92, kBlack, "UE-energy correlation");
    canvas->SaveAs("UE_energy_correlation.png");
    
    canvas->Close();
}
