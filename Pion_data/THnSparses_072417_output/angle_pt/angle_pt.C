// This macro makes an angle vs. photon pt chart of THnSparses_072417.root data
// Author: Ivan Chernyshev; Date created: 7/24/2017

#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <iostream>

//variables of hPion
const int axis_pion_Cen           = 0;
const int axis_pion_Zvtx          = 1;
const int axis_pionMass           = 2;
const int axis_pionPt             = 3;
const int axis_asymmetry          = 5;
const int axis_pionAngle          = 8;
const int axis_pionLambda1        = 9;
const int axis_pionLambda2        = 10;
const int axis_pionDisToCharged1  = 11;
const int axis_pionDisToCharged2  = 12;

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
void angle_pt() {
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Apply baseline cuts and get the data
    TFile* fIn = new TFile("THnSparses_072417.root","READ");
    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion);
    SetCut(h_Pion, axis_pionPt, 6.0, 20.0);
    SetCut(h_Pion, axis_pionDisToCharged1, 0.02, 0.14);
    SetCut(h_Pion, axis_pionDisToCharged2, 0.02, 0.14);
    SetCut(h_Pion, axis_pionLambda1, 0.1, 0.4);
    SetCut(h_Pion, axis_pionLambda2, 0.1, 0.4);
    SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
    TH2D* angle_pt_hist = h_Pion->Projection(axis_pionPt, axis_pionAngle);
    
    // Graph the histogram
    TCanvas* canvas = new TCanvas();
    angle_pt_hist->SetTitle("Opening angle-momentum relation; opening angle (mrad); transverse momentum (GeV/c)");
    angle_pt_hist->Draw("COLZ");
    myText(.35, .95, kBlack, "Opening angle-momentum relation");
    canvas->SaveAs("angle_pion_relation.png");
    
    canvas->Close();
}
