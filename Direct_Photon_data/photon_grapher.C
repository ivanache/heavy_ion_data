/**
   This macro creates momentum-photon graphs from THnSparses_070517
*/
// Author: Ivan Chernyshev

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TCanvas.h>
#include <iostream>

const int axis_photonPt           =  3;
const int axis_photonLambda       =  6;
const int axis_photonNcells       =  7;
const int axis_photonDisToBorder  = 10;
const int axis_photonDisToBadCell = 11;
const int axis_photonDisToCharged = 12;

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
void photon_grapher() {
    // Get the momentum-photon data file
    TFile* fIn = new TFile("THnSparses_070517.root", "READ");
    THnSparse* hPhoton = 0;
    fIn->GetObject("h_Cluster", hPhoton);
    TCanvas* canvas = new TCanvas();
    TCanvas* logcanvas = new TCanvas();
    //logcanvas->SetLogy();
    
    // Do the baseline cuts
    SetCut(hPhoton, axis_photonDisToCharged, 0.02, 0.14);
    SetCut(hPhoton, axis_photonNcells, 1.9, 30.0);
    SetCut(hPhoton, axis_photonDisToBorder, 2.0, 6.0);
    SetCut(hPhoton, axis_photonLambda, 0.1, 0.4);
    
    // Graph the whole sample
    logcanvas->cd();
    TH1D* hMomentum = hPhoton->Projection(axis_photonPt);
    hMomentum->SetAxisRange(0.0, 1400.0, "Y");
    hMomentum->SetTitle("Pt-Photon Plot (all data); Pt (GeV/c); Number of Entries");
    hMomentum->Draw();
    logcanvas->SaveAs("pt_photon.png");
    logcanvas->Clear();
    
    // For reference, also do a distance to bad cells-photon graph
    TH1D* hParam = hPhoton->Projection(axis_photonDisToBadCell);
    hParam->Draw();
    logcanvas->SaveAs("DisToBadCell_photon.png");
    logcanvas->Clear();
    
    //Cut Distance to bad cells to less than 3, then graph the momentum photon data again
    SetCut(hPhoton, axis_photonDisToBadCell, 0, 3);
    TH1D* hMomentum_DisToBadCells_lessthan_3 = hPhoton->Projection(axis_photonPt);
    hMomentum_DisToBadCells_lessthan_3->SetAxisRange(0.0, 1400.0, "Y");
    hMomentum_DisToBadCells_lessthan_3->SetTitle("Pt-Photon Plot (DisToBadCell < 3); Pt (GeV/c); Number of Entries");
    hMomentum_DisToBadCells_lessthan_3->Draw();
    logcanvas->SaveAs("pt_photon_distobadcell<3.png");
    logcanvas->Clear();
    
    //Cut Distance to bad cells to at least 3, then graph the momentum photon data again
    SetCut(hPhoton, axis_photonDisToBadCell, 3, 10);
    TH1D* hMomentum_DisToBadCells_atleast_3 = hPhoton->Projection(axis_photonPt);
    hMomentum_DisToBadCells_atleast_3->SetAxisRange(0.0, 1400.0, "Y");
    hMomentum_DisToBadCells_atleast_3->SetTitle("Pt-Photon Plot (DisToBadCell >= 3); Pt (GeV/c); Number of Entries");
    hMomentum_DisToBadCells_atleast_3->Draw();
    logcanvas->SaveAs("pt_photon_distobadcell>=3.png");
    logcanvas->Clear();
    
    // Get a histogram of the ratio of the photon count with a distance to bad cells of at least three
    // to the photon count with a distance to bad cells of less than three
    canvas->cd();
    TH1D* hRatio = (TH1D*) hMomentum->Clone("ratio");
    for (int i = 0; i < hMomentum->GetSize(); i++) { // Set the hRatio histogram's values to the photon count ratios
        if ((hMomentum_DisToBadCells_lessthan_3->GetBinContent(i)) != 0) {
            hRatio->SetBinContent(i, (hMomentum_DisToBadCells_atleast_3->GetBinContent(i))/(hMomentum_DisToBadCells_lessthan_3->GetBinContent(i)));
        }
        hRatio->SetBinError(i, 0.05); // Not a real error, just a placeholder
    }
    // Draw the histogram
    hRatio->SetTitle("Ratios of DisToBadCell>=3 to DisToBadCell<3 photons; Pt (GeV/c); Ratio");
    hRatio->SetAxisRange(0.0, 8.0, "Y");
    hRatio->Draw();
    canvas->SaveAs("photon_ratio.png");
    
    logcanvas->Close();
    canvas->Close();
}
