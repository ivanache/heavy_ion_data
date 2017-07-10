/**
   This macro creates momentum-photon graphs from THnSparses_070517
*/
// Author: Ivan Chernyshev

#include <unistd.h>
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
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Get the momentum-photon data file
    TFile* fIn = new TFile("THnSparses_070517.root", "READ");
    THnSparse* hPhoton = 0;
    fIn->GetObject("h_Cluster", hPhoton);
    TCanvas* canvas = new TCanvas();
    TCanvas* logcanvas = new TCanvas();
    //logcanvas->SetLogy();
    
    // Do the baseline cuts
    SetCut(hPhoton, axis_photonPt, 4.0, 20.0);
    SetCut(hPhoton, axis_photonDisToCharged, 0.02, 0.14);
    SetCut(hPhoton, axis_photonNcells, 1.9, 30.0);
    SetCut(hPhoton, axis_photonLambda, 0.1, 0.4);
    
    // Graph the whole sample
    logcanvas->cd();
    TH1D* hMomentum = hPhoton->Projection(axis_photonPt);
    hMomentum->SetAxisRange(0.0, 1400.0, "Y");
    hMomentum->SetTitle("Pt-Photon Plot (all data); Pt (GeV/c); Number of Entries");
    hMomentum->Draw();
    myText(.20,.95, kBlack, "Pt-Photon Plot (all data)");

    logcanvas->SaveAs("pt_photon.png");
    logcanvas->Clear();
    
    // For reference, also do a distance to bad cells-photon graph and a distance to border graph
    TH1D* hParam = hPhoton->Projection(axis_photonDisToBadCell);
    hParam->Draw();
    logcanvas->SaveAs("DisToBadCell_photon.png");
    logcanvas->Clear();
    
    hParam = hPhoton->Projection(axis_photonDisToBorder);
    hParam->Draw();
    logcanvas->SaveAs("DisToBorder_photon.png");
    logcanvas->Clear();
    
    //Cut Distance to bad cells to greater than 3, then graph the momentum photon data again
    SetCut(hPhoton, axis_photonDisToBadCell, 4, 10);
    TH1D* hMomentum_DisToBadCells_upper = hPhoton->Projection(axis_photonPt);
    hMomentum_DisToBadCells_upper->SetTitle("Pt-Photon Plot (DisToBadCell > 3); Pt (GeV/c); Fraction of total");
    // Normalize the plot
    hMomentum_DisToBadCells_upper->Scale(1/hMomentum_DisToBadCells_upper->Integral());
    hMomentum_DisToBadCells_upper->GetYaxis()->SetRangeUser(0.0, 0.09);
    hMomentum_DisToBadCells_upper->Draw();
    myText(.30,.95, kBlack, "Pt-Photon Plot for distance to bad cell > 3");
    logcanvas->SaveAs("pt_photon_distobadcell>3.png");
    
    const int num_of_distances = 3;
    int distances_to_bad_cells[num_of_distances] = {1, 2, 3};
    Color_t colors[num_of_distances] = {kRed, kBlue, kGreen};
    //Cut Distance to bad cells to 1, 2, and 3, and draw the resulting histogram for each amount of cuts
    logcanvas->cd();
    TH1D* hMomentum_DisToBadCells_lower[num_of_distances];
    for (int i = 0; i < num_of_distances; i++) {
        SetCut(hPhoton, axis_photonDisToBadCell, distances_to_bad_cells[i], distances_to_bad_cells[i] + 1);
        hMomentum_DisToBadCells_lower[i] = hPhoton->Projection(axis_photonPt);
        hMomentum_DisToBadCells_lower[i]->SetMarkerColor(colors[i]);
        hMomentum_DisToBadCells_lower[i]->SetTitle("Pt-Photon Plot (DisToBadCell = 1, 2, 3); Pt (GeV/c); Fraction of total");
        // Normalize the plot for the lower distance to bad cells
        hMomentum_DisToBadCells_lower[i]->Scale(1/hMomentum_DisToBadCells_lower[i]->Integral());
        //hMomentum_DisToBadCells_lower[i]->GetYaxis()->SetRangeUser(0.0, 0.5);
        hMomentum_DisToBadCells_lower[i]->Draw("same");
    }
    myMarkerText(0.42, 0.87, colors[0], 20, "Fraction of total photons, DisToBadCell=1", 1);
    myMarkerText(0.42, 0.82, colors[1], 20, "Fraction of total photons, DisToBadCell=2", 1);
    myMarkerText(0.42, 0.77, colors[2], 20, "Fraction of total photons, DisToBadCell=3", 1);
    myMarkerText(0.42, 0.72, kBlack, 20, "Fraction of total photons. DisToBadCell>3", 1);
    logcanvas->SaveAs("pt_photon_distobadcell.png");
    logcanvas->Clear();
    
    // Get a histogram of the ratio of the normalized photon count with a distance to bad cells over three
    // to the normalized photon count with a distance to bad cells of of 1, 2, and 3
    canvas->cd();
    for (int i = 0; i < num_of_distances; i++) {
        TH1D* hRatio = (TH1D*) hMomentum_DisToBadCells_upper->Clone("ratio");
        for (int j = 0; j < hMomentum_DisToBadCells_upper->GetSize(); j++) { // Set the hRatio histogram's values to the photon count ratios
            double new_bin_content = (hMomentum_DisToBadCells_lower[i]->GetBinContent(j))/(hMomentum_DisToBadCells_upper->GetBinContent(j));
            hRatio->SetBinContent(j, new_bin_content);
            hRatio->SetBinError(j, new_bin_content*TMath::Sqrt( ( ((hMomentum_DisToBadCells_lower[i]->GetBinError(j))/hMomentum_DisToBadCells_lower[i]->GetBinContent(j)) * ((hMomentum_DisToBadCells_lower[i]->GetBinError(j))/hMomentum_DisToBadCells_lower[i]->GetBinContent(j)) ) + ( ((hMomentum_DisToBadCells_upper->GetBinError(j))/hMomentum_DisToBadCells_upper->GetBinContent(j)) * ((hMomentum_DisToBadCells_upper->GetBinError(j))/hMomentum_DisToBadCells_upper->GetBinContent(j)) ) ));
        }
        hRatio->SetMarkerColor(colors[i]);
        // Draw the histogram
        hRatio->SetTitle(Form("Normalized Ratios of DisToBadCell>3 to DisToBadCell=%i photons; Pt (GeV/c); Ratio", distances_to_bad_cells[i]));
        hRatio->SetAxisRange(0.0, 2.0, "Y");
        if(i == 0)
            hRatio->Draw();
        else
            hRatio->Draw("same");
    }
    myText(.05,.95, kBlack, "Normalized Ratios of photons in DisToBadCell>3 to photons in DisToBadCell 1,2,3");
    myMarkerText(0.19, 0.28, colors[0], 20, "Normalized photon ratio of DisToBadCell=1 to DisToBadCell>3", 1);
    myMarkerText(0.19, 0.23, colors[1], 20, "Normalized photon ratio of DisToBadCell=2 to DisToBadCell>3", 1);
    myMarkerText(0.19, 0.18, colors[2], 20, "Normalized photon ratio of DisToBadCell=3 to DisToBadCell>3", 1);
    canvas->SaveAs("photon_ratios_distobadcells.png");
    canvas->Clear();
    
    // Cut distance to border to greater than 2, and graph the result
    logcanvas->cd();
    SetCut(hPhoton, axis_photonDisToBadCell, 0, 10); // Reset distance to bad cell cut to original value
    SetCut(hPhoton, axis_photonDisToBorder, 3, 5);
    TH1D* hMomentum_DisToBorder_upper = hPhoton->Projection(axis_photonPt);
    hMomentum_DisToBorder_upper->SetTitle("Pt-Photon Plot (DisToBorder > 2); Pt (GeV/c); Fraction of total");
    // Normalize the plot
    hMomentum_DisToBorder_upper->Scale(1/hMomentum_DisToBorder_upper->Integral());
    hMomentum_DisToBorder_upper->GetYaxis()->SetRangeUser(0.0, 0.09);
    hMomentum_DisToBorder_upper->Draw();
    myText(.30,.95, kBlack, "Pt-Photon Plot for distance to border > 2");
    logcanvas->SaveAs("pt_photon_distoborder>2.png");

    //Cut Distance to border to 0, 1, and 2, and draw the resulting histogram for each amount of cuts
    TH1D* hMomentum_DisToBorder_lower[num_of_distances];
    int distances_to_border[num_of_distances] = {0, 1, 2};
    for (int i = 0; i < num_of_distances; i++) {
        SetCut(hPhoton, axis_photonDisToBorder, distances_to_border[i], distances_to_border[i] + 1);
        hMomentum_DisToBorder_lower[i] = hPhoton->Projection(axis_photonPt);
        hMomentum_DisToBorder_lower[i]->SetMarkerColor(colors[i]);
        hMomentum_DisToBorder_lower[i]->SetTitle("Pt-Photon Plot (DisToBorder = 0, 1, 2); Pt (GeV/c); Fraction of total");
        // Normalize the plot for the lower distance to bad cells
        hMomentum_DisToBorder_lower[i]->Scale(1/hMomentum_DisToBorder_lower[i]->Integral());
        //hMomentum_DisToBorder_lower[i]->GetYaxis()->SetRangeUser(0.0, 0.5);
        hMomentum_DisToBorder_lower[i]->Draw("same");
    }
    myMarkerText(0.42, 0.87, colors[0], 20, "Fraction of total photons, DisToBorder=0", 1);
    myMarkerText(0.42, 0.82, colors[1], 20, "Fraction of total photons, DisToBorder=1", 1);
    myMarkerText(0.42, 0.77, colors[2], 20, "Fraction of total photons, DisToBorder=2", 1);
    myMarkerText(0.42, 0.72, kBlack, 20, "Fraction of total photons. DisToBorder>2", 1);
    logcanvas->SaveAs("pt_photon_distoborder.png");
    logcanvas->Clear();
    
    // Get a histogram of the ratio of the normalized photon count with a distance to border over two
    // to the normalized photon count with a distance to bad cells of of 0, 1, and 2
    canvas->cd();
    for (int i = 0; i < num_of_distances; i++) {
        TH1D* hRatio = (TH1D*) hMomentum_DisToBorder_upper->Clone("ratio");
        for (int j = 0; j < hMomentum_DisToBorder_upper->GetSize(); j++) { // Set the hRatio histogram's values to the photon count ratios
            double new_bin_content = (hMomentum_DisToBorder_lower[i]->GetBinContent(j))/(hMomentum_DisToBorder_upper->GetBinContent(j));
            hRatio->SetBinContent(j, new_bin_content);
            hRatio->SetBinError(j, new_bin_content*TMath::Sqrt( ( ((hMomentum_DisToBorder_lower[i]->GetBinError(j))/hMomentum_DisToBorder_lower[i]->GetBinContent(j)) * ((hMomentum_DisToBorder_lower[i]->GetBinError(j))/hMomentum_DisToBorder_lower[i]->GetBinContent(j)) ) + ( ((hMomentum_DisToBorder_upper->GetBinError(j))/hMomentum_DisToBorder_upper->GetBinContent(j)) * ((hMomentum_DisToBorder_upper->GetBinError(j))/hMomentum_DisToBorder_upper->GetBinContent(j)) ) ));
        }
        hRatio->SetMarkerColor(colors[i]);
        // Draw the histogram
        hRatio->SetTitle(Form("Normalized Ratios of DisToBorder>2 to DisToBorder=%i photons; Pt (GeV/c); Ratio", distances_to_border[i]));
        hRatio->SetAxisRange(0.0, 2.5, "Y");
        if(i == 0)
            hRatio->Draw();
        else
            hRatio->Draw("same");
    }
    myText(.05,.95, kBlack, "Normalized Ratios of photons in DisToBorder>2 to photons in DisToBorder 0,1,2");
    myMarkerText(0.19, 0.28, colors[0], 20, "Normalized photon ratio of DisToBorder=0 to DisToBorder>2", 1);
    myMarkerText(0.19, 0.23, colors[1], 20, "Normalized photon ratio of DisToBorder=1 to DisToBorder>2", 1);
    myMarkerText(0.19, 0.18, colors[2], 20, "Normalized photon ratio of DisToBorder=2 to DisToBorder>2", 1);
    canvas->SaveAs("photon_ratios_distoborder.png");
    canvas->Clear();
    
    logcanvas->Close();
    canvas->Close();
}
