// This file's main job is to remove the obstructive headers in front of the 2D efficiency maps

#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

#include "TH1F.h"


void twoD_efficiency_plotter() {
    // Open the NTuple output files
    TFile* efficiency_sourcev1 = new TFile("fout_efficiency_17g6a1_pthat2_clusterv1_small.root", "READ");
    TFile* efficiency_sourcev2 = new TFile("fout_efficiency_17g6a1_pthat2_clusterv2_small.root", "READ");
    
    // Get Atlas Style and define the Canvas
    TCanvas* c = new TCanvas();
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Get the files for each source
    TH2D* measured_v1 = 0;
    efficiency_sourcev1->GetObject("phi_eta_map_measured", measured_v1);
    TH2D* generated_v1 = 0;
    efficiency_sourcev1->GetObject("phi_eta_map_generated", generated_v1);
    TH2D* ratio_v1 = 0;
    efficiency_sourcev1->GetObject("phi_eta_map_ratio", ratio_v1);
    TH2D* measured_v2 = 0;
    efficiency_sourcev2->GetObject("phi_eta_map_measured", measured_v2);
    TH2D* generated_v2 = 0;
    efficiency_sourcev2->GetObject("phi_eta_map_generated", generated_v2);
    TH2D* ratio_v2 = 0;
    efficiency_sourcev2->GetObject("phi_eta_map_ratio", ratio_v2);
    
    // Plot
    measured_v1->SetTitle("Measured Photons; #phi (#frac{rad}{#pi}); #eta");
    measured_v1->GetXaxis()->SetTitleOffset(1.2);
    measured_v1->GetYaxis()->SetTitleOffset(1.2);
    measured_v1->Draw("COLZ");
    myText(0.4, 0.92, 1, "Measured Photons");
    c->SaveAs("measured_photons_phietamap__17g6a1_pthat2_clusterv1_small.png");
    c->Clear();
    
    generated_v1->SetTitle("Generated Photons; #phi (#frac{rad}{#pi}); #eta");
    generated_v1->GetXaxis()->SetTitleOffset(1.2);
    generated_v1->GetYaxis()->SetTitleOffset(1.2);
    generated_v1->Draw("COLZ");
    myText(0.4, 0.92, 1, "Generated Photons");
    c->SaveAs("generated_photons_phietamap__17g6a1_pthat2_clusterv1_small.png");
    c->Clear();
    
    ratio_v1->SetTitle("Efficiency; #phi (#frac{rad}{#pi}); #eta");
    ratio_v1->GetXaxis()->SetTitleOffset(1.2);
    ratio_v1->GetYaxis()->SetTitleOffset(1.2);
    ratio_v1->Draw("COLZ");
    myText(0.4, 0.92, 1, "Efficiency");
    c->SaveAs("efficiency_phietamap__17g6a1_pthat2_clusterv1_small.png");
    c->Clear();
    
    measured_v2->SetTitle("Measured Photons; #phi (#frac{rad}{#pi}); #eta");
    measured_v2->GetXaxis()->SetTitleOffset(1.2);
    measured_v2->GetYaxis()->SetTitleOffset(1.2);
    measured_v2->Draw("COLZ");
    myText(0.4, 0.92, 1, "Measured Photons");
    c->SaveAs("measured_photons_phietamap__17g6a1_pthat2_clusterv2_small.png");
    c->Clear();
    
    generated_v2->SetTitle("Generated Photons; #phi (#frac{rad}{#pi}); #eta");
    generated_v2->GetXaxis()->SetTitleOffset(1.2);
    generated_v2->GetYaxis()->SetTitleOffset(1.2);
    generated_v2->Draw("COLZ");
    myText(0.4, 0.92, 1, "Generated Photons");
    c->SaveAs("generated_photons_phietamap__17g6a1_pthat2_clusterv2_small.png");
    c->Clear();
    
    ratio_v2->SetTitle("Efficiency; #phi (#frac{rad}{#pi}); #eta");
    ratio_v2->GetXaxis()->SetTitleOffset(1.2);
    ratio_v2->GetYaxis()->SetTitleOffset(1.2);
    ratio_v2->Draw("COLZ");
    myText(0.4, 0.92, 1, "Efficiency");
    c->SaveAs("efficiency_phietamap__17g6a1_pthat2_clusterv2_small.png");
    c->Clear();
    
    c->Close();
}
