// This macro combines the efficiency charts from Cluster v1 and Cluster v2
// Author: Ivan Chernyshev; Date: February 1, 2018
#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

void efficiency_combiner() {
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    TCanvas* c = new TCanvas();
    
    // Get the data
    TFile* clusv1 = new TFile("fout_efficiency_17g6a1_pthat2_clusterv1_small.root", "READ");
    TFile* clusv2 = new TFile("fout_efficiency_17g6a1_pthat2_clusterv2_small.root", "READ");
    
    /**
     Efficiency
     */
    
    TH1D* v1_ratios = 0;
    clusv1->GetObject("ratio_photons", v1_ratios);
    TH1D* v2_ratios = 0;
    clusv2->GetObject("ratio_photons", v2_ratios);
    TH1D* v1_numerator = 0;
    clusv1->GetObject("measured_photons", v1_numerator);
    TH1D* v2_numerator = 0;
    clusv2->GetObject("measured_photons", v2_numerator);
    TH1D* v1_denominator = 0;
    clusv1->GetObject("generated_photons", v1_denominator);
    TH1D* v2_denominator = 0;
    clusv2->GetObject("generated_photons", v2_denominator);
    
    // Cut the histograms to a maximum pT of 20 GeV
    v1_ratios->GetXaxis()->SetRangeUser(0.0, 20.0);
    v2_ratios->GetXaxis()->SetRangeUser(0.0, 20.0);
    
    // Create the graph and the legend
    v1_ratios->SetMarkerColor(2);
    v1_ratios->SetMarkerStyle(33);
    v1_ratios->SetLineColor(2);
    v1_ratios->SetTitle("Efficiencies; Cluster pT (GeV); Efficiency");
    v1_ratios->GetXaxis()->SetTitleOffset(1.0);
    v1_ratios->GetYaxis()->SetTitleOffset(0.9);
    v1_ratios->Draw();
    v2_ratios->SetMarkerColor(4);
    v2_ratios->SetMarkerStyle(43);
    v2_ratios->SetLineColor(4);
    v2_ratios->Draw("same");
    
    myText(0.4, 0.92, 1, "Efficiencies");
    myMarkerText(0.5,0.8, 2, 33, "#scale[0.75]{cluster_v1}", 1);
    myMarkerText(0.5,0.7, 4, 43, "#scale[0.75]{cluster_v2}", 1);
    
    c->SaveAs("efficiencies.png");
    c->Clear();
    
    /**
     Numerator
     */
    
    // Cut the histograms to a maximum pT of 20 GeV
    v1_numerator->GetXaxis()->SetRangeUser(0.0, 20.0);
    v2_numerator->GetXaxis()->SetRangeUser(0.0, 20.0);
    
    // Create the graph and the legend
    v1_numerator->SetMarkerColor(2);
    v1_numerator->SetMarkerStyle(33);
    v1_numerator->SetLineColor(2);
    v1_numerator->SetTitle("Efficiency Numerator; Cluster pT (GeV); Detected Photons");
    v1_numerator->GetXaxis()->SetTitleOffset(1.0);
    v1_numerator->GetYaxis()->SetTitleOffset(1.1);
    v1_numerator->GetYaxis()->SetRangeUser(0, 100000);
    v1_numerator->Draw("E0");
    
    v2_numerator->SetMarkerColor(4);
    v2_numerator->SetMarkerStyle(43);
    v2_numerator->SetLineColor(4);
    v2_numerator->Draw("E0 SAME");
    
    myText(0.3, 0.92, 1, "Efficiency Numerator");
    myMarkerText(0.7,0.8, 2, 33, "#scale[0.75]{cluster_v1}", 1);
    myMarkerText(0.7,0.7, 4, 43, "#scale[0.75]{cluster_v2}", 1);
    
    c->SaveAs("efficiency_numerator.png");
    c->Clear();
    
    c->SetLogy();
    
    /**
     Denominator
     */
    
    // Cut the histograms to a maximum pT of 20 GeV
    v1_denominator->GetXaxis()->SetRangeUser(0.0, 20.0);
    v2_denominator->GetXaxis()->SetRangeUser(0.0, 20.0);
    
    // Create the graph and the legend
    v1_denominator->SetMarkerColor(2);
    v1_denominator->SetMarkerStyle(33);
    v1_denominator->SetLineColor(2);
    v1_denominator->SetTitle("Efficiency Denominator; Cluster pT (GeV); Generated Photons");
    v1_denominator->GetXaxis()->SetTitleOffset(1.0);
    v1_denominator->GetYaxis()->SetTitleOffset(1.1);
    v1_denominator->GetYaxis()->SetRangeUser(1000, 200000000);
    v1_denominator->Draw("E0");
    
    v2_denominator->SetMarkerColor(4);
    v2_denominator->SetMarkerStyle(43);
    v2_denominator->SetLineColor(4);
    v2_denominator->Draw("E0 SAME");
    
    myText(0.3, 0.92, 1, "Efficiency Denominator");
    myMarkerText(0.75,0.8, 2, 33, "#scale[0.75]{cluster_v1}", 1);
    myMarkerText(0.75,0.7, 4, 43, "#scale[0.75]{cluster_v2}", 1);
    
    c->SaveAs("efficiency_denominator.png");
    c->Clear();
    
    c->Close();
}
