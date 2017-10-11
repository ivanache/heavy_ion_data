// This macro takes the d, e, and f graphs from all three runs of the LHC: d, e, and f on one canvas\
// Intended to be used after all iterations of my_code are done
// Author: Ivan Chernyshev; Date: 10/10/17

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

void sumplots() {
    // Load files from runs
    TFile* d = new TFile("PionDataOutput_LHC13d.root", "READ");
    TFile* e = new TFile("PionDataOutput_LHC13e.root", "READ");
    TFile* f = new TFile("PionDataOutput_LHC13f.root", "READ");
    
    // Load data from each file
    // Num of pions
    TGraphErrors* d_num_pions = 0;
    d->GetObject("pion-integrals", d_num_pions);
    TGraphErrors* e_num_pions = 0;
    e->GetObject("pion-integrals", e_num_pions);
    TGraphErrors* f_num_pions = 0;
    f->GetObject("pion-integrals", f_num_pions);
    //Mean masses
    TGraphErrors* d_mean_masses = 0;
    d->GetObject("mean-masses", d_mean_masses);
    TGraphErrors* e_mean_masses = 0;
    e->GetObject("mean-masses", e_mean_masses);
    TGraphErrors* f_mean_masses = 0;
    f->GetObject("mean-masses", f_mean_masses);
    //Mass peak widths
    TGraphErrors* d_mass_widths = 0;
    d->GetObject("standard-dev-masses", d_mass_widths);
    TGraphErrors* e_mass_widths = 0;
    e->GetObject("standard-dev-masses", e_mass_widths);
    TGraphErrors* f_mass_widths = 0;
    f->GetObject("standard-dev-masses", f_mass_widths);
    
    // Create canvas and set ATLAS style
    TCanvas* canvas = new TCanvas();
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Draw the number of pion curves for all 3 runs. Be sure to include a legend
    d_num_pions->GetYaxis()->SetTitle("Number of Entries");
    d_num_pions->GetYaxis()->SetRangeUser(0, 9000);
    d_num_pions->GetYaxis()->SetTitleSize(.07);
    d_num_pions->GetYaxis()->SetTitleOffset(1.1);
    d_num_pions->GetXaxis()->SetTitleOffset(1.0);
    d_num_pions->SetMarkerColor(kBlack);
    d_num_pions->SetLineColor(kBlack);
    d_num_pions->SetMarkerStyle(1);
    d_num_pions->Draw();
    
    e_num_pions->GetYaxis()->SetTitle("Number of Entries");
    e_num_pions->GetYaxis()->SetTitleSize(.07);
    e_num_pions->GetYaxis()->SetTitleOffset(1.1);
    e_num_pions->GetXaxis()->SetTitleOffset(1.0);
    e_num_pions->SetMarkerColor(kRed);
    e_num_pions->SetLineColor(kRed);
    e_num_pions->Draw("same");
    
    f_num_pions->GetYaxis()->SetTitle("Number of Entries");
    f_num_pions->GetYaxis()->SetTitleSize(.07);
    f_num_pions->GetYaxis()->SetTitleOffset(1.1);
    f_num_pions->GetXaxis()->SetTitleOffset(1.0);
    f_num_pions->SetMarkerColor(kBlue);
    f_num_pions->SetLineColor(kBlue);
    f_num_pions->Draw("same");
    
    myText(0.20, 0.92, kBlack, "Number of Pions over momentum");
    myBoxText(0.65, 0.80, 0.05, 10, kBlack, "LHC13d run");
    myBoxText(0.65, 0.75, 0.05, 10, kRed, "LHC13e run");
    myBoxText(0.65, 0.70, 0.05, 10, kBlue, "LHC13f run");
    canvas->SaveAs("Num_of_entries.png");
    canvas->Clear();
    
    // Draw the mean mass curves for all 3 runs. Be sure to include a legend
    d_mean_masses->GetYaxis()->SetTitle("Mean Masses");
    d_mean_masses->GetYaxis()->SetRangeUser(0.13, 0.155);
    d_mean_masses->GetXaxis()->SetRangeUser(6, 16);
    d_mean_masses->GetYaxis()->SetTitleSize(.06);
    d_mean_masses->GetYaxis()->SetTitleOffset(1.3);
    d_mean_masses->GetXaxis()->SetTitleOffset(1.1);
    d_mean_masses->SetMarkerColor(kBlack);
    d_mean_masses->SetLineColor(kBlack);
    d_mean_masses->SetMarkerStyle(1);
    d_mean_masses->Draw();
    
    e_mean_masses->GetYaxis()->SetTitle("Mean Masses");
    e_mean_masses->GetYaxis()->SetTitleSize(.06);
    e_mean_masses->GetYaxis()->SetTitleOffset(1.3);
    e_mean_masses->GetXaxis()->SetTitleOffset(1.1);
    e_mean_masses->SetMarkerColor(kRed);
    e_mean_masses->SetLineColor(kRed);
    e_mean_masses->Draw("same");
    
    f_mean_masses->GetYaxis()->SetTitle("Mean Masses");
    f_mean_masses->GetYaxis()->SetTitleSize(.06);
    f_mean_masses->GetYaxis()->SetTitleOffset(1.3);
    f_mean_masses->GetXaxis()->SetTitleOffset(1.1);
    f_mean_masses->SetMarkerColor(kBlue);
    f_mean_masses->SetLineColor(kBlue);
    f_mean_masses->Draw("same");
    
    myText(0.20, 0.92, kBlack, "Mean Masses over momentum");
    myBoxText(0.65, 0.30, 0.05, 10, kBlack, "LHC13d run");
    myBoxText(0.65, 0.25, 0.05, 10, kRed, "LHC13e run");
    myBoxText(0.65, 0.20, 0.05, 10, kBlue, "LHC13f run");
    canvas->SaveAs("Mean_masses.png");
    canvas->Clear();
    
    // Draw the mass widths for all 3 runs. Be sure to include a legend
    d_mass_widths->GetYaxis()->SetTitle("Mass Peak Std. Dev");
    d_mass_widths->GetYaxis()->SetTitleSize(.07);
    d_mass_widths->GetYaxis()->SetTitleOffset(0.7);
    d_mass_widths->GetXaxis()->SetTitleOffset(1.1);
    d_mass_widths->SetMarkerColor(kBlack);
    d_mass_widths->SetLineColor(kBlack);
    d_mass_widths->SetMarkerStyle(1);
    d_mass_widths->Draw();
    
    e_mass_widths->GetYaxis()->SetTitle("Mass Peak Std. Dev");
    e_mass_widths->GetYaxis()->SetTitleSize(.07);
    e_mass_widths->GetYaxis()->SetTitleOffset(0.7);
    e_mass_widths->GetXaxis()->SetTitleOffset(1.1);
    e_mass_widths->SetMarkerColor(kRed);
    e_mass_widths->SetLineColor(kRed);
    e_mass_widths->Draw("same");
    
    f_mass_widths->GetYaxis()->SetTitle("Mass Peak Std. Dev");
    f_mass_widths->GetYaxis()->SetTitleSize(.07);
    f_mass_widths->GetYaxis()->SetTitleOffset(0.7);
    f_mass_widths->GetXaxis()->SetTitleOffset(1.1);
    f_mass_widths->SetMarkerColor(kBlue);
    f_mass_widths->SetLineColor(kBlue);
    f_mass_widths->Draw("same");
    
    myText(0.20, 0.92, kBlack, "Mass Peak Widths over momentum");
    myBoxText(0.65, 0.40, 0.05, 10, kBlack, "LHC13d run");
    myBoxText(0.65, 0.35, 0.05, 10, kRed, "LHC13e run");
    myBoxText(0.65, 0.30, 0.05, 10, kBlue, "LHC13f run");
    canvas->SaveAs("mass_widths.png");
    canvas->Clear();

    canvas->Close();
}
