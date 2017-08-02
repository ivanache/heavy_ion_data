/**
   This macro takes pion fit parameter data from the Monte Carlo simulation and from the actual ALICE data and graphs them on the same graph for comparison
*/
// Author: Ivan Chernyshev
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

void pion_combiner() {
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    TCanvas* canvas = new TCanvas();
    
    // Get the data
    TFile* realIn = new TFile("PionSparsesOutput_angle_17mrad.root", "READ");
    TFile* MCIn = new TFile("BackgroundMC_PionSparsesOutput_angle_17mrad.root", "READ");
    
    TGraph* realPionIntegrals = 0;
    TGraph* realPionMasses = 0;
    TGraph* realPionMassWidths = 0;
    realIn->GetObject("pion-integrals", realPionIntegrals);
    realIn->GetObject("mean-masses", realPionMasses);
    realIn->GetObject("standard-dev-masses", realPionMassWidths);
    
    TGraph* MCPionIntegrals = 0;
    TGraph* MCPionMasses = 0;
    TGraph* MCPionMassWidths = 0;
    MCIn->GetObject("pion-integrals", MCPionIntegrals);
    MCIn->GetObject("mean-masses", MCPionMasses);
    MCIn->GetObject("standard-dev-masses", MCPionMassWidths);
    
    // Assign colors based on the source, then graph, insert a legend, and save
    // Num of Pions
    MCPionIntegrals->SetMarkerColor(kRed);
    MCPionIntegrals->SetLineColor(kRed);
    MCPionIntegrals->GetYaxis()->SetRangeUser(0, 30000);
    realPionIntegrals->GetYaxis()->SetRangeUser(0, 30000);
    realPionIntegrals->SetTitle("Number of Pions with respect to Pt; Transverse Momentum (GeV); Number of pions");
    realPionIntegrals->Draw();
    MCPionIntegrals->Draw("same");
    myText(.3, .95, kBlack, "Number of Pions with respect to Pt");
    myBoxText(0.7, 0.8, 0.05, 10, kBlack, "Real data");
    myBoxText(0.7, 0.75, 0.05, 10, kRed, "Simulation data");
    canvas->SaveAs("NumOfPionsGraph.png");

    // Mean Masses
    canvas->Clear();
    MCPionMasses->SetMarkerColor(kRed);
    MCPionMasses->SetLineColor(kRed);
    
    // Expected pion mass, 134.98 MeV
    TF1* expectedMass = new TF1("mass_pdg", "[0]", 6, 16);
    expectedMass->SetParameter(0, 0.13498);
    expectedMass->SetLineWidth(2);
    expectedMass->SetLineColor(kBlue);
    expectedMass->SetTitle("Measured Pion Mass with respect to Pt; Transverse Momentum (GeV); Mass (GeV)");
    expectedMass->GetYaxis()->SetRangeUser(0.125, 0.155);
    expectedMass->Draw();
    
    realPionMasses->Draw("same");
    MCPionMasses->Draw("same");
    myText(.3, .95, kBlack, "Measured Pion Mass with respect to Pt");
    myBoxText(0.25, 0.8, 0.05, 10, kBlack, "Real data");
    myBoxText(0.25, 0.75, 0.05, 10, kRed, "Simulation data");
    myBoxText(0.25, 0.7, 0.05, 10, kBlue, "Expected Mass, 134.98 MeV");
    canvas->SaveAs("MeanMasses.png");
    
    // Mass Widths
    canvas->Clear();
    MCPionMassWidths->SetMarkerColor(kRed);
    MCPionMassWidths->SetLineColor(kRed);
    realPionMassWidths->Draw();
    MCPionMassWidths->Draw("same");
    myText(.3, .95, kBlack, "Pion Mass Standard Deviation with respect to Pt");
    myBoxText(0.25, 0.8, 0.05, 10, kBlack, "Real data");
    myBoxText(0.25, 0.75, 0.05, 10, kRed, "Simulation data");
    canvas->SaveAs("MassWidths.png");
    
    canvas->Close();
}

