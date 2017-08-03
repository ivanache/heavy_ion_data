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

const int num_of_bins = 6;

/**
 Error propogation functions
*/
double addition_error(double error_x, double error_y) {
    return TMath::Sqrt((error_x * error_x) + (error_y * error_y));
}
double mult_div_error(double x, double error_x, double y, double error_y, double total) {
    return total*TMath::Sqrt((error_x * error_x)/(x * x) + (error_y * error_y)/(y * y));
}
/**
 Main Function
*/
void pion_combiner() {
    // General-use storage variables
    double x_real[num_of_bins];
    double x_MC[num_of_bins];
    double y_real[num_of_bins];
    double y_MC[num_of_bins];
    
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    TCanvas* canvas = new TCanvas();
    
    // Get the data
    TFile* realIn = new TFile("PionSparsesOutput_angle_17mrad.root", "READ");
    TFile* MCIn = new TFile("BackgroundMC_PionSparsesOutput_angle_17mrad.root", "READ");
    
    TGraphErrors* realPionIntegrals = 0;
    TGraphErrors* realPionMasses = 0;
    TGraphErrors* realPionMassWidths = 0;
    realIn->GetObject("pion-integrals", realPionIntegrals);
    realIn->GetObject("mean-masses", realPionMasses);
    realIn->GetObject("standard-dev-masses", realPionMassWidths);
    
    TGraphErrors* MCPionIntegrals = 0;
    TGraphErrors* MCPionMasses = 0;
    TGraphErrors* MCPionMassWidths = 0;
    MCIn->GetObject("pion-integrals", MCPionIntegrals);
    MCIn->GetObject("mean-masses", MCPionMasses);
    MCIn->GetObject("standard-dev-masses", MCPionMassWidths);
    
    // Assign colors based on the source, then graph, insert a legend, and save
    // Num of Pions (normalized before graphing)
    // Normalization
    double total_MCPions = 0;
    double total_MCPions_error = 0;
    double total_realPions = 0;
    double total_realPions_error = 0;
    for (int i = 0; i < num_of_bins; i++) {
        MCPionIntegrals->GetPoint(i, x_MC[i], y_MC[i]);
        realPionIntegrals->GetPoint(i, x_real[i], y_real[i]);
        total_MCPions += y_MC[i];
        total_realPions += y_real[i];
        total_MCPions_error = addition_error(total_MCPions_error, MCPionIntegrals->GetErrorY(i));
        total_realPions_error = addition_error(total_realPions_error, realPionIntegrals->GetErrorY(i));
    }
    double MCPions_errors[num_of_bins];
    double realPions_errors[num_of_bins];
    for (int i = 0; i < num_of_bins; i++) {
        MCPionIntegrals->GetPoint(i, x_MC[i], y_MC[i]);
        realPionIntegrals->GetPoint(i, x_real[i], y_real[i]);
        MCPionIntegrals->SetPoint(i, x_MC[i], y_MC[i]/total_MCPions);
        realPionIntegrals->SetPoint(i, x_real[i], y_real[i]/total_realPions);
        MCPionIntegrals->SetPointError(i, MCPionIntegrals->GetErrorX(i), mult_div_error(y_MC[i], MCPionIntegrals->GetErrorY(i), total_MCPions, total_MCPions_error, y_MC[i]/total_MCPions));
        realPionIntegrals->SetPointError(i, realPionIntegrals->GetErrorX(i), mult_div_error(y_real[i], realPionIntegrals->GetErrorY(i), total_realPions, total_realPions_error, y_real[i]/total_realPions));
    }
    // Rest of procedure
    MCPionIntegrals->SetMarkerColor(kRed);
    MCPionIntegrals->SetLineColor(kRed);
    //MCPionIntegrals->GetYaxis()->SetRangeUser(0, 30000);
    //realPionIntegrals->GetYaxis()->SetRangeUser(0, 30000);
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
    
    
    // Do the same with the data to simulation ratios for the masses and mass widths
    canvas->Clear();
    TGraphErrors* PionMassRatios = new TGraphErrors(num_of_bins);
    for (int i = 0; i < num_of_bins; i++) {
        realPionMasses->GetPoint(i, x_real[i], y_real[i]);
        MCPionMasses->GetPoint(i, x_MC[i], y_MC[i]);
        PionMassRatios->SetPoint(i, x_MC[i], y_MC[i]/y_real[i]);
        PionMassRatios->SetPointError(i, realPionMasses->GetErrorX(i), (y_MC[i]/y_real[i])*TMath::Sqrt((realPionMasses->GetErrorY(i) * realPionMasses->GetErrorY(i))/(y_real[i] * y_real[i]) + (MCPionMasses->GetErrorY(i) * MCPionMasses->GetErrorY(i))/(y_MC[i] * y_MC[i])));
    }
    PionMassRatios->GetXaxis()->SetRangeUser(6,16);
    PionMassRatios->GetYaxis()->SetRangeUser(0.9, 1);
    PionMassRatios->SetTitle("Simulation/real data pion mass ratios; Transverse momentum (GeV); Simulation mass/real mass");
    PionMassRatios->Draw();
    myText(.3, .95, kBlack, "Simulation/real data pion mass ratios");
    canvas->SaveAs("PionMassRatios.png");
    
    
    // Do the same with the data to simulation ratios for the masses and mass widths
    canvas->Clear();
    TGraphErrors* PionMassWidthRatios = new TGraphErrors(num_of_bins);
    for (int i = 0; i < num_of_bins; i++) {
        realPionMassWidths->GetPoint(i, x_real[i], y_real[i]);
        MCPionMassWidths->GetPoint(i, x_MC[i], y_MC[i]);
        PionMassWidthRatios->SetPoint(i, x_MC[i], y_MC[i]/y_real[i]);
        PionMassWidthRatios->SetPointError(i, realPionMassWidths->GetErrorX(i), (y_MC[i]/y_real[i])*TMath::Sqrt((realPionMassWidths->GetErrorY(i) * realPionMassWidths->GetErrorY(i))/(y_real[i] * y_real[i]) + (MCPionMassWidths->GetErrorY(i) * MCPionMassWidths->GetErrorY(i))/(y_MC[i] * y_MC[i])));
    }
    PionMassWidthRatios->GetXaxis()->SetRangeUser(6,16);
    //PionMassWidthRatios->GetYaxis()->SetRangeUser(0.9, 1);
    PionMassWidthRatios->SetTitle("Simulation/real data pion mass ratios; Transverse momentum (GeV); Simulation mass/real mass");
    PionMassWidthRatios->Draw();
    myText(.3, .95, kBlack, "Simulation/real data pion mass ratios");
    canvas->SaveAs("PionMassWidthRatios.png");

    canvas->Clear();
    canvas->Close();
}

