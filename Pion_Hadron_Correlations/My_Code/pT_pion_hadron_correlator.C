// Creates a correlation function and projection out of a .root data file and a pion pT range that is one of the following: (6 GeV, 8 GeV), (8 GeV, 10 GeV), (10 GeV, 12 GeV), (12 GeV, 14 GeV), (14 GeV, 16 GeV), and a track pT range over which the data in the .root file is defined
// Does so by using the results obtained from mass_pion_modeler (which should have been called prior) to get the desired pion-hadron correlation from pion_hadron_corr
// Programmer: Ivan Chernyshev; Date 10/25/17

// Precodition: PionDataOutput.root exists in the same directory as this macro and comes from mass_pion_modeller.C

#include "TFile.h"
#include <TGraphErrors.h>

// Main function; all variables are in GeV
void pT_pion_hadron_correlator(double pion_pT_min, double pion_pT_max, double track_pT_min, double track_pT_max) {
    
    // Get the mass range from the .root file that mass_pion_modeller produced
    TFile* read_data = new TFile("PionDataOutput.root", "READ");
    TGraphErrors* masses_over_pT = 0;
    read_data->GetObject("mean-masses", masses_over_pT);
    TGraphErrors* masswidths_over_pT = 0;
    read_data->GetObject("standard-dev-masses", masswidths_over_pT);
    double pion_PT_center = (pion_pT_min + pion_pT_max)/2;
    double mass_center = masses_over_pT->Eval(pion_PT_center);
    double mass_width = masswidths_over_pT->Eval(pion_PT_center);
    double mass_min = mass_center - 2*mass_width;
    double mass_max = mass_center + 2*mass_width;
    
    // Get the directory name, form it from min and max pT bounds
    string filename = Form("TrackpT_%1.0f-%1.0fGeV/PionpT_%2.0f-%2.0fGeV", track_pT_min, track_pT_max, pion_pT_min, pion_pT_max);
    std::cout << "Filename: " << filename << std::endl;
    
    // Run pion_hadron_corr with the results
    string processcommand = std::string(Form(".x pion_hadron_corr.C(%5.5f, %5.5f, %5.5f, %5.5f, %5.5f, %5.5f, ", pion_pT_min, pion_pT_max, mass_min, mass_max, track_pT_min, track_pT_max)) + "\"" + filename + "\")";
    gROOT->ProcessLine(processcommand.c_str());
}
