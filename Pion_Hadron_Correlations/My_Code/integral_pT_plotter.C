// This macro takes pion-hadron correlations from Pi0_Hadron_Corr_Output.root (Created by pion_hadron_corr) and plots the integrals of the two Gaussian peaks over track pT bin
// Author: Ivan Chernyshev; Date: 11/10/17

#include "TFile.h"

// Main Function
void integral_pT_plotter() {
    // Get the data
    TFile* correlations = new TFile("Pi0_Hadron_Corr_Output.root", "READ");
}
