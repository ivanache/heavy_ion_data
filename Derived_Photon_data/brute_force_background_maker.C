/** 
    This macro takes the background lambda0 vs Number of Photons projection for the background, given a THnSparses and a pion data root file with signal-to-total ratio and fit data
    Afterwards combines the projection for the background with the one for the whole signal, given by photon_projector.C
    Author: Ivan Chernyshev; Date created: 6/27/17
*/

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


//variables of hPion
const int axis_pion_Cen    = 0;
const int axis_pion_Zvtx   = 1;
const int axis_pionMass    = 2;
const int axis_pionPt      = 3;
const int axis_photon1E    = 7;
const int axis_photon2E    = 8;
const int axis_asymmetry   = 9;
const int axis_pionAngle   = 16;
const int axis_pionLambda1 = 17;
const int axis_pionLambda2 = 18;
const int axis_pionNcells1 = 19;
const int axis_pionNcells2 = 20;
const int axis_pionMatchedTracks1 = 21;
const int axis_pionMatchedTracks2 = 22;
const int axis_pionDisToBorder1 = 27;
const int axis_pionDisToBorder2 = 28;

// Concatenates two strings and gives a char array
char* str_concat_converter(string str1, string str2){
    string sumstring = str1 + str2;
    char* output = new char[sumstring.length() + 1];
    strcpy(output, sumstring.c_str());
    return output;
}

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
 Main Function
*/
// Precondition: NumOfSigmasFromMean = 1, 2, 3, or 4
void brute_force_background_maker(int NumOfSigmasFromMean) {
    const double lambdamin = 0;
    const double lambdamax = 2;
    const int E_min = 6;
    const int E_max = 8;
    string directory_name = Form("%isigma/projections/background", NumOfSigmasFromMean);
    TCanvas* canvas = new TCanvas();

    
    // Set ATLAS Style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Take in the THnSparses, establish the link to the PhotonsOtput root file for output
    TFile* mainIn = new TFile("THnSparses_062717.root", "READ");
    TFile* pionIn = new TFile("Pion5CutsSparsesOutput.root", "READ");
    TFile* fOut = new TFile(Form("%isigmaPhotonsOutput.root", "UPDATE"));
    
    // Extract the data
    THnSparse* hPion = 0;
    TH1D* hFittedMass = 0;
    TGraph* sigTotal = 0;
    mainIn->GetObject("h_Pion", hPion);
    TH3D* hEvsLambdavsMass = hPion->Projection(axis_photon1E, axis_pionLambda1, axis_pionMass);
    pionIn->GetObject("mass_pion", hFittedMass);
    TF1* data_fit = (TF1*) hFittedMass->GetListOfFunctions()->FindObject("mass peak");
    pionIn->GetObject("Total_Sig_To_Total", sigTotal);
    double mean = data_fit->GetParameter(1);
    double sigma = data_fit->GetParameter(2);
    
    // Define mass bound constants, cut the mass data, momentum data, matched tracks data
    const double massmin = mean - NumOfSigmasFromMean*sigma;
    const double massmax = mean + NumOfSigmasFromMean*sigma;
    SetCut(h_Pion, axis_pionMatchedTracks1, -1.5, -0.75);
    SetCut(h_Pion, axis_pionMatchedTracks2, -1.5, -0.75);
    SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
    SetCut(h_Pion, axis_pionAngle, 0.015, 0.5);
    SetCut(h_Pion, axis_pionNcells1, 1.0, 30.0);
    SetCut(h_Pion, axis_pionNcells2, 1.0, 30.0);
    SetCut(h_Pion, axis_pionDisToBorder1, 2.0, 6.0);
    SetCut(h_Pion, axis_pionDisToBorder1, 2.0, 6.0);
    SetCut(h_photon, axis_pionPt, 8.0, 15.0);
    SetCut(h_photon, axis_pionMass, massmin, massmax);
    
    // Create a template based on a Mass and Lambda vs Entries projection in the THnSparses
    TH2D* MLambdaProjection = hPion->Projection(axis_pionLambda1, axis_pionMass);
    // Use that template to integrate the hEvsLambdavsMass data
    double y_bin_min = MLambdaProjection->GetYaxis()->FindBin(lambdamin);
    double y_bin_max = MLambdaProjection->GetYaxis()->FindBin(lambdamax);
    double x_bin_min = MLambdaProjection->GetXaxis()->FindBin(massmin);
    double x_bin_max = MLambdaProjection->GetXaxis()->FindBin(massmax);
    double actual_bin;
    for (int i = y_bin_min; i <= y_bin_max; i++)
        for (int j = x_bin_min; j <= x_bin_max; j++) {
            actual_bin = hEvsLambdavsMass->Integral(x_bin_min, x_bin_max, i, i, j, j);
            MLambdaProjection->SetBinContent(j, i, actual_bin);
        }
    
    // Save the graph as a .png
    
    // Now, for each mass in MLambdaProjection, starting from the bin that is right in the 
}
