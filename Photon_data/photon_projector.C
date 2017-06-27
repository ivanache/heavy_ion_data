// This macro is designed to make lambda vs E projections from the leading photon output of photon_analyzer
// Programmer: Ivan Chernyshev; Date Started :June 26, 2017

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

// Concatenates two strings and gives a char array
char* str_concat_converter(string str1, string str2){
    string sumstring = str1 + str2;
    char* output = new char[sumstring.length() + 1];
    strcpy(output, sumstring.c_str());
    return output;
}

/**
 Main function
 */
// Precondition: NumOfSigmasFromMeanMax = 1, 2, 3, or 4
void photon_projector(string NumOfSigmasFromMeanMax) {
    const double lambdamin = 0;
    const double lambdamax = 2;
    const int axis_pionLambda1 = 17;
    const int y_min = 6;
    const int y_max = 8;
    const string directory_name = Form("%ssigma/projections/", NumOfSigmasFromMeanMax.c_str());
    TCanvas* canvas = new TCanvas();
    
    // Set ATLAS Style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Open the root file, get the data
    TFile* fIn = new TFile(Form("%ssigmaPhotonsOutput.root", NumOfSigmasFromMeanMax.c_str()), "READ");
    fIn->Print();
    TH2D* hist = 0;
    fIn->GetObject("leading_Evslambda", hist);
    int yminbin = hist->GetYaxis()->FindBin(y_min);
    int ymaxbin = hist->GetYaxis()->FindBin(y_max);
    
    // Create the TH1D object, using the lambda vs entries graph from the THnSparses as a template
    TFile* TemplateSource = new TFile("THnSparses_062017.root", "READ");
    THnSparse* temp = 0;
    TemplateSource->GetObject("h_Pion", temp);
    TH1D* projection = temp->Projection(axis_pionLambda1);
    // Copy the energy data from the E vs lambda histogram, integrate it over Y, put it into the corresponding bin of the projection
    int maxbin = projection->GetXaxis()->FindBin(lambdamax);
    double actual_bin;
    
    for (int i = projection->GetXaxis()->FindBin(lambdamin); i <= maxbin; i++) {
        actual_bin = hist->Integral(i, i, yminbin, ymaxbin);
        projection->SetBinContent(i, actual_bin);
    }
    // Before graphing, create an object to write to the root file being read from
    TFile* fOut = new TFile(Form("%ssigmaPhotonsOutput.root", NumOfSigmasFromMeanMax.c_str()), "UPDATE");

    // Graph and update the root file
    projection->SetTitle("Photons from Pi0 decay; lambda0; Number of Leading Photons with energy 6-8 GeV");
    projection->GetXaxis()->SetTitleOffset(1.0);
    projection->GetYaxis()->SetTitleOffset(1.0);
    projection->Draw();
    myText(.35,.95, kBlack, "Photons from Pi0 decay");
    projection->Write("Lambda_projection");
    canvas->SaveAs(str_concat_converter(directory_name, "Lambda_vs_E_projection.png"));
    
    
    // Repeat for all intervals of momentum
    const int num_of_intervals = 5;
    double intervals[num_of_intervals][2] = {{8.0, 10.0}, {10.0, 11.0}, {11.0, 12.0}, {12.0, 13.0}, {13.0, 15.0}};
    for(int i = 0; i < num_of_intervals; i++) {
        double ptmin = intervals[i][0];
        double ptmax = intervals[i][1];
        canvas->Clear();
        
        // Open the root file, get the data
        fIn->GetObject(Form("leading_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax), hist);
        yminbin = hist->GetYaxis()->FindBin(y_min);
        ymaxbin = hist->GetYaxis()->FindBin(y_max);
        
        // Copy the energy data from the E vs lambda histogram, integrate it over Y, put it into the corresponding bin of the projection
        maxbin = projection->GetXaxis()->FindBin(lambdamax);
        for (int i = projection->GetXaxis()->FindBin(lambdamin); i <= maxbin; i++) {
            actual_bin = hist->Integral(i, i, yminbin, ymaxbin);
            projection->SetBinContent(i, actual_bin);
        }
        // Graph
        projection->SetTitle(Form("Leading Pi0 Photon Lambda vs. Leading Pi0 Photon Energy Projection (Momentum %2.2f-%2.2f Gev); lambda0; Number of Leading Photons with energy 6-8 GeV", ptmin, ptmax));
        projection->Draw();
        myText(.10,.95, kBlack, "Leading Pi0 Photon Lambda vs. Leading Pi0 Photon Energy Projection");
        projection->Write(Form("Lambda_projection_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax));
        canvas->SaveAs(str_concat_converter(directory_name, Form("Lambda_vs_E_projection_ptmin_%2.2fGeV_ptmax_%2.2fGeV.png", ptmin, ptmax)));

    }
    
    canvas->Close();
    
}
