// photon_analyzer is a program for analyzing photon data from a THnSparses
// Programmer: Ivan Chernyshev; Date Started :June 23, 2017

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
 Main function
*/
// Precondition: NumOfSigmasFromMeanMax = 1, 2, 3, or 4
void photon_analyzer(int NumOfSigmasFromMeanMax) {
    const string directory_name = Form("%isigma/", NumOfSigmasFromMeanMax);
    
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Load the THnSparses file, print its content, and get the data from it
    TFile* fIn = new TFile("THnSparses_062717.root","READ");
    fIn->Print();
    THnSparse* h_photon = 0;
    fIn->GetObject("h_Pion", h_photon);
    TCanvas* canvas = new TCanvas();
    
    // Set maximum cut to the momentum data
    SetCut(h_photon, axis_pionPt, 8.0, 15.0);
    
    // Cut the matched tracks, the asymmetry, the angle, and the Ncells
    SetCut(h_photon, axis_pionMatchedTracks1, -1.5, -0.75);
    SetCut(h_photon, axis_pionMatchedTracks2, -1.5, -0.75);
    SetCut(h_photon, axis_asymmetry, 0.0, 0.7);
    SetCut(h_photon, axis_pionAngle, 0.015, 0.5);
    SetCut(h_photon, axis_pionNcells1, 1.0, 30.0);
    SetCut(h_photon, axis_pionNcells2, 1.0, 30.0);
    
    
    // Load the THnSparses with the pion data
    // Use it to cut the mass to within 1 sigma of the mean pion mass
    TFile* pionIn = new TFile("PionCutsSparsesOutput.root", "READ");
    TH1D* piondata = 0;
    pionIn->GetObject("mass_pion", piondata);
    TF1* peakfunct = (TF1*) piondata->GetListOfFunctions()->FindObject("mass peak");
    double mean = peakfunct->GetParameter(1);
    double sigma = peakfunct->GetParameter(2);
    SetCut(h_photon, axis_pionMass, mean - NumOfSigmasFromMeanMax*sigma, mean + NumOfSigmasFromMeanMax*sigma);

    // Set up the root output file
    string rootfilename = Form("%isigmaPhotonsOutput.root", NumOfSigmasFromMeanMax);
    TFile* fOut = new TFile(rootfilename.c_str(), "RECREATE");
    
    // Plot energy vs. eigenvalue for the leading photon
    TH2D* hEnergy = h_photon->Projection(axis_photon1E, axis_pionLambda1);
    hEnergy->SetTitle("Leading Pi0 Photon Energy vs lambda0; lambda0; Leading photon energy (GeV)");
    hEnergy->SetAxisRange(6, 15, "Y");
    hEnergy->Draw("COLZ");
    myText(.35,.9, kBlack, "Leading  Pi0 Photon Energy vs lambda0");
    hEnergy->Write("leading_Evslambda");
    canvas->SaveAs(str_concat_converter(directory_name, "LeadingEvsLambda.png"));
    canvas->Clear();
    
    // Do the same for the trailing photon
    hEnergy = h_photon->Projection(axis_photon2E, axis_pionLambda2);
    hEnergy->SetTitle("Trailing Pi0 Photon Energy vs lambda0; lambda0; Trailing photon energy (GeV)");
    hEnergy->SetAxisRange(3, 9, "Y");
    hEnergy->Draw("COLZ");
    myText(.35,.9, kBlack, "Trailing Pi0 Photon Energy vs lambda0");
    hEnergy->Write("trailing_Evslambda");
    canvas->SaveAs(str_concat_converter(directory_name, "TrailingEvsLambda.png"));
    canvas->Clear();
    
    // Do both of the above for various momentum intervals
    const int num_of_intervals = 5;
    double intervals[num_of_intervals][2] = {{8.0, 10.0}, {10.0, 11.0}, {11.0, 12.0}, {12.0, 13.0}, {13.0, 15.0}};
    for(int i = 0; i < num_of_intervals; i++) {
        double ptmin = intervals[i][0];
        double ptmax = intervals[i][1];
        
        // Cut the data
        SetCut(h_photon, axis_pionPt, ptmin, ptmax);
        
        // Use the THnSparses with the pion data to cut the mass to within 1 sigma of the mean pion mass
        pionIn->GetObject(Form("mass-pion-%2.2fGeV-%2.2fGeV", ptmin, ptmax), piondata);
        TF1* peakfunct = (TF1*) piondata->GetListOfFunctions()->FindObject("mass peak");
        double mean = peakfunct->GetParameter(1);
        double sigma = peakfunct->GetParameter(2);
        SetCut(h_photon, axis_pionMass, mean - sigma, mean + sigma);

        
        // Plot energy vs. eigenvalue for the leading photon
        TH2D* hEnergy = h_photon->Projection(axis_photon1E, axis_pionLambda1);
        hEnergy->SetTitle("Leading Pi0 Photon Energy vs lambda0; lambda0; Leading photon energy (GeV)");
        hEnergy->SetAxisRange(6, 15, "Y");
        hEnergy->Draw("COLZ");
        myText(.35,.9, kBlack, "Leading Pi0 Photon Energy vs lambda0");
        hEnergy->Write(Form("leading_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax));
        canvas->SaveAs(str_concat_converter(directory_name, Form("LeadingEvsLambda_ptmin_%2.2f_ptmax_%2.2f.png", ptmin, ptmax)));
        canvas->Clear();
        
        // Do the same for the trailing photon
        hEnergy = h_photon->Projection(axis_photon2E, axis_pionLambda2);
        hEnergy->SetTitle("Trailing Pi0 Photon Energy vs lambda0; lambda0; Trailing photon energy (GeV)");
        hEnergy->SetAxisRange(3, 9, "Y");
        hEnergy->Draw("COLZ");
        myText(.35,.9, kBlack, "Trailing Pi0 Photon Energy vs lambda0");
        hEnergy->Write(Form("trailing_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax));
        canvas->SaveAs(str_concat_converter(directory_name, Form("TrailingEvsLambda_ptmin_%2.2f_ptmax_%2.2f.png", ptmin, ptmax)));
        canvas->Clear();

    }
    canvas->Close();
    return;
}
