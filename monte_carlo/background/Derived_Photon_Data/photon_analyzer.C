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
const int axis_pion_Cen           = 0;
const int axis_pion_Zvtx          = 1;
const int axis_pionMass           = 2;
const int axis_pionPt             = 3;
const int axis_asymmetry          = 5;
const int axis_photon1Pt          = 6;
const int axis_photon2Pt          = 7;
const int axis_pionAngle          = 8;
const int axis_pionLambda1        = 9;
const int axis_pionLambda2        = 10;
const int axis_pionDisToCharged1  = 11;
const int axis_pionDisToCharged2  = 12;

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
    const int x_min = 0;
    const int x_max = 2;
    const int y_min = 2;
    const int y_max = 20;
    
    // Set ATLAS style
    //gROOT->LoadMacro("AtlasStyle.C");
    //SetAtlasStyle();
    
    // Load the THnSparses file, print its content, and get the data from it
    TFile* fIn = new TFile("THnSparses_080117_MC.root","READ");
    fIn->Print();
    THnSparse* h_photon = 0;
    fIn->GetObject("h_Pion", h_photon);
    TCanvas* canvas = new TCanvas();
    
    // Set maximum cut to the momentum data
    SetCut(h_photon, axis_pionPt, 6.0, 20.0);
    
    // Cut the matched tracks, the asymmetry, the angle, and the Ncells
    //SetCut(h_photon, axis_pionDisToCharged1, 0.02, 0.14);
    //SetCut(h_photon, axis_pionDisToCharged2, 0.02, 0.14);
    //SetCut(h_photon, axis_pionLambda1, 0.1, 0.4);
    //SetCut(h_photon, axis_pionLambda2, 0.1, 0.4);
    SetCut(h_photon, axis_pionAngle, 0.015, 0.6);
    
    // Load the THnSparses with the pion data
    // Use it to cut the mass to within NumOfSigmasFromMeanMax sigma of the mean pion mass
    TFile* pionIn = new TFile("BackgroundMC_PionSparsesOutput_angle_17mrad.root", "READ");
    TH1D* piondata = 0;
    pionIn->GetObject("mass_pion", piondata);
    TF1* peakfunct = (TF1*) piondata->GetListOfFunctions()->FindObject("mass peak");
    double mean = peakfunct->GetParameter(1);
    double sigma = peakfunct->GetParameter(2);
    SetCut(h_photon, axis_pionMass, mean - NumOfSigmasFromMeanMax*sigma, mean + NumOfSigmasFromMeanMax*sigma);

    // Set up the root output file
    string rootfilename = Form("%isigmaPhotonsOutput.root", NumOfSigmasFromMeanMax);
    TFile* fOut = new TFile(rootfilename.c_str(), "RECREATE");
    
    // Plot Pt vs. eigenvalue for the leading photon
    TH2D* hPt_leading = h_photon->Projection(axis_photon1Pt, axis_pionLambda1);
    hPt_leading->SetTitle("Leading Pi0 Photon Momentum vs lambda0; lambda0; Leading Photon Momentum (GeV/c)");
    hPt_leading->SetAxisRange(2, 20, "Y");
    hPt_leading->Draw("COLZ");
    myText(.20,.97, kBlack, "Leading  Pi0 Photon Momentum vs lambda0, Pion Momentum 6-20 GeVc");
    myText(.35,.92, kBlack, Form("Mass within %i sigma of the mean", NumOfSigmasFromMeanMax));
    hPt_leading->Write("leading_Evslambda");
    canvas->SaveAs(str_concat_converter(directory_name, "LeadingEvsLambda.png"));
    canvas->Clear();
    
    // Do the same for the trailing photon
    TH2D* hPt_trailing = h_photon->Projection(axis_photon2Pt, axis_pionLambda2);
    hPt_trailing->SetTitle("Trailing Pi0 Photon Momentum vs lambda0; lambda0; Trailing Photon Momentum (GeV/c)");
    hPt_trailing->SetAxisRange(2, 20, "Y");
    hPt_trailing->Draw("COLZ");
    myText(.20,.97, kBlack, "Trailing Pi0 Photon Momentum vs lambda0, Pion Momentum 6-20 GeV,");
    myText(.35,.92, kBlack, Form("mass within %i sigma of the mean", NumOfSigmasFromMeanMax));
    hPt_trailing->Write("trailing_Evslambda");
    canvas->SaveAs(str_concat_converter(directory_name, "TrailingEvsLambda.png"));
    canvas->Clear();
    
    // Combine the results
    TH2D* hPt_total = (TH2D*)hPt_trailing->Clone();
    hPt_total->SetTitle("Total Pi0 Photon Momentum vs lambda0; lambda0; Photon Momentum (GeV/c)");
    for (double x = x_min; x <= x_max; x += hPt_trailing->GetXaxis()->GetBinWidth(0))
        for (double y = y_min; y <= y_max; y += hPt_trailing->GetYaxis()->GetBinWidth(0)) {
            int x_leading_binnum = hPt_leading->GetXaxis()->FindBin(x);
            int x_trailing_binnum = hPt_trailing->GetXaxis()->FindBin(x);
            int y_leading_binnum = hPt_leading->GetYaxis()->FindBin(y);
            int y_trailing_binnum = hPt_trailing->GetYaxis()->FindBin(y);
            
            double sum = hPt_leading->GetBinContent(x_leading_binnum, y_leading_binnum) + hPt_trailing->GetBinContent(x_trailing_binnum, y_trailing_binnum);
            hPt_total->SetBinContent(x_leading_binnum, y_leading_binnum, sum);
        }
    hPt_total->Draw("COLZ");
    myText(.20,.97, kBlack, "Total  Pi0 Photon Momentum vs lambda0, Pion Momentum 6-20 GeVc");
    myText(.35,.92, kBlack, Form("Mass within %i sigma of the mean", NumOfSigmasFromMeanMax));
    hPt_total->Write("total_Evslambda");
    canvas->SaveAs(str_concat_converter(directory_name, "TotalEvsLambda.png"));
    canvas->Clear();
    
    // Do both of the above for various momentum intervals
    const int num_of_intervals = 7;
    double intervals[num_of_intervals][2] = {{6.0, 20.0}, {6.0, 8.0}, {8.0, 10.0}, {10.0, 12.0}, {12.0, 14.0}, {14.0, 16.0}, {16.0, 20.0}};
    for(int i = 0; i < num_of_intervals; i++) {
        
        double ptmin = intervals[i][0];
        double ptmax = intervals[i][1];
        
        // Cut the data
        SetCut(h_photon, axis_pionPt, ptmin, ptmax);
        
        // Plot Momentum vs. eigenvalue for the leading photon
        hPt_leading = h_photon->Projection(axis_photon1Pt, axis_pionLambda1);
        hPt_leading->SetTitle("Leading Pi0 Photon Momentum vs lambda0; lambda0; Leading photon momentum (GeV/c)");
        hPt_leading->SetAxisRange(6, 20, "Y");
        hPt_leading->Draw("COLZ");
        myText(.20,.97, kBlack, Form("Leading Pi0 Photon Momentum vs lambda0, Pt %2.2f-%2.2f GeV", ptmin, ptmax));
        myText(.35,.92, kBlack, Form("mass within %i sigma of the mean", NumOfSigmasFromMeanMax));
        hPt_leading->Write(Form("leading_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax));
        canvas->SaveAs(str_concat_converter(directory_name, Form("LeadingEvsLambda_ptmin_%2.2f_ptmax_%2.2f.png", ptmin, ptmax)));
        canvas->Clear();
        
        // Do the same for the trailing photon
        hPt_trailing = h_photon->Projection(axis_photon2Pt, axis_pionLambda2);
        hPt_trailing->SetTitle("Trailing Pi0 Photon Momentum vs lambda0; lambda0; Trailing photon momentum (GeV/c)");
        hPt_trailing->SetAxisRange(3, 15, "Y");
        hPt_trailing->Draw("COLZ");
        myText(.20,.97, kBlack, Form("Trailing Pi0 Photon Momentum vs lambda0, Pt %2.2f-%2.2f GeV", ptmin, ptmax));
        myText(.35,.92, kBlack, Form("mass within %i sigma of the mean", NumOfSigmasFromMeanMax));
        hPt_trailing->Write(Form("trailing_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax));
        canvas->SaveAs(str_concat_converter(directory_name, Form("TrailingEvsLambda_ptmin_%2.2f_ptmax_%2.2f.png", ptmin, ptmax)));
        canvas->Clear();
        
        // Combine the results
        for (double x = x_min; x <= x_max; x += hPt_trailing->GetXaxis()->GetBinWidth(1))
            for (double y = y_min; y <= y_max; y += hPt_trailing->GetYaxis()->GetBinWidth(1)) {
                int x_leading_binnum = hPt_leading->GetXaxis()->FindBin(x);
                int x_trailing_binnum = hPt_trailing->GetXaxis()->FindBin(x);
                int y_leading_binnum = hPt_leading->GetYaxis()->FindBin(y);
                int y_trailing_binnum = hPt_trailing->GetYaxis()->FindBin(y);
                
                double sum;
                sum = hPt_leading->GetBinContent(x_leading_binnum, y_leading_binnum) + hPt_trailing->GetBinContent(x_trailing_binnum, y_trailing_binnum);
                hPt_total->SetBinContent(x_leading_binnum, y_leading_binnum, sum);
            }
        hPt_total->Draw("COLZ");
        myText(.20,.97, kBlack, Form("Total Pi0 Photon Momentum vs lambda0, Pt %2.2f-%2.2f GeV", ptmin, ptmax));
        myText(.35,.92, kBlack, Form("mass within %i sigma of the mean", NumOfSigmasFromMeanMax));
        hPt_total->Write(Form("total_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax));
        canvas->SaveAs(str_concat_converter(directory_name, Form("TotalEvsLambda_ptmin_%2.2f_ptmax_%2.2f.png", ptmin, ptmax)));
        canvas->Clear();
    }
    canvas->Close();
    return;
}
