// This macro takes two distances from the mean and uses the graphs made by photon_analyzer and photon_projector to make versions of these graphs that only includes data between the two distances from the mean
// Programmer: Ivan Chernyshev; Date created: June 26, 2017

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
// Precondition: NumOfSigmasFromMeanMin = 1, 2, 3, or 4, same with NumOfSigmasFromMeanMax,
// NumOfSigmasFromMeanMin < NumOfSigmasFromMeanMax
void photon_subtractor(int NumOfSigmasFromMeanMin, int NumOfSigmasFromMeanMax) {
    TCanvas* canvas = new TCanvas();
    const double lambdamin = 0;
    const double lambdamax = 2;
    const int y_min = 6;
    const int y_max = 15;
    string directory_name = Form("%i-%isigma/", NumOfSigmasFromMeanMin, NumOfSigmasFromMeanMax);
    
    //Grab the input root files, get the data
    TFile* fInMin = new TFile(Form("%isigmaPhotonsOutput.root", NumOfSigmasFromMeanMin), "READ");
    TFile* fInMax = new TFile(Form("%isigmaPhotonsOutput.root", NumOfSigmasFromMeanMax), "READ");
    fInMin->Print();
    fInMax->Print();
    TH2D* histmin = 0;
    TH2D* histmax = 0;
    fInMin->GetObject("leading_Evslambda", histmin);
    fInMax->GetObject("leading_Evslambda", histmax);
    
    // Create a TH2D to store the data, use one of the input histograms as a template
    TH2D* histoutput = (TH2D*)histmin->Clone("leading_Evslambda_subtracted");
    // Loop through all bins in the two th2d histograms, take the difference, store in the output histogram
    double x_bin_min = histoutput->GetXaxis()->FindBin(lambdamin);
    double x_bin_max = histoutput->GetXaxis()->FindBin(lambdamax);
    double y_bin_min = histoutput->GetYaxis()->FindBin(y_min);
    double y_bin_max = histoutput->GetYaxis()->FindBin(y_max);
    double new_content;
    for (int i = x_bin_min; i <= x_bin_max; i++)
        for (int j = y_bin_min; j <= y_bin_max; j++) {
            new_content = histmax->GetBinContent(i, j) - histmin->GetBinContent(i, j);
            histoutput->SetBinContent(i, j, new_content);
        }
    
    // Create an output root file to store the output
    TFile* fOut = new TFile(Form("%i-%isigmaPhotonsOutput.root", NumOfSigmasFromMeanMin, NumOfSigmasFromMeanMax),"RECREATE");
    
    // Graph, save, and write
    histoutput->SetTitle("; lambda0; Leading photon energy (GeV)");
    histoutput->Draw("COLZ");
    myText(.20,.97, kBlack, "Leading Pi0 Photon Energy vs lambda0, Pt 8-15 GeV");
    myText(.35,.92, kBlack, Form("mass %i-%i sigma from mean", NumOfSigmasFromMeanMin, NumOfSigmasFromMeanMax));
    histoutput->Write("leading_Evslambda");
    canvas->SaveAs(str_concat_converter(directory_name, "LeadingEvsLambda.png"));
    
    // Repeat for all intervals of momentum
    const int num_of_intervals = 5;
    double intervals[num_of_intervals][2] = {{8.0, 10.0}, {10.0, 11.0}, {11.0, 12.0}, {12.0, 13.0}, {13.0, 15.0}};
    double ptmin;
    double ptmax;
    for(int i = 0; i < num_of_intervals; i++) {
        ptmin = intervals[i][0];
        ptmax = intervals[i][1];
        
        // get the data
        fInMin->GetObject(Form("leading_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax), histmin);
        fInMax->GetObject(Form("leading_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax), histmax);
        
        // Loop through all bins in the two th2d histograms, take the difference, store in the output histogram
        for (int i = x_bin_min; i <= x_bin_max; i++)
            for (int j = y_bin_min; j <= y_bin_max; j++) {
                new_content = histmax->GetBinContent(i, j) - histmin->GetBinContent(i, j);
                histoutput->SetBinContent(i, j, new_content);
            }
        
        // Graph, save, and write
        histoutput->SetTitle("; lambda0; Leading photon energy (GeV)");
        histoutput->Draw("COLZ");
        myText(.20,.97, kBlack, Form("Leading Pi0 Photon Energy vs lambda0, Pt %2.2f-%2.2f GeV", ptmin, ptmax));
        myText(.35,.92, kBlack, Form("mass %i-%i sigma from mean", NumOfSigmasFromMeanMin, NumOfSigmasFromMeanMax));
        myText(.35,.9, kBlack, Form("Leading Pi0 Photon Energy vs lambda0, mass %i-%i sigma from mean", NumOfSigmasFromMeanMin, NumOfSigmasFromMeanMax));
        histoutput->Write(Form("leading_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax));
        canvas->SaveAs(str_concat_converter(directory_name, Form("leading_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV.png", ptmin, ptmax)));

    }
    
    // Now, repeat with the projections
    directory_name = Form("%i-%isigma/projections/", NumOfSigmasFromMeanMin, NumOfSigmasFromMeanMax);
    // Get the data
    TH1D* projmin = 0;
    TH1D* projmax = 0;
    fInMin->GetObject("Lambda_projection", projmin);
    fInMax->GetObject("Lambda_projection", projmax);
    
    // Create a TH1D to store the data, use one of the input projections as a template
    TH1D* projoutput = (TH1D*)projmin->Clone("Lambda_projection_subtracted");
    // Loop through all bins in the two th1d projections, take the difference, store in the output projection
    x_bin_min = projoutput->GetXaxis()->FindBin(lambdamin);
    x_bin_max = projoutput->GetXaxis()->FindBin(lambdamax);
    for (int i = x_bin_min; i <= x_bin_max; i++) {
        new_content = projmax->GetBinContent(i) - projmin->GetBinContent(i);
        projoutput->SetBinContent(i, new_content);
        projoutput->SetBinError(i, TMath::Sqrt(projoutput->GetBinContent(i)));
    }

    // Graph, save, and write
    projoutput->SetTitle("; lambda0; Leading photon energy (GeV)");
    projoutput->Draw();
    myText(.10, .97, kBlack, "Leading Pi0 Photon Energy vs lambda0, Pt 8-15 GeV");
    myText(.35,.92, kBlack, Form("mass %i-%i sigma from mean", NumOfSigmasFromMeanMin, NumOfSigmasFromMeanMax));
    projoutput->Write("Lambda_projection");
    canvas->SaveAs(str_concat_converter(directory_name, "Lambda_vs_E_projection.png"));
    
    // Do the same for each momentum interval
    for(int i = 0; i < num_of_intervals; i++)
    {
        ptmin = intervals[i][0];
        ptmax = intervals[i][1];
        
        fInMin->GetObject(Form("Lambda_projection_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax), projmin);
        fInMax->GetObject(Form("Lambda_projection_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax), projmax);
        
        // Loop through all bins in the two th1d projections, take the difference, store in the output projection
        x_bin_min = projoutput->GetXaxis()->FindBin(lambdamin);
        x_bin_max = projoutput->GetXaxis()->FindBin(lambdamax);
        double new_content;
        for (int i = x_bin_min; i <= x_bin_max; i++) {
            new_content = projmax->GetBinContent(i) - projmin->GetBinContent(i);
            projoutput->SetBinContent(i, new_content);
            projoutput->SetBinError(i, TMath::Sqrt(projoutput->GetBinContent(i)));
        }
        
        // Graph, save, and write
        projoutput->SetTitle("; lambda0; Leading photon energy (GeV)");
        projoutput->Draw();
        myText(.10,.97, kBlack, Form("Leading Pi0 Photon Energy vs lambda0, Pt %2.2f-%2.2f GeV", ptmin, ptmax));
        myText(.35,.92, kBlack, Form("mass %i-%i sigma from mean", NumOfSigmasFromMeanMin, NumOfSigmasFromMeanMax));
        projoutput->Write(Form("Lambda_projection_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax));
        canvas->SaveAs(str_concat_converter(directory_name, Form("Lambda_vs_E_projection.png_%2.2fGeV_ptmax_%2.2fGeV.png", ptmin, ptmax)));

    }
}
