/**
   This macro approximates the background on a given interval based on the signal between 3 and 4 sigma, given by photon_subtractor, as well as a signal to total ratio and the signal plus background data given
   Programmer: Ivan Chernyshev; Date created: 6/28/2017
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

// Concatenates two strings and gives a char array
char* str_concat_converter(string str1, string str2){
    string sumstring = str1 + str2;
    char* output = new char[sumstring.length() + 1];
    strcpy(output, sumstring.c_str());
    return output;
}

/**
 Main Function
*/
void background_maker(int NumOfSigmasFromMean) {
    const double lambda_min = 0;
    const double lambda_max = 2;
    const double energy_min = 6;
    const double energy_max = 15;
    string background_directory_name = Form("%isigma/background/", NumOfSigmasFromMean);
    string total_directory_name = Form("%isigma/background+total/", NumOfSigmasFromMean);
    TCanvas* canvas = new TCanvas();
    
    // Set ATLAS Style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Load the file with the 3-5 background sample, signal-to-total ratio, and data
    TFile* backIn = new TFile("3-5sigmaPhotonsOutput.root", "READ");
    TFile* ratioIn = new TFile("Pion3CutsSparsesOutput.root", "READ");
    TFile* totalIn = new TFile(Form("%isigmaPhotonsOutput.root", NumOfSigmasFromMean), "READ");
    
    // Load data
    TH2D* total_lambda_E_graph = 0;
    TH1D* total_lambda_projection = 0;
    TGraph* sig_tot_graph = 0;
    TH2D* back_lambda_E_graph = 0;
    TH1D* back_lambda_projection = 0;
    totalIn->GetObject("leading_Evslambda", total_lambda_E_graph);
    totalIn->GetObject("Lambda_projection", total_lambda_projection);
    ratioIn->GetObject("Total_Sig_To_Total", sig_tot_graph);
    backIn->GetObject("leading_Evslambda", back_lambda_E_graph);
    backIn->GetObject("Lambda_projection", back_lambda_projection);
    
    // Take the integrals of the two backgrounds and the two totals
    int x_bin_min = total_lambda_E_graph->GetXaxis()->FindBin(lambda_min);
    int x_bin_max = total_lambda_E_graph->GetXaxis()->FindBin(lambda_max);
    int y_bin_min = total_lambda_E_graph->GetYaxis()->FindBin(energy_min);
    int y_bin_max = total_lambda_E_graph->GetYaxis()->FindBin(energy_max);
    double h2d_total_int = total_lambda_E_graph->Integral(x_bin_min, x_bin_max, y_bin_min, y_bin_max);
    double h1d_total_int = total_lambda_projection->Integral(x_bin_min, x_bin_max);
    double h2d_back_int = back_lambda_E_graph->Integral(x_bin_min, x_bin_max, y_bin_min, y_bin_max);
    double h1d_back_int = back_lambda_projection->Integral(x_bin_min, x_bin_max);
    // Print it out, so it can be checked if it makes sense
    std::cout << "\n3-5 sigma background integrals:\n2D histogram: " << h2d_back_int << "\n1D histogram: " << h1d_back_int << Form("\n%i sigma total integrals:\n2Dhistogram: ", NumOfSigmasFromMean) << h2d_total_int << "\n1D histogram: " << h1d_total_int << std::endl;
    
    // Now take the background-to-total ratio using 1 - signal to total
    // Find the ratio by which to multiply all of the points using (b-t ratio)(total integral)/(experimental background integral)
    double bt_ratio = 1.0 - (sig_tot_graph->Eval(NumOfSigmasFromMean));
    double h2d_adjustment_factor = (bt_ratio * h2d_total_int)/(h2d_back_int);
    double h1d_adjustment_factor = (bt_ratio * h1d_total_int)/(h1d_back_int);
    // Print these out, so that they can be checked if they make sense
    std::cout << "\nCumulative background to total ratio: " << bt_ratio << "\nAdjustment factor, calculated from 2D histograms: " << h2d_adjustment_factor << "\nAdjustment factor, calculated from 1D histograms: " << h1d_adjustment_factor << std::endl;
    
    // For all points in both background graphs, multiply each point by the adjustment factor
    double value_3_4_sigma;
    for (int i = x_bin_min; i <= x_bin_max; i++) {
        for (int j = y_bin_min; j <= y_bin_max; j++) {
            value_3_4_sigma = back_lambda_E_graph->GetBinContent(i, j);
            back_lambda_E_graph->SetBinContent(i, j, value_3_4_sigma * h2d_adjustment_factor);
        }
        value_3_4_sigma = back_lambda_projection->GetBinContent(i);
        back_lambda_projection->SetBinContent(i, value_3_4_sigma * h1d_adjustment_factor);
    }
    
    // Graph out the two results and save them
    back_lambda_E_graph->GetXaxis()->SetTitleOffset(0.9);
    back_lambda_E_graph->GetYaxis()->SetTitleOffset(0.9);
    back_lambda_E_graph->Draw("COLZ");
    myText(0.10, 0.95, kBlack, Form("Leading Pi0 Photon Energy vs lambda0, Pt 8-15 GeV, mass within %i sigma of mean", NumOfSigmasFromMean));
    canvas->SaveAs(str_concat_converter(background_directory_name, "LeadingEvsLambda.png"));
    canvas->Clear();
    back_lambda_projection->GetXaxis()->SetTitleOffset(1.0);
    back_lambda_projection->GetYaxis()->SetTitleOffset(1.0);
    back_lambda_projection->Draw();
    myText(0.10, 0.95, kBlack, Form("Photons from Pi0 decay, Pt 8-15 GeV, mass within %i sigma of mean", NumOfSigmasFromMean));
    canvas->SaveAs(str_concat_converter(background_directory_name, "Lambda_vs_E_projection.png"));
    total_lambda_projection->SetMarkerColor(kRed);
    total_lambda_projection->Draw();
    total_lambda_projection->GetXaxis()->SetTitleOffset(1.0);
    total_lambda_projection->GetYaxis()->SetTitleOffset(1.0);
    myText(0.10, 0.95, kBlack, Form("Photons from Pi0 decay, Pt 8-15 GeV, mass within %i sigma of mean", NumOfSigmasFromMean));
    back_lambda_projection->Draw("same");
    myMarkerText(0.30, 0.85, kBlack, 20, Form("Background data, from within %i sigma of the mean", NumOfSigmasFromMean), 1);
    myMarkerText(0.30, 0.80, kRed, 20, Form("Total data, from within %i sigma of the mean", NumOfSigmasFromMean), 1);
    canvas->SaveAs(str_concat_converter(total_directory_name, "Lambda_vs_E_projection.png"));
    
    // Get the signal-to-total ratios for the individual momentum intervals from the signal to toatal ratio source file (Pion3CutsSparsesOutput.root)
    TMultiGraph* interval_ratios = 0;
    ratioIn->GetObject("signal-over-total", interval_ratios);
    TList* listofgraphs = interval_ratios->GetListOfGraphs();
    
    // Repeat for all intervals of momentum
    const int num_of_intervals = 6;
    double intervals[num_of_intervals][2] = {{6.0, 8.0}, {8.0, 10.0}, {10.0, 12.0}, {12.0, 14.0}, {14.0, 16.0}, {16.0, 20.0}};
    //TH1D* interval_back_lambda_projection
    for(int i = 0; i < num_of_intervals; i++) {
        double ptmin = intervals[i][0];
        double ptmax = intervals[i][1];
        canvas->Clear();
        
        // Get Data
        totalIn->GetObject(Form("leading_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax), total_lambda_E_graph);
        totalIn->GetObject(Form("Lambda_projection_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax), total_lambda_projection);
        sig_tot_graph = (TGraph*)listofgraphs->FindObject(Form("sig_to_tot_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax));
        backIn->GetObject(Form("leading_Evslambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax), back_lambda_E_graph);
        backIn->GetObject(Form("Lambda_projection_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmax), back_lambda_projection);
        //backIn->GetObject(Form("Lambda_projection_ptmin_%2.2fGeV_ptmax_%2.2fGeV", ptmin, ptmin), back_lambda_projection);
        
        // Take the integrals of the two backgrounds and the two totals
        x_bin_min = total_lambda_E_graph->GetXaxis()->FindBin(lambda_min);
        x_bin_max = total_lambda_E_graph->GetXaxis()->FindBin(lambda_max);
        y_bin_min = total_lambda_E_graph->GetYaxis()->FindBin(energy_min);
        y_bin_max = total_lambda_E_graph->GetYaxis()->FindBin(energy_max);
        h2d_total_int = total_lambda_E_graph->Integral(x_bin_min, x_bin_max, y_bin_min, y_bin_max);
        h1d_total_int = total_lambda_projection->Integral(x_bin_min, x_bin_max);
        h2d_back_int = back_lambda_E_graph->Integral(x_bin_min, x_bin_max, y_bin_min, y_bin_max);
        h1d_back_int = back_lambda_projection->Integral(x_bin_min, x_bin_max);
        // Print it out, so it can be checked if it makes sense
        std::cout << "\n3-5 sigma background integrals:\n2D histogram: " << h2d_back_int << "\n1D histogram: " << h1d_back_int << Form("\n%i sigma total integrals:\n2Dhistogram: ", NumOfSigmasFromMean) << h2d_total_int << "\n1D histogram: " << h1d_total_int << std::endl;
        
        // Now take the background-to-total ratio using 1 - signal to total
        // Find the ratio by which to multiply all of the points using (b-t ratio)(total integral)/(experimental background integral)
        bt_ratio = 1.0 - (sig_tot_graph->Eval(NumOfSigmasFromMean));
        h2d_adjustment_factor = (bt_ratio * h2d_total_int)/(h2d_back_int);
        h1d_adjustment_factor = (bt_ratio * h1d_total_int)/(h1d_back_int);
        // Print these out, so that they can be checked if they make sense
        std::cout << "\nCumulative background to total ratio: " << bt_ratio << "\nAdjustment factor, calculated from 2D histograms: " << h2d_adjustment_factor << "\nAdjustment factor, calculated from 1D histograms: " << h1d_adjustment_factor << std::endl;
        
        // For all points in both background graphs, multiply each point by the adjustment factor
        for (int i = x_bin_min; i <= x_bin_max; i++) {
            for (int j = y_bin_min; j <= y_bin_max; j++) {
                value_3_4_sigma = back_lambda_E_graph->GetBinContent(i, j);
                back_lambda_E_graph->SetBinContent(i, j, value_3_4_sigma * h2d_adjustment_factor);
            }
            value_3_4_sigma = back_lambda_projection->GetBinContent(i);
            back_lambda_projection->SetBinContent(i, value_3_4_sigma * h1d_adjustment_factor);
        }
        
        // Graph out the two results and save them
        back_lambda_E_graph->GetXaxis()->SetTitleOffset(0.9);
        back_lambda_E_graph->GetYaxis()->SetTitleOffset(0.9);
        back_lambda_E_graph->Draw("COLZ");
        myText(0.10, 0.95, kBlack, Form("Leading Pi0 Photon Energy vs lambda0, Pt 8-15 GeV, mass within %i sigma of mean", NumOfSigmasFromMean));
        canvas->SaveAs(str_concat_converter(background_directory_name, Form("LeadingEvsLambda_ptmin_%2.2fGeV_ptmax_%2.2fGeV.png", ptmin, ptmax)));
        canvas->Clear();
        back_lambda_projection->GetXaxis()->SetTitleOffset(1.0);
        back_lambda_projection->GetYaxis()->SetTitleOffset(1.0);
        back_lambda_projection->Draw();
        myText(0.10, 0.95, kBlack, Form("Photons from Pi0 decay, Pt 8-15 GeV, mass within %i sigma of mean", NumOfSigmasFromMean));
        canvas->SaveAs(str_concat_converter(background_directory_name, Form("Lambda_vs_E_projection_ptmin_%2.2fGeV_ptmax_%2.2fGeV.png", ptmin, ptmax)));
        total_lambda_projection->SetMarkerColor(kRed);
        total_lambda_projection->Draw();
        total_lambda_projection->GetXaxis()->SetTitleOffset(1.0);
        total_lambda_projection->GetYaxis()->SetTitleOffset(1.0);
        myText(0.10, 0.95, kBlack, Form("Photons from Pi0 decay, Pt 8-15 GeV, mass within %i sigma of mean", NumOfSigmasFromMean));
        back_lambda_projection->Draw("same");
        myMarkerText(0.30, 0.85, kBlack, 20, Form("Background data, from within %i sigma of the mean", NumOfSigmasFromMean), 1);
        myMarkerText(0.30, 0.80, kRed, 20, Form("Total data, from within %i sigma of the mean", NumOfSigmasFromMean), 1);
        canvas->SaveAs(str_concat_converter(total_directory_name, Form("Lambda_vs_E_projection.png_ptmin_%2.2fGeV_ptmax_%2.2fGeV.png", ptmin, ptmax)));
    }
    canvas->Close();
}
