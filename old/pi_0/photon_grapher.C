/**
   This macro creates momentum-photon graphs from THnSparses_070517
*/
// Author: Ivan Chernyshev

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

const int axis_photonPt           =  3;
const int axis_photonLambda       =  6;
const int axis_photonNcells       =  7;
const int axis_photonDisToBorder  = 9;
const int axis_photonDisToBadCell = 10;
const int axis_photonDisToCharged = 11;

// The cutting function
void SetCut(THnSparse* h, const int axis, double min, double max){
    //make a selection on the chosen variable
    double width = h->GetAxis(axis)->GetBinWidth(1);
    int binmin = h->GetAxis(axis)->FindBin(min);
    int binmax = h->GetAxis(axis)->FindBin(max);
    h->GetAxis(axis)->SetRange(binmin, binmax - 1);
    return;
}

// Concatenates two strings and gives a char array
char* str_concat_converter(string str1, string str2){
    string sumstring = str1 + str2;
    char* output = new char[sumstring.length() + 1];
    strcpy(output, sumstring.c_str());
    return output;
}

/**
 Main function
 lambda_option is used to denote which portion of lambda0 is to be used:
    "0.1-0.4" for all lambda between 0.1 and 0.4
    ">0.4" for all lambda greater than 0.4
 pT_option is used to denote which pT interval to use
    "4-20" for all pT between 4 and 20 GeV
    "10-20" for all pT between 10 and 20 GeV
    All option values other than the above will result in the program shutting down
*/
void photon_grapher(string lambda_option, string pT_option) {
    // Initialize the lambda bounds, pT lower bound, and the directory to be used based on the options, quit the program if an option is invalid
    string directory_name;
    double lambdamin;
    double lambdamax;
    if (lambda_option == "0.1-0.4"){
        directory_name = "0.1<lambda<0.4/";
        lambdamin = 0.1;
        lambdamax = 0.4;
    }
    else if (lambda_option == ">0.4"){
        directory_name = "0.4<lambda/";
        lambdamin = 0.4;
        lambdamax = 2.1;
    }
    else {
        std::cout << "ERROR: input lambda_option must be \"0.1-0.4\" or \">0.4\"" << std::endl;
        return;
    }
    
    double pTmin;
    double pTmax;
    if (pT_option == "4-20") {
        directory_name += "4<pT<20/";
        pTmin = 4;
        pTmax = 20;
    }
    else if (pT_option == "10-20") {
        directory_name += "10<pT<20/";
        pTmin = 10;
        pTmax = 20;
    }
    else {
        std::cout << "ERROR: input pT_option must be \"4-20\" or \"10-20\"" << std::endl;
        return;
    }
    
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Define the constant ratio = 1 function, to be used in every ratio plot.
    TF1* ratio1 = new TF1("ratio1", "[0]", 4, 20);
    ratio1->SetParameter(0, 1);
    ratio1->SetLineColor(kBlack);
    ratio1->SetLineStyle(kDashed);
    
    // Get the momentum-photon data file
    TFile* fIn = new TFile("THnSparses_LHC13d.root", "READ");
    THnSparse* hPhoton = 0;
    fIn->GetObject("h_Cluster", hPhoton);
    //TCanvas* canvas = new TCanvas();
    TCanvas* logcanvas = new TCanvas();
    //logcanvas->SetLogy();
    
    // Do the baseline cuts
    SetCut(hPhoton, axis_photonPt, pTmin, pTmax);
    SetCut(hPhoton, axis_photonDisToCharged, 0.02, 0.14);
    SetCut(hPhoton, axis_photonNcells, 1.9, 30.0);
    SetCut(hPhoton, axis_photonLambda, lambdamin, lambdamax);
    
    // Graph the whole sample
    logcanvas->cd();
    TH1D* hMomentum = hPhoton->Projection(axis_photonPt);
    hMomentum->SetAxisRange(0.0, 1400.0, "Y");
    hMomentum->SetTitle("Pt-Photon Plot (all data); Pt (GeV/c); Number of Entries");
    hMomentum->Draw();
    myText(.20,.95, kBlack, "Pt-Photon Plot (all data)");

    logcanvas->SaveAs(str_concat_converter(directory_name, "pt_photon.png"));
    logcanvas->Clear();
    
    // For reference, also do a distance to bad cells-photon graph and a distance to border graph
    TH1D* hParam = hPhoton->Projection(axis_photonDisToBadCell);
    hParam->Draw();
    logcanvas->SaveAs(str_concat_converter(directory_name, "DisToBadCell_photon.png"));
    logcanvas->Clear();
    
    hParam = hPhoton->Projection(axis_photonDisToBorder);
    hParam->Draw();
    logcanvas->SaveAs(str_concat_converter(directory_name, "DisToBorder_photon.png"));
    logcanvas->Clear();
    
    //Cut Distance to bad cells to greater than 3, then graph the momentum photon data again
    SetCut(hPhoton, axis_photonDisToBadCell, 4, 10);
    TH1D* hMomentum_DisToBadCells_upper = hPhoton->Projection(axis_photonPt);
    hMomentum_DisToBadCells_upper->Rebin(2);
    hMomentum_DisToBadCells_upper->SetTitle("Pt-Photon Plot; Pt (GeV/c); Fraction of total");
    // Normalize the plot
    hMomentum_DisToBadCells_upper->Scale(1/hMomentum_DisToBadCells_upper->Integral());
    if (pT_option == "4-20")
        hMomentum_DisToBadCells_upper->GetYaxis()->SetRangeUser(0.0, 0.3);
    if (pT_option == "10-20")
        hMomentum_DisToBadCells_upper->GetYaxis()->SetRangeUser(0.0, 0.5);
    hMomentum_DisToBadCells_upper->SetMarkerStyle(21);
    hMomentum_DisToBadCells_upper->Draw();
    myText(.30,.95, kBlack, "Pt-Photon Plot");
    logcanvas->SaveAs(str_concat_converter(directory_name, "pt_photon_distobadcell>3.png"));
    
    // Some constants to use for the rest of the analysis
    const int num_of_distances = 3;
    int distances_to_bad_cells[num_of_distances] = {1, 2, 3};
    Color_t colors[num_of_distances] = {kRed, kBlue, kGreen};
    int marker_styles[num_of_distances] = {20, 22, 34};
    double bin_offsets[num_of_distances] = {-0.05, 0.05, 0.10};
    
    //Cut Distance to bad cells to 1, 2, and 3, and draw the resulting histogram for each amount of cuts
    logcanvas->cd();
    TH1D* hMomentum_DisToBadCells_lower[num_of_distances];
    for (int i = 0; i < num_of_distances; i++) {
        SetCut(hPhoton, axis_photonDisToBadCell, distances_to_bad_cells[i], distances_to_bad_cells[i] + 1);
        hMomentum_DisToBadCells_lower[i] = hPhoton->Projection(axis_photonPt);
        hMomentum_DisToBadCells_lower[i]->Rebin(2);
        hMomentum_DisToBadCells_lower[i]->SetMarkerColor(colors[i]);
        hMomentum_DisToBadCells_lower[i]->SetMarkerStyle(marker_styles[i]);
        hMomentum_DisToBadCells_lower[i]->SetLineColor(colors[i]);
        hMomentum_DisToBadCells_lower[i]->SetTitle("Pt-Photon Plot; Pt (GeV/c); Fraction of total");
        // Normalize the plot for the lower distance to bad cells
        hMomentum_DisToBadCells_lower[i]->Scale(1/hMomentum_DisToBadCells_lower[i]->Integral());
        hMomentum_DisToBadCells_lower[i]->GetXaxis()->SetLimits(pTmin + bin_offsets[i], pTmax + bin_offsets[i]);
        hMomentum_DisToBadCells_lower[i]->Draw("same");
    }
    myMarkerText(0.42, 0.87, colors[0], marker_styles[0], "DisToBadCell=1", 1);
    myMarkerText(0.42, 0.82, colors[1], marker_styles[1], "DisToBadCell=2", 1);
    myMarkerText(0.42, 0.77, colors[2], marker_styles[2], "DisToBadCell=3", 1);
    myMarkerText(0.42, 0.72, kBlack, 21, "DisToBadCell>3", 1);
    logcanvas->SaveAs(str_concat_converter(directory_name, "pt_photon_distobadcell.png"));
    logcanvas->Clear();
    
    // Get a histogram of the ratio of the normalized photon count with a distance to bad cells over three
    // to the normalized photon count with a distance to bad cells of of 1, 2, and 3
    logcanvas->cd();
    double chisquares[num_of_distances] = {0, 0, 0};
    int num_of_bins[num_of_distances] = {hMomentum_DisToBadCells_upper->GetSize(), hMomentum_DisToBadCells_upper->GetSize(), hMomentum_DisToBadCells_upper->GetSize()};
    for (int i = 0; i < num_of_distances; i++) {
        TH1D* hRatio = (TH1D*) hMomentum_DisToBadCells_upper->Clone("ratio");
        for (int j = 0; j < hMomentum_DisToBadCells_upper->GetSize(); j++) {      // Set the hRatio histogram's values to the photon count ratios
            
            if ((hMomentum_DisToBadCells_upper->GetBinContent(j)) == 0) // AVOID DIVIDING BY ZERO
                continue;
            
            double new_bin_content = (hMomentum_DisToBadCells_lower[i]->GetBinContent(j))/(hMomentum_DisToBadCells_upper->GetBinContent(j));
            //std::cout << "new ratio is " << new_bin_content << std::endl;
            hRatio->SetBinContent(j, new_bin_content);
            hRatio->SetBinError(j, new_bin_content*TMath::Sqrt( ( ((hMomentum_DisToBadCells_lower[i]->GetBinError(j))/hMomentum_DisToBadCells_lower[i]->GetBinContent(j)) * ((hMomentum_DisToBadCells_lower[i]->GetBinError(j))/hMomentum_DisToBadCells_lower[i]->GetBinContent(j)) ) + ( ((hMomentum_DisToBadCells_upper->GetBinError(j))/hMomentum_DisToBadCells_upper->GetBinContent(j)) * ((hMomentum_DisToBadCells_upper->GetBinError(j))/hMomentum_DisToBadCells_upper->GetBinContent(j)) ) ));
            
            // Now, add this bin's contribution to the chi-square, provided that the bin isn't a division by zero, in which case decrement the bin count by 1 to show that the bin is not countede
            //std::cout << "This bin chi2: " << (new_bin_content - ratio1->GetParameter(0)) * (new_bin_content - ratio1->GetParameter(0))/(ratio1->GetParameter(0)) << std::endl;
            if (hMomentum_DisToBadCells_upper->GetBinContent(j) == 0)
                num_of_bins[i]--;
            else
                chisquares[i] += (new_bin_content - ratio1->GetParameter(0)) * (new_bin_content - ratio1->GetParameter(0))/((hRatio->GetBinError(j)) * (hRatio->GetBinError(j)));
        }
        hRatio->SetMarkerColor(colors[i]);
        hRatio->SetMarkerStyle(marker_styles[i]);
        hRatio->SetLineColor(colors[i]);
        // Draw the histogram
        hRatio->SetTitle(Form("Normalized Ratios of DisToBadCell>3 to DisToBadCell=%i photons; Pt (GeV/c); Ratio", distances_to_bad_cells[i]));
        hRatio->SetAxisRange(0.5, 1.4, "Y");
        hRatio->GetXaxis()->SetLimits(pTmin + bin_offsets[i], pTmax + bin_offsets[i]);
        if(i == 0)
            hRatio->Draw();
        else
            hRatio->Draw("same");
    }
    // Insert the ratio = 1 line, and save the combined graph as a file with a legend (include reduced chi-square in the legend)
    ratio1->Draw("same");
    myText(.05,.95, kBlack, "Normalized dBadCell>3/dBadCell=1,2,3 Photon Ratios");
    myMarkerText(0.19, 0.43, colors[0], marker_styles[0], "dBadCell=1/dBadCell>3 ratio", 1);
    if (lambda_option == "0.1-0.4")
        myText(0.19, 0.38, kBlack, Form("chi2_v = %3.1f; p-val = %2.2f", chisquares[0]/(num_of_bins[0] - 1), TMath::Prob(chisquares[0], (num_of_bins[0] - 1))));
    else
        myText(0.19, 0.38, kBlack, Form("chi2_v = %3.1f; p-val = %2.6f", chisquares[0]/(num_of_bins[0] - 1), TMath::Prob(chisquares[0], (num_of_bins[0] - 1))));
    myMarkerText(0.19, 0.34, colors[1], marker_styles[1], "dBadCell=2/dBadCell>3 ratio;", 1);
    myText(0.19, 0.29, kBlack, Form("chi2_v = %3.1f; p-val = %2.2f", chisquares[1]/(num_of_bins[1] - 1), TMath::Prob(chisquares[1], (num_of_bins[1] - 1))));
    myMarkerText(0.19, 0.25, colors[2], marker_styles[2], "dBadCell=3/dBadCell>3 ratio;", 1);
    myText(0.19, 0.20, kBlack, Form("chi2_v = %3.1f; p-val = %2.2f", chisquares[2]/(num_of_bins[2] - 1), TMath::Prob(chisquares[2], (num_of_bins[2] - 1))));
    myText(0.15, 0.16, kBlack, "#scale[1.5]{...} Ratio = 1 fit");
    logcanvas->SaveAs(str_concat_converter(directory_name, "photon_ratios_distobadcells.png"));
    logcanvas->Clear();
    
    // Cut distance to border to greater than 2, and graph the result
    logcanvas->Clear();
    logcanvas->cd();
    SetCut(hPhoton, axis_photonDisToBadCell, 0, 10); // Reset distance to bad cell cut to original value
    SetCut(hPhoton, axis_photonDisToBorder, 3, 5);
    TH1D* hMomentum_DisToBorder_upper = hPhoton->Projection(axis_photonPt);
    hMomentum_DisToBorder_upper->Rebin(2);
    hMomentum_DisToBorder_upper->SetTitle("Pt-Photon Plot; Pt (GeV/c); Fraction of total");
    // Normalize the plot
    hMomentum_DisToBorder_upper->Scale(1/hMomentum_DisToBorder_upper->Integral());
    if (pT_option == "4-20")
        hMomentum_DisToBorder_upper->GetYaxis()->SetRangeUser(0.0, 0.3);
    if (pT_option == "10-20")
        hMomentum_DisToBorder_upper->GetYaxis()->SetRangeUser(0.0, 0.5);
    hMomentum_DisToBorder_upper->SetMarkerStyle(21);
    hMomentum_DisToBorder_upper->Draw();
    myText(.30,.95, kBlack, "Pt-Photon Plot");
    logcanvas->SaveAs(str_concat_converter(directory_name, "pt_photon_distoborder>2.png"));

    //Cut Distance to border to 0, 1, and 2, and draw the resulting histogram for each amount of cuts
    TH1D* hMomentum_DisToBorder_lower[num_of_distances];
    int distances_to_border[num_of_distances] = {0, 1, 2};
    for (int i = 0; i < num_of_distances; i++) {
        SetCut(hPhoton, axis_photonDisToBorder, distances_to_border[i], distances_to_border[i] + 1);
        hMomentum_DisToBorder_lower[i] = hPhoton->Projection(axis_photonPt);
        hMomentum_DisToBorder_lower[i]->Rebin(2);
        hMomentum_DisToBorder_lower[i]->SetMarkerColor(colors[i]);
        hMomentum_DisToBorder_lower[i]->SetMarkerStyle(marker_styles[i]);
        hMomentum_DisToBorder_lower[i]->SetLineColor(colors[i]);
        hMomentum_DisToBorder_lower[i]->SetTitle("Pt-Photon Plot (DisToBorder = 0, 1, 2); Pt (GeV/c); Fraction of total");
        // Normalize the plot for the lower distance to bad cells
        hMomentum_DisToBorder_lower[i]->Scale(1/hMomentum_DisToBorder_lower[i]->Integral());
        hMomentum_DisToBorder_lower[i]->GetXaxis()->SetLimits(pTmin + bin_offsets[i], pTmax + bin_offsets[i]);
        hMomentum_DisToBorder_lower[i]->Draw("same");
    }
    myMarkerText(0.42, 0.87, colors[0], marker_styles[0], "DisToBorder=0", 1);
    myMarkerText(0.42, 0.82, colors[1], marker_styles[1], "DisToBorder=1", 1);
    myMarkerText(0.42, 0.77, colors[2], marker_styles[2], "DisToBorder=2", 1);
    myMarkerText(0.42, 0.72, kBlack, 21, "DisToBorder>2", 1);
    logcanvas->SaveAs(str_concat_converter(directory_name, "pt_photon_distoborder.png"));
    logcanvas->Clear();
    
    // Get a histogram of the ratio of the normalized photon count with a distance to border over two
    // to the normalized photon count with a distance to bad cells of of 0, 1, and 2
    logcanvas->cd();
    // Reassign chisquares and num_of_bins parameters to default values
    for (int k = 0; k < num_of_distances; k++)
        chisquares[k] = 0;
    for (int k = 0; k < num_of_distances; k++)
        num_of_bins[k] = hMomentum_DisToBorder_upper->GetSize();
    for (int i = 0; i < num_of_distances; i++) {
        TH1D* hRatio = (TH1D*) hMomentum_DisToBorder_upper->Clone("ratio");
        //cout << "\nDistance to border = " << distances_to_border[i] << std::endl;
        for (int j = 0; j < hMomentum_DisToBorder_upper->GetSize(); j++) { // Set the hRatio histogram's values to the photon count ratios
            double new_bin_content = (hMomentum_DisToBorder_lower[i]->GetBinContent(j))/(hMomentum_DisToBorder_upper->GetBinContent(j));
            
            if ((hMomentum_DisToBorder_upper->GetBinContent(j)) == 0) // AVOID DIVIDING BY ZERO
                continue;
            
            //std::cout << "Upper distance is " << hMomentum_DisToBorder_upper->GetBinContent(j) << std::endl;
            //std::cout << "Lower distance is " << hMomentum_DisToBorder_lower[i]->GetBinContent(j) << std::endl;
            //std::cout << "new ratio is " << new_bin_content << std::endl;
            hRatio->SetBinContent(j, new_bin_content);
            hRatio->SetBinError(j, new_bin_content*TMath::Sqrt( ( ((hMomentum_DisToBorder_lower[i]->GetBinError(j))/hMomentum_DisToBorder_lower[i]->GetBinContent(j)) * ((hMomentum_DisToBorder_lower[i]->GetBinError(j))/hMomentum_DisToBorder_lower[i]->GetBinContent(j)) ) + ( ((hMomentum_DisToBorder_upper->GetBinError(j))/hMomentum_DisToBorder_upper->GetBinContent(j)) * ((hMomentum_DisToBorder_upper->GetBinError(j))/hMomentum_DisToBorder_upper->GetBinContent(j)) ) ));
            
            // Now, add this bin's contribution to the chi-square, provided that the bin isn't a division by zero, in which case decrement the bin count by 1 to show that the bin is not counted
            //std::cout << "This bin chi2: " << (new_bin_content - ratio1->GetParameter(0)) * (new_bin_content - ratio1->GetParameter(0))/(ratio1->GetParameter(0)) << std::endl;
            if (hMomentum_DisToBorder_upper->GetBinContent(j) == 0)
                num_of_bins[i]--;
            else
                chisquares[i] += (new_bin_content - ratio1->GetParameter(0)) * (new_bin_content - ratio1->GetParameter(0))/((hRatio->GetBinError(j)) * (hRatio->GetBinError(j)));
        }
        // Draw the histogram
        hRatio->SetMarkerColor(colors[i]);
        hRatio->SetMarkerStyle(marker_styles[i]);
        hRatio->SetLineColor(colors[i]);
        hRatio->SetTitle(Form("Normalized Ratios of DisToBorder>2 to DisToBorder=%i photons; Pt (GeV/c); Ratio", distances_to_border[i]));
        hRatio->SetAxisRange(0.5, 1.4, "Y");
        hRatio->GetXaxis()->SetLimits(pTmin + bin_offsets[i], pTmax + bin_offsets[i]);
        if(i == 0)
            hRatio->Draw();
        else
            hRatio->Draw("same");
    }
    // Insert the ratio = 1 line, and save the combined graph as a file with a legend (include reduced chi-square in the legend)
    ratio1->Draw("same");
    myText(.05,.95, kBlack, "Normalized dBorder>2/dBorder=0,1,2 Photon Ratios");
    myMarkerText(0.19, 0.43, colors[0], marker_styles[0], "dBorder=0/dBorder>2 ratio", 1);
    myText(0.19, 0.38, kBlack, Form("chi2_v = %3.1f; p-val = %2.2f", chisquares[0]/(num_of_bins[0] - 1), TMath::Prob(chisquares[0], (num_of_bins[0] - 1))));
    myMarkerText(0.19, 0.34, colors[1], marker_styles[1], "dBorder=1/dBorder>2 ratio:", 1);
    myText(0.19, 0.29, kBlack, Form("chi2_v = %3.1f; p-val = %2.2f", chisquares[1]/(num_of_bins[1] - 1), TMath::Prob(chisquares[1], (num_of_bins[1] - 1))));
    myMarkerText(0.19, 0.25, colors[2], marker_styles[2], "dBorder=2/dBorder>2 ratio;", 1);
    myText(0.19, 0.20, kBlack, Form("chi2_v = %3.1f; p-val = %2.2f", chisquares[2]/(num_of_bins[2] - 1), TMath::Prob(chisquares[2], (num_of_bins[2] - 1))));
    myText(0.15, 0.16, kBlack, "#scale[1.5]{...} Ratio = 1 fit");
    logcanvas->SaveAs(str_concat_converter(directory_name, "photon_ratios_distoborder.png"));
    logcanvas->Clear();
    
    logcanvas->Close();
    //canvas->Close();
}
