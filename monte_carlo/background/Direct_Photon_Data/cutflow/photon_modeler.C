/**
   This macro is responsible for graphing a spectrum of the pion data from THnSparses_080117_MC.root
   It also outputs data on the integrals to a .tex file for representation on a table
*/
// Programmer: Ivan Chernyshev

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

const int axis_photonPt           =  3;
const int axis_photonPseudoRapidity= 4;
const int axis_photonLambda       =  6;
const int axis_photonNcells       =  7;
const int axis_photonDisToBorder  =  9;
const int axis_photonDisToBadCell = 10;
const int axis_photonDisToCharged = 11;
const int axis_photonExoticity    = 14;
const int axis_photonTime         = 15;
const int axis_photonIsolation    = 18;

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
void photon_modeler() {
    // Initialize constants and the .tex table file headers and footers, both for the test .tex table, which includes both the number of photons and the percentage of the previous amount, and the main file, which includes only the percentage of the previous amount.
    const int num_of_cuts = 8;
    const string parameter_names[num_of_cuts] = {"lambda", "distance to charged", "distance to border", "distance to bad cells", "number of cells", "exoticity", "time", "isolation"};
    const string title_cut_names_line1[num_of_cuts] = {"Cuts: .1<lambda<.4", "Cuts: .1<lambda<.4, dR>.02", "Cuts: .1<lambda<.4, dR>.02, dBorder>0", "Cuts: .1<lambda<.4, dR>.02, dBorder>0, dBadCells>1", "Cuts: .1<lambda<.4, dR>.02, dBorder>0, dBadCells>1, Ncells>1", "Cuts: .1<lambda<.4, dR>.02, dBorder>0, dBadCells>1, Ncells>1, exo<.97", "Cuts: .1<lambda<.4, dR>.02, dBorder>0, dBadCells>1, Ncells>1, exo<.97, |t|<30ns", "Cuts: .1<lambda<.4, dR>.02, dBorder>0, dBadCells>1, Ncells>1, exo<.97, |t|<30ns,"};
    const string title_cut_names_line2[num_of_cuts] = {"", "", "", "", "", "", "", "isolation<4GeV" };
    const string table_line_headers[num_of_cuts] = {"\n+$0.1 < \\lambda <$ 0.4", "\n+dR $> 20$ mrad", "\n+DisToBorder $>$ 0", "\n+DisToBadCell $>$ 1", "\n+Ncells $> 1$", "\n+Exoticity $< 0.97$", "\n+$|Time| < 30$ ns", "\n+isolation $< 4 GeV/c$"};
    const string file_header = "\\documentclass{beamer} \n\\usepackage{graphicx} \n\n\\title{Cut Data Tables} \n\\author{Ivan Chernyshev} \n\\date{\\today} \n\n\\begin{document} \n\n\\frame \n{ \n\\frametitle{Analysis of Cuts: Pt 10-50 GeV, Jet-Jet Monte-Carlo simulation} \n\\begin{table} \n\\caption{How cuts affect number of photons} \n\\centering \n\\begin{tabular}{c c c c} \n\\hline\\hline ";
    const string table_header_test = "\nCuts & Num of Photons & Percentage of Previous\\\\ [0.5ex] \n\\hline";
    const string table_header_main = "\nCuts & Percentage of Previous\\\\ [0.5ex] \n\\hline";
    const string table_footer = "\n[1ex] \n\\hline \n\\end{tabular} \n\\label{table:nonlin} \n\\end{table} \n Final value's portion of the original: ";
    const string file_footer = "\n } \n\\end{document}";
    
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Get the momentum-photon data file
    TFile* fIn = new TFile("THnSparses_080117_MC.root", "READ");
    THnSparse* hPhoton = 0;
    fIn->GetObject("h_Cluster", hPhoton);
    TCanvas* canvas = new TCanvas();
    canvas->SetLogy();
    
    // Set up the output file, where all of the cut parameter plots and cutflow plots will go
    TFile* fOut = new TFile("background_gamma_cutflow.root", "RECREATE");

    // Do the baseline cut, then make an array of parameters for the other cuts
    SetCut(hPhoton, axis_photonPt, 10, 50);
    double cut_params[num_of_cuts][3] = {{axis_photonLambda, 0.1, 0.4}, {axis_photonDisToCharged, 0.02, 0.15}, {axis_photonDisToBorder, 0.9, 6}, {axis_photonDisToBadCell, 1.9, 10}, {axis_photonNcells, 1.9, 30.0}, {axis_photonExoticity, 0, 0.97}, {axis_photonTime, -30, 30}, {axis_photonIsolation, 0, 5}};
    
    // Graph the pT-photon data, calculate the total number of photons in the set, and store this data in a string intended for use in a .tex file with a table
    TH1D* hpT = hPhoton->Projection(axis_photonPt);
    hpT->SetTitle("pT-Photon Spectrum (No cuts); pT (GeV/c); Number of Entries");
    hpT->Draw();
    myText(.45, .92, kBlack, "pT-Photon Spectrum (No cuts)");
    canvas->SaveAs("pT_photon_0_cuts.png");
    hpT->Write("0_cumulativecuts_spectrum");
    string table_lines_main = "\nNone & n/a\\\\";
    string table_lines_test = Form("\nNone & %4.0f & n/a\\\\", hpT->Integral());
    
    // Store the input for this table entry in an array for comparison to the values in the next data entries
    double previous_table_entries[num_of_cuts + 1] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    previous_table_entries[0] = hpT->Integral();
    
    // Repeat for every consequent cut added to the baseline, in the following order: 0.1<lambda<0.4,  DisToCharged>0.02, DisToBorder>0, DisToBadCells>1, Ncells>1
    // For each table entry, also include the percentage of the previous number of photons that the new number of photons is
    for (int i = 0; i < num_of_cuts; i++) {
        canvas->Clear();
        // Before doing any non-baseline cut, plot the parameter to be cut vs photons
        TH1D* hParam = hPhoton->Projection(cut_params[i][0]);
        hParam->SetTitle(Form("%s-photon Spectrum; %s; Number of Entries", parameter_names[i].c_str(), parameter_names[i].c_str()));
        if (i != 6) { // If the parameter is not time, make a separate canvas and plot on normal scale
            TCanvas* normalcanvas = new TCanvas();
            normalcanvas->cd();
            //hParam->GetYaxis()->SetRangeUser(1, 30000);
            hParam->Draw();
            hParam->Write(Form("cut_parameter_%i", i + 1));
            myText(.20, .97, kBlack, Form("%s-photon Spectrum", parameter_names[i].c_str()));
            if (i == 0)
                myText(.20, .92, kBlack, Form("Cuts: None"));
            else {
                myText(.02, .92, kBlack, Form("%s", title_cut_names_line1[i - 1].c_str()));
                myText(.5, .87, kBlack, Form("%s", title_cut_names_line2[i - 1].c_str()));
            }
            normalcanvas->SaveAs(Form("%s_photon.png", parameter_names[i].c_str()));
            normalcanvas->Close();
            canvas->cd();
        }
        else {
            //hParam->GetYaxis()->SetRangeUser(1, 30000);
            hParam->Draw();
            hParam->Write(Form("cut_parameter_%i", i + 1));
            myText(.20, .97, kBlack, Form("%s-photon Spectrum", parameter_names[i].c_str()));
            if (i == 0)
                myText(.20, .92, kBlack, Form("Cuts: None"));
            else {
                myText(.02, .92, kBlack, Form("%s", title_cut_names_line1[i - 1].c_str()));
                myText(.5, .87, kBlack, Form("%s", title_cut_names_line2[i - 1].c_str()));
            }
            canvas->SaveAs(Form("%s_photon.png", parameter_names[i].c_str()));
            canvas->Clear();
        }
        //fIn->GetObject("h_Cluster", hPhoton);
        // Cuts
        SetCut(hPhoton, cut_params[i][0], cut_params[i][1], cut_params[i][2]);
        /**if (i == 6) {
            hParam = hPhoton->Projection(axis_photonTime);
            hParam->Draw();
            canvas->SaveAs("after-cut-time.png");
        }*/
        
        hpT = hPhoton->Projection(axis_photonPt);
        hpT->SetTitle("pT-Photon Spectrum; pT (GeV/c); Number of Entries");
        //hpT->GetYaxis()->SetRangeUser(1, 1000000);
        hpT->Draw();
        hpT->Write(Form("%i_cumulativecuts_spectrum", i + 1));
        myText(.20, .97, kBlack, "pT-Photon Spectrum");
        myText(.02, .92, kBlack, Form("%s", title_cut_names_line1[i].c_str()));
        myText(.5, .87, kBlack, Form("%s", title_cut_names_line2[i].c_str()));
        canvas->SaveAs(Form("pT_photon_%i_cuts.png", i + 1));
        // Put a "not applicable" message at the time cut, as time isn't well-modeled in this simulation.
        if (i == 6) {
            table_lines_test += Form("%s & n/a & (time is not well-modeled) ", table_line_headers[i].c_str());
            table_lines_test += " \\\\";
            table_lines_main += Form("%s & n/a (time is not well-modeled) ", table_line_headers[i].c_str());
            table_lines_main += " \\\\";
        }
        else {
            table_lines_test += Form("%s & %4.0f & %2.1f ", table_line_headers[i].c_str(), hpT->Integral(), (100*hpT->Integral())/previous_table_entries[i]);
            table_lines_test += "$\\%$ \\\\";
            table_lines_main += Form("%s & %2.1f ", table_line_headers[i].c_str(), (100*hpT->Integral())/previous_table_entries[i]);
            table_lines_main += "$\\%$ \\\\";
        }
        
        // Store the input for this table entry in an array for comparison to the values in the next data entries (unless the cut parameter is time, in which case just transfer the previous vaue
        if (i != 6)
            previous_table_entries[i + 1] = hpT->Integral();
        else
            previous_table_entries[i + 1] = previous_table_entries[i];
    }
    
    // Output to the .tex files, include the ratio between the initial and final values
    ofstream table_file_output_test;
    table_file_output_test.open("table_file_test.tex");
    table_file_output_test << file_header << table_header_test << table_lines_test << table_footer << Form("%2.1f",100*previous_table_entries[num_of_cuts]/previous_table_entries[0]) << "\\%" << file_footer;
    table_file_output_test.close();
    
    ofstream table_file_output_main;
    table_file_output_test.open("table_file_main.tex");
    table_file_output_test << file_header << table_header_main << table_lines_main << table_footer << Form("%2.1f",100*previous_table_entries[num_of_cuts]/previous_table_entries[0]) << "\\%" << file_footer;
    table_file_output_main.close();
    
    canvas->Close();
    TH1D* hParam = hPhoton->Projection(axis_photonTime);
    hParam->Draw();
}
