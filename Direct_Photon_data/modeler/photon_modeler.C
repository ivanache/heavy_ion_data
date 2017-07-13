/**
   This macro is responsible for graphing a spectrum of the pion data from THnSparses_071217.root
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
const int axis_photonLambda       =  6;
const int axis_photonNcells       =  7;
const int axis_photonDisToBorder  =  9;
const int axis_photonDisToBadCell = 10;
const int axis_photonDisToCharged = 11;
const int axis_photonExoticity    = 12;
const int axis_photonTime         = 13;

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
    // Initialize constants and the .tex table file header and footer
    const int num_of_cuts = 7;
    const string parameter_names[num_of_cuts] = {"lambda", "distance to charged", "distance to border", "distance to bad cells", "number of cells", "exoticity", "time"};
    const string title_cut_names[num_of_cuts] = {"Cuts: .1<lambda<.4", "Cuts: .1<lambda<.4, dR>.02", "Cuts: .1<lambda<.4, dR>.02, dBorder>0", "Cuts: .1<lambda<.4, dR>.02, dBorder>0, dBadCells>1", "Cuts: .1<lambda<.4, dR>.02, dBorder>0, dBadCells>1, Ncells>1", "Cuts: .1<lambda<.4, dR>.02, dBorder>0, dBadCells>1, Ncells>1, exo<.97", "Cuts: .1<lambda<.4, dR>.02, dBorder>0, dBadCells>1, Ncells>1, exo<.97, |t|<30ns"};
    const string table_line_headers[num_of_cuts] = {"\n+$0.1 < \\lambda <$ 0.4", "\n+dR $> 20$ mrad", "\n+DisToBorder $>$ 0", "\n+DisToBadCell $>$ 1", "\n+Ncells $> 1$", "\n+Exoticity $< 0.97$", "\n+$|Time| < 30$ ns"};
    const string table_header = "\\documentclass{beamer} \n\\usepackage{graphicx} \n\n\\title{Cut Data Tables} \n\\author{Ivan Chernyshev} \n\\date{\\today} \n\n\\begin{document} \n\n\\frame \n{ \n\\frametitle{Analysis of Cuts: Pt 10-20 GeV} \n\\begin{table} \n\\caption{How cuts affect number of photons} \n\\centering \n\\begin{tabular}{c c c c} \n\\hline\\hline \nCuts & Num of Photons\\\\ [0.5ex] \n\\hline";
    const string table_footer = "\n[1ex] \n\\hline \n\\end{tabular} \n\\label{table:nonlin} \n\\end{table} \n} \n\\end{document}";
    
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    // Get the momentum-photon data file
    TFile* fIn = new TFile("THnSparses_071217.root", "READ");
    THnSparse* hPhoton = 0;
    fIn->GetObject("h_Cluster", hPhoton);
    TCanvas* canvas = new TCanvas();

    // Do the baseline cut, then make an array of parameters for the other cuts
    SetCut(hPhoton, axis_photonPt, 10, 20);
    double cut_params[num_of_cuts][3] = {{axis_photonLambda, 0.1, 0.4}, {axis_photonDisToCharged, 0.02, 0.15}, {axis_photonDisToBorder, 0.9, 6}, {axis_photonDisToBadCell, 1.9, 10}, {axis_photonNcells, 1.9, 30.0}, {axis_photonExoticity, 0, 0.97}, {axis_photonTime, -30, 30}};
    
    // Graph the pT-photon data, calculate the total number of photons in the set, and store this data in a string intended for use in a .tex file with a table
    TH1D* hpT = hPhoton->Projection(axis_photonPt);
    hpT->SetTitle("pT-Photon Spectrum (No cuts); pT (GeV/c); Number of Entries");
    hpT->Draw();
    myText(.45, .92, kBlack, "pT-Photon Spectrum (No cuts)");
    canvas->SaveAs("pT_photon_0_cuts.png");
    string table_lines = Form("\nNone & %4.0f\\\\", hpT->Integral());
    
    // Repeat for every consequent cut added to the baseline, in the following order: 0.1<lambda<0.4,  DisToCharged>0.02, DisToBorder>0, DisToBadCells>1, Ncells>1
    for (int i = 0; i < num_of_cuts; i++) {
        canvas->Clear();
        // Before doing any non-baseline cut, plot the parameter to be cut vs photons
        TH1D* hParam = hPhoton->Projection(cut_params[i][0]);
        hParam->SetTitle(Form("%s-photon Spectrum; %s; Number of Entries", parameter_names[i].c_str(), parameter_names[i].c_str()));
        if (i == 6) { // If the parameter is time, make a separate canvas and plot on log scale
            TCanvas* logcanvas = new TCanvas();
            logcanvas->SetLogy();
            logcanvas->cd();
            hParam->Draw();
            myText(.20, .97, kBlack, Form("%s-photon Spectrum", parameter_names[i].c_str()));
            if (i == 0)
                myText(.20, .92, kBlack, Form("Cuts: None"));
            else
                myText(.02, .92, kBlack, Form("%s", title_cut_names[i - 1].c_str()));
            logcanvas->SaveAs(Form("%s_photon.png", parameter_names[i].c_str()));
            logcanvas->Clear();
            canvas->cd();
        }
        else {
            hParam->Draw();
            myText(.20, .97, kBlack, Form("%s-photon Spectrum", parameter_names[i].c_str()));
            if (i == 0)
                myText(.20, .92, kBlack, Form("Cuts: None"));
            else
                myText(.02, .92, kBlack, Form("%s", title_cut_names[i - 1].c_str()));
            canvas->SaveAs(Form("%s_photon.png", parameter_names[i].c_str()));
            canvas->Clear();
        }
        SetCut(hPhoton, cut_params[i][0], cut_params[i][1], cut_params[i][2]);
        
        hpT = hPhoton->Projection(axis_photonPt);
        hpT->SetTitle(Form("pT-Photon Spectrum %s; pT (GeV/c); Number of Entries", title_cut_names[i].c_str()));
        hpT->Draw();
        myText(.20, .97, kBlack, "pT-Photon Spectrum");
        myText(.02, .92, kBlack, Form("%s", title_cut_names[i].c_str()));
        canvas->SaveAs(Form("pT_photon_%i_cuts.png", i + 1));
        table_lines += Form("%s & %4.0f\\\\", table_line_headers[i].c_str(), hpT->Integral());
    }
    
    // Output to a .tex file
    ofstream table_file_output;
    table_file_output.open("table_file.tex");
    table_file_output << table_header << table_lines << table_footer;
    table_file_output.close();
    
    canvas->Close();
}
