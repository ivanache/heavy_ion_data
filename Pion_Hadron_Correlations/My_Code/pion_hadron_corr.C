// This macro is my first attempt at pion-hadron correlations, inspired by the sample that Miguel sent me, Plotting.C in the Pion_hadron_correlations/Miguel_sample_code
//Author: Ivan Chernyhsev
//Date: 10/19/17

#include "TList.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TVirtualFitter.h"
#include <iostream>

#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

//variables of hPion
const int axis_pion_Cen         = 0;
const int axis_pion_Zvtx        = 1;
const int axis_pionMass         = 2;
const int axis_pionPt           = 3;
const int axis_pionRapidity     = 4;
const int axis_pion_asymmetry   = 5;
const int axis_pion_Ph1_Pt      = 6;
const int axis_pion_Ph2_Pt      = 7;
const int axis_pionOpeningAngle = 8;
const int axis_pion_Ph1_lambda02= 9;
const int axis_pion_Ph2_lambda02= 10;
const int axis_pion_Ph1_dR      = 11;
const int axis_pion_Ph2_dR      = 12;

//variables of hPionTrack
const int axis_corr_centrality  = 0;
const int axis_corr_zvertex     = 1;
const int axis_corr_triggerpT   = 2;
const int axis_corr_trackpT     = 3;
const int axis_corr_dphi        = 4;
const int axis_corr_deta        = 5;
const int axis_corr_zt          = 6;
const int axis_corr_xi          = 7;
const int axis_corr_mass        = 8;
const int axis_corr_photon1Pt   = 9;
const int axis_corr_photon2Pt   = 10;
const int axis_corr_photon1M02  = 11;
const int axis_corr_photon2M02  = 12;

//variables of h_Cluster
const int axis_Cluster_RunNumber = 0;
const int axis_Cluster_Cen  = 1;
const int axis_Cluster_Zvtx = 2;
const int axis_Cluster_Pt   = 3;
const int axis_Cluster_Eta  = 4;
const int axis_Cluster_Phi  = 5;
const int axis_Cluster_M02  = 6;
const int axis_Cluster_Ncells = 7;
const int axis_Cluster_Nmaxima = 8;
const int axis_Cluster_DisToBorder = 9;
const int axis_Cluster_DisToBad = 10;
const int axis_Cluster_dR = 11;
const int axis_Cluster_deta = 12;
const int axis_Cluster_dphi = 13;
const int axis_Cluster_Exoticity = 14;
const int axis_Cluster_time = 15;
const int axis_Cluster_nTracks =16;

// Concatenates two strings and gives a char array
char* str_concat_converter(string str1, string str2){
    string sumstring = str1 + str2;
    char* output = new char[sumstring.length() + 1];
    strcpy(output, sumstring.c_str());
    return output;
}

// Cutting function
void SetCut(THnSparse* h, const int axis, double min, double max){
    
    int binmin = h->GetAxis(axis)->FindBin(min);
    int binmax = h->GetAxis(axis)->FindBin(max);
    h->GetAxis(axis)->SetRange(binmin, binmax-1);
    return;
}
/**
// Graphing and labeling function for TH2D histograms
void graph2D(TH2D* graph, string title, string ylabel, string xlabel, double xseparation, double yseparation, TCanvas* canvas, string graph_option) {
    // Set all titles, set all title separations
    graph->SetTitle(Form("%s; %s; %s", title.c_str(), xlabel.c_str(), ylabel.c_str()));
    graph->GetYaxis()->SetTitleOffset(yseparation);
    graph->GetXaxis()->SetTitleOffset(xseparation);
    
    // Draw the graph
    graph->Draw(graph_option.c_str());
}
*/
// Graphing and labeling function for histograms
void graph(TH1* graph, string title, string ylabel, string xlabel, double xseparation, double yseparation, TCanvas* canvas, string graph_option= "") {
    // Set all titles, set all title separations
    graph->SetTitle(Form("%s; %s; %s", title.c_str(), xlabel.c_str(), ylabel.c_str()));
    graph->GetYaxis()->SetTitleOffset(yseparation);
    graph->GetXaxis()->SetTitleOffset(xseparation);
    
    // Draw the graph
    graph->Draw(graph_option.c_str());
}


// 2D histogram dividing function
// Precondition: the two histograms must have the same y- and x-dimensions
TH2D* divide_histograms2D(TH2D* graph1, TH2D* graph2){
    // Make the 2D histogram to contain the quotient and find minimum and maximum bins along both axes
    TH2D* quotient = new TH2D(*graph1);
    double x_bin_min = quotient->GetXaxis()->FindBin(quotient->GetXaxis()->GetXmin());
    double x_bin_max = quotient->GetXaxis()->FindBin(quotient->GetXaxis()->GetXmax());
    double y_bin_min = quotient->GetYaxis()->FindBin(quotient->GetYaxis()->GetXmin());
    double y_bin_max = quotient->GetYaxis()->FindBin(quotient->GetYaxis()->GetXmax());
    
    // Loop over all bins, divide the element from graph 1 by its counterpart on graph 2. Manually propogate the error
    for(int i = x_bin_min; i <= x_bin_max; i++)
        for(int j = y_bin_min; j <= y_bin_max; j++) {
            if (graph2->GetBinContent(i, j) != 0)
                quotient->SetBinContent(i, j, (graph1->GetBinContent(i, j))/(graph2->GetBinContent(i, j)));
            // Failsafe
            else
                quotient->SetBinContent(i, j, 0);
            
            // Error propogation
            if (graph2->GetBinContent(i, j) == 0) {
                if (graph1->GetBinContent(i, j) == 0)
                    quotient->SetBinError(i, j, 0);
                else
                    quotient->SetBinError(i, j, (quotient->GetBinContent(i, j)*TMath::Sqrt(( (graph1->GetBinError(i, j)/graph1->GetBinContent(i, j)) * (graph1->GetBinError(i, j)/graph1->GetBinContent(i, j)) ) )));
            }
            else if (graph1->GetBinContent(i, j) == 0)
                quotient->SetBinError(i, j, (quotient->GetBinContent(i, j)*TMath::Sqrt(( (graph2->GetBinError(i, j)/graph2->GetBinContent(i, j)) * (graph2->GetBinError(i, j)/graph2->GetBinContent(i, j)) ) )));
            else
                quotient->SetBinError(i, j, (quotient->GetBinContent(i, j)*TMath::Sqrt(( (graph1->GetBinError(i, j)/graph1->GetBinContent(i, j)) * (graph1->GetBinError(i, j)/graph1->GetBinContent(i, j)) ) + ( (graph2->GetBinError(i, j)/graph2->GetBinContent(i, j)) * (graph2->GetBinError(i, j)/graph2->GetBinContent(i, j)) ) )));
    }
    
    
    return quotient;
}

// Projects a 2D histogram onto a 1D histogram by integrating its y-axis from a minimum to a maximum
TH1D* project_2Dhistogram(TH2D* hist2D, double ymin, double ymax) {
    // Make the 1D histogram by calling the project function of the TH2D
    TH1D* projection = hist2D->ProjectionX();
    
    //Get the x-min and x-max bin numbers
    double x_bin_min = hist2D->GetXaxis()->FindBin(hist2D->GetXaxis()->GetXmin());
    double x_bin_max = hist2D->GetXaxis()->FindBin(hist2D->GetXaxis()->GetXmax());
    //Get the y-min and y-max bin numbers
    double y_bin_min = hist2D->GetYaxis()->FindBin(ymin);
    double y_bin_max = hist2D->GetYaxis()->FindBin(ymax);
    double* error;
    
    // Loop over the x-bins from x_bin_min to x_bin_max, replace each term in projection with the within-bin integral from
    for(int i = x_bin_min; i <= x_bin_max; i++) {
        projection->SetBinContent(i, hist2D->IntegralAndError(i, i, y_bin_min, y_bin_max, *error));
        projection->SetBinError(i, *error);
        
    }
    
    return projection;
}

// The 2Gaussians + constant fit
// Must be called with one independent variable x and 7 parameters par
double Two_Gaussian_fit(Double_t* x, Double_t* par) {
    // Define which parameter is which
    double A = par[0];
    
    double B = par[1];
    double mean_1 = par[2];
    double sigma_1 = par[3];
    
    double C = par[4];
    double mean_2 = par[5];
    double sigma_2 = par[6];
    
    //Set the exponents of the two Gaussian functios
    double arg1;
    if (sigma_1 != 0)
        arg1 = -0.5*((x[0] - mean_1)/sigma_1)*((x[0] - mean_1)/sigma_1);
    
    double arg2;
    if (sigma_2 != 0)
        arg2 = -0.5*((x[0] - mean_2)/sigma_2)*((x[0] - mean_2)/sigma_2);
    
    double fitval = A + B*TMath::Exp(arg1)/TMath::Sqrt(2*TMath::Pi()*sigma_1*sigma_1) + C*TMath::Exp(arg2)/TMath::Sqrt(2*TMath::Pi()*sigma_2*sigma_2);
    return fitval;
}
const int Two_Gaussian_params = 7;

// Main function
// Must be called with the minimum and maximum values for the following parameters: trigger pT, mass, track pT. All are in GeV except the mass-related terms
// Optional parameter of directory name
void pion_hadron_corr(double triggerpT_min, double triggerpT_max, double mass_min, double mass_max, double trackpT_min, double trackpT_max, string directory_name = "") {
    TCanvas* canvas = new TCanvas();

    
    // Import data
    TFile* input = new TFile("THnSparses_LHC13d_101517.root", "READ");
    input->Print();
    
    // Get the THnSparses
    THnSparse* hPionTrack = 0;
    input->GetObject("h_PionTrack", hPionTrack);
    THnSparse* hPionTrack_Mixed = 0;
    input->GetObject("h_PionTrack_Mixed", hPionTrack_Mixed);
    
    //Graph the trackPt curve for both THnSparses, then set atlas style and the sumw2 error evaluator and format the directory name string
    TH1D* trackPt_curve = hPionTrack->Projection(axis_corr_trackpT);
    trackPt_curve->SetTitle("h_PionTrack Track Spectrum; Track Pt (GeV); Counts");
    trackPt_curve->Draw();
    canvas->SaveAs("h_PionTrack_trackspectrum.png");
    canvas->Clear();
    
    TH1D* trackPtMixed_curve = hPionTrack_Mixed->Projection(axis_corr_trackpT);
    trackPtMixed_curve->SetTitle("h_PionTrack_Mixed Track Spectrum; Track Pt (GeV); Counts");
    trackPtMixed_curve->Draw();
    canvas->SaveAs("h_PionTrackMixed_trackspectrum.png");
    canvas->Clear();
    
    SetAtlasStyle();
    
    if (directory_name != "")
        directory_name += "/";
    
    // Cut the pion pT of both THnSparses to 10-12 GeV, the mass to 110-150 MeV, and track pT to 1-2 GeV
    SetCut(hPionTrack, axis_corr_triggerpT, triggerpT_min, triggerpT_max);
    SetCut(hPionTrack, axis_corr_mass, mass_min/1000, mass_max/1000);
    SetCut(hPionTrack, axis_corr_trackpT, trackpT_min, trackpT_max);
    
    SetCut(hPionTrack_Mixed, axis_corr_triggerpT, triggerpT_min, triggerpT_max);
    SetCut(hPionTrack_Mixed, axis_corr_mass, mass_min/1000, mass_max/1000);
    SetCut(hPionTrack_Mixed, axis_corr_trackpT, trackpT_min, trackpT_max);
    
    // Make a 2D projection over both delta-phi and delta-eta for both THnSparses
    TH2D* Pion_Track_Projection = hPionTrack->Projection(axis_corr_deta, axis_corr_dphi);
    TH2D* Pion_Track_Mixed_Projection = hPionTrack_Mixed->Projection(axis_corr_deta, axis_corr_dphi);
    
    // Output the 2D projections
    graph(Pion_Track_Projection, "Pion Track", "#Delta #eta", "#Delta #phi [rad]", 1.0, 1.0, canvas, "COLZ");
    myText(.40,.92, kBlack, "Pion Track");
    canvas->SaveAs(str_concat_converter(directory_name, "pion_track_graph.png"));
    canvas->Clear();
    
    graph(Pion_Track_Mixed_Projection, "Mixed Pion Track", "#Delta #eta", "#Delta #phi [rad]", 1.0, 1.0, canvas, "COLZ");
    myText(.40,.92, kBlack, "Mixed Pion Track");
    canvas->SaveAs(str_concat_converter(directory_name, "mixed_pion_track_graph.png"));
    canvas->Clear();
    
    // Get and graph the quotient of all bins from the pion track divided by all bins from the mixed pion track
    // This is the correlation function
    // Use both a surface plot and a 2D intensity chart
    TH2D* correlation_function = divide_histograms2D(Pion_Track_Projection, Pion_Track_Mixed_Projection);
    
    graph(correlation_function, "Correlation Function", "#Delta #eta", "#Delta #phi [rad]", 1.0, 1.0, canvas, "COLZ");
    myText(.40,.92, kBlack, "Correlation Function");
    canvas->SaveAs(str_concat_converter(directory_name, "correlation_function_intensitychart.png"));
    
    graph(correlation_function, "Correlation Function", "#Delta #eta", "#Delta #phi [rad]", 1.0, 1.0, canvas, "SURF2");
    myText(.40,.92, kBlack, "Correlation Function");
    canvas->SaveAs(str_concat_converter(directory_name, "correlation_function_surfaceplot.png"));
    
    // Get the projection of the correlation function over |delta eta| < 0.8
    TH1D* correlation_projection = project_2Dhistogram(correlation_function, -0.8, 0.8);
    
    // Fit the projection to a 2Gaussians + constant curve
    TF1* fitfunc = new TF1("fit", Two_Gaussian_fit, -0.4, 1.5, Two_Gaussian_params);
    fitfunc->SetParNames("Constant", "Magnitude 1", "Mean 1", "Sigma 1", "Magnitude 2", "Mean 2", "Sigma 2");
    fitfunc->SetParameters(150, 225, 0, 0.1, 75, 1, 0.3);
    correlation_projection->Fit(fitfunc);
    
    // Graph the projection
    graph(correlation_projection, "Correlation Function: Projection over |#Delta #eta| < 0.8", "Correlation ratio", "#Delta #phi [rad]", 1.0, 1.0, canvas);
    myText(.20,.92, kBlack, "Correlation Function: Projection over |#Delta #eta| < 0.8");
    // Label regarding pion pT, mass, track pT cuts
    myText(.6, 0.7, kBlack, "#scale[0.7]{Param Borders}");
    myText(.6, 0.67, kBlack, Form("#scale[0.5]{%2.0f GeV < #pi^{0} pT < %2.0f GeV}", triggerpT_min, triggerpT_max));
    myText(.6, 0.64, kBlack, Form("#scale[0.5]{%3.0f MeV < #pi^{0} mass < %3.0f MeV}", mass_min, mass_max));
    myText(.6, 0.61, kBlack, Form("#scale[0.5]{%2.0f GeV < Track pT < %2.0f GeV}", trackpT_min, trackpT_max));
    //fitfunc->Print();
    // myText(.35, 0.81, kBlack, "#scale[0.75]{Fit Function: A + #frac{B}{#sqrt{2 #pi #sigma_{1}^{2}}} e^{#frac{(x-#bar{x}_{1})^{2}}{2 #sigma_{1}^{2}}} + #frac{C}{#sqrt{2 #pi #sigma_{2}^{2}}} e^{#frac{(x-#bar{x}_{2})^{2}}{2 #sigma_{2}^{2}}}}");
    //myText(.18, 0.81, kBlack, Form("#scale[0.7]{Fit Function: %4.0f + #frac{%4.1f}{#sqrt{2 #pi %0.4f^{2}}} e^{#frac{(x-%2.2f)^{2}}{2*%0.4f^{2}}} + #frac{%4.1f}{#sqrt{2 #pi %0.3f^{2}}} e^{#frac{(x-%2.2f)^{2}}{2 * %0.3f^{2}}}}", fitfunc->GetParameter(0), fitfunc->GetParameter(1), fitfunc->GetParameter(3), fitfunc->GetParameter(2), fitfunc->GetParameter(3), fitfunc->GetParameter(4), fitfunc->GetParameter(6), fitfunc->GetParameter(5), fitfunc->GetParameter(6)));
    canvas->SaveAs(str_concat_converter(directory_name, "correlation_function_projection.png"));
    
    canvas->Close();
}
