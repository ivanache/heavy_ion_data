/**
 my_code is a data-processing macro meant for processing THnSparses root files, specifically the h_Pion one in THnSparses_LHC13d_090717.root, THnSparses_LHC13e_090717.root, THnSparses_LHC13f_090717.root
 Programmer: Ivan Chernyshev
 */
#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"
#include "TAxis.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

// The output file
TFile* fOut;

// The general filepath string
string directory_name;

//variables of hPion
const int axis_pion_Cen           = 0;
const int axis_pion_Zvtx          = 1;
const int axis_pionMass           = 2;
const int axis_pionPt             = 3;
const int axis_asymmetry          = 5;
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

// The Gaussian fit-function, for the residual distributions
double gaussian_peak(Double_t *x, Double_t *par) {
    double A = par[0];
    double mean = par[1];
    double sigma = par[2];
    
    double arg;
    if (sigma != 0)
        arg = (x[0] - mean)/sigma;
    
    double fitval = A*TMath::Exp(-0.5*arg*arg)/TMath::Sqrt(2*TMath::Pi()*sigma*sigma);
    return fitval;
}

// The Gaussian fit-function, for simplifying fitting for symmetric peaks
double gaussian_model(Double_t *x, Double_t *par) {
    double A = par[0];
    double mean = par[1];
    double sigma = par[2];
    double B = par[3];
    double C = par[4];
    double D = par[5];
    
    double arg;
    if (sigma != 0)
        arg = (x[0] - mean)/sigma;
    
    double fitval = A*TMath::Exp(-0.5*arg*arg)/TMath::Sqrt(2*TMath::Pi()*sigma*sigma) + B*x[0]*x[0] + C*x[0] + D;
    return fitval;
}

/**
// The crystal ball fitting function
// The peak for the crystal ball function (because it is piecewise)
double crystal_ball_function_peak(Double_t *x, Double_t *par) {
    double A = par[0];
    double mean = par[1];
    double sigma = par[2];
    double alpha = par[3];
    double n = par[4];
    
    double fitval = A*ROOT::Math::crystalball_pdf(x[0], alpha, n, sigma, mean);
    return fitval;
}
// Crystal ball function, quadric background model
double crystal_ball_model_quadric(Double_t *x, Double_t *par) {
    double A = par[0];
    double mean = par[1];
    double sigma = par[2];
    double alpha = par[3];
    double n = par[4];
    double B = par[5];
    double C = par[6];
    double D = par[7];
    double E = par[8];
    double F = par[9];
    
    double fitval = A*ROOT::Math::crystalball_pdf(x[0], alpha, n, sigma, mean) + B*TMath::Power(x[0], 4.0) + C*TMath::Power(x[0], 3.0) + D*x[0]*x[0] + E*x[0] + F;
    return fitval;
}

// Crystal ball function, quadratic background model
double crystal_ball_model_quadratic(Double_t *x, Double_t *par) {
    double A = par[0];
    double mean = par[1];
    double sigma = par[2];
    double alpha = par[3];
    double n = par[4];
    double B = par[5];
    double C = par[6];
    double D = par[7];
    
    double fitval = A*ROOT::Math::crystalball_pdf(x[0], alpha, n, sigma, mean) + B*x[0]*x[0] + C*x[0] + D;
    return fitval;
}
*/

// The three functions used in this program: the first for the entire data, second for the peak alone, third for the background alone
TF1 *func;
TF1 *peak;
TF1 *background;

// Takes the seven parameters that would be passed to func, along with the number of sigmas away from the mean ("x"), and returns the ratio of the integral of
// peak within the interval that is within x sigmas from the mean and the integral of the whole fit over the same integral
double signal_over_total(Double_t *x, Double_t *par) {
    // Extract the parameters
    //func->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]);
    //peak->SetParameters(par[0], par[1], par[2]);
    double Nsigma = x[0];
    // Based on the peak model, derive the mean and sigma
    double mean;
    double sigma;
    mean = par[1];
    sigma = par[2];
    
    // Evaluate the integrals and return their quotient
    double signal = peak->Integral(mean-Nsigma*sigma, mean+Nsigma*sigma);
    double total = func->Integral(mean-Nsigma*sigma, mean+Nsigma*sigma);
    
    return signal/total;
}

// Takes a THnSparse pointer (almost always h_Pion), an hPion variable number, a Canvas object pointer, and a file name
// Graphs the data under the hPion variable number and saves the graph as a .png
// SetCuts done to the THnSparse object before the function call do apply
// Not advisable if you want to save any variable from the data or the TH1D object used to graph it
void graph_raw_data(THnSparse* data, const int hPion_var, TCanvas* can, char* filename, string rootname){
    can->Clear();
    TH1D* data_subset = data->Projection(hPion_var);
    data_subset->Draw();
    myText(.20,.95, kBlack, str_concat_converter("#scale[1]{" + rootname, " vs entries}"));
    can->SaveAs(filename);
    data_subset->Write(rootname.c_str());
    can->Clear();
}

/**
 Main function
 Precondition: runletter, the one-letter string that denotes the run letter is either "d", "e", or "f"
 */
void my_code(string runletter) {
    // Leave the file if runletter is not "d", "e", or "f"
    if (runletter != "d" && runletter != "e" && runletter != "f") {
        cout << "ERROR: run letter must be \"d\", \"e\", or \"f\"" << endl;
    }
    
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    // Headers for various cut nums, for use both for table output (another version is defined in the to-table input section and for titles
    const int cutnum_option_quantity = 5;
    directory_name = Form("data/%s/", runletter.c_str());
    //string headers[cutnum_option_quantity] = {"None ", "angle > 17 mrad ", "0.1 < lambda < 0.4 ", "dR > 20 mrad ", "asymmetry < 0.7 "};
    //if (option == "TESTSENSITIVITY")
        //headers[1] = "angle > 18 mrad";
    
    //Open the files
    TFile* fIn = new TFile(Form("THnSparses_LHC13%s_090717.root", runletter.c_str()),"READ"); //get file
    string rootfilename;
    rootfilename = Form("PionDataOutput_LHC13%s.root", runletter.c_str());
    fOut = new TFile(rootfilename.c_str(), "RECREATE"); // Create an output file
    fIn->Print(); //print file content
    
    // Get the data
    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion); //get array
    
    //The TCanvases and TPads (for sub-canvassing)
    TCanvas* graphcanvas = new TCanvas();
    TPad *pad[2] = {new TPad("pad0","",0,0.36,1,1), new TPad("pad1","",0,0.1,1,0.45)};
    TCanvas* logcanvas = new TCanvas();
    logcanvas->SetLogy();
    
    //For the mass plot, restrict to mass between 0.08 and 0.25 and the momentum to between 6 GeV and 20 GeV
    //Apply all other baseline cuts, adjust the angle cut based on the my_code option
    // Graph the data with respect to all of the parameters while you're at it
    SetCut(h_Pion, axis_pionPt, 6.0, 20.0);
    SetCut(h_Pion, axis_pionMass, 0.1, 0.22);
    
    graph_raw_data(h_Pion, axis_pionDisToCharged1, logcanvas, str_concat_converter(directory_name, "DisToChargedParticle1_pion_plot.png"), "DisToChargedParticle1");
    graph_raw_data(h_Pion, axis_pionDisToCharged2, logcanvas, str_concat_converter(directory_name, "DisToChargedParticle2_pion_plot.png"), "DisToChargdParticle2");
    //SetCut(h_Pion, axis_pionDisToCharged1, 0.02, 0.14);
    //SetCut(h_Pion, axis_pionDisToCharged2, 0.02, 0.14);
    
    graph_raw_data(h_Pion, axis_pionLambda1, graphcanvas, str_concat_converter(directory_name, "lambda1_pion_plot.png"), "lambda1");
    graph_raw_data(h_Pion, axis_pionLambda2, graphcanvas, str_concat_converter(directory_name, "lambda2_pion_plot.png"), "lambda2");
    SetCut(h_Pion, axis_pionLambda1, 0.1, 0.4);
    SetCut(h_Pion, axis_pionLambda2, 0.1, 0.4);
    
    graph_raw_data(h_Pion, axis_pionAngle, graphcanvas, str_concat_converter(directory_name, "angle_pion_plot.png"), "angle");
    //if (option == "DEFAULT")
    SetCut(h_Pion, axis_pionAngle, 17, 60);
    /**else if (option == "TESTSENSITIVITY")
        SetCut(h_Pion, axis_pionAngle, 18.5, 60);
    else
        SetCut(h_Pion, axis_pionAngle, 19.5, 60);*/
    
    // The asymmetry needs a separate graph per momentum interval
    // The 14-16 GeV interval and above are the ones whose asymmetry is not cut at all
    double min;
    double max;
    const int num_of_intervals = 6;
    double intervals[num_of_intervals][2] = {{6.0, 8.0}, {8.0, 10.0}, {10.0, 12.0}, {12.0, 14.0}, {14.0, 16.0}, {16.0, 20.0}};
    graph_raw_data(h_Pion, axis_asymmetry, graphcanvas, str_concat_converter(directory_name,"asymmetry_pion_plot.png"), "asymmetry");
    for (int i = 0; i < num_of_intervals; i++) {
        min = intervals[i][0];
        max = intervals[i][1];
        
        // Cut the momentum to within the interval
        SetCut(h_Pion, axis_pionPt, min, max);
        
        // Create a function, draw a graph, format it, and label it
        TH1D* asymmetry_curve;
        asymmetry_curve = h_Pion->Projection(axis_asymmetry);
        asymmetry_curve->SetTitle(Form("Counts at momenta %2.2f-%2.2f GeV over asymmetry vs entries; Asymmetry; Number of Entries", min, max));
        asymmetry_curve->GetYaxis()->SetTitleOffset(1.0);
        asymmetry_curve->GetXaxis()->SetTitleOffset(1);
        asymmetry_curve->Draw();
        myText(.20,.95, kBlack, Form("#scale[0.75]{Counts at momenta %2.2f-%2.2f GeV over asymmetry vs entries}", min, max));
        
        //Save and clear the canvas for the next use
        graphcanvas->SaveAs(str_concat_converter(directory_name,Form("asymmetry_%2.2f-%2.2f_GeV_pion_plot.png", min, max)));
        asymmetry_curve->Write(Form("asymmetry_%2.2fGeV-%2.2fGeV", min, max));
        graphcanvas->Clear();
    }
    // Set pt cut to what it originally was
    SetCut(h_Pion, axis_pionPt, 6, 20.0);
    // Asymmetry cut will be applied in the individual pt bins
    
    //std::cout << "Cuts: distance to charged particles, lambda, angle, and asymmetry\n\n";
    
    // plot mass data, write it into the root file, and set up the fit function
    TH1D* hMass = h_Pion->Projection(axis_pionMass);
    //hMass->Rebin(2);
    TH1D* residual = (TH1D*)hMass->Clone("residual");
    double MASSWIDTH = hMass->GetBinWidth(1);
    //hMass->SetAxisRange(0., 5000., "Y");
    hMass->GetYaxis()->SetTitle("Number of Entries");
    hMass->GetYaxis()->SetTitleSize(.05);
    hMass->GetYaxis()->SetTitleOffset(0.8);
    graphcanvas->cd();
    pad[0]->Draw();
    pad[0]->cd();
    hMass->Draw();
    hMass->Write("unfitted_mass_pion");
    
    //Start making the fit, restrict the parameters to reasonable ranges, insert guess values, and give understandable names
    const int num_of_params = 6;
    int num_of_peak_params = 3;
    //background = new TF1("background curve", "[0]*TMath::Power(x, 4.0) + [1]*TMath::Power(x, 3.0) + [2]*x*x + [3]*x + [4]", 0.08, 0.26);
    background = new TF1("background curve", "[0]*x*x + [1]*x + [2]", 0.08, 0.26);
    //num_of_peak_params = 5;
    func = new TF1("fit", gaussian_model,0.05,0.5,num_of_params);
    peak = new TF1("mass peak", gaussian_peak, 0.05, 0.5, num_of_peak_params);
    func->SetParNames("Integral", "Mean", "Sigma", "Quadratic coeff", "Linear coeff", "Constant");
    func->SetParameters(40,  0.13, 0.012, 25, 19, -120000, 28000, -1600);
    func->SetParLimits(0, 0.001, 400.0);//integral
    func->SetParLimits(1, 0.11, 0.16); //mean
    func->SetParLimits(2, 0.008, 0.016); // width
    func->SetParLimits(3, -1000000.0, 0.0); // Quadric and quadratic factors
    //func->SetParLimits(7, -10000000.0, 0.0);
    //func->SetParLimits(8, 0.0, 1000000.0); // Linear
    //func->SetParLimits(9, -1000000.0, 0.0); //Constant
    
    // Plot the fit for the mass and (separately) the Gaussian and background components of it, print out the number of pions and the other parameters, write them on the graph
    hMass->Fit(func, "0");
    hMass->Fit(func);
    //hMass->Fit(func);
    func->SetLineColor(kRed);
    std::cout << "\n\nNumber of pions:" << (func->GetParameter(0))/MASSWIDTH << std::endl << "Error:" << (func->GetParError(0))/MASSWIDTH << std::endl;
    std::cout << "\nMean:" << func->GetParameter(1) << std::endl << "Error:" << func->GetParError(1) << std::endl;
    std::cout << "\nSigma:" << func->GetParameter(2) << std::endl << "Error:" << func->GetParError(2) << std::endl;
    //std::cout << "\nAlpha:" << func->GetParameter(3) << std::endl << "Error:" << func->GetParError(3) << std::endl;
    //std::cout << "\nN:" << func->GetParameter(4) << std::endl << "Error:" << func->GetParError(4) << std::endl;
    //std::cout << "\nQuadric:" << func->GetParameter(5) << std::endl << "Error:" << func->GetParError(5) << std::endl;
    //std::cout << "\nCubic:" << func->GetParameter(6) << std::endl << "Error:" << func->GetParError(6) << std::endl;
    std::cout << "\nQuadratic:" << func->GetParameter(3) << std::endl << "Error:" << func->GetParError(5) << std::endl;
    std::cout << "\nLinear:" << func->GetParameter(4) << std::endl << "Error:" << func->GetParError(6) << std::endl;
    std::cout << "\nConstant:" << func->GetParameter(5) << std::endl << "Error:" << func->GetParError(7) << std::endl;
    func->SetLineColor(kRed);
    func->Draw("same");
    
    int i = 0;
    for(; i < num_of_peak_params; i++) {
        peak->SetParameter(i, func->GetParameter(i));
        peak->SetParError(i, func->GetParError(i));
    }
    for(; i < num_of_params; i++) {
        background->SetParameter(i - num_of_peak_params, func->GetParameter(i));
        background->SetParError(i - num_of_peak_params, func->GetParError(i));
    }
    peak->SetLineColor(kBlue);
    //peak->SetFillColor(kBlue);
    peak->Draw("same");
    background->SetLineColor(kGreen);
    //background->SetFillColor(kGreen);
    background->Draw("same");
    hMass->Draw("same");
    hMass->GetListOfFunctions()->Add(peak);
    hMass->GetListOfFunctions()->Add(background);
    hMass->Write("mass_pion");// Load into the ROOT file
    myText(.20,.95, kBlack, "#scale[1]{Mass vs. Entries, Pt 6-20 GeV/c, cuts:}");
    myText(.20,.90, kBlack, "#scale[1]{0.1 < lambda < 0.4, asymmetry < 0.7}");
    myText(.6, .85, kBlack, Form("#scale[0.75]{Mean Mass: %3.3f +/- %3.4f GeV}", func->GetParameter(1), func->GetParError(1)));
    myText(.6, .80, kBlack, Form("#scale[0.75]{Mass Width: %2.3f +/- %2.4f GeV}", func->GetParameter(2), func->GetParError(2)));
    
    // Create a new TH1D to measure the distribution of residuals, set all bin values to zero
    TH1D* residual_dist = new TH1D("residual_distribution", "Residual_Distribution", 10, -5, 5);
    residual_dist->Sumw2();
    residual_dist->SetTitle("Residual Distribution; Residual value; Num of residuals");
    // Load the individual residuals into the residual histogram (the one that goes below the mass-pion graph), count the number of residuals in each bin of the residual histogram
    // In the meantime, sum up the squares of those residuals to make a chi-square
    double chisquare = 0;
    std::cout << "Differences:\n";
    for (int i = 0; i < hMass->GetSize(); i++) {
        
        double residualvalue = (hMass->GetBinContent(i) - func->Eval(hMass->GetBinCenter(i)))/hMass->GetBinError(i);
        if ((hMass->GetBinError(i)) != 0) {
            residual->SetBinContent(i, residualvalue);
        }
        else
            residual->SetBinContent(i, 0);
        residual->SetBinError(i, 0); //Residuals don't have errors
        
        residual_dist->Fill(residualvalue);
        chisquare += residualvalue*residualvalue;
    }
    std::cout << std::endl;
    for(int i = 0; i < residual_dist->GetSize(); i++) {
        std::cout << "Final residual count for " << residual_dist->GetBinCenter(i) << " is " << residual_dist->GetBinContent(i) << std::endl;
    }
    myText(.6, .75, kBlack, Form("#scale[0.75]{Reduced Chi-square: %2.1f/%i}", func->GetChisquare(), (hMass->GetSize() - num_of_params - 1)));
    myText(.6, .70, kBlack, Form("#scale[0.75]{P-val: %2.2f}", TMath::Prob(func->GetChisquare(), (hMass->GetSize() - (num_of_params) - 1))));
    std::cout << "Reduced Chi Square " << func->GetChisquare()/(hMass->GetSize() - (num_of_params) - 1) << std::endl;
    
    // Plot the residual; save as a PDF, print out the individual residuals
    graphcanvas->cd();
    pad[1]->Draw("p");
    pad[1]->cd();
    residual->SetTitle("; Pion Mass (GeV); Residuals");
    residual->SetAxisRange(-4., 4., "Y");
    residual->GetXaxis()->SetTitle("Pion Mass (GeV)");
    residual->GetXaxis()->SetTitleSize(.1);
    residual->GetYaxis()->SetTitleSize(.1);
    residual->GetYaxis()->SetTitleOffset(0.3);
    residual->GetXaxis()->SetTitleOffset(0.3);
    residual->SetMarkerStyle(20);
    residual->SetMarkerSize(1);
    residual->Draw("p");
    residual->Write("residual"); // Load into the ROOT file
    graphcanvas->SaveAs(str_concat_converter(directory_name, "MyFit_.png"));
    
    // Fit and plot the residuals distribution
    graphcanvas->Clear();
    TF1* residual_dist_fit = new TF1("fit", gaussian_peak, -5, 5, 3);
    residual_dist_fit->SetLineColor(2);
    residual_dist_fit->SetParameters(50,  0, 1);
    residual_dist_fit->SetParLimits(0, 0, 10000);
    residual_dist_fit->SetParLimits(1, -10, 10);
    residual_dist_fit->SetParLimits(2, 0, 10);
    residual_dist->Fit(residual_dist_fit);
    residual_dist->Draw();
    residual_dist_fit->Draw("same");
    myText(.20, .92, kBlack, "#scale[1]{Residual Distribution, Pt 6-20 GeV/c}");
    myText(.2, .75, kBlack, Form("#scale[0.75]{Mean: %2.1f+/-%1.1f}", residual_dist_fit->GetParameter(1), residual_dist_fit->GetParError(1)));
    myText(.2, .7, kBlack, Form("#scale[0.75]{Sigma: %2.1f+/-%1.1f}", residual_dist_fit->GetParameter(2), residual_dist_fit->GetParError(2)));
    graphcanvas->SaveAs(str_concat_converter(directory_name, "MyFit_residual_dist.png"));
    
    // Plot signal to noise ratio for the entire sample over various distances from the mean, print out values for 1 sigma and 2 sigmas
    graphcanvas->Clear();
    TF1* sig_over_tot_funct = new TF1("Signal over Total", signal_over_total, 0, 4, num_of_params);
    sig_over_tot_funct->SetParameters(func->GetParameters());
    sig_over_tot_funct->SetTitle("Signal over Total vs Distance From Mean; Num of Standard Deviations From Mean; Signal to Total Ratio");
    sig_over_tot_funct->GetYaxis()->SetRangeUser(0.5, 1.0);
    sig_over_tot_funct->Draw();
    std::cout << "\n\nS/T for 0 sigma: " << sig_over_tot_funct->Eval(0.000001) << "\n\nS/T for 1 sigma: " << sig_over_tot_funct->Eval(1) << "\nS/T for 2 sigma: " << sig_over_tot_funct->Eval(2) << std::endl;
    myText(.2,.95, kBlack, "#scale[1]{8-20 GeV/c Signal to Total, cuts:}");
    myText(.20,.90, kBlack, "#scale[1]{angle > %i mrad, 0.1 < lambda < 0.4, asymmetry < 0.7}");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "WholeSample_Signal_Over_Total.png"));
    TGraph* total_sigtot = new TGraph(sig_over_tot_funct);
    total_sigtot->Write("Total_Sig_To_Total");
     
    /**
    // Plot the energies of the two photons against each other
    graphcanvas->Clear();
    auto h2D = h_Pion->Projection(axis_photon2E, axis_photon1E);
    h2D->SetTitle("Photon momentum correlation; Leading Photon Energy (GeV); Trailing Photon Energy (GeV)");
    h2D->GetXaxis()->SetRangeUser(6.0, 15.0);
    h2D->GetYaxis()->SetRangeUser(3.0, 10.0);
    h2D->Draw("COLZ");
    //myText(.20,.95, kBlack, Form("#scale[1.5]{Photon energies, latest cut: %s}", headers[NumOfCuts].c_str()));
    graphcanvas->SaveAs(str_concat_converter(directory_name, "PhotonEs.png"));
    h2D->Write("Photon_Energies"); 
     */

    // Plot the data for the momentum, save into the Root file
    graphcanvas->Clear();
    TH1D* hPt = h_Pion->Projection(axis_pionPt);
    hPt->Draw();
    hPt->Write("Momentum-entries_chart");
    myText(.20,.95, kBlack, "#scale[1]{Momentum vs. Entries, cuts:}");
    myText(.20,.90, kBlack, "#scale[1]{angle > %i mrad, 0.1 < lambda < 0.4, asymmetry < 0.7}");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "momentum_pion_plot.png"));
    
    // A collection of variables that is needed for the next steps
    // num_of_intervals is already declaresd and initialized
    TMultiGraph* peaks_over_totals = new TMultiGraph();
    Color_t graph_colors[num_of_intervals] = {kRed, kBlue, kGreen, kYellow, kCyan, kBlack};
    
    // Intervals is already declared and initialized, along with min and max
    double chisquares[num_of_intervals];
    double means[num_of_intervals];
    double mean_errors[num_of_intervals];
    double sigmas[num_of_intervals];
    double sigma_errors[num_of_intervals];
    double gaussian_integrals[num_of_intervals];
    double integral_errors[num_of_intervals];
    
    double center[num_of_intervals];
    double widths[num_of_intervals];
    graphcanvas->Clear();
    
    //func->SetParameters(60,  0.14, 0.011, 1.3, 4.0,  -100000, 30000, -60000, 100, 10000);
    //func->SetParLimits(2, 0.01, 0.016); // width
    // start cutting the data up; plot the mass data for momenta of 6-8, 8-10, 10-12, 12-14, 14-16, 16-18, 18-20
    func->SetParameter(0, 30);
    for(int i = 0; i < num_of_intervals; i++) {
         
        pad[0] = new TPad("pad0","",0,0.36,1,1);
        pad[1] = new TPad("pad1","",0,0.05,1,0.45);
        min = intervals[i][0];
        max = intervals[i][1]; // Interval bounds
        
        // Asymmetry cut
        if (min >= 12)
            SetCut(h_Pion, axis_asymmetry, 0.0, 0.9);
        else
            SetCut(h_Pion, axis_asymmetry, 0.0, 0.6);
        
        // Cut the momentum to within the interval
        SetCut(h_Pion, axis_pionPt, min, max);
        // Plot the data and load it into the root file, after rebinning it properly and applying adaptive cuts
        hMass = h_Pion->Projection(axis_pionMass);
        TH1D* hMass = h_Pion->Projection(axis_pionMass);
        //hMass->Rebin(2);
        TH1D* residual = (TH1D*)hMass->Clone("residual");
        double MASSWIDTH = hMass->GetBinWidth(1);
        // Rebin if the interval is 16-20 GeV
        if (i == 5) {
            func->SetParLimits(1, 0.13, 0.19); //mean
            func->SetParLimits(2, 0.008, 0.018); // width
            
            hMass->Rebin(2);
            MASSWIDTH = hMass->GetBinWidth(1);
            residual = (TH1D*)hMass->Clone("residual");
        }
        //hMass->SetAxisRange(0.0, 1400.0, "Y");
        graphcanvas->cd();
        pad[0]->Draw();
        pad[0]->cd();
        hMass->GetYaxis()->SetTitle("Number of Entries");
        hMass->GetYaxis()->SetTitleOffset(0.8);
        hMass->Draw();
        hMass->Write(Form("unfitted_mass_pion-%2.2fGeV-%2.2fGeV", min, max));
        myText(.20,.95, kBlack, Form("#scale[1]{Mass vs Entries, Pt %2.2f-%2.2f, cuts:}", min, max));
        myText(.20,.90, kBlack, "#scale[1]{angle > %i mrad, 0.1 < lambda < 0.4, asymmetry < 0.7}");
        
        // Graph the fit and (separately) the Gaussian component of it, write the measured mean mass and the mass width (standard deviation) on the graph
        hMass->Fit(func);
        //hMass->Fit(func);
        if (i == 3 || i == 4)
            hMass->Fit(func);
        func->SetLineColor(kRed);
        func->Draw("same");
        int j = 0;
        for(; j < num_of_peak_params; j++) {
            peak->SetParameter(j, func->GetParameter(j));
            peak->SetParError(j, func->GetParError(j));
        }
        for(; j < num_of_params; j++) {
            background->SetParameter(j - num_of_peak_params, func->GetParameter(j));
            peak->SetParError(j - num_of_peak_params, func->GetParError(j));
        }
        peak->SetLineColor(kBlue);
        peak->Draw("same");
        background->SetLineColor(kGreen);
        background->Draw("same");
        //hMass->Draw("same");
        myText(.6, .85, kBlack, Form("#scale[0.75]{Mean Mass: %3.3f +/- %3.4f GeV}", func->GetParameter(1), func->GetParError(1)));
        myText(.6, .80, kBlack, Form("#scale[0.75]{Mass Width: %2.3f +/- %2.4f GeV}", func->GetParameter(2), func->GetParError(2)));
        
        /**
        // Print out all parameters
        std::cout << "\n\nNumber of pions:" << (func->GetParameter(0))/MASSWIDTH << std::endl << "Error:" << (func->GetParError(0))/MASSWIDTH << std::endl;
        std::cout << "\nMean:" << func->GetParameter(1) << std::endl << "Error:" << func->GetParError(1) << std::endl;
        std::cout << "\nSigma:" << func->GetParameter(2) << std::endl << "Error:" << func->GetParError(2) << std::endl;
        //std::cout << "\nAlpha:" << func->GetParameter(3) << std::endl << "Error:" << func->GetParError(3) << std::endl;
        //std::cout << "\nN:" << func->GetParameter(4) << std::endl << "Error:" << func->GetParError(4) << std::endl;
        //std::cout << "\nQuadric:" << func->GetParameter(5) << std::endl << "Error:" << func->GetParError(5) << std::endl;
        //std::cout << "\nCubic:" << func->GetParameter(6) << std::endl << "Error:" << func->GetParError(6) << std::endl;
        std::cout << "\nQuadratic:" << func->GetParameter(3) << std::endl << "Error:" << func->GetParError(5) << std::endl;
        std::cout << "\nLinear:" << func->GetParameter(4) << std::endl << "Error:" << func->GetParError(6) << std::endl;
        std::cout << "\nConstant:" << func->GetParameter(5) << std::endl << "Error:" << func->GetParError(7) << std::endl;
        */
        
        // Reinitialize residual_dist
        delete residual_dist;
        residual_dist = new TH1D("residual_distribution", "Residual_Distribution", 10, -5, 5);
        residual_dist->Sumw2();
        residual_dist->SetTitle("Residual Distribution; Residual value; Num of residuals");
        chisquare = 0;
        // Now add the residuals, creating values for residual_dist and counting up the chi-square in the process
        for (int i = 0; i < hMass->GetSize(); i++) {
            
            double residualvalue = (hMass->GetBinContent(i) - func->Eval(hMass->GetBinCenter(i)))/hMass->GetBinError(i);
            if ((hMass->GetBinError(i)) != 0) {
                residual->SetBinContent(i, residualvalue);
            }
            else
                residual->SetBinContent(i, 0);
            residual->SetBinError(i, 0); //Residuals don't have errors
            
            residual_dist->Fill(residualvalue);
            /**
            int dist_bin_num = residual_dist->FindBin(residual->GetBinContent(i));
            residual_dist->SetBinContent(dist_bin_num, residual_dist->GetBinContent(dist_bin_num) + 1);
            */
            chisquare += residualvalue*residualvalue;
            //chisquare += ((residual->GetBinContent(i))*(residual->GetBinContent(i)));
        }
        chisquares[i] = func->GetChisquare()/(hMass->GetSize() - (num_of_params) - 1); //Reduced Chi Square
        std::cout << Form("Reduced Chi Square: %2.2f", chisquares[i]) << std::endl;
        myText(.6, .75, kBlack, Form("#scale[0.75]{Reduced Chi-square: %2.1f/%i}", func->GetChisquare(), (hMass->GetSize() - num_of_params - 1)));
        myText(.6, .70, kBlack, Form("#scale[0.75]{P-val: %2.2f}", TMath::Prob(func->GetChisquare(), (hMass->GetSize() - num_of_params - 1))));
        graphcanvas->cd();
        residual->SetAxisRange(-4., 4., "Y");
        residual->GetXaxis()->SetTitleSize(.08);
        residual->GetYaxis()->SetTitleSize(.1);
        residual->GetYaxis()->SetTitleOffset(0.3);
        residual->GetXaxis()->SetTitleOffset(0.8);
        residual->SetTitle("; Pion Mass (GeV); Residuals");
        pad[1]->Draw("p");
        pad[1]->cd();
        residual->Draw("p");
        graphcanvas->SaveAs(str_concat_converter(directory_name, Form("MyFit_Ptmin_%2.2f_Ptmax_%2.2f.png", min, max)));
        
        // Add the mean mass parameter and its error to means and mean_errors, respectively; the standard deviation and
        // its error to sigmas and sigma_errors, the integral of the Gaussian peak and its error to gaussian_integrals
        // and integral_errors, and the center point of the interval into center (and its error into widths)
        means[i] = func->GetParameter(1);
        mean_errors[i] = func->GetParError(1);
        sigmas[i] = func->GetParameter(2) * 1000;
        sigma_errors[i] = func->GetParError(2) * 1000;
        
        gaussian_integrals[i] = (func->GetParameter(0))/MASSWIDTH;
        integral_errors[i] = (func->GetParError(0))/MASSWIDTH;
        
        center[i] = (min+max)/2.0;
        widths[i] = max - center[i];
        
        // Print out the number of pions and the signal to noise ratio
        std::cout << "\nNumber of pions: " << func->GetParameter(0)/MASSWIDTH << "\nError: " << func->GetParError(0)/MASSWIDTH << std::endl;
        std::cout << "\n\nS/T for 0 sigma: " << sig_over_tot_funct->Eval(0.000001) << "\n\nS/T for 1 sigma: " << sig_over_tot_funct->Eval(1) << "\nS/T for 2 sigma: " << sig_over_tot_funct->Eval(2) << std::endl;
        
        /**
        // print out the mean, standard deviation, and their corresponding errors; and save the graph in a new file
        std::cout << "Mean: " << means[i] << std::endl << "Standard Deviation: " << sigmas[i] << std::endl << "Standard Deviation Error: " << sigma_errors[i] << std::endl;
        graphcanvas->SaveAs(Form(str_concat_converter(directory_name, "MyFit_Ptmin_%2.2f_Ptmax_%2.2f.png"), min, max));
         */
        
        // Fit and plot the residuals distribution
        graphcanvas->Clear();
        residual_dist_fit->SetLineColor(2);
        residual_dist->Fit(residual_dist_fit);
        residual_dist_fit->SetParLimits(0, 0, 10000);
        residual_dist_fit->SetParLimits(1, -10, 10);
        residual_dist->GetYaxis()->SetTitleOffset(0.8);
        residual_dist->GetXaxis()->SetTitleOffset(0.9);
        residual_dist->Draw();
        residual_dist_fit->Draw("same");
        myText(.20, .92, kBlack, Form("#scale[1]{Residual Distribution, Pt %2.2f-%2.2f GeV/c}", min, max));
        myText(.2, .85, kBlack, Form("#scale[0.5]{Mean: %2.1f+/-%1.1f}", residual_dist_fit->GetParameter(1), residual_dist_fit->GetParError(1)));
        myText(.2, .8, kBlack, Form("#scale[0.5]{Sigma: %2.1f+/-%1.1f}", residual_dist_fit->GetParameter(2), residual_dist_fit->GetParError(2)));
        graphcanvas->SaveAs(str_concat_converter(directory_name, Form("MyFit_residual_dist_Ptmin_%2.2f_Ptmax_%2.2f.png", min, max)));

        
        //Now use the signal_over_total function to get a graph of sigma vs. signal/total
        graphcanvas->Clear();
        sig_over_tot_funct = new TF1("Signal over Total", signal_over_total, 0, 4, num_of_params);
        sig_over_tot_funct->SetParameters(func->GetParameters());
        graphcanvas->Clear();
        
        sig_over_tot_funct->SetTitle(Form("Signal over Total vs Distance From Mean: %2.2f to %2.2f GeV; Num of Standard Deviations From Mean; Signal to Total Ratio", min, max));
        TGraph* g_sig_over_tot = new TGraph(sig_over_tot_funct);
        g_sig_over_tot->SetLineColor(graph_colors[i]);
        g_sig_over_tot->SetName(Form("sig_to_tot_ptmin_%2.2fGeV_ptmax_%2.2fGeV", min, max));
        peaks_over_totals->Add(g_sig_over_tot);
        
        
        /**
        // Write data on number of pions and signal to total to the table file
        // Edit the common header to place the sub-interval's momentum bounds in the frame title
        string common_header = Form("\\frame\n{\n\\frametitle{Analysis of Cuts: Pt %2.0f-%2.0f GeV}\nNote: asymmetry without cuts is always less than 1\n\\begin{table}\n\\caption{How asymmetry cuts affect data quality}\n\\centering\n\\begin{tabular}{c c c c c}\n\\hline\\hline\nAsymmetry Cut & \\# Pions1 & S/T at 1 $\\sigma$ & S/T at 2 $\\sigma$ & p-val \\\\ [0.5ex]\n\\hline\n", min, max);
        string common_footer = "[1ex]\n\\hline\n\\end{tabular}\n\\label{table:nonlin}\n\\end{table}\n}\n";
        // First take in what is already in the file, put the string intended for input in its proper order, then write to the file
        ifstream table_file_input;
        table_file_input.open(Form("data/table_file_Ptmin_%2.2f_Ptmax_%2.2f.tex", min, max));
        // Print an error message if file cannot be opened
        if (!table_file_input)
            cout << "WARNING: TABLE FILE NOT FOUND" << std::endl;
        string before_in = "";
        string after_in = "";
        string current_line;
        string line_component;
        double order_num;
        while (!table_file_input.eof()) {
            std::getline(table_file_input, current_line);
            std::stringstream linestream(current_line);
            //cout << current_line << " is current line" << std::endl;
            std::getline(linestream, line_component, '%');
            //cout << line_component << " is line component" << std::endl;
            std::getline(linestream, line_component, '%');
            //cout << line_component << " is line component" << std::endl;
            // Continue to the next line if the line component's first character is not an integer
            // If the table file is not found, then break out of the loop if the line component is improper
            if (!table_file_input) {
                if (!isdigit(line_component[0]))
                    break;
            }
            else {
                if (!isdigit(line_component[0]))
                    continue;
            }

            order_num = std::stod(line_component);
            //cout << "Line originally in file: order number " << order_num << std::endl << std::endl;
            if (order_num > asymmetry_cut_num) {
                before_in = before_in + current_line + "\n";
                //cout << "Line goes before input\n";
            }
            else if (order_num < asymmetry_cut_num) {
                after_in = after_in + current_line + "\n";
                //cout << "Line goes after input\n";
            }
            //else
              //  cout << "Line overwritten\n";
            
        }
        table_file_input.close();
        ofstream table_file_output;
        table_file_output.open(Form("data/table_file_Ptmin_%2.2f_Ptmax_%2.2f.tex", min, max));
        if (!table_file_output)
            cout << "WARNING: TABLE FILE NOT FOUND" << std::endl;
        table_file_output << common_header << before_in << "asym $<$ " << asymmetry_cut_num << " & " << Form("%4.0f",(func->GetParameter(0))/MASSWIDTH) << " +/- " << Form("%4.0f",func->GetParError(0)/MASSWIDTH) << " & " << Form("%2.2f", sig_over_tot_funct->Eval(1)) <<  " & " << Form("%2.2f", sig_over_tot_funct->Eval(2)) << " & " << Form("%2.3f", TMath::Prob(func->GetChisquare(), (hMass->GetSize() - num_of_params - 1))) << " \\\\ %" << asymmetry_cut_num << "%" << std::endl << after_in << common_footer;
        table_file_output.close();
        */
        
        //Load onto the ROOT file
        hMass->GetListOfFunctions()->Add(peak);
        hMass->GetListOfFunctions()->Add(background);
        hMass->Write(Form("mass-pion-%2.2fGeV-%2.2fGeV", min, max));
        residual->Write(Form("residual-%2.2fGeV-%2.2fGeV", min, max));
        
        /**
        // Graph the photon energy ratio
        graphcanvas->Clear();
        auto h2D = h_Pion->Projection(axis_photon2E, axis_photon1E);
        h2D->SetTitle("Photon energy correlation; Leading Photon Energy (GeV); Trailing Photon Energy (GeV)");
        h2D->GetXaxis()->SetRangeUser(6.0, 15.0);
        h2D->GetYaxis()->SetRangeUser(3.0, 10.0);
        h2D->Draw("COLZ");
        //myText(.20,.95, kBlack, Form("#scale[1.5]{Photon energies: Pion momenta %2.2f to %2.2f, latest cut: %s}", min, max, headers[NumOfCuts].c_str()));
        graphcanvas->SaveAs(Form(str_concat_converter(directory_name, "PhotonEs_Ptmin_%2.2f_Ptmax_%2.2f.png"), min, max));
        h2D->Write(Form("Photon_Energies_Ptmin_%2.2f_Ptmax_%2.2f", min, max));
         */
    }// end of loop over pt
    
    // Graph the signal/total curves for each momentum increment
    graphcanvas->Clear();
    peaks_over_totals->SetTitle("Signal over Total vs Distance From Mean; #scale[1]{Num of Standard Deviations From Mean}; #scale[1]{Signal to Total Ratio}");
    //peaks_over_totals->GetYaxis()->SetRangeUser(0.4, 1.0);
    peaks_over_totals->SetMaximum(1.0);
    peaks_over_totals->SetMinimum(0.2);
    peaks_over_totals->Draw("Al");
    myBoxText(0.25, 0.45, 0.05, 10, graph_colors[0], "6-8 GeV");
    myBoxText(0.25, 0.40, 0.05, 10, graph_colors[1], "8-10 GeV");
    myBoxText(0.25, 0.35, 0.05, 10, graph_colors[2], "10-12 GeV");
    myBoxText(0.25, 0.30, 0.05, 10, graph_colors[3], "12-14 GeV");
    myBoxText(0.25, 0.25, 0.05, 10, graph_colors[4], "14-16 GeV");
    myBoxText(0.25, 0.20, 0.05, 10, graph_colors[5], "16-20 GeV");
    peaks_over_totals->Write("signal-over-total");
    myText(.30,.93, kBlack, "#scale[1]{Signal-to-total ratios}");
    //myText(.20,.92, kBlack, "#scale[1]{angle > %i mrad, 0.1 < lambda < 0.4, asymmetry < 0.7}");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "Overall_Signal_Over_Total.png"));
    graphcanvas->Clear();
    
    // Add a constant "expected mass" function to be graphed alongside the mean mass data
    TF1* mass_pdg = new TF1("mass_pdg", "[0]", 0, 20);
    mass_pdg->SetParameter(0, 0.13498);
    mass_pdg->SetLineWidth(2);
    mass_pdg->SetLineColor(kRed);
    
    // Graph mean masses with error bars
    TGraphErrors* g_mean = new TGraphErrors(num_of_intervals, center, means, widths, mean_errors);
    graphcanvas->Clear();
    g_mean->Print();
    g_mean->SetTitle("Mean Masses for Various Momenta; Momentum (GeV); Mass (GeV/c^2)");
    g_mean->GetYaxis()->SetRangeUser(0.125, 0.157);
    //g_mean->SetMarkerSize(2);
    //g_mean->SetMarkerStyle(20);
    g_mean->Draw("AP");
    mass_pdg->Draw("same");
    g_mean->Write("mean-masses");
    myText(.20,.97, kBlack, "#scale[1]{Mean masses over Pt, cuts:}");
    myText(.20,.92, kBlack, "#scale[1]{angle > %i mrad, 0.1 < lambda < 0.4, asymmetry < 0.7}");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "meanMass_v_pT.png"));
    
    // Graph mass standard deviations with error bars
    graphcanvas->Clear();
    TGraphErrors* g_sigma = new TGraphErrors(num_of_intervals, center, sigmas, widths, sigma_errors);
    g_sigma->Print();
    g_sigma->SetTitle("Mass Peak Widths for Various Momenta; Momentum (GeV); Mass Width (MeV/c^2)");
    g_sigma->GetYaxis()->SetTitleOffset(.7);
    g_sigma->GetXaxis()->SetTitleOffset(.9);
    g_sigma->GetYaxis()->SetRangeUser(3.0, 15.7);
    g_sigma->GetXaxis()->SetRangeUser(6.0, 16.0);
    //g_sigma->SetMarkerSize(2);
    //g_sigma->SetMarkerStyle(20);
    g_sigma->Write("standard-dev-masses");
    g_sigma->Draw("AP");
    myText(.20,.97, kBlack, "#scale[1]{Mass widths over Pt, cuts:}");
    myText(.20,.92, kBlack, "#scale[1]{angle > %i mrad, 0.1 < lambda < 0.4, asymmetry < 0.7}");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "massWidths_v_pT.png"));
    
    // Graph reduced chi squares over the momentum interval
    graphcanvas->Clear();
    TGraph* g_chisquare = new TGraphErrors(num_of_intervals, center, chisquares);
    g_chisquare->Print();
    g_chisquare->SetTitle("ReducedChi-Squares for Various Momenta; Momentum (GeV); Reduced Chi-Squares");
    g_chisquare->SetMarkerSize(2);
    g_chisquare->SetMarkerStyle(20);
    g_chisquare->GetYaxis()->SetRangeUser(2.0, 20.0);
    g_chisquare->Write("chi-square");
    g_chisquare->Draw("AP");
    myText(.20,.97, kBlack, "#scale[1]{Reduced chi-squares over Pt, cuts:}");
    myText(.20,.92, kBlack, "#scale[1]{angle > %i mrad, 0.1 < lambda < 0.4, asymmetry < 0.7}");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "reduced_chisquare_v_pT.png"));
    
    // Graph the distribution integrals over momentum
    graphcanvas->Clear();
    TGraphErrors* g_integral = new TGraphErrors(num_of_intervals, center, gaussian_integrals, widths, integral_errors);
    g_integral->Print();
    g_integral->SetTitle("Peak integrals for Various Momenta; Momentum (GeV); Number of Pions");
    g_integral->GetXaxis()->SetRangeUser(6.0, 16.0);
    //g_integral->GetYaxis()->SetRangeUser(0.0, 6000.0);
    g_integral->Draw("AP");
    g_integral->Write("pion-integrals");
    myText(.20,.97, kBlack, "#scale[1]{Integrals over Pt, cuts:}");
    myText(.20,.92, kBlack, "#scale[1]{angle > %i mrad, 0.1 < lambda < 0.4, asymmetry < 0.7}");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "peakIntegrals_v_pT.png"));
    
    graphcanvas->Close();
    logcanvas->Close();
    return;
}
