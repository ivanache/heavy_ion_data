/**
 my_code is a data-processing macro meant for processing THnSparses root files, specifically the h_Pion one in THnSparses_070517.root
 the macro must be called with two bools, the first to tell if the lambda must be cut, the second to tell if asymmetry must be cut
 */

#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TGraphErrors.h>
#include <TCanvas.h>
#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"
#include <iostream>
#include <fstream>

// The output file
TFile* fOut;

// The general filepath string
string directory_name;

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
const int axis_pionDisToCharged1 = 25;
const int axis_pionDisToCharged2 = 26;
const int axis_pionDisToBorder1 = 27;
const int axis_pionDisToBorder2 = 28;
const int axis_pionDisToBadChannel1 = 29;
const int axis_pionDisToBadChannel2 = 30;

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

// The fitting function
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
// Crystal ball function
double crystal_ball_model(Double_t *x, Double_t *par) {
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
    can->SaveAs(filename);
    data_subset->Write(rootname.c_str());
    can->Clear();
}

/**
 Main function
 */
// Precondition: NumOfCuts is 0, 1, 2, 3, 4, 5, 6, or 7
// Precondition: option is either "1" (for the cut on the distance to bad channel to be > 1) or "2" (for the cut on the distance to bad channel to be > 2)
void my_code(int NumOfCuts, string option = "1") {
    // If NumOfCuts or option is not one of the accepted values, print an error message
    if (NumOfCuts < 0 || NumOfCuts > 7){
        std::cout << "ERROR: NumOfCuts must be:\n0, for no cuts,\n1, for distance to charged particles cut,\n2, for all of the above cuts plus asymmetry cut,\n3, for all of the above cuts plus angle cut,\n4, for all of the above cuts plus Ncells cut,\n5, for all of the above cuts plus distance to border cut, or\n6, for all of the above cuts plus lambda cut" << std::endl;
        return;
    }
    if (option != "1" && option != "2") {
        std::cout << "ERROR: option must be:\n\"1\", for distance to bad channel > 1, or\n\"2\", for distance to bad channel > 2" << std::endl;
        return;
    }
    
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    directory_name = Form("data/%icuts/", NumOfCuts);
    
    //Open the files
    TFile* fIn = new TFile("THnSparses_070517.root","READ"); //get file
    string rootfilename = Form("data/Pion%iCutsSparsesOutput.root", NumOfCuts);
    fOut = new TFile(rootfilename.c_str(), "RECREATE"); // Create an output file
    fIn->Print(); //print file content
    
    // Get the data
    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion); //get array
    
    //The TCanvases and TPads (for sub-canvassing)
    TCanvas* graphcanvas = new TCanvas();
    TPad *pad[2] = {new TPad("pad0","",0,0.38,1,1), new TPad("pad1","",0,0,1,0.46)};
    TCanvas* logcanvas = new TCanvas();
    logcanvas->SetLogy();
    
    //For the mass plot, restrict to mass between 0.08 and 0.25 and the momentum to between 5 GeV and 18 GeV
    //plot the data for the asymmetry, both lambdas, the angle, and number of cells
    SetCut(h_Pion, axis_pionPt, 8.0, 20.0);
    SetCut(h_Pion, axis_pionMass, 0.08, 0.25);
    
    
    //Apply cuts, as described in the cout output
    if(NumOfCuts == 0) {
        cout << "No cuts done\n\n";
        //fit_y_max = 1600.0;
    }
    else if (NumOfCuts == 1) {
        SetCut(h_Pion, axis_pionDisToCharged1, 0.02, 0.14);
        SetCut(h_Pion, axis_pionDisToCharged2, 0.02, 0.14);
        std::cout << "Cuts: distance to charged particles\n\n";
        //fit_y_max = 1000.0;
    }
    else if (NumOfCuts == 2) {
        SetCut(h_Pion, axis_pionDisToCharged1, 0.02, 0.14);
        SetCut(h_Pion, axis_pionDisToCharged2, 0.02, 0.14);
        SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
        std::cout << "Cuts: distance to charged particles and asymmetry\n\n";
        //fit_y_max = 1000.0;
    }
    else if (NumOfCuts == 3) {
        SetCut(h_Pion, axis_pionDisToCharged1, 0.02, 0.14);
        SetCut(h_Pion, axis_pionDisToCharged2, 0.02, 0.14);
        SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
        SetCut(h_Pion, axis_pionAngle, 0.015, 0.5);
        std::cout << "Cuts: distance to charged particles, asymmetry, and angle\n\n";
        //fit_y_max = 600.0;
    }
    else if (NumOfCuts == 4) {
        SetCut(h_Pion, axis_pionDisToCharged1, 0.02, 0.14);
        SetCut(h_Pion, axis_pionDisToCharged2, 0.02, 0.14);
        SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
        SetCut(h_Pion, axis_pionAngle, 0.015, 0.5);
        SetCut(h_Pion, axis_pionNcells1, 1.9, 30.0);
        SetCut(h_Pion, axis_pionNcells2, 1.9, 30.0);
        std::cout << "Cuts: distance to charged particles, asymmetry, angle, and Ncells\n\n";
        //fit_y_max = 600.0;
    }
    else if (NumOfCuts == 5) {
        SetCut(h_Pion, axis_pionDisToCharged1, 0.02, 0.14);
        SetCut(h_Pion, axis_pionDisToCharged2, 0.02, 0.14);
        SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
        SetCut(h_Pion, axis_pionAngle, 0.015, 0.5);
        SetCut(h_Pion, axis_pionNcells1, 1.9, 30.0);
        SetCut(h_Pion, axis_pionNcells2, 1.9, 30.0);
        SetCut(h_Pion, axis_pionDisToBorder1, 2.0, 6.0);
        SetCut(h_Pion, axis_pionDisToBorder2, 2.0, 6.0);
        std::cout << "Cuts: distance to charged particles, asymmetry, angle, Ncells, and DisToBorder\n\n";
    }
    else if (NumOfCuts == 6) {
        SetCut(h_Pion, axis_pionDisToCharged1, 0.02, 0.14);
        SetCut(h_Pion, axis_pionDisToCharged2, 0.02, 0.14);
        SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
        SetCut(h_Pion, axis_pionAngle, 0.015, 0.5);
        SetCut(h_Pion, axis_pionNcells1, 1.9, 30.0);
        SetCut(h_Pion, axis_pionNcells2, 1.9, 30.0);
        SetCut(h_Pion, axis_pionDisToBorder1, 2.0, 6.0);
        SetCut(h_Pion, axis_pionDisToBorder2, 2.0, 6.0);
        if (option == "1") {
            SetCut(h_Pion, axis_pionDisToBadChannel1, 1.9, 10.0);
            SetCut(h_Pion, axis_pionDisToBadChannel2, 1.9, 10.0);
        }
        if (option == "2") {
            SetCut(h_Pion, axis_pionDisToBadChannel1, 2.9, 10.0);
            SetCut(h_Pion, axis_pionDisToBadChannel2, 2.9, 10.0);
        }
        std::cout << "Cuts: distance to charged particles, asymmetry, angle, Ncells, DisToBorder, and distance to bad channel\n\n";
    }
    
    else if (NumOfCuts == 7) {
        SetCut(h_Pion, axis_pionDisToCharged1, 0.02, 0.14);
        SetCut(h_Pion, axis_pionDisToCharged2, 0.02, 0.14);
        SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
        SetCut(h_Pion, axis_pionAngle, 0.015, 0.5);
        SetCut(h_Pion, axis_pionNcells1, 1.9, 30.0);
        SetCut(h_Pion, axis_pionNcells2, 1.9, 30.0);
        SetCut(h_Pion, axis_pionDisToBorder1, 2.0, 6.0);
        SetCut(h_Pion, axis_pionDisToBorder2, 2.0, 6.0);
        if (option == "1") {
            SetCut(h_Pion, axis_pionDisToBadChannel1, 1.9, 10.0);
            SetCut(h_Pion, axis_pionDisToBadChannel2, 1.9, 10.0);
        }
        if (option == "2") {
            SetCut(h_Pion, axis_pionDisToBadChannel1, 2.9, 10.0);
            SetCut(h_Pion, axis_pionDisToBadChannel2, 2.9, 10.0);
        }
        SetCut(h_Pion, axis_pionLambda1, 0.1, 0.4);
        SetCut(h_Pion, axis_pionLambda2, 0.1, 0.4);
        std::cout << "Cuts: distance to charged particles, asymmetry, angle, Ncells, DisToBorder, distance to bad channel, and lambda02\n\n";
    }
    
    // Graph the data with respect to all of the parameters
    graph_raw_data(h_Pion, axis_asymmetry, graphcanvas, str_concat_converter(directory_name,"asymmetry_pion_plot.png"), "asymmetry");
    graph_raw_data(h_Pion, axis_pionLambda1, graphcanvas, str_concat_converter(directory_name, "lambda1_pion_plot.png"), "lambda1");
    graph_raw_data(h_Pion, axis_pionLambda2, graphcanvas, str_concat_converter(directory_name, "lambda2_pion_plot.png"), "lambda2");
    graph_raw_data(h_Pion, axis_pionAngle, graphcanvas, str_concat_converter(directory_name, "angle_pion_plot.png"), "angle");
    graph_raw_data(h_Pion, axis_pionNcells1, graphcanvas, str_concat_converter(directory_name, "Ncells1_pion_plot.png"), "Ncells1");
    graph_raw_data(h_Pion, axis_pionNcells2, graphcanvas, str_concat_converter(directory_name, "Ncells2_pion_plot.png"), "Ncells2");
    graph_raw_data(h_Pion, axis_pionMatchedTracks1, graphcanvas, str_concat_converter(directory_name, "MatchedTracks1_pion_plot.png"), "MatchedTracks1");
    graph_raw_data(h_Pion, axis_pionMatchedTracks2, graphcanvas, str_concat_converter(directory_name, "MatchedTracks2_pion_plot.png"), "MatchedTracks2");
    graph_raw_data(h_Pion, axis_pionDisToBorder1, graphcanvas, str_concat_converter(directory_name, "DisToBorder1_pion_plot.png"), "DisToBorder1");
    graph_raw_data(h_Pion, axis_pionDisToBorder2, graphcanvas, str_concat_converter(directory_name, "DisToBorder2_pion_plot.png"), "DisToBorder2");
    graph_raw_data(h_Pion, axis_pionDisToCharged1, logcanvas, str_concat_converter(directory_name, "DisToChargedParticle1_pion_plot.png"), "DisToChargedParticle1");
    graph_raw_data(h_Pion, axis_pionDisToCharged2, logcanvas, str_concat_converter(directory_name, "DisToChargedParticle2_pion_plot.png"), "DisToChargdParticle2");
    graph_raw_data(h_Pion, axis_pionDisToBadChannel1, graphcanvas, str_concat_converter(directory_name, "DisToBadChannel1_pion_plot.png"), "DisToBadChannel1");
    graph_raw_data(h_Pion, axis_pionDisToBadChannel2, graphcanvas, str_concat_converter(directory_name, "DisToBadChannel2_pion_plot.png"), "DisToBadChannel2");
    
    // plot mass data, write it into the root file, and set up the fit function
    TH1D* hMass = h_Pion->Projection(axis_pionMass);
    hMass->Rebin(2);
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
    int num_of_params = 8;
    int num_of_peak_params = 3;
    background = new TF1("background curve", "[0]*TMath::Power(x, 4.0) + [1]*TMath::Power(x, 3.0) + [2]*x*x + [3]*x + [4]", 0.08, 0.26);
    num_of_params = 10;
    num_of_peak_params = 5;
    func = new TF1("fit", crystal_ball_model,0.05,0.5,num_of_params);
    peak = new TF1("mass peak", crystal_ball_function_peak, 0.05, 0.5, num_of_peak_params);
    func->SetParNames("Integral", "Mean", "Sigma", "Alpha", "N", "Quadric coeff", "Cubic coeff", "Quadratic coeff", "Linear coeff", "Constant");
    func->SetParameters(60,  0.14, 0.013, 1, 4.0,  -100000, 30000, -60000, 100, 10000);
    func->SetParLimits(0, 0.1, 400.0);//integral
    func->SetParLimits(1, 0.13, 0.151); //mean
    func->SetParLimits(2, 0.01, 0.016); // width
    func->SetParLimits(3, 0.5, 1000.0); // alpha
    func->SetParLimits(4, 1.2, 1000.0); // n
    func->SetParLimits(5, -1000000.0, 0.0); // Quadric and quadratic factors
    func->SetParLimits(7, -1000000.0, 0.0);
    //func->SetParLimits(8, 0.0, 1000000.0); // Linear
    //func->SetParLimits(9, -1000000.0, 0.0); //Constant
    
    // Plot the fit for the mass and (separately) the Gaussian and background components of it, print out the number of pions and the other parameters
    hMass->Fit(func, "0");
    std::cout << "\n\nNumber of pions:" << (func->GetParameter(0))/MASSWIDTH << std::endl << "Error:" << (func->GetParError(0))/MASSWIDTH << std::endl;
    std::cout << "\nMean:" << func->GetParameter(1) << std::endl << "Error:" << func->GetParError(1) << std::endl;
    std::cout << "\nSigma:" << func->GetParameter(2) << std::endl << "Error:" << func->GetParError(2) << std::endl;
    std::cout << "\nAlpha:" << func->GetParameter(3) << std::endl << "Error:" << func->GetParError(3) << std::endl;
    std::cout << "\nN:" << func->GetParameter(4) << std::endl << "Error:" << func->GetParError(4) << std::endl;
    std::cout << "\nQuadric:" << func->GetParameter(5) << std::endl << "Error:" << func->GetParError(5) << std::endl;
    std::cout << "\nCubic:" << func->GetParameter(6) << std::endl << "Error:" << func->GetParError(6) << std::endl;
    std::cout << "\nQuadratic:" << func->GetParameter(7) << std::endl << "Error:" << func->GetParError(7) << std::endl;
    std::cout << "\nLinear:" << func->GetParameter(8) << std::endl << "Error:" << func->GetParError(8) << std::endl;
    std::cout << "\nConstant:" << func->GetParameter(9) << std::endl << "Error:" << func->GetParError(9) << std::endl;
    func->SetLineColor(kRed);
    func->Draw("same");
    std::cout << "Reduced Chi Square " << (func->GetChisquare())/10 << std::endl; //Reduced Chi Square of the mass vs entries curve (function has 18 degrees of freedom, 7 parameters)
    int i = 0;
    for(; i < num_of_peak_params; i++)
        peak->SetParameter(i, func->GetParameter(i));
    for(; i < num_of_params; i++)
        background->SetParameter(i - num_of_peak_params, func->GetParameter(i));
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
    
    // Plot the residual; save as a PDF, print out the individual residuals
    std::cout << "Differences:\n";
    for (int i = 0; i < hMass->GetSize(); i++) {
        if (hMass->GetBinError(i))
            residual->SetBinContent(i, ((hMass->GetBinContent(i) - func->Eval(hMass->GetBinCenter(i)))/hMass->GetBinError(i)));
        residual->SetBinError(i, 0); //Residuals don't have errors
        std::cout << hMass->GetBinCenter(i) << " GeV momentum: " << (hMass->GetBinContent(i) - func->Eval(hMass->GetBinCenter(i))) << std::endl;
    }
    std::cout << std::endl;
    graphcanvas->cd();
    pad[1]->Draw("p");
    pad[1]->cd();
    residual->SetTitle("; Pion Mass (GeV); Residuals");
    residual->SetAxisRange(-4., 4., "Y");
    residual->GetXaxis()->SetTitleSize(.06);
    residual->GetYaxis()->SetTitleSize(.07);
    residual->GetYaxis()->SetTitleOffset(0.3);
    residual->GetXaxis()->SetTitleOffset(0.5);
    residual->SetMarkerStyle(20);
    residual->SetMarkerSize(1);
    residual->Draw("p");
    residual->Write("residual"); // Load into the ROOT file
    graphcanvas->SaveAs(str_concat_converter(directory_name, "mass_pion_plot.png"));
    
    
    // Plot signal to noise ratio for the entire sample over various distances from the mean, print out values for 1 sigma and 2 sigmas
    graphcanvas->Clear();
    TF1* sig_over_tot_funct = new TF1("Signal over Total", signal_over_total, 0, 4, num_of_params);
    sig_over_tot_funct->SetParameters(func->GetParameters());
    sig_over_tot_funct->SetTitle("Signal over Total vs Distance From Mean; Num of Standard Deviations From Mean; Signal to Total Ratio");
    sig_over_tot_funct->Draw();
    std::cout << "\n\nS/T for 0 sigma: " << sig_over_tot_funct->Eval(0.000001) << "\n\nS/T for 1 sigma: " << sig_over_tot_funct->Eval(1) << "\nS/T for 2 sigma: " << sig_over_tot_funct->Eval(2) << std::endl;
    graphcanvas->SaveAs(str_concat_converter(directory_name, "WholeSample_Signal_Over_Total.png"));
    TGraph* total_sigtot = new TGraph(sig_over_tot_funct);
    total_sigtot->Write("Total_Sig_To_Total");
    
    // Write to the table file
        // Common header: goes on top of the table content commands. Creates the frame and table, sets the frame title, table title, and frame and table formatting
    string common_header = "\\frame\n{\n\\frametitle{Analysis of Cuts: Pt 7-20 GeV}\n\\begin{table}\n\\caption{How cuts affect data quality}\n\\centering\n\\begin{tabular}{c c c c}\n\\hline\\hline\nCuts & \\# Pions1 & S/T at 1 $\\sigma$ & S/T at 2 $\\sigma$ \\\\ [0.5ex]\n\\hline\n";
        // Common footer: goes on the bottom of the table content commands. Closes up the frame and the table
    string common_footer = "[1ex]\n\\hline\n\\end{tabular}\n\\label{table:nonlin}\n\\end{table}\n}\n";
    const int cutnum_option_quantity = 8;
    string headers[cutnum_option_quantity] = {"None ", "+dR $> 20$ mrad ", "+asymmetry $< 0.7$ ", "+angle $> 15$ mrad ", "+Ncells $> 1$", "+DisToBorder $>$ 2 ", "+DisToBadChannel $>$ 1", "+$0.1 < \\lambda <$ 0.4 "};
    if (option == "2")
        headers[6] = "+DisToBadChannel $>$ 2";
    // First take in what is already in the file, put the string intended for input in its proper order, then write to the file
    ifstream table_file_input;
    table_file_input.open("data/table_file.tex");
    // Print an error message if file cannot be opened
    if (!table_file_input)
        cout << "WARNING: TABLE FILE NOT FOUND" << std::endl;
    string before_in = "";
    string after_in = "";
    string current_line;
    string line_component;
    int order_num;
    while (!table_file_input.eof()) {
        std::getline(table_file_input, current_line);
        std::stringstream linestream(current_line);
        //cout << current_line << " is current line" << std::endl;
        std::getline(linestream, line_component, '%');
        //cout << line_component << " is line component" << std::endl;
        std::getline(linestream, line_component, '%');
        //cout << line_component << " is line component" << std::endl;
        // Continue to the next line if the line component's first character is not an integer
        if (!isdigit(line_component[0]))
            continue;
        order_num = std::stoi(line_component);
        //cout << "Line originally in file: order number " << order_num << std::endl << std::endl;
        if (order_num < NumOfCuts) {
            before_in = before_in + current_line + "\n";
            //cout << "Line goes before input\n";
        }
        else if (order_num > NumOfCuts) {
            after_in = after_in + current_line + "\n";
            //cout << "Line goes after input\n";
        }
        //else
            //cout << "Line overwritten\n";
        
    }
    table_file_input.close();
    ofstream table_file_output;
    table_file_output.open("data/table_file.tex");
    if (!table_file_output)
        cout << "WARNING: TABLE FILE NOT FOUND" << std::endl;
    table_file_output << common_header << before_in << headers[NumOfCuts] << "& " << Form("%4.0f",(func->GetParameter(0))/MASSWIDTH) << " +/- " << Form("%4.0f",(func->GetParError(0))/MASSWIDTH) << " & " << Form("%2.2f", sig_over_tot_funct->Eval(1)) <<  " & " << Form("%2.2f", sig_over_tot_funct->Eval(2)) << " \\\\ %" << NumOfCuts << "%" << std::endl << after_in << common_footer;
    table_file_output.close();
    
    // Plot the energies of the two photons against each other
    graphcanvas->Clear();
    auto h2D = h_Pion->Projection(axis_photon2E, axis_photon1E);
    h2D->SetTitle("Photon momentum correlation; Leading Photon Energy (GeV); Trailing Photon Energy (GeV)");
    h2D->GetXaxis()->SetRangeUser(6.0, 15.0);
    h2D->GetYaxis()->SetRangeUser(3.0, 10.0);
    h2D->Draw("COLZ");
    myText(.35,.9, kBlack, "Photon energies");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "PhotonEs.png"));
    h2D->Write("Photon_Energies");

    // Plot the data for the momentum, save into the Root file
    graphcanvas->Clear();
    TH1D* hPt = h_Pion->Projection(axis_pionPt);
    hPt->Draw();
    hPt->Write("Momentum-entries_chart");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "momentum_pion_plot.png"));
    
    // A collection of variables that is needed for the next steps
    const int num_of_intervals = 6;
    TMultiGraph* peaks_over_totals = new TMultiGraph();
    Color_t graph_colors[num_of_intervals] = {kRed, kBlue, kGreen, kYellow, kCyan, kBlack};
    
    double min;
    double max;
    
    double intervals[num_of_intervals][2] = {{8.0, 10.0}, {10.0, 11.0}, {11.0, 12.0}, {12.0, 13.0}, {13.0, 15.0}, {15.0, 20.0}};
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
    
    // start cutting the data up; plot the mass data for momenta of 8-10, 10-11, 11-12, 12-13, 13-15, 15-20
    for(int i = 0; i < num_of_intervals; i++) {
        //Adaptive Cuts go here
         
        pad[0] = new TPad("pad0","",0,0.38,1,1);
        pad[1] = new TPad("pad1","",0,0,1,0.46);
        min = intervals[i][0];
        max = intervals[i][1]; // Interval bounds
        
        // Cut the momentum to withing the interval, and graph the parameter vs number of entries graphs for parameters whose cuts should potentially be momentum-dependent
        SetCut(h_Pion, axis_pionPt, min, max);
            // lambda
        graph_raw_data(h_Pion, axis_pionLambda1, graphcanvas, str_concat_converter(directory_name, Form("lambda1_pion_plot_Ptmin_%2.2f_Ptmax_%2.2f.png", min, max)), Form("lambda1_Ptmin_%2.2f_Ptmax_%2.2f", min, max));
        graph_raw_data(h_Pion, axis_pionLambda1, graphcanvas, str_concat_converter(directory_name, Form("lambda2_pion_plot_Ptmin_%2.2f_Ptmax_%2.2f.png", min, max)), Form("lambda2_Ptmin_%2.2f_Ptmax_%2.2f", min, max));
            // distance to bad channel
        graph_raw_data(h_Pion, axis_pionDisToBadChannel1, graphcanvas, str_concat_converter(directory_name, Form("DisToBadChannel1_pion_plot.png_Ptmin_%2.2f_Ptmax_%2.2f.png", min, max)), Form("DisToBadChannel1_Ptmin_%2.2f_Ptmax_%2.2f", min, max));
        graph_raw_data(h_Pion, axis_pionDisToBadChannel2, graphcanvas, str_concat_converter(directory_name, Form("DisToBadChannel2_pion_plot.png_Ptmin_%2.2f_Ptmax_%2.2f.png", min, max)), Form("DisToBadChannel2_Ptmin_%2.2f_Ptmax_%2.2f", min, max));

        
        // Plot the data and load it into the root file, after rebinning it properly and applying adaptive cuts
        hMass = h_Pion->Projection(axis_pionMass);
        if (i == num_of_intervals - 1) {
            hMass->Rebin(4);
            MASSWIDTH = hMass->GetBinWidth(1);
            residual = (TH1D*)hMass->Clone("residual");
            func->SetParLimits(1, 0.13, 0.165); //mean
            func->SetParLimits(5, -100000.0, 0.0); // Quadric
        }
        else if (i == 0 || i == 4) {
            hMass->Rebin(2);
            func->SetParLimits(5, -100000.0, 0.0); // Quadric
            func->SetParLimits(2, 0.01, 0.016); // width
            if (i == 0) {
                func->SetParLimits(1, 0.13, 0.145); //mean
                func->SetParLimits(2, 0.009, 0.016); // width
            }
        }
        else if (i == 3) {
            hMass->Rebin(2);
            func->SetParLimits(2, 0.009, 0.016); // width
            func->SetParLimits(5, -100000.0, 0.0); // Quadric
        }
        else {
            hMass->Rebin(2);
            func->SetParLimits(5, -1000000.0, 0.0); // Quadric
            func->SetParLimits(1, 0.13, 0.151); //mean
            func->SetParLimits(2, 0.01, 0.016); // width
        }

        //hMass->SetAxisRange(0.0, 1400.0, "Y");
        graphcanvas->cd();
        pad[0]->Draw();
        pad[0]->cd();
        hMass->Draw();
        hMass->Write(Form("unfitted_mass_pion-%2.2fGeV-%2.2fGeV", min, max));
        
        // Find a fit just as you did for the entire data set and the reduced chi square of the fit into its respective array
        // Graph the fit and (separately) the Gaussian component of it
        hMass->Fit(func);
        func->SetLineColor(kRed);
        chisquares[i] = (func->GetChisquare())/10; //Reduced Chi Square (function has 18 degrees of freedom, 7 parameters)
        std::cout << Form("Reduced Chi Square: %2.2f", chisquares[i]) << std::endl;
        func->Draw("same");
        int j = 0;
        for(; j < num_of_peak_params; j++)
            peak->SetParameter(j, func->GetParameter(j));
        for(; j < num_of_params; j++)
            background->SetParameter(j - num_of_peak_params, func->GetParameter(j));
        peak->SetLineColor(kBlue);
        peak->Draw("same");
        background->SetLineColor(kGreen);
        background->Draw("same");
        hMass->Draw("same");
        
        // Now add the residuals
        for (int i = 0; i < hMass->GetSize(); i++) {
            if ((hMass->GetBinError(i)) != 0) {
                residual->SetBinContent(i, (hMass->GetBinContent(i) - func->Eval(hMass->GetBinCenter(i)))/hMass->GetBinError(i));
            }
            else
                residual->SetBinContent(i, 0);
            residual->SetBinError(i, 0); //Residuals don't have errors
            //std::cout<< "Bin: " << i << std::endl;
            //std::cout<< (hMass->GetBinContent(i)) << std::endl;
            //std::cout<< (func->Eval(hMass->GetBinCenter(i))) << std::endl;
            //std::cout<< (hMass->GetBinError(i)) << std::endl << std::endl;
            
        }
        graphcanvas->cd();
        residual->SetAxisRange(-4., 4., "Y");
        //residual->SetTitle("; Pion Mass (GeV); Residuals");
        pad[1]->Draw("p");
        pad[1]->cd();
        residual->Draw("p");
        
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
        
        // print out the mean, standard deviation, and their corresponding errors; and save the graph in a new file
        std::cout << "Mean: " << means[i] << std::endl << "Standard Deviation: " << sigmas[i] << std::endl << "Standard Deviation Error: " << sigma_errors[i] << std::endl;
        graphcanvas->SaveAs(Form(str_concat_converter(directory_name, "MyFit_Ptmin_%2.2f_Ptmax_%2.2f.png"), min, max));
        
        //Now use the signal_over_total function to get a graph of sigma vs. signal/total
        sig_over_tot_funct = new TF1("Signal over Total", signal_over_total, 0, 4, num_of_params);
        sig_over_tot_funct->SetParameters(func->GetParameters());
        graphcanvas->Clear();
        
        sig_over_tot_funct->SetTitle(Form("Signal over Total vs Distance From Mean: %2.2f to %2.2f GeV; Num of Standard Deviations From Mean; Signal to Total Ratio", min, max));
        TGraph* g_sig_over_tot = new TGraph(sig_over_tot_funct);
        g_sig_over_tot->SetLineColor(graph_colors[i]);
        peaks_over_totals->Add(g_sig_over_tot);
        
        // Print out the number of pions and the signal to noise ratio
        std::cout << "\nNumber of pions: " << func->GetParameter(0)/MASSWIDTH << "\nError: " << func->GetParError(0)/MASSWIDTH << std::endl;
        std::cout << "\n\nS/T for 0 sigma: " << sig_over_tot_funct->Eval(0.000001) << "\n\nS/T for 1 sigma: " << sig_over_tot_funct->Eval(1) << "\nS/T for 2 sigma: " << sig_over_tot_funct->Eval(2) << std::endl;
        
        // Write data on number of pions and signal to total to the table file
        // Edit the common header to place the sub-interval's momentum bounds in the frame title
        common_header = Form("\\frame\n{\n\\frametitle{Analysis of Cuts: Pt %2.0f-%2.0f GeV}\n\\begin{table}\n\\caption{How cuts affect data quality}\n\\centering\n\\begin{tabular}{c c c c}\n\\hline\\hline\nCuts & \\# Pions1 & S/T at 1 $\\sigma$ & S/T at 2 $\\sigma$ \\\\ [0.5ex]\n\\hline\n", min, max);
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
        int order_num;
        while (!table_file_input.eof()) {
            std::getline(table_file_input, current_line);
            std::stringstream linestream(current_line);
            //cout << current_line << " is current line" << std::endl;
            std::getline(linestream, line_component, '%');
            //cout << line_component << " is line component" << std::endl;
            std::getline(linestream, line_component, '%');
            //cout << line_component << " is line component" << std::endl;
            // Continue to the next line if the line component's first character is not an integer
            if (!isdigit(line_component[0]))
                continue;
            order_num = std::stoi(line_component);
            //cout << "Line originally in file: order number " << order_num << std::endl << std::endl;
            if (order_num < NumOfCuts) {
                before_in = before_in + current_line + "\n";
                //cout << "Line goes before input\n";
            }
            else if (order_num > NumOfCuts) {
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
        table_file_output << common_header << before_in << headers[NumOfCuts] << "& " << Form("%4.0f",(func->GetParameter(0))/MASSWIDTH) << " +/- " << Form("%4.0f",(func->GetParError(0))/MASSWIDTH) << " & " << Form("%2.2f", sig_over_tot_funct->Eval(1)) <<  " & " << Form("%2.2f", sig_over_tot_funct->Eval(2)) << " \\\\ %" << NumOfCuts << "%" << std::endl << after_in << common_footer;
        table_file_output.close();

        
        //Load onto the ROOT file
        hMass->GetListOfFunctions()->Add(peak);
        hMass->GetListOfFunctions()->Add(background);
        hMass->Write(Form("mass-pion-%2.2fGeV-%2.2fGeV", min, max));
        residual->Write(Form("residual-%2.2fGeV-%2.2fGeV", min, max));
        
        // Graph the photon energy ratio
        graphcanvas->Clear();
        auto h2D = h_Pion->Projection(axis_photon2E, axis_photon1E);
        h2D->SetTitle("Photon energy correlation; Leading Photon Energy (GeV); Trailing Photon Energy (GeV)");
        h2D->GetXaxis()->SetRangeUser(6.0, 15.0);
        h2D->GetYaxis()->SetRangeUser(3.0, 10.0);
        h2D->Draw("COLZ");
        myText(.35,.9, kBlack, Form("Photon energies: Pion momenta %2.2f to %2.2f", min, max));
        graphcanvas->SaveAs(Form(str_concat_converter(directory_name, "PhotonEs_Ptmin_%2.2f_Ptmax_%2.2f.png"), min, max));
        h2D->Write(Form("Photon_Energies_Ptmin_%2.2f_Ptmax_%2.2f", min, max));
    }// end of loop over pt
    
    // Graph the signal/total curves for each momentum increment
    graphcanvas->Clear();
    peaks_over_totals->SetTitle("Signal over Total vs Distance From Mean; Num of Standard Deviations From Mean; Signal to Total Ratio");
    //peaks_over_totals->GetYaxis()->SetRangeUser(0.4, 1.0);
    peaks_over_totals->SetMaximum(1.0);
    peaks_over_totals->SetMinimum(0.2);
    peaks_over_totals->Draw("Al");
    myBoxText(0.25, 0.45, 0.05, 10, graph_colors[0], "8-10 GeV");
    myBoxText(0.25, 0.40, 0.05, 10, graph_colors[1], "10-11 GeV");
    myBoxText(0.25, 0.35, 0.05, 10, graph_colors[2], "11-12 GeV");
    myBoxText(0.25, 0.30, 0.05, 10, graph_colors[3], "12-13 GeV");
    myBoxText(0.25, 0.25, 0.05, 10, graph_colors[4], "13-15 GeV");
    myBoxText(0.25, 0.20, 0.05, 10, graph_colors[5], "15-20 GeV");
    peaks_over_totals->Write("signal-over-total");
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
    g_mean->GetYaxis()->SetRangeUser(0.125, 0.156);
    //g_mean->SetMarkerSize(2);
    //g_mean->SetMarkerStyle(20);
    g_mean->Draw("AP");
    mass_pdg->Draw("same");
    g_mean->Write("mean-masses");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "meanMass_v_pT.png"));
    
    // Graph mass standard deviations with error bars
    graphcanvas->Clear();
    TGraphErrors* g_sigma = new TGraphErrors(num_of_intervals, center, sigmas, widths, sigma_errors);
    g_sigma->Print();
    g_sigma->SetTitle("Mass Peak Widths for Various Momenta; Momentum (GeV); Mass Width (MeV/c^2)");
    g_sigma->GetYaxis()->SetRangeUser(3.0, 16.0);
    //g_sigma->SetMarkerSize(2);
    //g_sigma->SetMarkerStyle(20);
    g_sigma->Write("standard-dev-masses");
    g_sigma->Draw("AP");
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
    graphcanvas->SaveAs(str_concat_converter(directory_name, "reduced_chisquare_v_pT.png"));
    
    // Graph the distribution integrals over momentum
    graphcanvas->Clear();
    TGraphErrors* g_integral = new TGraphErrors(num_of_intervals, center, gaussian_integrals, widths, integral_errors);
    g_integral->Print();
    g_integral->SetTitle("Peak integrals for Various Momenta; Momentum (GeV); Number of Pions");
    //g_integral->GetYaxis()->SetRangeUser(0.0, 6000.0);
    g_integral->Draw("AP");
    g_integral->Write("pion-integrals");
    graphcanvas->SaveAs(str_concat_converter(directory_name, "peakIntegrals_v_pT.png"));
    
    graphcanvas->Close();
    logcanvas->Close();
    return;
}
