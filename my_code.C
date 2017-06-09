/**
 my_code is a data-processing macro meant for processing THnSparses root files, specifically the h_Pion one in THnSparses_060717.root
 the macro must be called with two bools, the first to tell if the lambda must be cut, the second to tell if asymmetry must be cut
 */

#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TGraphErrors.h>
#include <iostream>

// The general filepath string
string directory_name;

//variables of hPion
const int axis_pion_Cen    = 0;
const int axis_pion_Zvtx   = 1;
const int axis_pionMass    = 2;
const int axis_pionPt      = 3;
const int axis_asymmetry   = 9;
const int axis_pionAngle   = 16;
const int axis_pionLambda1 = 17;
const int axis_pionLambda2 = 18;
const int axis_pionNcells1 = 19;
const int axis_pionNcells2 = 20;

/**
// Concatenates two C-strings without changing the contents of the input parameters
char* output_strcat(char* input1, char* input2){
    char *sumstring = new char[strlen(input1) + strlen(input2)];
    strcpy(sumstring, input1);
    strcat(sumstring, input2);
    
    return sumstring;
}
*/

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
// Fitting models: Ae^(x-mean)^2/(2*sigma^2) + B*sin((x-C)/D)/x^E + F
// Ae^(x-mean)^2/(2*sigma^2) + B*(x - C)^D + E
// Currently chosen one: Ae^(x-mean)^2/(2*sigma^2) + B*x^4 + C*x^3 + D*x^2 + E*x + F
// Ae^(x-mean)^2/(2*sigma^2) + B/(1 + Ce^D(x-E)) + F
double compound_model(Double_t *x,Double_t *par) {
    double arg = 0;

    
    double A = par[0];
    double mean = par[1];
    double sigma = par[2];
    double B = par[3];
    double C = par[4];
    double D = par[5];
    double E = par[6];
    double F = par[7];
    
    if (sigma != 0)
        arg = (x[0] - mean)/sigma;

    
    //double fitval = A*TMath::Exp(-0.5*arg*arg) + B*TMath::Sin((x[0] - C)/D)/TMath::Power(x[0], E) + F;
    //double fitval = A*TMath::Exp(-0.5*arg*arg) + B*TMath::Power((x[0] - C), D) + E;
    double fitval = A*(1.0/(sigma*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*arg*arg) + B*TMath::Power(x[0], 4.0) + C*TMath::Power(x[0], 3.0) + D*x[0]*x[0] + E*x[0] + F;
    //double fitval = A*TMath::Exp(-0.5*arg*arg) + B/(1 + C*TMath::Exp(D*(x[0]-E))) + F;
    return fitval;
}

// The two modeling functions used in this program: the first for the entire data, second for the peak alone
TF1 *func = new TF1("fit", compound_model,0.05,0.5,8);
TF1 *peak = new TF1("mass peak", "[0]*(1.0/([2]*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*((x - [1])/[2])*((x - [1])/[2]))", 0.08, 0.26);

// Takes the seven parameters that would be passed to func, along with the number of sigmas away from the mean ("x"), and returns the ratio of the integral of
// peak within the interval that is within x sigmas from the mean and the integral of peak over the same integral
double signal_over_total(Double_t *x, Double_t *par) {
    // Extract the parameters
    //func->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]);
    //peak->SetParameters(par[0], par[1], par[2]);
    double Nsigma = x[0];
    double mean = par[1];
    double sigma = par[2];
    
    // Evaluate the integrals and return their quotient
    double signal = peak->Integral(mean-Nsigma*sigma, mean+Nsigma*sigma);
    double total = func->Integral(mean-Nsigma*sigma, mean+Nsigma*sigma);
    
    return signal/total;
}

// Takes a THnSparse pointer (almost always h_Pion), an hPion variable number, a Canvas object pointer, and a file name
// Graphs the data under the hPion variable number and saves the graph as a .png
// SetCuts done to the THnSparse object before the function call do apply
// Not advisable if you want to save any variable from the data or the TH1D object used to graph it
void graph_raw_data(THnSparse* data, const int hPion_var, TCanvas* can, char* filename){
    TH1D* data_subset = data->Projection(hPion_var);
    data_subset->Draw();
    
    can->SaveAs(filename);
}

/**
 Main function
 */
void my_code(int num_of_cuts){
    // Generate the directory name to store files in, using the given int
    // num_of_cuts = 0 means no cuts, 1 means only lambda_2 is cut, 2 means lambda_2 and asymmetry are cut,
    // 3 means lambda_2, asymmetry, and angle are cut, 4 means lambda_2, asymmetry, angle, and ncells are cut
    // Throw an exception if any other int is given, and re-request cut input
    const int int_out_of_bounds_exception = 100;
    bool all_clear;
    do {
        try {
            if (num_of_cuts < 0 || num_of_cuts > 4)
                throw int_out_of_bounds_exception;
            else
                all_clear = true;
        }
    
        catch (int exc) {
            if (exc == int_out_of_bounds_exception){
                cout << "ERROR: only acceptable numbers of cuts are 0, 1, 2, 3, and 4. Please enter another number of cuts below:\n";
                cin >> num_of_cuts;
                all_clear = false;
            }
        }
    }
    while(!all_clear);
    directory_name = std::to_string(num_of_cuts) + "cuts/";
    
    std::cout << "Directory chosen: " << directory_name << std::endl;
    
    //Open the file
    TFile* fIn = new TFile("THnSparses_060717.root","READ"); //get file
    fIn->Print(); //print file content
    
    // Get the data
    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion); //get array
    
    //For the mass plot, restrict to mass between 0.08 and 0.25, define the Canvas,
    //plot the data for the asymmetry, both lambdas, the angle, and number of cells
    SetCut(h_Pion, axis_pionMass, 0.08, 0.25);
    TCanvas* canvas = new TCanvas();
    
    graph_raw_data(h_Pion, axis_asymmetry, canvas, str_concat_converter(directory_name,"asymmetry_pion_plot.png"));
    graph_raw_data(h_Pion, axis_pionLambda1, canvas, str_concat_converter(directory_name, "lambda1_pion_plot.png"));
    graph_raw_data(h_Pion, axis_pionLambda2, canvas, str_concat_converter(directory_name, "lambda2_pion_plot.png"));
    graph_raw_data(h_Pion, axis_pionAngle, canvas, str_concat_converter(directory_name, "angle_pion_plot.png"));
    graph_raw_data(h_Pion, axis_pionNcells1, canvas, str_concat_converter(directory_name, "Ncells1_pion_plot.png"));
    graph_raw_data(h_Pion, axis_pionNcells2, canvas, str_concat_converter(directory_name, "Ncells2_pion_plot.png"));
    
    // restrict asymmetry to below 0.7, lambda 1 and 2 to below 0.4, the angle absolute value to above 0.015, and Ncells 1 and 2 to above 1.5
    // as requested by num_of_cuts according to the first block of comments in this function
    if(num_of_cuts == 0) {
        cout << "No cuts done\n\n";
    }
    else if(num_of_cuts == 1) {
        SetCut(h_Pion, axis_pionLambda1, 0.0, 0.4);
        SetCut(h_Pion, axis_pionLambda2, 0.0, 0.4);
        cout << "Cuts: lambda02\n\n";
    }
    else if(num_of_cuts == 2) {
        SetCut(h_Pion, axis_pionLambda1, 0.0, 0.4);
        SetCut(h_Pion, axis_pionLambda2, 0.0, 0.4);
        SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
        cout << "Cuts: lambda02 and asymmetry\n\n";
    }
    else if(num_of_cuts == 3) {
        SetCut(h_Pion, axis_pionLambda1, 0.0, 0.4);
        SetCut(h_Pion, axis_pionLambda2, 0.0, 0.4);
        SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
        SetCut(h_Pion, axis_pionAngle, 0.015, 0.5);
        cout << "Cuts: lambda02, asymmetry, and angle\n\n";
    }
    else if(num_of_cuts == 4) {
        SetCut(h_Pion, axis_pionLambda1, 0.0, 0.4);
        SetCut(h_Pion, axis_pionLambda2, 0.0, 0.4);
        SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
        SetCut(h_Pion, axis_pionAngle, 0.015, 0.5);
        SetCut(h_Pion, axis_pionNcells1, 1.0, 30.0);
        SetCut(h_Pion, axis_pionNcells2, 1.0, 30.0);
        cout << "Cuts: lambda02, asymmetry, angle, and Ncells\n\n";

    }
    
    // plot mass data and set up the fit function
    TH1D* hMass = h_Pion->Projection(axis_pionMass);
    const double MASSWIDTH = hMass->GetBinWidth(1);
    hMass->SetAxisRange(0., 7000., "Y");
    hMass->Draw();
    
    //Restrict the parameters to reasonable ranges, insert guess values, and give understandable names
    func->SetParameters(600,  0.14, 0.3,  1, 0.03, 0.6, 0.1, 1);
    func->SetParLimits(0, 1, 10000.0);//integral
    func->SetParLimits(1, 0.1, 0.2); //mean
    func->SetParLimits(2, 0.01, 0.05); // width
    func->SetParLimits(3, -100000.0, 0.0); // Quadric and quadratic factors
    func->SetParLimits(5, -100000.0, 0.0);
    //func->SetParameters(600,  0.14, 0.3,  1, 0.03, 0.6, 0.1, 1);
    //func->SetParLimits(0, 1, 10000.0);//integral
    //func->SetParLimits(1, 0.1, 0.16); //mean
    //func->SetParLimits(2, 0.005, 0.03); // width
    //func->SetParLimits(3, -10000000.0, 0.0); // Quadric and quadratic factors
    //func->SetParLimits(4, -100000.0, 1000000.0);
    //func->SetParLimits(5, -1000000.0, 0.0);
    //func->SetParNames("Amplitude", "Mean", "Sigma", "Damped Sine coeff", "Shift", "Period Factor", "Damping Exponent", "Constant");
    func->SetParNames("Integral", "Mean", "Sigma", "Quadric coeff", "Cubic coeff", "Quadratic coeff", "Linear coeff", "Constant");
    //func->SetParNames("Amplitude", "Mean", "Sigma", "Logistic asymptote", "e-coeff", "In-exponent coeff", "Shift", "Constant");
    
    // Plot the fit for the mass and (separately) the Gaussian component of it; save as a PDF
    hMass->Fit(func);
    func->Draw("same");
    std::cout << "Reduced Chi Square " << (func->GetChisquare())/10 << std::endl; //Reduced Chi Square of the mass vs entries curve (function has 18 degrees of freedom, 7 parameters)
    for(int i = 0; i < 3; i++) {
         peak->SetParameter(i, func->GetParameter(i));
    }
    peak->SetLineColor(kBlue);
    peak->Draw("same");
    canvas->SaveAs(str_concat_converter(directory_name, "mass_pion_plot.png"));
    
    // Plot the data for the momentum
    TH1D* hPt = h_Pion->Projection(axis_pionPt);
    hPt->Draw();
    canvas->SaveAs(str_concat_converter(directory_name, "momentum_pion_plot.png"));
    
    // A collection of variables that is needed for the next steps
    const int num_of_intervals = 7;
    TMultiGraph* peaks_over_totals = new TMultiGraph();
    Color_t graph_colors[num_of_intervals] = {kRed, kBlue, 8, kYellow, kCyan, kMagenta, kBlack};
    
    double intervals[num_of_intervals][2] = {{5.0, 7.5}, {7.5, 10.0}, {10.0, 11.0}, {11.0, 12.0}, {12.0, 13.0}, {13.0, 15.0}, {15.0, 18.0}};
    double chisquares[num_of_intervals];
    double means[num_of_intervals];
    double mean_errors[num_of_intervals];
    double sigmas[num_of_intervals];
    double sigma_errors[num_of_intervals];
    double gaussian_integrals[num_of_intervals];
    double integral_errors[num_of_intervals];
    
    double center[num_of_intervals];
    double widths[num_of_intervals];

    // start cutting the data up; plot the mass data for momenta of 5-10, 10-15, 15-20, and so forth
    for(int i = 0; i < num_of_intervals; i++) {
    //for(int i = 0; i < 1; i++) {
        double ptmin = intervals[i][0];
        double ptmax = intervals[i][1]; // Interval bounds
        
        // Plot the data
        SetCut(h_Pion, axis_pionPt, ptmin, ptmax);
        
        hMass = h_Pion->Projection(axis_pionMass);
        hMass->SetAxisRange(0.0, 1600.0, "Y");
        hMass->Draw();
        
        // Find a fit just as you did for the entire data set
        // Graph the fit and (separately) the Gaussian component of it
        hMass->Fit(func);
        chisquares[i] = (func->GetChisquare())/10; //Reduced Chi Square (function has 18 degrees of freedom, 7 parameters)
        std::cout << Form("Reduced Chi Square: %2.2f", chisquares[i]) << std::endl;
        func->Draw("same");
        for(int i = 0; i < 3; i++) {
            peak->SetParameter(i, func->GetParameter(i));
        }
        peak->SetLineColor(kBlue);
        peak->Draw("same");
        
        // Add the mean mass parameter and its error to means and mean_errors, respectively; the standard deviation and
        // its error to sigmas and sigma_errors, the integral of the Gaussian peak and its error to gaussian_integrals
        // and integral_errors, the center point of the interval into center (and its error into widths), and chi square of the fit into its respective array
        means[i] = func->GetParameter(1);
        mean_errors[i] = func->GetParError(1);
        sigmas[i] = func->GetParameter(2) * 1000;
        sigma_errors[i] = func->GetParError(2) * 1000;
        
        gaussian_integrals[i] = (func->GetParameter(0))/MASSWIDTH;
        integral_errors[i] = (func->GetParError(0))/MASSWIDTH;
        
        center[i] = (ptmin+ptmax)/2.0;
        widths[i] = ptmax - center[i];
        
        // print out the mean, standard deviation, and their corresponding errors; and save the graph in a new file
        std::cout << "Mean: " << means[i] << std::endl << "Standard Deviation: " << sigmas[i] << std::endl;
        canvas->SaveAs(Form(str_concat_converter(directory_name, "MyFit_Ptmin_%2.2f_Ptmax_%2.2f.png"), ptmin, ptmax));
        
        //Now use the signal_over_total function to get a graph of sigma vs. signal/total
        TF1* sig_over_tot_funct = new TF1("Signal over Total", signal_over_total, 0, 4, 7);
        sig_over_tot_funct->SetParameters(func->GetParameters());
        canvas->Clear();
        
        sig_over_tot_funct->SetTitle(Form("Signal over Total vs Distance From Mean: %2.2f to %2.2f GeV; Num of Standard Deviations From Mean; Signal to Total Ratio", ptmin, ptmax));
        ///sig_over_tot_funct->Draw();
        TGraph* g_sig_over_tot = new TGraph(sig_over_tot_funct);
        g_sig_over_tot->SetLineColor(graph_colors[i]);
        //g_sig_over_tot->GetYaxis()->SetRangeUser(0.4, 1.0);
        peaks_over_totals->Add(g_sig_over_tot);
        //canvas->SaveAs(Form(str_concat_converter(directory_name, "Signal_Over_Total_Ptmin_%2.2f_Ptmax_%2.2f.png"), ptmin, ptmax));
        
    }
    
    // Graph the signal/total curves for each momentum increment
    // Red = 5-7.5, Blue = 7.5-10, Green = 10-11, Yellow = 11-12, Brown = 12-13, Magenta = 13-15, Black = 15-18
    canvas->Clear();
    peaks_over_totals->SetTitle("Signal over Total vs Distance From Mean; Num of Standard Deviations From Mean; Signal to Total Ratio");
    //peaks_over_totals->GetYaxis()->SetRangeUser(0.4, 1.0);
    peaks_over_totals->SetMaximum(1.0);
    peaks_over_totals->SetMinimum(0.4);
    peaks_over_totals->Draw("Al");
    canvas->SaveAs(str_concat_converter(directory_name, "Overall_Signal_Over_Total.png"));
    canvas->Clear();

    // Add a constant "expected mass" function to be graphed alongside the mean mass data
    TF1* mass_pdg = new TF1("mass_pdg", "[0]", 0, 20);
    mass_pdg->SetParameter(0, 0.13498);
    mass_pdg->SetLineWidth(2);
    mass_pdg->SetLineColor(kRed);
    
    // Graph mean masses with error bars
    TGraphErrors* g_mean = new TGraphErrors(num_of_intervals, center, means, widths, mean_errors);
    canvas->Clear();
    g_mean->Print();
    g_mean->SetTitle("Mean Masses for Various Momenta; Momentum (GeV); Mass (GeV/c^2)");
    g_mean->GetYaxis()->SetRangeUser(0.133, 0.156);
    //g_mean->SetMarkerSize(2);
    //g_mean->SetMarkerStyle(20);
    g_mean->Draw("AP");
    mass_pdg->Draw("same");
    canvas->SaveAs(str_concat_converter(directory_name, "meanMass_v_pT.png"));
    
    // Graph mass standard deviations with error bars
    canvas->Clear();
    TGraphErrors* g_sigma = new TGraphErrors(num_of_intervals, center, sigmas, widths, sigma_errors);
    g_sigma->Print();
    g_sigma->SetTitle("Mass Peak Widths for Various Momenta; Momentum (GeV); Mass (MeV/c^2)");
    g_sigma->GetYaxis()->SetRangeUser(8.5, 16.5);
    //g_sigma->SetMarkerSize(2);
    //g_sigma->SetMarkerStyle(20);
    g_sigma->Draw("AP");
    canvas->SaveAs(str_concat_converter(directory_name, "massWidths_v_pT.png"));
    
    // Graph reduced chi squares over the momentum interval
    canvas->Clear();
    TGraph* g_chisquare = new TGraphErrors(num_of_intervals, center, chisquares);
    g_chisquare->Print();
    g_chisquare->SetTitle("ReducedChi-Squares for Various Momenta; Momentum (GeV); Reduced Chi-Squares");
    g_chisquare->SetMarkerSize(2);
    g_chisquare->SetMarkerStyle(20);
    g_chisquare->GetYaxis()->SetRangeUser(2.0, 20.0);
    g_chisquare->Draw("AP");
    canvas->SaveAs(str_concat_converter(directory_name, "reduced_chisquare_v_pT.png"));
    
    // Graph the Gaussian distribution integrals over momentum
    canvas->Clear();
    TGraphErrors* g_integral = new TGraphErrors(num_of_intervals, center, gaussian_integrals, widths, integral_errors);
    g_integral->Print();
    g_integral->SetTitle("Gaussian Peak integrals for Various Momenta; Momentum (GeV); Number of Pions");
    g_integral->GetYaxis()->SetRangeUser(0.0, 6000.0);
    g_integral->Draw("AP");
    canvas->SaveAs(str_concat_converter(directory_name, "peakIntegrals_v_pT.png"));


    canvas->Close();
}
