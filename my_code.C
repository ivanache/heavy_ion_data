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

// The output file
TFile* fOut;

// The general filepath string
string directory_name;

// The string with the model name
string model;

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

// The fitting functions
// Fitting models:
// Primary one: Ae^(x-mean)^2/(2*sigma^2) + B*x^4 + C*x^3 + D*x^2 + E*x + F
//Ae^(x-mean)^2/(2*sigma^2) + B*sin((x-C)/D)/x^E + F
// Ae^(x-mean)^2/(2*sigma^2) + B*(x - C)^D + E
// Ae^(x-mean)^2/(2*sigma^2) + B/(1 + Ce^D(x-E)) + F
// Ae^(x-mean)^2/(2*sigma^2) + B*ln(x[0])/ln(C) + E
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
    
    double fitval = A*(1.0/(sigma*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*arg*arg) + B*TMath::Power(x[0], 4.0) + C*TMath::Power(x[0], 3.0) + D*x[0]*x[0] + E*x[0] + F;
    return fitval;
}
// Non-bell curve peaks
// Gaussian plus exponential
double exp_gaussian_model(Double_t *x, Double_t* par) {
    double arg = 0;
    
    double A = par[0];
    double mean = par[1];
    double sigma = par[2];
    double lambda = par[3];
    double B = par[4];
    double C = par[5];
    double D = par[6];
    double E = par[7];
    double F = par[8];
    
    if (sigma != 0)
        arg = (mean + lambda*sigma*sigma - x[0])/(TMath::Sqrt(2)*sigma);
    
    double lambda_over_two = lambda/2.0;
    
    double fitval = A*lambda_over_two*TMath::Exp(lambda_over_two*(2*mean + lambda*sigma*sigma - 2*x[0]))*TMath::Erfc(arg) + B*TMath::Power(x[0], 4.0) + C*TMath::Power(x[0], 3.0) + D*x[0]*x[0] + E*x[0] + F;
    return fitval;
}

// Modified Rayleigh distribution
/**double modified_rayleigh_model(Double_t *x, Double_t* par) {
    double arg = 0;
    
    double A = par[0];
    double sigma = par[1];
    double horizontal_shift = par[2];
    double B = par[3];
    double C = par[4];
    double D = par[5];
    double E = par[6];
    double F = par[7];
    
    double modified_x = x[0] - horizontal_shift;
    if (sigma != 0)
        arg = modified_x/sigma;
    
    double fitval = A*(modified_x*TMath::Exp(-0.5*arg*arg))/(sigma*sigma) + B*TMath::Power(x[0], 4.0) + C*TMath::Power(x[0], 3.0) + D*x[0]*x[0] + E*x[0] + F;
    return fitval;
} */

double skew_normal_dist(Double_t *x, Double_t* par) {
    double arg = 0;
    
    double A = par[0];
    double position = par[1];
    double scale = par[2];
    double shape = par[3];
    double B = par[4];
    double C = par[5];
    double D = par[6];
    double E = par[7];
    double F = par[8];
    
    if (scale != 0)
        arg = ((x[0] - position)/scale);
    
    double sqrt_pi_2 = TMath::Sqrt(TMath::Pi())/TMath::Sqrt(2);
    
    double fitval = (A/(scale*TMath::Pi()))*TMath::Exp(-0.5*arg*arg)*sqrt_pi_2*(TMath::Erf(shape*(x[0] - position)/(scale*TMath::Sqrt(2))) + 1.0) + B*TMath::Power(x[0], 4.0) + C*TMath::Power(x[0], 3.0) + D*x[0]*x[0] + E*x[0] + F;
    return fitval;
}

/**
 double sin_model(Double_t *x,Double_t *par) {
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
 
 double fitval = A*(1.0/(sigma*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*arg*arg) + B*TMath::Sin((x[0] - C)/D)/TMath::Power(x[0], E) + F;
 return fitval;
 }
 double logistic_model(Double_t *x,Double_t *par){
 double arg = 0;
 
 double A = par[0];
 double mean = par[1];
 double sigma = par[2];
 double B = par[3];
 double C = par[4];
 double D = par[5];
 double E = par[6];
 
 if (sigma != 0)
 arg = (x[0] - mean)/sigma;
 
 double fitval = A*(1.0/(sigma*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*arg*arg) + B*TMath::Log(x[0] - C)/TMath::Log(D) + E;
 return fitval;
 
 }
 double logarithmic_model(Double_t *x,Double_t *par) {
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
 
 double fitval = A*(1.0/(sigma*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*arg*arg) + B/(1 + C*TMath::Exp(D*(x[0]-E))) + F;
 return fitval;
 }
 double power_model(Double_t *x,Double_t *par) {
 double arg = 0;
 
 double A = par[0];
 double mean = par[1];
 double sigma = par[2];
 double B = par[3];
 double C = par[4];
 double D = par[5];
 double E = par[6];
 
 if (sigma != 0)
 arg = (x[0] - mean)/sigma;
 
 double fitval = A*(1.0/(sigma*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*arg*arg) + B*TMath::Power((x[0] - C), D) + E;
 return fitval;
 }
 */

// The three functions used in this program: the first for the entire data, second for the peak alone, third for the background alone
TF1 *func;
TF1 *peak;
TF1 *background;

// Takes the seven parameters that would be passed to func, along with the number of sigmas away from the mean ("x"), and returns the ratio of the integral of
// peak within the interval that is within x sigmas from the mean and the integral of peak over the same integral
// Needs to be modernized for the new plots
double signal_over_total(Double_t *x, Double_t *par) {
    // Extract the parameters
    //func->SetParameters(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]);
    //peak->SetParameters(par[0], par[1], par[2]);
    double Nsigma = x[0];
    // Based on the peak model, derive the mean and sigma
    double mean;
    double sigma;
    if (model == "Gaussian"){
        mean = par[1];
        sigma = par[2];
    }
    if (model == "Exponential_Gaussian") {
        mean = par[1] + (1.0/par[3]);
        sigma = TMath::Sqrt((par[2]*par[2]) + (1.0/(par[3]*par[3])));
    }
    if (model == "Skew_Normal") {
        double delta = par[3]/TMath::Sqrt(1 + par[3]*par[3]);
        mean = par[1] + (par[2]*delta*TMath::Sqrt(2.0/TMath::Pi()));
        sigma = par[2]*TMath::Sqrt(1 - (2.0*delta*delta/TMath::Pi()));
    }
    
    // Evaluate the integrals and return their quotient
    double signal = peak->Integral(mean-Nsigma*sigma, mean+Nsigma*sigma);
    double total = func->Integral(mean-Nsigma*sigma, mean+Nsigma*sigma);
    
    return signal/total;
    //return 0.0;
}

// Takes a THnSparse pointer (almost always h_Pion), an hPion variable number, a Canvas object pointer, and a file name
// Graphs the data under the hPion variable number and saves the graph as a .png
// SetCuts done to the THnSparse object before the function call do apply
// Not advisable if you want to save any variable from the data or the TH1D object used to graph it
void graph_raw_data(THnSparse* data, const int hPion_var, TCanvas* can, char* filename, string rootname){
    TH1D* data_subset = data->Projection(hPion_var);
    data_subset->Draw();
    can->SaveAs(filename);
    data_subset->Write(rootname.c_str());
}

/**
 Main function
 */
// Precondition: model_name is "Gaussian", "Exponential_Gaussian", "Skew_Normal"
void my_code(string model_name) {
    model = model_name;
    directory_name = "data/" + model_name + "/";
    
    //Open the files
    TFile* fIn = new TFile("THnSparses_060717.root","READ"); //get file
    string rootfilename = "data/Pion" + model_name + "SparsesOutput.root";
    fOut = new TFile(rootfilename.c_str(), "RECREATE"); // Create an output file
    fIn->Print(); //print file content
    
    // Get the data
    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion); //get array
    
    //The TCanvas and TPads (for sub-canvassing)
    TCanvas* canvas = new TCanvas();
    TPad *pad[2] = {new TPad("pad0","",0,0.42,1,1), new TPad("pad1","",0,0,1,0.42)};
    
    //For the mass plot, restrict to mass between 0.08 and 0.25 and the momentum to between 5 GeV and 18 GeV
    //plot the data for the asymmetry, both lambdas, the angle, and number of cells
    SetCut(h_Pion, axis_pionPt, 8.0, 15.0);
    //SetCut(h_Pion, axis_pionPt, 5.0, 18.0);
    SetCut(h_Pion, axis_pionMass, 0.08, 0.25);
    
    graph_raw_data(h_Pion, axis_asymmetry, canvas, str_concat_converter(directory_name,"asymmetry_pion_plot.png"), "asymmetry");
    graph_raw_data(h_Pion, axis_pionLambda1, canvas, str_concat_converter(directory_name, "lambda1_pion_plot.png"), "lambda1");
    graph_raw_data(h_Pion, axis_pionLambda2, canvas, str_concat_converter(directory_name, "lambda2_pion_plot.png"), "lambda2");
    graph_raw_data(h_Pion, axis_pionAngle, canvas, str_concat_converter(directory_name, "angle_pion_plot.png"), "angle");
    graph_raw_data(h_Pion, axis_pionNcells1, canvas, str_concat_converter(directory_name, "Ncells1_pion_plot.png"), "Ncells1");
    graph_raw_data(h_Pion, axis_pionNcells2, canvas, str_concat_converter(directory_name, "Ncells2_pion_plot.png"), "Ncells2");
    
    // restrict asymmetry to below 0.7, lambda 1 and 2 to below 0.4, the angle absolute value to above 0.015, and Ncells 1 and 2 to above 1.5
    // Also set the maximum y value of the pion entries vs mass to 600
    SetCut(h_Pion, axis_pionLambda1, 0.0, 0.4);
    SetCut(h_Pion, axis_pionLambda2, 0.0, 0.4);
    SetCut(h_Pion, axis_asymmetry, 0.0, 0.7);
    SetCut(h_Pion, axis_pionAngle, 0.015, 0.5);
    SetCut(h_Pion, axis_pionNcells1, 1.0, 30.0);
    SetCut(h_Pion, axis_pionNcells2, 1.0, 30.0);
    std::cout << "Cuts: lambda02, asymmetry, angle, and Ncells\n\n";
    double fit_y_max = 600.0;
    
    // plot mass data and set up the fit function
    TH1D* hMass = h_Pion->Projection(axis_pionMass);
    TH1D* residual = (TH1D*)hMass->Clone("residual");
    const double MASSWIDTH = hMass->GetBinWidth(1);
    hMass->SetAxisRange(0., 2500., "Y");
    hMass->GetYaxis()->SetTitle("Number of Entries");
    hMass->GetYaxis()->SetTitleSize(.05);
    hMass->GetYaxis()->SetTitleOffset(0.5);
    canvas->cd();
    pad[0]->Draw();
    pad[0]->cd();
    hMass->Draw();
    
    
    //Start making the fit, restrict the parameters to reasonable ranges, insert guess values, and give understandable names
    int num_of_params = 8;
    int num_of_peak_params = 3;
    background = new TF1("background curve", "[0]*TMath::Power(x, 4.0) + [1]*TMath::Power(x, 3.0) + [2]*x*x + [3]*x + [4]", 0.08, 0.26);
    if (model_name == "Gaussian") {
        func = new TF1("fit", compound_model,0.05,0.5,num_of_params);
        peak = new TF1("mass peak", "[0]*(1.0/([2]*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*((x - [1])/[2])*((x - [1])/[2]))", 0.08, 0.26);
        func->SetParNames("Integral", "Mean", "Sigma", "Quadric coeff", "Cubic coeff", "Quadratic coeff", "Linear coeff", "Constant");
        func->SetParameters(60,  0.14, 0.3,  -100000, 30000, -60000, 0, 10000);
        func->SetParLimits(0, 1, 10000.0);//integral
        func->SetParLimits(1, 0.1, 0.2); //mean
        func->SetParLimits(2, 0.005, 0.05); // width
        func->SetParLimits(3, -1000000.0, 0.0); // Quadric and quadratic factors
        func->SetParLimits(5, -1000000.0, 0.0);
        
    }
    if (model_name == "Exponential_Gaussian") {
        num_of_params = 9;
        num_of_peak_params = 4;
        func = new TF1("fit", exp_gaussian_model,0.05,0.5,num_of_params);
        peak = new TF1("mass peak", "[0]*([3]/2.0)*TMath::Exp(([3]/2.0)*(2*[1] + [3]*[2]*[2] - 2*x[0]))*TMath::Erfc(([1] + [3]*[2]*[2] - x[0])/(TMath::Sqrt(2)*[2]))", 0.08, 0.26);
        func->SetParNames("Integral", "Mean", "Sigma", "Lambda", "Quadric coeff", "Cubic coeff", "Quadratic coeff", "Linear coeff", "Constant");
        func->SetParameters(60,  0.14, 0.3, 100, -100000, 1000000, -60000, 10, 10000);
        func->SetParLimits(0, 1.0, 50.0);//integral
        func->SetParLimits(1, 0.1, 0.2); //Mean
        func->SetParLimits(2, 0.00005, 0.02); // width
        func->SetParLimits(3, 0.0, 1000.0); //lambda
        func->SetParLimits(4, -1000000.0, 0.0); // Quadric factor
        func->SetParLimits(5, 0.0, 1000000.0); // Cubic factor
        func->SetParLimits(6, -1000000.0, 0.0); // Quadratic factor
        func->SetParLimits(7, 0.0, 10000000.0); // Linear factor
        func->SetParLimits(8, -10000000.0, 0.0); // Constant
    }
    
    if (model_name == "Skew_Normal") {
        num_of_params = 9;
        num_of_peak_params = 4;
        func = new TF1("fit", skew_normal_dist,0.05,0.5,num_of_params);
        peak = new TF1("mass peak", "([0]/([2]*TMath::Pi()))*TMath::Exp(-0.5*((x[0] - [1])/[2])*((x - [1])/[2]))*TMath::Sqrt(TMath::Pi())/TMath::Sqrt(2)*(TMath::Erf([3]*(x - [1])/([2]*TMath::Sqrt(2))) + 1.0)");
        func->SetParNames("Integral", "Location", "Scale", "Shape", "Quadric coeff", "Cubic coeff", "Quadratic coeff", "Linear coeff", "Constant");
        func->SetParameters(60,  0.14, 0.3, 100, -100000, 1000000, -60000, 10, 10000);
        func->SetParLimits(0, 1.0, 100.0); // Integral
        func->SetParLimits(1, 0, 0.3); // Position
        func->SetParLimits(2, 0.005, 0.05); // width
        func->SetParLimits(3, -10, 10); //Shape
        func->SetParLimits(4, -1000000.0, 0.0); // Quadric and quadratic factors
        func->SetParLimits(6, -1000000.0, 0.0);
    }
    
    /**
     if (model_name == "Logistic") {
     func = new TF1("fit", logistic_model,0.05,0.5,num_of_params);
     background = new TF1("background curve", "[0]/(1 + [1]*TMath::Exp([2]*(x-[3]))) + [4]", 0.08, 0.26);
     }
     if (model_name == "Sine") {
     func = new TF1("fit", sin_model,0.05,0.5,num_of_params);
     background = new TF1("background curve", "[0]*TMath::Sin((x - [1])/[2])/TMath::Power(x, [3]) + [4]", 0.08, 0.26);
     }
     if (model_name == "Logarithmic") {
     num_of_params = 7;
     func = new TF1("fit", logarithmic_model,0.05,0.5,num_of_params);
     background = new TF1("background curve", "[0]*TMath::Log(x - [1])/TMath::Log([2]) + [3]", 0.08, 0.26);
     }
     if (model_name == "Power") {
     num_of_params = 7;
     func = new TF1("fit", power_model,0.05,0.5,num_of_params);
     background = new TF1("background curve", "[0]*TMath::Power((x - [1]), [2]) + [3]", 0.08, 0.26);
     }
     if (model_name == "Logarithmic")
     func->SetParNames("Integral", "Mean", "Sigma", "Logarithmic coeff", "Horiz shift", "Base", "Constant");
     func->SetParameters(60,  0.14, 0.3,  3.0, 0.13, 9, 1, 10);
     func->SetParLimits(0, 1, 10000.0);//integral
     func->SetParLimits(1, 0.1, 0.2); //mean
     func->SetParLimits(2, 0.005, 0.05); // width
     if (model_name == "Sine") {
     func->SetParNames("Integral", "Mean", "Sigma", "Damped Sine coeff", "Shift", "Period Factor", "Damping Exponent", "Constant");
     func->SetParameters(60,  0.14, 0.3,  3.0, 0.13, 9, 1, 10);
     func->SetParLimits(0, 1, 10000.0);//integral
     func->SetParLimits(1, 0.1, 0.2); //mean
     func->SetParLimits(2, 0.005, 0.05); // width
     }
     if (model_name == "Logistic") {
     func->SetParNames("Integral", "Mean", "Sigma", "Logistic asymptote", "e-coeff", "In-exponent coeff", "Shift", "Constant");
     func->SetParameters(6,  0.14, 0.3,  3.0, 0.13, 9, 1, 10);
     func->SetParLimits(0, 1, 10000.0);//integral
     func->SetParLimits(1, 0.1, 0.2); //mean
     func->SetParLimits(2, 0.005, 0.05); // width
     }
     if (model_name == "Power") {
     func->SetParNames("Integral", "Mean", "Sigma", "Power Coeff", "Horiz shft", "Exponent", "Constant");
     func->SetParameters(60,  0.14, 0.3,  3.0, 0.13, 9, 1, 10);
     func->SetParLimits(0, 1, 10000.0);//integral
     func->SetParLimits(1, 0.1, 0.2); //mean
     func->SetParLimits(2, 0.005, 0.05); // width
     }
     */
    //func->SetParNames("Amplitude", "Mean", "Sigma", "Logistic asymptote", "e-coeff", "In-exponent coeff", "Shift", "Constant");
    
    // Plot the fit for the mass and (separately) the Gaussian and background components of it
    hMass->Fit(func);
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
    hMass->GetListOfFunctions()->Add(peak);
    hMass->GetListOfFunctions()->Add(background);
    hMass->Write("mass_pion");// Load into the ROOT file
    
    // Plot the residual; save as a PDF
    for (int i = 0; i < hMass->GetSize(); i++) {
        if (hMass->GetBinError(i))
            residual->SetBinContent(i, ((hMass->GetBinContent(i) - func->Eval(hMass->GetBinCenter(i)))/hMass->GetBinError(i)));
        residual->SetBinError(i, 0); //Residuals don't have errors
        //std::cout<< "Bin: " << i << std::endl;
        //std::cout<< (hMass->GetBinContent(i)) << std::endl;
        //std::cout<< (func->Eval(hMass->GetBinCenter(i))) << std::endl;
        //std::cout<< (hMass->GetBinError(i)) << std::endl << std::endl;
    }
    cout << std::endl;
    canvas->cd();
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
    canvas->SaveAs(str_concat_converter(directory_name, "mass_pion_plot.png"));
    
    // Plot the data for the momentum
    canvas->Clear();
    TH1D* hPt = h_Pion->Projection(axis_pionPt);
    hPt->Draw();
    canvas->SaveAs(str_concat_converter(directory_name, "momentum_pion_plot.png"));
    
    // A collection of variables that is needed for the next steps
    const int num_of_intervals = 5;
    TMultiGraph* peaks_over_totals = new TMultiGraph();
    Color_t graph_colors[num_of_intervals] = {kBlack, kRed, kBlue, kGreen, kYellow};
    
    double intervals[num_of_intervals][2] = {{8.0, 10.0}, {10.0, 11.0}, {11.0, 12.0}, {12.0, 13.0}, {13.0, 15.0}};
    double chisquares[num_of_intervals];
    double means[num_of_intervals];
    double mean_errors[num_of_intervals];
    double sigmas[num_of_intervals];
    double sigma_errors[num_of_intervals];
    double gaussian_integrals[num_of_intervals];
    double integral_errors[num_of_intervals];
    
    double center[num_of_intervals];
    double widths[num_of_intervals];
    canvas->Clear();
    
    // start cutting the data up; plot the mass data for momenta of 8-10, 10-11, 11-12, 12-13, 13-15
    for(int i = 0; i < num_of_intervals; i++) {
        //Adaptive Cuts
        if (model_name == "Skew_Normal") {
            if (i == 0) {
                //func->SetParLimits(0, 1.0, 10.0); // Integral
                //func->SetParLimits(3, 0.0, 10); //Shape
                //func->SetParLimits(2, 0, 0.01); // Scale
            }
            else {
                //func->SetParLimits(0, 1.0, 100.0); // Integral
                //func->SetParLimits(3, -10, 10); //Shape
                //func->SetParLimits(2, 0, 100.0); // Scale
            }
        }
        //for(int i = 0; i < 1; i++) {
        pad[0] = new TPad("pad0","",0,0.42,1,1);
        pad[1] = new TPad("pad1","",0,0.02,1,0.42);
        double ptmin = intervals[i][0];
        double ptmax = intervals[i][1]; // Interval bounds
        
        // Plot the data
        SetCut(h_Pion, axis_pionPt, ptmin, ptmax);
        
        hMass = h_Pion->Projection(axis_pionMass);
        hMass->SetAxisRange(0.0, fit_y_max, "Y");
        canvas->cd();
        pad[0]->Draw();
        pad[0]->cd();
        hMass->Draw();
        
        // Find a fit just as you did for the entire data set
        // Graph the fit and (separately) the Gaussian component of it
        hMass->Fit(func);
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
        canvas->cd();
        residual->SetAxisRange(-4., 4., "Y");
        //residual->SetTitle("; Pion Mass (GeV); Residuals");
        pad[1]->Draw("p");
        pad[1]->cd();
        residual->Draw("p");
        
        // Add the mean mass parameter and its error to means and mean_errors, respectively; the standard deviation and
        // its error to sigmas and sigma_errors, the integral of the Gaussian peak and its error to gaussian_integrals
        // and integral_errors, the center point of the interval into center (and its error into widths), and chi square of the fit into its respective array
        if (model_name == "Gaussian") {
            means[i] = func->GetParameter(1);
            mean_errors[i] = func->GetParError(1);
            sigmas[i] = func->GetParameter(2) * 1000;
            sigma_errors[i] = func->GetParError(2) * 1000;
        }
        if (model_name == "Exponential_Gaussian") {
            double mu = func->GetParameter(1);
            double mu_error = func->GetParError(1);
            double lambda = func->GetParameter(3);
            double lambda_error = func->GetParError(3);
            double sigma = func->GetParameter(2);
            double sigma_error = func->GetParError(2);
            
            means[i] = mu + (1.0/lambda);
            mean_errors[i] = TMath::Sqrt(mu_error*mu_error + (lambda_error*lambda_error)/(lambda*lambda*lambda*lambda));
            sigmas[i] = 1000 * TMath::Sqrt((sigma*sigma) + (1.0/(lambda*lambda)));
            sigma_errors[i] = 1000 * 0.5 * TMath::Sqrt((4.0 * sigma_error*sigma_error*(sigma*sigma)) + (4.0*lambda_error*lambda_error/TMath::Power(lambda, 6.0)))/(sigmas[i]);
        }
        
        if (model_name == "Skew_Normal") {
            double delta = (func->GetParameter(3))/TMath::Sqrt(1 + (func->GetParameter(3))*(func->GetParameter(3)));
            double delta_error = (func->GetParError(3))/TMath::Power((1 + (func->GetParameter(3))*(func->GetParameter(3))), 1.5);
            double sigma_coeff = 1 - 2*delta*delta/TMath::Pi();
            double sigma_coeff_error = (4.0/TMath::Pi())*delta*delta_error;
            
            means[i] = (func->GetParameter(1)) + (func->GetParameter(2))*delta*TMath::Sqrt(2.0/TMath::Pi());
            mean_errors[i] = TMath::Sqrt((func->GetParError(1))*(func->GetParError(1)) + (2*(func->GetParameter(2))*(func->GetParameter(2))*delta*delta/TMath::Pi())*(((func->GetParError(2))/(func->GetParameter(2)))*((func->GetParError(2))/(func->GetParameter(2))) + (delta_error/delta)*(delta_error/delta)));
            sigmas[i] = ((func->GetParameter(2))*TMath::Sqrt(sigma_coeff))*1000;
            sigma_errors[i] = (TMath::Sqrt(4*(func->GetParError(2))*(func->GetParError(2))/((func->GetParameter(2))*(func->GetParameter(2))) + (sigma_coeff_error*sigma_coeff_error)/(sigma_coeff*sigma_coeff))/(2.0*TMath::Power(sigmas[i], 3.0)))*1000;
        }
        
        gaussian_integrals[i] = (func->GetParameter(0))/MASSWIDTH;
        integral_errors[i] = (func->GetParError(0))/MASSWIDTH;
        
        center[i] = (ptmin+ptmax)/2.0;
        widths[i] = ptmax - center[i];
        
        // print out the mean, standard deviation, and their corresponding errors; and save the graph in a new file
        std::cout << "Mean: " << means[i] << std::endl << "Standard Deviation: " << sigmas[i] << std::endl << "Standard Deviation Error: " << sigma_errors[i] << std::endl;
        canvas->SaveAs(Form(str_concat_converter(directory_name, "MyFit_Ptmin_%2.2f_Ptmax_%2.2f.png"), ptmin, ptmax));
        
        //Now use the signal_over_total function to get a graph of sigma vs. signal/total
        TF1* sig_over_tot_funct = new TF1("Signal over Total", signal_over_total, 0, 4, num_of_params);
        sig_over_tot_funct->SetParameters(func->GetParameters());
        canvas->Clear();
        
        sig_over_tot_funct->SetTitle(Form("Signal over Total vs Distance From Mean: %2.2f to %2.2f GeV; Num of Standard Deviations From Mean; Signal to Total Ratio", ptmin, ptmax));
        ///sig_over_tot_funct->Draw();
        TGraph* g_sig_over_tot = new TGraph(sig_over_tot_funct);
        g_sig_over_tot->SetLineColor(graph_colors[i]);
        //g_sig_over_tot->GetYaxis()->SetRangeUser(0.4, 1.0);
        peaks_over_totals->Add(g_sig_over_tot);
        //canvas->SaveAs(Form(str_concat_converter(directory_name, "Signal_Over_Total_Ptmin_%2.2f_Ptmax_%2.2f.png"), ptmin, ptmax));
        
        //Load onto the ROOT file
        hMass->GetListOfFunctions()->Add(peak);
        hMass->GetListOfFunctions()->Add(background);
        hMass->Write(Form("mass-pion-%2.2fGeV-%2.2fGeV", ptmin, ptmax));
        residual->Write(Form("residual-%2.2fGeV-%2.2fGeV", ptmin, ptmax));
    }
    
    // Graph the signal/total curves for each momentum increment
    // Black = 7.5-10, Red = 10-11, Blue = 11-12, Green = 12-13, Yellow = 13-15
    canvas->Clear();
    peaks_over_totals->SetTitle("Signal over Total vs Distance From Mean; Num of Standard Deviations From Mean; Signal to Total Ratio");
    //peaks_over_totals->GetYaxis()->SetRangeUser(0.4, 1.0);
    peaks_over_totals->SetMaximum(1.0);
    peaks_over_totals->SetMinimum(0.4);
    peaks_over_totals->Draw("Al");
    peaks_over_totals->Write("signal-over-total");
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
    g_mean->Write("mean-masses");
    canvas->SaveAs(str_concat_converter(directory_name, "meanMass_v_pT.png"));
    
    // Graph mass standard deviations with error bars
    canvas->Clear();
    TGraphErrors* g_sigma = new TGraphErrors(num_of_intervals, center, sigmas, widths, sigma_errors);
    g_sigma->Print();
    g_sigma->SetTitle("Mass Peak Widths for Various Momenta; Momentum (GeV); Mass (MeV/c^2)");
    g_sigma->GetYaxis()->SetRangeUser(7.5, 16.0);
    //g_sigma->SetMarkerSize(2);
    //g_sigma->SetMarkerStyle(20);
    g_sigma->Write("standard-dev-masses");
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
    g_chisquare->Write("chi-square");
    g_chisquare->Draw("AP");
    canvas->SaveAs(str_concat_converter(directory_name, "reduced_chisquare_v_pT.png"));
    
    // Graph the Gaussian distribution integrals over momentum
    canvas->Clear();
    TGraphErrors* g_integral = new TGraphErrors(num_of_intervals, center, gaussian_integrals, widths, integral_errors);
    g_integral->Print();
    g_integral->SetTitle("Gaussian Peak integrals for Various Momenta; Momentum (GeV); Number of Pions");
    //g_integral->GetYaxis()->SetRangeUser(0.0, 6000.0);
    g_integral->Draw("AP");
    g_integral->Write("pion-integrals");
    canvas->SaveAs(str_concat_converter(directory_name, "peakIntegrals_v_pT.png"));
    
    canvas->Close();
    return;
}
