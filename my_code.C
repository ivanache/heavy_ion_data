#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TGraphErrors.h>
#include <iostream>

//variables of hPion
const int axis_pion_Cen  = 0;
const int axis_pion_Zvtx = 1;
const int axis_pionMass  = 2;
const int axis_pionPt    = 3;

// The cutting function
void SetCut(THnSparse* h, const int axis, double min, double max){
    //make a selection on the chosen variable
    double width = h->GetAxis(axis)->GetBinWidth(1);
    int binmin = h->GetAxis(axis)->FindBin(min);
    int binmax = h->GetAxis(axis)->FindBin(max);
    h->GetAxis(axis)->SetRange(binmin, binmax);
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

void my_code(){
    //Open the file
    TFile* fIn = new TFile("Ntuple.root","READ"); //get file
    fIn->Print(); //print file content
    
    // Get the data
    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion); //get array
    
    //For the mass plot, restrict to mass between 0.08 and 0.25, plot data, and set up the fit function
    SetCut(h_Pion, axis_pionMass, 0.08, 0.25);
    TH1D* hMass = h_Pion->Projection(axis_pionMass);
    const double MASSWIDTH = hMass->GetBinWidth(1);
    TCanvas* canvas = new TCanvas();
    hMass->Draw();
    TF1 *func = new TF1("fit", compound_model,0.05,0.5,8);
    
    //Restrict the parameters to reasonable ranges, insert guess values, and give understandable names
    func->SetParameters(600,  0.14, 0.3,  1, 0.03, 0.6, 0.1, 1);
    func->SetParLimits(0, 1, 10000.0);//integral
    func->SetParLimits(1, 0.1, 0.2); //mean
    func->SetParLimits(2, 0.01, 0.05); // width
    func->SetParLimits(3, -10000.0, 0.0); // Quadric and quadratic factors
    func->SetParLimits(5, -100000.0, 0.0);
    //func->SetParLimits(7, 0.0, 0.1);
    //func->SetParLimits(3, 100.0, 5000.0);
    //func->SetParLimits(6, 100.0, 5000.0);
    //func->SetParNames("Amplitude", "Mean", "Sigma", "Damped Sine coeff", "Shift", "Period Factor", "Damping Exponent", "Constant");
    func->SetParNames("Integral", "Mean", "Sigma", "Quadric coeff", "Cubic coeff", "Quadratic coeff", "Linear coeff", "Constant");
    //func->SetParNames("Amplitude", "Mean", "Sigma", "Logistic asymptote", "e-coeff", "In-exponent coeff", "Shift", "Constant");
    
    // Plot the fit for the mass and (separately) the Gaussian component of it; save as a PDF
    hMass->Fit(func);
    func->Draw("same");
    TF1* peak = new TF1("mass peak", "[0]*(1.0/([2]*TMath::Sqrt(2*TMath::Pi())))*TMath::Exp(-0.5*((x - [1])/[2])*((x - [1])/[2]))", 0.08, 0.26);
    for(int i = 0; i < 3; i++) {
         peak->SetParameter(i, func->GetParameter(i));
    }
    peak->SetLineColor(kBlue);
    peak->Draw("same");
    canvas->SaveAs("mass_pion_plot.png");
    
    // Plot the data for the momentum
    TH1D* hPt = h_Pion->Projection(axis_pionPt);
    hPt->Draw();
    canvas->SaveAs("momentum_pion_plot.png");
    
    // A collection of variables that is needed for the next steps
    const int num_of_intervals = 7;
    double intervals[num_of_intervals][2] = {{5.0, 7.5}, {7.5, 10.0}, {10.0, 11.0}, {11.0, 12.0}, {12.0, 13.0}, {13.0, 15.0}, {15.0, 18.0}};
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
        hMass->Draw();
        
        // Find a fit just as you did for the entire data set
        // Graph the fit and (separately) the Gaussian component of it
        hMass->Fit(func);
        func->Draw("same");
        for(int i = 0; i < 3; i++) {
            peak->SetParameter(i, func->GetParameter(i));
        }
        peak->SetLineColor(kBlue);
        peak->Draw("same");
        
        // Add the mean mass parameter and its error to means and mean_errors, respectively; the standard deviation and
        // its error to sigmas and sigma_errors, the integral of the Gaussian peak and its error to gaussian_integrals
        // and integral_errors, and the center point of the interval into center (and its error into widths), print out
        // the mean, standard deviation, and their corresponding errors, and save the graph in a new file
        means[i] = func->GetParameter(1);
        mean_errors[i] = func->GetParError(1);
        sigmas[i] = func->GetParameter(2) * 1000;
        sigma_errors[i] = func->GetParError(2) * 1000;
        
        center[i] = (ptmin+ptmax)/2.0;
        widths[i] = ptmax - center[i];
        
        gaussian_integrals[i] = (func->GetParameter(0))/MASSWIDTH;
        integral_errors[i] = (func->GetParError(0))/MASSWIDTH;
        
        std::cout << "Mean: " << means[i] << std::endl << "Standard Deviation: " << sigmas[i] << std::endl;
        canvas->SaveAs(Form("MyFit_Ptmin_%2.2f_Ptmax_%2.2f.png", ptmin, ptmax));
    }
    
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
    //g_mean->SetMarkerSize(2);
    //g_mean->SetMarkerStyle(20);
    g_mean->Draw("AP");
    mass_pdg->Draw("same");
    canvas->SaveAs("meanMass_v_pT.png");
    
    // Graph mass standard deviations with error bars
    canvas->Clear();
    TGraphErrors* g_sigma = new TGraphErrors(num_of_intervals, center, sigmas, widths, sigma_errors);
    g_sigma->Print();
    g_sigma->SetTitle("Mass Peak Widths for Various Momenta; Momentum (GeV); Mass (MeV/c^2)");
    //g_sigma->SetMarkerSize(2);
    //g_sigma->SetMarkerStyle(20);
    g_sigma->Draw("AP");
    canvas->SaveAs("massWidths_v_pT.png");
    
    // Graph the Gaussian distribution integrals over momentum
    canvas->Clear();
    TGraphErrors* g_integral = new TGraphErrors(num_of_intervals, center, gaussian_integrals, widths, integral_errors);
    g_integral->Print();
    g_integral->SetTitle("Gaussian Peak integrals for Various Momenta; Momentum (GeV); Gaussian Peak Integral");
    g_integral->Draw("AP");
    canvas->SaveAs("peakIntegrals_v_pT.png");


    canvas->Close();
}
