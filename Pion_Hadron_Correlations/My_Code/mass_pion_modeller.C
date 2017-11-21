// This macro creates a mass-pion model and outputs resultant parameters of it per pT bin in the form of a .root file
// Is based on my_code from much of the Pion_data folder
// Programmer: Ivan Chernyshev; Date: 10/25/17

#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TGraphErrors.h>
#include <TCanvas.h>

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
    Main function
*/
void mass_pion_modeller() {
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    TH1::SetDefaultSumw2();
    
    //Open the files
    TFile* fIn = new TFile("THnSparses_LHC13d_101517.root","READ"); //get file
    TFile* fOut = new TFile(Form("PionDataOutput.root"), "RECREATE"); // Create an output file
    fIn->Print(); //print file content
    
    // Get the data
    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion); //get array
    
    //Create the TCanvas and TPads
    TCanvas* graphcanvas = new TCanvas();
    //TPad *pad[2] = {new TPad("pad0","",0,0.36,1,1), new TPad("pad1","",0,0.1,1,0.45)};
    
    // Apply cuts
    SetCut(h_Pion, axis_pionLambda1, 0.1, 0.4);
    SetCut(h_Pion, axis_pionLambda2, 0.1, 0.4);
    SetCut(h_Pion, axis_pionAngle, 17, 60);
    
    // Create the fit-functions: total, background, and peak, create initial parameters, restrict functions to certain values
    const int num_of_params = 6;
    const int num_of_peak_params = 3;
    TF1* background = new TF1("background curve", "[0]*x*x + [1]*x + [2]", 0.08, 0.26);
    TF1* func = new TF1("fit", gaussian_model,0.05,0.5,num_of_params);
    TF1* peak = new TF1("mass peak", gaussian_peak, 0.05, 0.5, num_of_peak_params);
    func->SetParNames("Integral", "Mean", "Sigma", "Quadratic coeff", "Linear coeff", "Constant");
    func->SetParameters(40,  0.13, 0.012, -120000, 28000, -1600);
    func->SetParLimits(0, 0.001, 400.0);//integral
    func->SetParLimits(1, 0.11, 0.16); //mean
    func->SetParLimits(2, 0.008, 0.016); // standard deviation
    func->SetParLimits(3, -1000000.0, 0.0); // Quadratic factor
    
    // Prepare to evaluate the data per pT interval: define the specific intervals, and arrays for all of the quantities that are later to be graphed against pT and stored in the .root file
        // Which pT intervals there are (they're 2 GeV wide, and range from 6 GeV to 16 GeV
    const int num_of_intervals = 6;
    double intervals[num_of_intervals][2] = {{6.0, 8.0}, {8.0, 10.0}, {10.0, 12.0}, {12.0, 14.0}, {14.0, 16.0}, {16.0, 20.0}};
    
        // Inependent variables, that represent pT when anything is graphed against pT, as well as its error
    double center[num_of_intervals];
    double widths[num_of_intervals];
    
        // Dependent variables, the variables to be graphed against pT, and their errors
    double means[num_of_intervals];
    double mean_errors[num_of_intervals];
    double sigmas[num_of_intervals];
    double sigma_errors[num_of_intervals];
    double gaussian_integrals[num_of_intervals];
    double integral_errors[num_of_intervals];
    
    // Model the data for each pT interval
    for(int i = 0; i < num_of_intervals; i++) {
        
        //pad[0] = new TPad("pad0","",0,0.36,1,1);
        //pad[1] = new TPad("pad1","",0,0.05,1,0.45);
        double min = intervals[i][0];
        double max = intervals[i][1]; // Interval bounds
        
        // pT-dependent cuts
        if (min >= 12)
            SetCut(h_Pion, axis_asymmetry, 0.0, 0.9);
        else
            SetCut(h_Pion, axis_asymmetry, 0.0, 0.6);
        SetCut(h_Pion, axis_pionPt, min, max);
        
        // Declare the variables that will store: the mass-pion data and its fit (hMass), the residual of each point in hMass, and the width of each mass bin in hMass
        // This data will later be used in creating the output for the .root file
        TH1D* hMass = h_Pion->Projection(axis_pionMass);
        TH1D* residual = (TH1D*)hMass->Clone("residual");
        double MASSWIDTH = hMass->GetBinWidth(1);
        
        // Graph the raw data
        //graphcanvas->cd();
        //pad[0]->Draw();
        //pad[0]->cd();
        hMass->GetYaxis()->SetTitle("Number of Entries");
        hMass->GetYaxis()->SetTitleOffset(0.8);
        hMass->Draw();
        
        // Find the fits, graph them on the same graph as the raw data, write to the .root file
        hMass->Fit(func);
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
        hMass->Write(Form("mass-pion-%2.2fGeV-%2.2fGeV", min, max));
        graphcanvas->Clear();
        
        // Now, fill residual with the residuals of each point
        for (int i = 0; i < hMass->GetSize(); i++) {
            
            double residualvalue = (hMass->GetBinContent(i) - func->Eval(hMass->GetBinCenter(i)))/hMass->GetBinError(i);
            if ((hMass->GetBinError(i)) != 0) {
                residual->SetBinContent(i, residualvalue);
            }
            else
                residual->SetBinContent(i, 0);
            residual->SetBinError(i, 0); //Residuals don't have errors
            
        }
        
        // Graph residual and write it to the root file
        residual->SetAxisRange(-4., 4., "Y");
        residual->GetXaxis()->SetTitleSize(.08);
        residual->GetYaxis()->SetTitleSize(.1);
        residual->GetYaxis()->SetTitleOffset(0.3);
        residual->GetXaxis()->SetTitleOffset(0.8);
        residual->SetTitle("; Pion Mass (GeV); Residuals");
        residual->Draw("p");
        residual->Write(Form("residual-%2.2fGeV-%2.2fGeV", min, max));
        graphcanvas->Clear();
        
        // Fill the independent and dependent variable array elements corresponding to the pT bin currently being studied
        means[i] = func->GetParameter(1) * 1000;
        mean_errors[i] = func->GetParError(1) * 1000;
        sigmas[i] = func->GetParameter(2) * 1000;
        sigma_errors[i] = func->GetParError(2) * 1000;
        
        gaussian_integrals[i] = (func->GetParameter(0))/MASSWIDTH;
        integral_errors[i] = (func->GetParError(0))/MASSWIDTH;
        
        center[i] = (min+max)/2.0;
        widths[i] = max - center[i];
        
        graphcanvas->Clear();
    }
    
    // Now, output the pT-dependent graphs
    
    // Mean mass data
        // Expected value
    TF1* mass_pdg = new TF1("mass_pdg", "[0]", 0, 20);
    mass_pdg->SetParameter(0, 0.13498);
    mass_pdg->SetLineWidth(2);
    mass_pdg->SetLineColor(kRed);
        // Mean masse
    TGraphErrors* g_mean = new TGraphErrors(num_of_intervals, center, means, widths, mean_errors);
    graphcanvas->Clear();
    g_mean->Print();
    g_mean->SetTitle("Mean Masses for Various Momenta; Momentum (GeV); Mass (MeV)");
    g_mean->Draw("AP");
    mass_pdg->Draw("same");
    g_mean->Write("mean-masses");
    graphcanvas->Clear();
    
    //Mass standard deviation data
    TGraphErrors* g_sigma = new TGraphErrors(num_of_intervals, center, sigmas, widths, sigma_errors);
    g_sigma->Print();
    g_sigma->SetTitle("Mass Peak Widths for Various Momenta; Momentum (GeV); Mass Width (MeV)");
    g_sigma->GetYaxis()->SetTitleOffset(.7);
    g_sigma->GetXaxis()->SetTitleOffset(.9);
    g_sigma->GetYaxis()->SetRangeUser(3.0, 15.7);
    g_sigma->GetXaxis()->SetRangeUser(6.0, 16.0);
    g_sigma->Draw("AP");
    g_sigma->Write("standard-dev-masses");
    graphcanvas->Clear();

    // Graph the distribution integrals over momentum
    TGraphErrors* g_integral = new TGraphErrors(num_of_intervals, center, gaussian_integrals, widths, integral_errors);
    g_integral->Print();
    g_integral->SetTitle("Peak integrals for Various Momenta; Momentum (GeV); Number of Pions");
    g_integral->GetXaxis()->SetRangeUser(6.0, 16.0);
    //g_integral->GetYaxis()->SetRangeUser(0.0, 6000.0);
    g_integral->Draw("AP");
    g_integral->Write("pion-integrals");
    
    graphcanvas->Close();
}
