#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TGraph.h>
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
// Ae^(x-mean)^2/(2*sigma^2) + B*x^4 + C*x^3 + D*x^2 + E*x + F
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
    double fitval = A*TMath::Exp(-0.5*arg*arg) + B*TMath::Power(x[0], 4.0) + C*TMath::Power(x[0], 3.0) + D*x[0]*x[0] + E*x[0] + F;
    // +
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
    TCanvas* canvas = new TCanvas();
    hMass->Draw();
    TF1 *func = new TF1("fit", compound_model,0.05,0.5,8);
    
    //Restrict the parameters to reasonable ranges, insert guess values, and give understandable names
    func->SetParameters(6000,  0.14, 0.3,  1, 0.03, 0.6, 0.1, 1);
    func->SetParLimits(0, 0.0, 10000.0);
    func->SetParLimits(1, 0.1, 0.2);
    func->SetParLimits(2, 0.0001, 0.1); // widths
    func->SetParLimits(3, -10000.0, 0.0);
    func->SetParLimits(5, -100000.0, 0.0);
    //func->SetParLimits(7, 0.0, 0.1);
    //func->SetParLimits(3, 100.0, 5000.0);
    //func->SetParLimits(6, 100.0, 5000.0);
    //func->SetParNames("Amplitude", "Mean", "Sigma", "Damped Sine coeff", "Shift", "Period Factor", "Damping Exponent", "Constant");
    func->SetParNames("Amplitude", "Mean", "Sigma", "Quadric coeff", "Cubic coeff", "Quadratic coeff", "Linear coeff", "Constant");
    //func->SetParNames("Amplitude", "Mean", "Sigma", "Logistic asymptote", "e-coeff", "In-exponent coeff", "Shift", "Constant");
    
    // Plot the fit for the mass, save as a PDF
    hMass->Fit(func);
    func->Draw("same");
    canvas->SaveAs("mass_pion_plot.png");
    
    // Plot the data for the momentum
    TH1D* hPt = h_Pion->Projection(axis_pionPt);
    hPt->Draw();
    canvas->SaveAs("momentum_pion_plot.png");
    
    // start cutting the data up; plot the mass data for momenta of 5-10, 10-15, 15-20, and so forth
    const int num_of_intervals = 4;
    int intervals[num_of_intervals][2] = {{5, 10}, {10, 12}, {12, 15}, {15, 20}};
    double means[num_of_intervals];
    double sigmas[num_of_intervals];
    double center[num_of_intervals];

    for(int i = 0; i < num_of_intervals; i++) {
        double ptmin = intervals[i][0];
        double ptmax = intervals[i][1]; // Interval bounds
        
        // Plot the data
        SetCut(h_Pion, axis_pionPt, ptmin, ptmax);
        hMass = h_Pion->Projection(axis_pionMass);
        hMass->Draw();
        
        // Find a fit just as above
        hMass->Fit(func);
        func->Draw("same");
        
        means[i] = func->GetParameter(1);
        sigmas[i] = func->GetParameter(2);
        center[i] = (ptmin+ptmax)/2.0;
        std::cout << "Mean: " << means[i] << std::endl << "Standard Deviation: " << sigmas[i] << std::endl;
        canvas->SaveAs(Form("MyFit_Ptmin_%2.f_Ptmax_%2.f.png", ptmin, ptmax));
    }
    
    // Graph means
    TGraph* g = new TGraph(num_of_intervals, center, means);
    canvas->Clear();
    g->Print();
    g->SetMarkerSize(2);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    canvas->SaveAs("meanMass_v_pT.png");
    
    // Graph standard deviations
    TGraph* g_sigma = new TGraph(num_of_intervals, center, sigmas);
    canvas->Clear();
    g_sigma->Print();
    g_sigma->SetMarkerSize(2);
    g_sigma->SetMarkerStyle(20);
    g_sigma->Draw("AP");
    canvas->SaveAs("massWidths_v_pT.png");

    //plot the mass data for momenta of 15-20
    //SetCut(h_Pion, axis_pionPt, 15.0, 20.0);
    //hMass = h_Pion->Projection(axis_pionMass);
    //hMass->Draw();
    
    // Find a fit just as above
    //hMass->Fit(func);
    //func->Draw("same");
    //canvas->SaveAs("mass_pion_plot_15-20.png");

    canvas->Close();
}
