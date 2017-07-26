// This macro graphs photon momentum/energy vs isolation and photon momentum/energy vs. UE
// Author: Ivan Chernyshev; Date: 7/24/2017

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

const int axis_photonPt                  =  3;
const int axis_photonPseudoRapidity      =  4;
const int axis_photonLambda              =  6;
const int axis_photonNcells              =  7;
const int axis_photonDisToBorder         =  9;
const int axis_photonDisToBadCell        = 10;
const int axis_photonDisToCharged        = 11;
const int axis_photonExoticity           = 14;
const int axis_photonTime                = 15;
const int axis_photonIsolation_R_track   = 18;
const int axis_photonIsolation_R_cluster = 22;
const int axis_photonUE_track            = 20;
const int axis_photonUE_cluster          = 24;

// The cutting function
void SetCut(THnSparse* h, const int axis, double min, double max){
    //make a selection on the chosen variable
    double width = h->GetAxis(axis)->GetBinWidth(1);
    int binmin = h->GetAxis(axis)->FindBin(min);
    int binmax = h->GetAxis(axis)->FindBin(max);
    h->GetAxis(axis)->SetRange(binmin, binmax - 1);
    return;
}

void momentum_correlator() {
    // Set ATLAS style
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    TCanvas* canvas = new TCanvas();
    
    // Load isolation vs. momentum data from THnSparses_072417.root
    TFile* fIn = new TFile("THnSparses_072417.root", "READ");
    THnSparse* hPhoton = 0;
    fIn->GetObject("h_Cluster", hPhoton);
    
    // Cut the data
    SetCut(hPhoton, axis_photonPt, 6, 50);
    SetCut(hPhoton, axis_photonLambda, 0.1, 0.4);
    SetCut(hPhoton, axis_photonDisToCharged, 0.02, 0.15);
    SetCut(hPhoton, axis_photonDisToBorder, 0.9, 6);
    SetCut(hPhoton, axis_photonDisToBadCell, 1.9, 10);
    SetCut(hPhoton, axis_photonNcells, 1.9, 30.0);
    SetCut(hPhoton, axis_photonExoticity, 0, 0.97);
    SetCut(hPhoton, axis_photonTime, -30, 30);
    SetCut(hPhoton, axis_photonPseudoRapidity, -0.27, 0.27);
    
    // Graph the data
    TH2D* pt_isolation_hist = hPhoton->Projection(axis_photonIsolation_R_track, axis_photonPt);
    pt_isolation_hist->SetTitle("Isolation-energy correlation; Photon energy (GeV); #Sigma E^{cone}_{T} [GeV], track");
    gPad->SetLogz(kTRUE);
    pt_isolation_hist->Draw("COLZ");
    myText(.20, .92, kBlack, "Isolation-energy correlation");
    canvas->SaveAs("isolation_energy_correlation.png");
    
    // Repeat for UE vs. momentum for tracks, but scale each bin by 1/0.937.
    double scale_factor = 1/0.937;
    SetCut(hPhoton, axis_photonPt, 6, 25);
    SetCut(hPhoton, axis_photonUE_track, 0, 25);
    canvas->Clear();
    TH2D* pt_UE_hist_track = hPhoton->Projection(axis_photonUE_track, axis_photonPt);
    TAxis* yaxis = pt_UE_hist_track->GetYaxis();
    yaxis->Set(yaxis->GetNbins(), scale_factor*yaxis->GetXmin(), scale_factor*yaxis->GetXmax());
    /**int xbinmin = pt_UE_hist_track->GetXaxis()->FindBin(6);
    int xbinmax = pt_UE_hist_track->GetXaxis()->FindBin(25);
    int ybinmin = pt_UE_hist_track->GetYaxis()->FindBin(0.0);
    int ybinmax = pt_UE_hist_track->GetYaxis()->FindBin(25);
    for (int i = xbinmin; i <= xbinmax; i++)
        for (int j = ybinmin; j <= ybinmax; j++) {
            pt_UE_hist_track->SetBinContent(i, j, (pt_UE_hist->GetBinContent(i, j))/0.937);
        }*/
    pt_UE_hist_track->SetTitle("Track UE density-photon energy correlation; Photon energy (GeV); Track UE Density (GeV)");
    pt_UE_hist_track->Draw("COLZ");
    myText(.20, .92, kBlack, "Track UE density-photon energy correlation");
    canvas->SaveAs("Track_UE_energy_density_correlation.png");
    
    // Now, take a projection of the mean UE in every photon energy bin
    canvas->Clear();
    TH1D* pt_UEmeans_hist_track = new TH1D("mean_UE_histogram", "Mean track UE densities over photon energy; Photon energy (GeV); Mean Track UE Density (GeV)", pt_UE_hist_track->GetXaxis()->GetNbins(), pt_UE_hist_track->GetXaxis()->GetXmin(), pt_UE_hist_track->GetXaxis()->GetXmax());
    for (int i = 0; i < pt_UEmeans_hist_track->GetSize(); i++) {
        // Create another histogram, for storing UE-entry data for each photon energy bin so that its mean could be put into the mean histogram
        TH1D* ptbin_UE_hist = pt_UE_hist_track->ProjectionY("transitional_storage", i, i + 1);
        //for (int j = 0; j < ptbin_UE_hist->GetSize(); j++)
            //ptbin_UE_hist->SetBinContent(j, pt_UEmeans_hist_track->GetBinContent(i, j));
        pt_UEmeans_hist_track->SetBinContent(i, ptbin_UE_hist->GetMean());
        std::cout << "Mean UE for bin " << i << " is " << ptbin_UE_hist->GetMean() << std::endl;
        
    }
    pt_UEmeans_hist_track->Draw();
    myText(.20, .92, kBlack, "Mean track UE densities over photon energy");
    canvas->SaveAs("Mean_track_UEs_over_photonE.png");
    
    // Repeat for UE vs. momentum for clusters, but scale UE by 1/0.57.
    scale_factor = 1/0.57;
    SetCut(hPhoton, axis_photonUE_cluster, 0, 25);
    canvas->Clear();
    TH2D* pt_UE_hist_cluster = hPhoton->Projection(axis_photonUE_cluster, axis_photonPt);
    yaxis = pt_UE_hist_cluster->GetYaxis();
    yaxis->Set(yaxis->GetNbins(), scale_factor*yaxis->GetXmin(), scale_factor*yaxis->GetXmax());
    /**xbinmin = pt_UE_hist_cluster->GetXaxis()->FindBin(6);
    xbinmax = pt_UE_hist_cluster->GetXaxis()->FindBin(25);
    ybinmin = pt_UE_hist_cluster->GetYaxis()->FindBin(0.0);
    ybinmax = pt_UE_hist_cluster->GetYaxis()->FindBin(25);
    for (int i = xbinmin; i <= xbinmax; i++)
        for (int j = ybinmin; j <= ybinmax; j++) {
            pt_UE_hist_cluster->SetBinContent(i, j, (pt_UE_hist_cluster->GetBinContent(i, j))/0.57);
        }*/
    pt_UE_hist_cluster->SetTitle("Cluster UE-energy correlation; Photon energy (GeV); Cluster UE Density (GeV)");
    pt_UE_hist_cluster->Draw("COLZ");
    myText(.20, .92, kBlack, "Cluster UE density-photon energy correlation");
    canvas->SaveAs("Cluster_UE_energy_density_correlation.png");
    
    // Now, take a projection of the mean UE in every photon energy bin
    canvas->Clear();
    TH1D* pt_UEmeans_hist_cluster = new TH1D("mean_UE_histogram", "Mean cluster UE densities over photon energy; Photon energy (GeV); Mean Cluster UE Density (GeV)", pt_UE_hist_cluster->GetXaxis()->GetNbins(), pt_UE_hist_cluster->GetXaxis()->GetXmin(), pt_UE_hist_cluster->GetXaxis()->GetXmax());
    for (int i = 0; i < pt_UEmeans_hist_cluster->GetSize(); i++) {
        // Create another histogram, for storing UE-entry data for each photon energy bin so that its mean could be put into the mean histogram
        TH1D* ptbin_UE_hist = pt_UE_hist_cluster->ProjectionY("transitional_storage", i, i + 1);
        //TH1D* ptbin_UE_hist = new TH1D("transitional_storage", "transitional_storage", pt_UE_hist_cluster->GetYaxis()->GetNbins(), pt_UE_hist_cluster->GetYaxis()->GetXmin(), pt_UE_hist_cluster->GetYaxis()->GetXmin());
        //for (int j = 0; j < ptbin_UE_hist->GetSize(); j++)
            //ptbin_UE_hist->SetBinContent(j, pt_UEmeans_hist_cluster->GetBinContent(i, j));
        pt_UEmeans_hist_cluster->SetBinContent(i, ptbin_UE_hist->GetMean());
        std::cout << "Mean UE for bin " << i << " is " << ptbin_UE_hist->GetMean() << std::endl;
    }
    pt_UEmeans_hist_cluster->Draw();
    myText(.20, .92, kBlack, "Mean cluster UE densities over photon energy");
    canvas->SaveAs("Mean_cluster_UEs_over_photonE.png");

    
    canvas->Close();
}
