// This file exists for the sole purpose of removing the count tags from the energy response charts
// Author: Ivan Chernyshev; Date: 3/6/2018

#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

#include "TH1F.h"

void twoD_energy_response_plotter() {
// Open the NTuple output files
TFile* energyresponsev1 = new TFile("fout_energyresponse_17g6a1_pthat2_clusterv1_small.root", "READ");
TFile* energyresponsev2 = new TFile("fout_energyresponse_17g6a1_pthat2_clusterv2_small.root", "READ");

// Get Atlas Style and define the Canvas
TCanvas* c = new TCanvas();
gROOT->LoadMacro("AtlasStyle.C");
SetAtlasStyle();

// Get the files for each source
TH2D* ptresponse_v1 = 0;
energyresponsev1->GetObject("energy_response", ptresponse_v1);
TH2D* ptresponse_v2 = 0;
energyresponsev2->GetObject("energy_response", ptresponse_v2);

// Plot
ptresponse_v1->SetTitle("p_{T} response; True p_{T} (GeV); Measured p_{T} (GeV)");
ptresponse_v1->GetXaxis()->SetTitleOffset(1.2);
ptresponse_v1->GetYaxis()->SetTitleOffset(1.2);
ptresponse_v1->Draw("COLZ");
myText(0.4, 0.92, 1, "p_{T} response");
c->SaveAs("energy_response_scatter__17g6a1_pthat2_clusterv1_small.png");
c->Clear();
    
ptresponse_v2->SetTitle("p_{T} response; True p_{T} (GeV); Measured p_{T} (GeV)");
ptresponse_v2->GetXaxis()->SetTitleOffset(1.2);
ptresponse_v2->GetYaxis()->SetTitleOffset(1.2);
ptresponse_v2->Draw("COLZ");
myText(0.4, 0.92, 1, "p_{T} response");
c->SaveAs("energy_response_scatter__17g6a1_pthat2_clusterv2_small.png");
c->Clear();


c->Close();
}
