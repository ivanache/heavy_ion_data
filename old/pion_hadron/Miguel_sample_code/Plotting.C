// This macro is a template of all future pion-hardon correlations
// Author: Miguel Arratia
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
const int axis_pion_PhPt_1      = 6;
const int axis_pion_PhPt_2      = 7;
const int axis_pionOpeningAngle = 8;
const int axis_pion_PhM02_1     = 9;
const int axis_pion_PhM02_2     = 10;
const int axis_pion_PhdR_1      = 11;
const int axis_pion_PhdR_2      = 12;

//variables common to hPionTrack and hpionCluster
const int axis_corr_centrality  = 0;
const int axis_corr_zvertex     = 1;
const int axis_corr_triggerpT   = 2;
const int axis_corr_trackpT     = 3;
const int axis_corr_dphi        = 4;
const int axis_corr_deta        = 5;
const int axis_corr_zt          = 6;
const int axis_corr_xi          = 7;
//non-shared, (Pion)
const int axis_corr_mass        = 8;
const int axis_corr_photon1Pt   = 9;
const int axis_corr_photon2Pt   = 10;
const int axis_corr_photon1M02  = 11;
const int axis_corr_photon2M02  = 12;
//non-shared (Cluster)
const int axis_corr_M02     = 8;
const int axis_corr_dR      = 9;

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


//fit function
double fitf(Double_t *x,Double_t *par) {
    double arg1 = 0;
    if (par[1]!=0) arg1 = (x[0] - 0.0)/par[1];
    double arg2 = 0;
    if (par[3]!=0) arg2 = (x[0] - 1.0)/par[3];
    
    double normfactor_1 = 1.0/(sqrt(2*TMath::Pi())*par[1]);
    double normfactor_2 = 1.0/(sqrt(2*TMath::Pi())*par[3]);
    
    double fitval = par[0]*normfactor_1*TMath::Exp(-0.5*arg1*arg1)+par[2]*normfactor_2*TMath::Exp(-0.5*arg2*arg2)+par[4];
    return fitval;
}

double GaussP0(Double_t *x, Double_t *par){
    double arg1 = 0;
    if (par[1]!=0) arg1 = (x[0] - 0.0)/par[1];
    double fitval = par[0]*TMath::Exp(-0.5*arg1*arg1) + par[2];
    return fitval;
}

void SetCut(THnSparse* h, const int axis, double min, double max){
    
    int binmin = h->GetAxis(axis)->FindBin(min);
    int binmax = h->GetAxis(axis)->FindBin(max);
    h->GetAxis(axis)->SetRange(binmin, binmax-1);
    return;
}

void Guardar(TCanvas* c, TString name){
    
    c->SaveAs("PDFOUTPUT/"+name+".pdf");
    c->SaveAs("PNGOUTPUT/"+name+".png");
    
}


void Looping(THnSparse* h, THnSparse* h_mixed, const int axisNumber, std::vector<double> bins, const char* name, TMultiGraph& multi)
{
    
    TH2F* h_dphi_deta = 0x0;
    TH2F* h_dphi_deta_Mixed = 0x0;
    TH1F* h_dphi_near = 0x0;
    TH1F* h_dphi_far = 0x0;
    TH1F* h_trackpt = 0x0;
    TH1F* h_deta = 0x0;
    TH1F* h_deta_Mixed = 0x0;
    
    auto c = new TCanvas("c","c",800,800);
    
    std::vector<double> Yield_near;
    std::vector<double> Yield_near_err;
    std::vector<double> Yield_far;
    std::vector<double> Yield_far_err;
    
    for(int n=0; n<bins.size()-1; n++){    //loop over bins of hadron pT/zt or xi
        
        double min = bins.at(n);
        double max = bins.at(n+1);
        
        const char* tag =  Form("%s_%1.f_%1.f", name, 10*min, 10*max);
        
        SetCut(h, axisNumber, min, max);
        SetCut(h_mixed, axisNumber, min, max);
        
        ///////////////////////////////
        c->cd();
        h_trackpt = (TH1F*)h->Projection(axisNumber);
        h_trackpt->Draw("hist");
        //Guardar(c, Form("%s_testing_%s", name, tag ));
        c->Clear();
        
        h_dphi_deta = (TH2F*)h->Projection(axis_corr_dphi,  axis_corr_deta);
        h_dphi_deta->Sumw2();
        
        auto c2 = new TCanvas("c2","c2",1000,600);
        c2->Divide(2);
        c2->cd(1);
        gPad->SetPhi(-55);
        gPad->SetTheta(+30);
        
        c2->cd(1);
        h_dphi_deta->Draw("SURF1 FB BB");
        h_dphi_deta->SetLineWidth(1.0);
        myText(0.4,0.82,kBlack, "Same Events");
        h_dphi_deta->GetXaxis()->SetNdivisions(1);
        h_dphi_deta->GetYaxis()->SetNdivisions(4);
        h_dphi_deta->GetZaxis()->SetNdivisions(1);
        
        h_dphi_deta->GetXaxis()->CenterTitle(kTRUE);
        h_dphi_deta->GetYaxis()->CenterTitle(kTRUE);
        h_dphi_deta->SetTitle("; |#Delta#eta|; #Delta#phi/#pi");
        
        h_dphi_deta_Mixed = (TH2F*)h_mixed->Projection(axis_corr_dphi,  axis_corr_deta);
        h_dphi_deta_Mixed->Sumw2();
        
        c->cd();
        h_deta =  (TH1F*)h_dphi_deta_Mixed->ProjectionX("h_deta",0.0,200.0);
        
        h_deta->Draw();
        double scale = h_deta->GetMaximum()/h_dphi_deta_Mixed->GetNbinsY();
        std::cout << "SCALE " << h_deta->GetMaximum() << " scale:" << scale << std::endl;
        Guardar(c, Form("%s_Deta_%s", name, tag ));
        c->Clear();
        
        h_dphi_deta_Mixed->Scale(1.0/scale);
        
        c2->cd(2);
        gPad->SetPhi(-55);
        gPad->SetTheta(+20);
        h_dphi_deta_Mixed->SetTitle("Mixed");
        h_dphi_deta_Mixed->Draw("SURF1 FB BB");
        h_dphi_deta_Mixed->SetLineWidth(1.0);
        //h_dphi_deta_Mixed->SetLineColor(kGray);
        myText(0.4,0.82,kBlack, "Mixed Events");
        h_dphi_deta_Mixed->GetXaxis()->SetNdivisions(1);
        h_dphi_deta_Mixed->GetYaxis()->SetNdivisions(4);
        h_dphi_deta_Mixed->GetZaxis()->SetNdivisions(1);
        h_dphi_deta_Mixed->GetZaxis()->SetRangeUser(0,1.0);
        h_dphi_deta_Mixed->SetMinimum(0);
        h_dphi_deta_Mixed->GetXaxis()->CenterTitle(kTRUE);
        h_dphi_deta_Mixed->GetYaxis()->CenterTitle(kTRUE);
        h_dphi_deta_Mixed->SetTitle("; |#Delta#eta|; #Delta#phi/#pi");
        
        Guardar(c2, Form("%s_2D_%s", name, tag ));
        c2->Clear();
        c2->Close();
        
        c->cd();
        h_deta =  (TH1F*)h_dphi_deta_Mixed->ProjectionX("h_deta",0,200);
        h_deta->Scale(1.0/h_dphi_deta_Mixed->GetNbinsY());
        h_deta->Draw();
        c->Clear();
        
        //Divide
        c->cd();
        h_dphi_deta->Divide(h_dphi_deta_Mixed);
        h_dphi_deta->Draw("SURF1 FB BB");
        myText(0.4,0.76,kBlack, "Corrected correlation");
        gPad->SetPhi(-55);
        gPad->SetTheta(+30);
        h_dphi_deta->GetXaxis()->SetNdivisions(1);
        h_dphi_deta->GetYaxis()->SetNdivisions(4);
        //h_dphi_deta->SetMinimum(0);
        h_dphi_deta->GetZaxis()->SetNdivisions(1);
        h_dphi_deta->GetXaxis()->CenterTitle(kTRUE);
        h_dphi_deta->GetYaxis()->CenterTitle(kTRUE);
        h_dphi_deta->SetLineWidth(1.0);
        
        h_dphi_deta->SetTitle("; |#Delta#eta|; #Delta#phi/#pi");
        Guardar(c, Form("%s_2DCorr_%s", name, tag ));
        c->Clear();
        //Projecting dphi
        std::cout << "Projecting into eta bins" << std::endl;
        for(int i=1; i< h_dphi_deta->GetXaxis()->GetNbins()+1 ; i++){
            double width =  h_dphi_deta->GetXaxis()->GetBinWidth(i) ;
            std::cout << " GetBinEdge" << h_dphi_deta->GetXaxis()->GetBinLowEdge(i) << " " <<  h_dphi_deta->GetXaxis()->GetBinLowEdge(i) + width << " i " << i << std::endl;
        }
        //h_dphi_near = (TH1F*)h_dphi_deta->ProjectionY("h_dphi_near",7, 14);
        int binmin = h_dphi_deta->GetXaxis()->FindBin(-0.6);
        int binmax = h_dphi_deta->GetXaxis()->FindBin(+0.6)-1;
        std::cout <<" Minimum and maximum" << binmin << " " << binmax << std::endl;
        std::cout << " binmin " << binmin << " binmax " << binmax << std::endl;
        h_dphi_near = (TH1F*)h_dphi_deta->ProjectionY("h_dphi_near", binmin, binmax);
        
        binmin =  h_dphi_deta->GetXaxis()->FindBin(-1.4);
        binmax = h_dphi_deta->GetXaxis()->FindBin(-0.8)-1;
        std::cout <<" Minimum and maximum" << binmin << " " << binmax << std::endl;
        auto h_dphi_far1 = (TH1F*)h_dphi_deta->ProjectionY("h_dphi_far1",binmin, binmax);
        h_dphi_far1->Sumw2();
        
        binmin =  h_dphi_deta->GetXaxis()->FindBin(+0.8);
        binmax = h_dphi_deta->GetXaxis()->FindBin(+1.4)-1;
        std::cout <<" Minimum and maximum" << binmin << " " << binmax << std::endl;
        h_dphi_far = (TH1F*)h_dphi_deta->ProjectionY("h_dphi_far", binmin, binmax);
        h_dphi_far->Sumw2();
        
        h_dphi_far->Add(h_dphi_far1);
        
        h_dphi_near->SetLineColor(kAzure-3);
        h_dphi_near->SetMarkerColor(kAzure-3);
        h_dphi_far->SetMarkerColorAlpha(kRed,0.75);
        h_dphi_far->SetLineColorAlpha(kRed, 0.75);
        h_dphi_far->SetMarkerStyle(20);
        h_dphi_far->SetLineStyle(kDashed);
        
        h_dphi_near->Sumw2();
        h_dphi_far->Sumw2();
        h_dphi_near->Draw();
        h_dphi_near->SetMinimum(0);
        h_dphi_near->GetYaxis()->SetNdivisions(5);
        h_dphi_far->Draw("same");
        h_dphi_near->GetYaxis()->SetTitle("arb units");
        
        
        TString label;
        if(name=="pt"){
            label = Form("%2.1f< p_{T}^{h} <%2.1f GeV", min, max);
        }
        else if(name=="zt"){
            label = Form("%2.1f<Z_{T}<%2.1f", min, max);
        }
        else if(name=="xi"){
            label = Form("%2.1f< #xi# <%2.1f", min, max);
        }
        ///myText(.58,.85,kBlack, Form("%2.1f < p_{T}^{#pi^{0}} < %2.0f GeV",8.0,20.0));
        myText(.58,.85,kBlack, label);
        myText(.60,.79,kAzure-3, "|#Delta#eta|<0.6");
        myText(.60,.73,kRed,     "0.8<|#Delta#eta|<1.4");
        
        //Call fit function:
        TF1 *func = new TF1("fit",fitf,-0.5,1.5,6);
        func->SetLineColorAlpha(kAzure,0.25);
        func->SetLineWidth(2.0);
        func->SetParameters(100,  0.2, 100,  0.5, 3000);
        func->SetParLimits(1, 0.05, 1.0); // widths
        func->SetParLimits(3, 0.05, 1.0); // widths
        func->SetParLimits(0, 1, 10000);  //yield
        func->SetParLimits(2, 1, 10000);  //yield
        
        h_dphi_near->Fit(func);
        func->Draw("same");
        h_dphi_near->GetYaxis()->CenterTitle(kTRUE);
        Guardar(c, Form("%s_hdphi_withFit_%s", name, tag ));
        
        //fill vectors with fit results
        Yield_near.push_back(func->GetParameter(0));
        Yield_near_err.push_back(func->GetParError(0));
        Yield_far.push_back(func->GetParameter(2));
        Yield_far_err.push_back(func->GetParError(2));
        ///Do |eta| projection
        
        c->Clear();
        h_deta =  (TH1F*)h_dphi_deta->ProjectionX("h_deta",1,10); // project the deta bu only in the near side
        h_deta->Draw();
        h_deta->GetXaxis()->SetNdivisions(9);
        h_deta->GetYaxis()->SetNdivisions(4);
        myText(.20,.85,kBlack, "|#Delta#phi| < #pi/2");
        
        std::cout << "phi bins" << std::endl;
        
        for(int i=1; i< h_dphi_deta->GetYaxis()->GetNbins()+1 ; i++){
            double width =  h_dphi_deta->GetYaxis()->GetBinWidth(i) ;
            std::cout << " GetBinEdge" << h_dphi_deta->GetYaxis()->GetBinLowEdge(i) << " " <<  h_dphi_deta->GetYaxis()->GetBinLowEdge(i) + width << " i " << i << std::endl;
        }
        
        std::cout << " FIT GAUSS + p0" << std::endl;
        TF1 *f = new TF1("f",GaussP0,-1.0,1.0,3);
        f->SetLineColor(kRed);
        f->SetLineWidth(2.0);
        f->SetParameters(100, 0.5, 3000);
        f->SetParLimits(0, 0.0, 1e5);
        f->SetParLimits(1, 0.05, 2.0 );
        
        h_deta->Fit(f);
        myText(.20,.79,kRed, Form("#sigma = %2.2f", f->GetParameter(1)));
        myText(.55,.85,kBlack, label);
        c->SaveAs(Form("PNGOUTPUT/hdeta_withFit_%s_%1.f_%1.f.png", name, 10*min, 10*max));
        c->SaveAs(Form("PDFOUTPUT/hdeta_withFit_%s_%1.f_%1.f.pdf", name, 10*min, 10*max));
        //undo the cuts on track pt for next iteration
        SetCut(h, axisNumber, 0, 100);
        SetCut(h_mixed, axisNumber, 0, 100);
    }//end loop
    
    
    TGraphErrors* g_yieldnear = new TGraphErrors();
    TGraphErrors* g_yieldfar = new TGraphErrors();
    double midpoint;
    double widthbin;
    double error;
    for(int n=0; n < Yield_near.size(); n++){
        std::cout<< Form("%2.2f--%2.2f = %2.f ", bins.at(n), bins.at(n+1), Yield_near.at(n)) << std::endl;
        midpoint = (bins.at(n)+bins.at(n+1))/2.0;
        widthbin = (bins.at(n+1)-bins.at(n));
        error = widthbin/2.0;
        g_yieldnear->SetPoint(n, midpoint, Yield_near.at(n)/widthbin);
        g_yieldnear->SetPointError(n, error, Yield_near_err.at(n)/widthbin);
        
        g_yieldfar->SetPoint(n, midpoint, Yield_far.at(n)/widthbin);
        g_yieldfar->SetPointError(n, error, Yield_far_err.at(n)/widthbin);
    }
    
    g_yieldfar->SetMarkerColor(kRed);
    g_yieldfar->SetLineColor(kRed);
    g_yieldnear->SetMarkerColor(kAzure-3);
    g_yieldnear->SetLineColor(kAzure-3);
    
    
    //TMultiGraph* multi = new TMultiGraph();
    multi.Add(g_yieldnear);
    multi.Add(g_yieldfar);
    
    multi.SetMinimum(0);
    
    TString title;
    if(name=="xi"){
        title = "#xi";
    }
    
    c->Clear();
    multi.Draw("AP");
    multi.GetXaxis()->SetNdivisions(6);
    multi.GetYaxis()->SetNdivisions(5);
    if(name=="xi"){
        multi.GetXaxis()->SetTitle("#xi");
        multi.GetYaxis()->SetTitle("dN/d#xi");
    }
    else if(name=="pt"){
        multi.GetXaxis()->SetTitle("p^{h}_{T} [GeV]");
        multi.GetYaxis()->SetTitle("dN/dp_{T}");
    }
    else if(name=="zt"){
        multi.GetXaxis()->SetTitle("Z_{T}");
        multi.GetYaxis()->SetTitle("dN/dz_{T}");
    }
    
    myText(.65,.85,kAzure-3, "Near side");
    myText(.65,.80,kRed, "Away side");
    //multi->SetTitle("; " + title + " ; Yield");
    Guardar(c, Form("%s_Yields", name));
    
    
    c->Close();
    delete c;
    
    return;
}

void PlotCorrelation(THnSparse* h, THnSparse* h_mixed, double axisAssoc, std::vector<double> bins){
    //Write this in such a way that it is agnostic on whether you work on pion , cluster , merged cluster.
    auto c = new TCanvas();
    
    std::cout<< "Line 404" << std::endl;
    
    //Assoc variable (pT, xi or zt)
    auto h_1D =  h->Projection(axisAssoc);
    gPad->SetLogy(kTRUE);
    //    gPad->SetLogx(kTRUE);
    h_1D->Draw("hist");
    h_1D->SetTitle("; Zt; entries");
    h_1D = h_mixed->Projection(axisAssoc);
    h_1D->SetLineColor(2);
    h_1D->Draw("histsame");
    Guardar(c, "Track_zt");
    c->Clear();
    
    std::cout<< "Line 418" << std::endl;
    
    //cumulative track pt
    h_1D =  h->Projection(axisAssoc);
    h_1D->Scale(1.0/h_1D->Integral());
    auto h_cumulative = h_1D->GetCumulative();
    h_cumulative->Draw("hist");
    gPad->SetLogy(kFALSE);
    h_cumulative->SetTitle("; zt; cumulative prob");
    h_1D =  h_mixed->Projection(axisAssoc);
    h_1D->Scale(1.0/h_1D->Integral());
    h_cumulative = h_1D->GetCumulative();
    h_cumulative->SetLineColor(kRed);
    h_cumulative->Draw("hist same");
    Guardar(c, "Track_zt_cumulative");
    //Loop over pT, xi or Zt bins and perform analysis.
    
    std::cout<< "Line 435" << std::endl;
    
    TMultiGraph multi;
    Looping(h, h_mixed, axisAssoc, bins, "zt", multi);
    c->Close();
    
    std::cout<< "Line 441" << std::endl;
    
    multi.Print();
    delete c;
    return;
}





void PlotPionHistograms(THnSparse* h_Pion){
    ////
    
    auto c = new TCanvas();
    h_Pion->GetAxis(axis_pionMass)->SetTitle("m_{#gamma#gamma} [GeV]");
    h_Pion->GetAxis(axis_pionPt)->SetTitle("p^{#gamma#gamma}_{T} [GeV]");
    
    auto h_1D = h_Pion->Projection(axis_pionMass);
    h_1D->Draw("hist");
    Guardar(c, "Mass_All");
    c->Clear();
    
    
    SetCut(h_Pion, axis_pionMass, 0.100, 0.180); //Invariant mass selection
    
    h_1D = h_Pion->Projection(axis_pionPt);
    gPad->SetLogy(1);
    h_1D->Draw("hist");
    Guardar(c, "PionPt_All");
    gPad->SetLogy(0);
    c->Clear();
    
    SetCut(h_Pion, axis_pionPt, 8.0, 12.0);  // Pion pT selection
    
    h_1D = h_Pion->Projection(axis_pionMass);
    h_1D->Draw("hist");
    Guardar(c, "Mass_withpTCut");
    c->Clear();
    
    //Mass vs pT
    auto h_2D = h_Pion->Projection(axis_pionMass,axis_pionPt);
    h_2D->Draw("colz");
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(6);
    h_2D->GetXaxis()->SetNdivisions(4);
    gPad->SetLogz(1);
    auto h_proj = h_2D->ProfileX("h_proj",0,1e3);
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same");
    Guardar(c, "2D_MassPT");
    c->Clear();
    
    //lambda vs pT photon 1
    SetCut(h_Pion, axis_pion_PhPt_1, 0.0, 15.0);
    h_2D = h_Pion->Projection( axis_pion_PhPt_1, axis_pion_PhM02_1);
    h_2D->Draw("colz");
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(6);
    h_2D->GetXaxis()->SetNdivisions(4);
    Guardar(c, "2D_Ph1Lambda_pT");
    c->Clear();
    
    SetCut(h_Pion,axis_pion_PhPt_2, 0.0, 15.0);
    auto temp = h_Pion->Projection( axis_pion_PhPt_2, axis_pion_PhM02_2);
    temp->Draw("colz");
    temp->GetZaxis()->SetNdivisions(3);
    temp->GetYaxis()->SetNdivisions(6);
    temp->GetXaxis()->SetNdivisions(4);
    Guardar(c, "2D_Ph2Lambda_pT");
    c->Clear();
    
    //Get The Summed histogram.
    h_2D->Add(temp);
    h_2D->Draw("colz");
    h_2D->SetTitle("; Cluster #lambda_{02}; #gamma Energy [GeV]");
    h_2D->GetXaxis()->SetNdivisions(7);
    gPad->SetLogz(kTRUE);
    myText(0.35,0.92,kBlack, "Photons from #pi^{0}#rightarrow#gamma#gamma decay");
    myText(0.6,0.22,kBlack, "120 < m_{#gamma#gamma} < 160 MeV");
    myText(0.6,0.17,kBlack, "8 < p_{T}^{#gamma#gamma}< 20 GeV");
    Guardar(c,"2D_PhLambda_pT");
    c->Clear();
    //project the lambda 02.
    
    auto lambda_1 = h_2D->ProjectionX("lambda02_1", h_2D->GetYaxis()->FindBin(8.0), h_2D->GetYaxis()->FindBin(15.0)-1);
    auto lambda_2 = h_2D->ProjectionX("lambda02_2", h_2D->GetYaxis()->FindBin(10.0), h_2D->GetYaxis()->FindBin(12.0)-1);
    
    
    lambda_1->DrawNormalized("hist");
    //    lambda_1->GetYaxis()->SetRangeUser(0,0.25);
    lambda_2->SetLineColor(kRed);
    //lambda_2->DrawNormalized("hist same");
    myText(0.35,0.92,kBlack, "Photons from #pi^{0}#rightarrow#gamma#gamma decay");
    myText(0.6,0.82,kBlack, "8 < p_{T}^{#gamma}< 20 GeV");
    Guardar(c, "FromPion_Lambda_02");
    
    
    //Mass vs Centrality
    h_2D = h_Pion->Projection(axis_pionMass,axis_pion_Cen);
    h_2D->Draw("colz");
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(6);
    h_2D->GetXaxis()->SetNdivisions(4);
    gPad->SetLogz(0);
    h_proj = h_2D->ProfileX("h_massproj",0,100);
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same");
    Guardar(c, "2D_MassCen");
    
    //Mass vs Zvertex
    h_2D = h_Pion->Projection(axis_pionMass,    axis_pion_Zvtx);
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(6);
    h_2D->GetXaxis()->SetNdivisions(4);
    h_2D->Draw("colz");
    gPad->SetLogz(0);
    h_proj = h_2D->ProfileX("h_proj",0,100);
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same");
    Guardar(c, "2D_MassZvtx");
    
    //Mass vs asymmetry
    SetCut(h_Pion, axis_pionPt, 12.0, 20.0);
    h_2D = h_Pion->Projection(axis_pionMass, axis_pion_asymmetry);
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(6);
    h_2D->GetXaxis()->SetNdivisions(4);
    h_2D->Draw("colz");
    gPad->SetLogz(0);
    h_proj = h_2D->ProfileX("h_proj",0,100);
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same");
    myText(0.6,0.92,kBlack, "12 < p_{T}^{#gamma#gamma}< 20 GeV");
    Guardar(c, "2D_MassAsym");
    c->Clear();
    
    
    
    //Mass vs opening angle
    SetCut(h_Pion, axis_pionOpeningAngle, -0.05,0.05);
    h_2D = h_Pion->Projection(axis_pionMass, axis_pionOpeningAngle);
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(6);
    h_2D->GetXaxis()->SetNdivisions(6);
    h_2D->Draw("colz");
    h_2D->SetTitle("; #Delta#phi; m_{#gamma#gamma} [GeV]");
    h_proj = h_2D->ProfileX("h_proj",0,100);
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same");
    myText(0.6,0.92,kBlack, "12 < p_{T}^{#gamma#gamma}< 20 GeV");
    Guardar(c, "2D_MassOpeningangle");
    c->Clear();
    
    //Get back to 8.0--20 Range.
    SetCut(h_Pion, axis_pionPt, 8.0, 20.0);
    //Asymmetry vs pT
    h_2D =  h_Pion->Projection(axis_pion_asymmetry,axis_pionPt);
    h_2D->Draw("COLZ");
    h_2D->GetZaxis()->SetNdivisions(3);
    h_2D->GetYaxis()->SetNdivisions(5);
    h_2D->GetXaxis()->SetNdivisions(4);
    h_proj = h_2D->ProfileX("h_proj",0,100);
    h_proj->SetMarkerColor(kRed);
    h_proj->SetLineColor(kRed);
    h_proj->Draw("same");
    Guardar(c, "2D_AsymPT");
    c->Clear();
    
    //Lambda0 photon 1 vs photon 2
    h_2D = h_Pion->Projection( axis_pion_PhM02_1 ,  axis_pion_PhM02_2 );
    h_2D->Draw("colz");
    gPad->SetLogz(0);
    Guardar(c, "2D_PhotonM02");
    c->Clear();
    
    //Energy photon 1 vs photon2;
    h_2D = h_Pion->Projection( axis_pion_PhPt_1   , axis_pion_PhPt_2);
    h_2D->Draw("colz");
    h_2D->SetTitle("; #gamma_{1} Energy [GeV]; #gamma_{2} Energy [GeV]");
    h_2D->GetXaxis()->SetRangeUser(0,15.0);
    h_2D->GetYaxis()->SetRangeUser(0,15.0);
    gPad->SetLogz(0);
    Guardar(c, "2D_PhotonPt");
    c->Clear();
    
    //Pion eta
    h_1D = h_Pion->Projection(axis_pionRapidity);
    h_1D->Draw("PL");
    h_1D->Sumw2();
    Guardar(c, "PionRapidity");
    c->Clear();
    
    
    c->Close();
    delete c;
}

void PlotClusterHistograms(THnSparse* h){
    ////
    //SetCut(h, axis_Cluster_Ncells, 3.0, 100);
    
    auto c = new TCanvas();
    h->GetAxis(axis_Cluster_Pt)->SetTitle("Cluster E_{T} [GeV]");
    h->GetAxis(axis_Cluster_Phi)->SetTitle("Cluster #phi [rad]");
    h->GetAxis(axis_Cluster_Eta)->SetTitle("Cluster #eta");
    h->GetAxis(axis_Cluster_Ncells)->SetTitle("Number of cells in cluster");
    
    //Plotting ET distribution
    auto h_1D = h->Projection(axis_Cluster_Pt);
    h_1D->Draw("hist");
    gPad->SetLogy(kTRUE);
    h_1D->GetYaxis()->SetTitle("entries");
    Guardar(c, "Clusters_All" );
    gPad->SetLogy(kFALSE);
    c->Clear();
    
    SetCut(h,  axis_Cluster_Pt, 8, 100); // Cut Et
    SetCut(h , axis_Cluster_Phi, 1.2, TMath::Pi());
    auto h_2D  = h->Projection(axis_Cluster_Eta, axis_Cluster_Phi);
    
    h_2D->RebinY(3);
    h_2D->Draw("COLZ");
    myText(0.25,0.85,kBlack, "E_{T}> 8 GeV");
    Guardar(c, "Cluster_EtaPhi_NoCut");
    c->Clear();
    
    //Cutting on distance-to-border
    for (int n=1; n<4; n++){
        SetCut(h, axis_Cluster_DisToBorder, 1.0*n, 6.0);
        h_2D  = h->Projection(axis_Cluster_Eta, axis_Cluster_Phi);
        h_2D->RebinY(3);
        h_2D->Draw("COLZ");
        myText(0.25,0.85,kBlack, "E_{T}> 8 GeV");
        myText(0.35,0.92,kBlack, Form("%i cell away from border", n));
        Guardar(c, Form("Cluster_EtaPhi_%i_fromBorder", n));
        c->Clear();
    }
    
    
    
    
    //Projecting phi
    h_1D = h->Projection(axis_Cluster_Phi);
    h_1D->Draw();
    std::cout << " BinWidth Phi " << h_1D->GetBinWidth(1) << "rads" << std::endl;
    Guardar(c, "Cluster_Phi");
    c->Clear();
    
    h_1D = h->Projection(axis_Cluster_Eta);
    h_1D->Draw();
    h_1D->Rebin(3);
    std::cout << " BinWidth Eta " << h_1D->GetBinWidth(1) << std::endl;
    Guardar(c, "Cluster_Eta");
    c->Clear();
    
    //Number of cells in the cluster:
    h_1D = h->Projection(axis_Cluster_Ncells);
    h_1D->Draw("hist");
    Guardar(c, "Ncells");
    c->Clear();
    
    //Number of cells vs Et
    h_2D  = h->Projection(axis_Cluster_Ncells, axis_Cluster_Pt);
    h_2D->Draw("COLZ");
    //gPad->SetLogx(kTRUE);
    Guardar(c, "NcellsVsEt");
    gPad->SetLogx(kFALSE);
    c->Clear();
    
    //Cutting on Ncells
    SetCut(h, axis_Cluster_Ncells, 3.0, 100);
    
    //Distance to bad channel
    h_1D = h->Projection(axis_Cluster_DisToBad);
    h_1D->Draw("hist");
    Guardar(c, "DistanceToBad");
    c->Clear();
    //Cutting on distance-to-bad-channel, plotting 2D eta--phi distribution
    for (int n=1; n<4; n++){
        SetCut(h, axis_Cluster_DisToBad, 1.0*n, 9.0);
        h_2D  = h->Projection(axis_Cluster_Eta, axis_Cluster_Phi);
        h_2D->RebinY(3);
        h_2D->Draw("COLZ");
        myText(0.25,0.85,kBlack, "E_{T}> 8 GeV");
        myText(0.35,0.92,kBlack, Form("%i cell away from bad channel", n));
        Guardar(c, Form("Cluster_EtaPhi_%i_fromBadChannel", n));
        c->Clear();
    }
    
    c->Close();
    delete c;
    delete h_1D;
}




void Plotting(){
    SetAtlasStyle();
    
    TFile* fIn = new TFile("THnSparses_LHC13d_101517.root","READ");
    fIn->Print();
    
    //Getting the different THnSparses
    THnSparse* h_Pion = 0;
    fIn->GetObject("h_Pion",h_Pion);
    THnSparse *h_PionTrack=0;
    fIn->GetObject("h_PionTrack", h_PionTrack);
    THnSparse *  h_PionTrack_Mixed = 0;
    fIn->GetObject("h_PionTrack_Mixed", h_PionTrack_Mixed);
    
    THnSparse* h_Cluster = 0;
    fIn->GetObject("h_Cluster",h_Cluster);
    //THnSparse *h_ClusterTrack=0;
    //fIn->GetObject("h_ClusterTrack", h_ClusterTrack);
    //THnSparse *  h_ClusterTrack_Mixed = 0;
    //fIn->GetObject("h_ClusterTrack_Mixed", h_ClusterTrack_Mixed);
    
    //THnSparse* h_Track = 0;
    //fIn->GetObject("h_Track", h_Track);
    
    std::vector<double> bins_trackpt;
    bins_trackpt.push_back(1.0);
    bins_trackpt.push_back(2.0);
    bins_trackpt.push_back(4.0);
    bins_trackpt.push_back(10.0);
    
    std::vector<double> bins_trackzt;
    bins_trackzt.push_back(0.0);
    bins_trackzt.push_back(0.2);
    bins_trackzt.push_back(0.5);
    
    //Cutting the pion-h THnSparse
    SetCut(h_PionTrack, axis_corr_photon1M02, 0.0, 0.4);
    SetCut(h_PionTrack, axis_corr_photon2M02, 0.0, 0.4);
    SetCut(h_PionTrack_Mixed, axis_corr_photon1M02, 0.0, 0.4);
    SetCut(h_PionTrack_Mixed, axis_corr_photon2M02, 0.0, 0.4);
    
    //Cut on pT of pion and mass.
    SetCut(h_PionTrack,  axis_corr_triggerpT, 8.0 ,15.0);
    SetCut(h_PionTrack_Mixed,  axis_corr_triggerpT, 8.0 ,15.0);
    SetCut(h_PionTrack, axis_corr_mass, 0.100,0.180);
    SetCut(h_PionTrack_Mixed, axis_corr_mass, 0.100, 0.180);
    
    //Cutting Cluster-track:
    
    //SetCut(h_ClusterTrack, axis_corr_triggerpT, 15.0, 50.0);
    //SetCut(h_ClusterTrack_Mixed, axis_corr_triggerpT, 15.0, 50.0);
    //SetCut(h_ClusterTrack, axis_corr_M02, 0.0, 0.4);
    //SetCut(h_ClusterTrack_Mixed, axis_corr_M02, 0.0, 0.4);
    //SetCut(h_ClusterTrack, axis_corr_dR, 0.02, 100);
    //SetCut(h_ClusterTrack_Mixed, axis_corr_dR, 0.02, 100);
    
    PlotCorrelation(h_PionTrack, h_PionTrack_Mixed,axis_corr_zt, bins_trackzt);
    //PlotCorrelation(h_ClusterTrack, h_ClusterTrack_Mixed, axis_corr_zt, bins_trackzt);
    
    //PlotClusterHistograms(h_Cluster);
}
/*
 
 //Define bin edges for pT, xi and Zt intervals.
 std::vector<double> bins_trackpt;
 bins_trackpt.push_back(0.5);
 bins_trackpt.push_back(1.0);
 bins_trackpt.push_back(2.0);
 //bins_pt.push_back(10.0);
 std::vector<double> bins_trackxi;
 bins_trackxi.push_back(0.0);
 bins_trackxi.push_back(1.5);
 bins_trackxi.push_back(2.5);
 bins_trackxi.push_back(3.0);
 //bins_xi.push_back(10.0);
 std::vector<double> bins_trackzt;
 bins_trackzt.push_back(0.0);
 bins_trackzt.push_back(0.5);
 bins_trackzt.push_back(1.0);
 bins_trackzt.push_back(2.0);
 
 auto c = new TCanvas();
 auto h = h_PionTrack->Projection(axis_corr_mass);
 h->SetTitle("; m_{#gamma#gamma} [GeV]; entries");
 [6~    h->Draw("hist");
 Guardar(c, "MassCorr");
 c->Clear();
 
 //Set limit on pion variables before running the correlation analysis.
 //Invariant mass cut
 SetCut(h_PionTrack, axis_corr_mass, 0.100,0.180);
 SetCut(h_PionTrack_Mixed, axis_corr_mass, 0.100, 0.180);
 //Pt cuts:
 SetCut(h_PionTrack,  axis_corr_triggerpT, 8.0 ,16.0);
 SetCut(h_PionTrack_Mixed,  axis_corr_triggerpT, 8.0 ,16.0);
 //PlotCorrelation(h_PionTrack, h_PionTrack_Mixed,axis_corr_trackpT, bins_trackpt);
 
 //Perform correlation for pions with pT between 8 and 16 GeV.
 
 //h = h_ClusterTrack->Projection(axis_Clustercorr_M02);
 //h->SetTitle("; Cluster #lambda_{02}; entries");
 //h->Draw("hist");
 //Guardar(c, "Lambda02Corr");
 //c->Clear();
 //c->Close();
 
 SetCut(h_ClusterTrack,  axis_Clustercorr_M02, 0.1, 0.4 );
 SetCut(h_ClusterTrack_Mixed,  axis_Clustercorr_M02, 0.1, 0.4 );
 // PlotCorrelation(h_ClusterTrack, h_ClusterTrack_Mixed, axis_corr_trackpT, bins_trackpt);
 
 PlotClusterHistograms(h_Cluster);
 //PlotPionHistograms(h_Pion);
 }*/
