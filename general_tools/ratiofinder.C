// This macro finds the ratios of various graphs from the two inputs submitted to it and outputs a ROOT file with the quotients

#include <iostream>

// 1D histogram dividing function
// Precondition: the two histograms must have the same x-dimensions
TH1D* divide_histograms1D(TH1D* graph1, TH1D* graph2, std::string name, std::string title){
    // Make the 1D histogram to contain the quotient and find minimum and maximum bins along both axes
    TH1D* quotient = new TH1D(*graph2);
    quotient->SetNameTitle(name.c_str(), title.c_str());
    double x_bin_min = quotient->GetXaxis()->FindBin(quotient->GetXaxis()->GetXmin());
    double x_bin_max = quotient->GetXaxis()->FindBin(quotient->GetXaxis()->GetXmax());
    
    // Loop over all bins, divide the element from graph 1 by its counterpart on graph 2. Manually propogate the error
    for(int i = x_bin_min; i <= x_bin_max; i++){
        if (graph2->GetBinContent(i) != 0)
            quotient->SetBinContent(i, (graph1->GetBinContent(i))/(graph2->GetBinContent(i)));
        // Failsafe
        else
            quotient->SetBinContent(i, 0);
        
        // Error propogation
        if (graph2->GetBinContent(i) == 0) {
            if (graph1->GetBinContent(i) == 0)
                quotient->SetBinError(i, 0);
            else
                quotient->SetBinError(i, (quotient->GetBinContent(i)*TMath::Sqrt(( (graph1->GetBinError(i)/graph1->GetBinContent(i)) * (graph1->GetBinError(i)/graph1->GetBinContent(i)) ) )));
        }
        else if (graph1->GetBinContent(i) == 0)
            quotient->SetBinError(i, (quotient->GetBinContent(i)*TMath::Sqrt(( (graph2->GetBinError(i)/graph2->GetBinContent(i)) * (graph2->GetBinError(i)/graph2->GetBinContent(i)) ) )));
        else
            quotient->SetBinError(i, (quotient->GetBinContent(i)*TMath::Sqrt(( (graph1->GetBinError(i)/graph1->GetBinContent(i)) * (graph1->GetBinError(i)/graph1->GetBinContent(i)) ) + ( (graph2->GetBinError(i)/graph2->GetBinContent(i)) * (graph2->GetBinError(i)/graph2->GetBinContent(i)) ) )));
    }
    
    
    return quotient;
}

void ratiofinder(std::string dividenddataname, std::string divisordataname, std::string outfilename) {
    TFile *divisordata = new TFile(divisordataname.c_str(), "READ");
    TFile *dividenddata = new TFile(dividenddataname.c_str(), "READ");
    
    TH1D* XobsPb_sig_dividend = 0;
    TH1D* XobsPb_bkg_dividend = 0;
    TH1D* XobsPb_adj_dividend = 0;
    TH1D* XobsPb_sig_divisor = 0;
    TH1D* XobsPb_bkg_divisor = 0;
    TH1D* XobsPb_adj_divisor = 0;
    
    TH1D* dPhi_sig_dividend = 0;
    TH1D* dPhi_bkg_dividend = 0;
    TH1D* dPhi_adj_dividend = 0;
    TH1D* dPhi_sig_divisor = 0;
    TH1D* dPhi_bkg_divisor = 0;
    TH1D* dPhi_adj_divisor = 0;
    
    TH1D* Xj_sig_dividend = 0;
    TH1D* Xj_bkg_dividend = 0;
    TH1D* Xj_adj_dividend = 0;
    TH1D* Xj_sig_divisor = 0;
    TH1D* Xj_bkg_divisor = 0;
    TH1D* Xj_adj_divisor = 0;
    
    TH1D* pTD_sig_dividend = 0;
    TH1D* pTD_bkg_dividend = 0;
    TH1D* pTD_adj_dividend = 0;
    TH1D* pTD_sig_divisor = 0;
    TH1D* pTD_bkg_divisor = 0;
    TH1D* pTD_adj_divisor = 0;
    
    TH1D* Multiplicity_sig_dividend = 0;
    TH1D* Multiplicity_bkg_dividend = 0;
    TH1D* Multiplicity_adj_dividend = 0;
    TH1D* Multiplicity_sig_divisor = 0;
    TH1D* Multiplicity_bkg_divisor = 0;
    TH1D* Multiplicity_adj_divisor = 0;
    
    dividenddata->GetObject("sig_XobsPb", XobsPb_sig_dividend);
    dividenddata->GetObject("bkg_XobsPb", XobsPb_bkg_dividend);
    dividenddata->GetObject("hadj_XobsPb", XobsPb_adj_dividend);
    divisordata->GetObject("sig_XobsPb", XobsPb_sig_divisor);
    divisordata->GetObject("bkg_XobsPb", XobsPb_bkg_divisor);
    divisordata->GetObject("hadj_XobsPb", XobsPb_adj_divisor);
    
    dividenddata->GetObject("sig_dPhi", dPhi_sig_dividend);
    dividenddata->GetObject("bkg_dPhi", dPhi_bkg_dividend);
    dividenddata->GetObject("hadj_dPhi", dPhi_adj_dividend);
    divisordata->GetObject("sig_dPhi", dPhi_sig_divisor);
    divisordata->GetObject("bkg_dPhi", dPhi_bkg_divisor);
    divisordata->GetObject("hadj_dPhi", dPhi_adj_divisor);
    
    dividenddata->GetObject("sig_Xj", Xj_sig_dividend);
    dividenddata->GetObject("bkg_Xj", Xj_bkg_dividend);
    dividenddata->GetObject("hadj_Xj", Xj_adj_dividend);
    divisordata->GetObject("sig_Xj", Xj_sig_divisor);
    divisordata->GetObject("bkg_Xj", Xj_bkg_divisor);
    divisordata->GetObject("hadj_Xj", Xj_adj_divisor);
    
    dividenddata->GetObject("sig_pTD", pTD_sig_dividend);
    dividenddata->GetObject("bkg_pTD", pTD_bkg_dividend);
    dividenddata->GetObject("hadj_pTD", pTD_adj_dividend);
    divisordata->GetObject("sig_pTD", pTD_sig_divisor);
    divisordata->GetObject("bkg_pTD", pTD_bkg_divisor);
    divisordata->GetObject("hadj_pTD", pTD_adj_divisor);
    
    dividenddata->GetObject("sig_Multiplicity", Multiplicity_sig_dividend);
    dividenddata->GetObject("bkg_Multiplicity", Multiplicity_bkg_dividend);
    dividenddata->GetObject("hadj_Multiplicity", Multiplicity_adj_dividend);
    divisordata->GetObject("sig_Multiplicity", Multiplicity_sig_divisor);
    divisordata->GetObject("bkg_Multiplicity", Multiplicity_bkg_divisor);
    divisordata->GetObject("hadj_Multiplicity", Multiplicity_adj_divisor);
    
    TH1D* XobsPb_sig_quotient = divide_histograms1D(XobsPb_sig_dividend, XobsPb_sig_divisor, "sig_XobsPb_ratio", "x_{obs}^{Pb} Signal Region #frac{pPb}{pp} ratio; x_{obs}^{Pb} ratio;");
    TH1D* XobsPb_bkg_quotient = divide_histograms1D(XobsPb_bkg_dividend, XobsPb_bkg_divisor, "bkg_XobsPb_ratio", "x_{obs}^{Pb} Background Region #frac{pPb}{pp} ratio; x_{obs}^{Pb} ratio;");
    TH1D* XobsPb_adj_quotient = divide_histograms1D(XobsPb_adj_dividend, XobsPb_adj_divisor, "hadj_XobsPb_ratio", "x_{obs}^{Pb} #frac{pPb}{pp} ratio; x_{obs}^{Pb} ratio;");
    TH1D* dPhi_sig_quotient = divide_histograms1D(dPhi_sig_dividend, dPhi_sig_divisor, "sig_dPhi_ratio", "#Delta #phi Signal Region #frac{pPb}{pp} ratio; #Delta #phi ratio;");
    TH1D* dPhi_bkg_quotient = divide_histograms1D(dPhi_bkg_dividend, dPhi_bkg_divisor, "bkg_dPhi_ratio", "#Delta #phi Background Region  #frac{pPb}{pp}ratio; #Delta #phi ratio;");
    TH1D* dPhi_adj_quotient = divide_histograms1D(dPhi_adj_dividend, dPhi_adj_divisor, "hadj_dPhi_ratio", "#Delta #phi #frac{pPb}{pp} ratio; #Delta #phi ratio;");
    TH1D* Xj_sig_quotient = divide_histograms1D(Xj_sig_dividend, Xj_sig_divisor, "sig_Xj_ratio", "Xj Signal Region #frac{pPb}{pp} ratio; Xj ratio;");
    TH1D* Xj_bkg_quotient = divide_histograms1D(Xj_bkg_dividend, Xj_bkg_divisor, "bkg_Xj_ratio", "Xj Background #frac{pPb}{pp} Region ratio; Xj ratio;");
    TH1D* Xj_adj_quotient = divide_histograms1D(Xj_adj_dividend, Xj_adj_divisor, "hadj_ratio", "Xj #frac{pPb}{pp} ratio; Xj ratio;");
    TH1D* pTD_sig_quotient = divide_histograms1D(pTD_sig_dividend, pTD_sig_divisor, "sig_pTD_ratio", "pTD Signal Region #frac{pPb}{pp} ratio; p_{T}D ratio;");
    TH1D* pTD_bkg_quotient = divide_histograms1D(pTD_bkg_dividend, pTD_bkg_divisor, "bkg_pTD_ratio", "pTD Background Region #frac{pPb}{pp} ratio; p_{T}D ratio;");
    TH1D* pTD_adj_quotient = divide_histograms1D(pTD_adj_dividend, pTD_adj_divisor, "hadj_pTD_ratio", "pTD #frac{pPb}{pp} ratio; p_{T}D ratio;");
    TH1D* Multiplicity_sig_quotient = divide_histograms1D(Multiplicity_sig_dividend, Multiplicity_sig_divisor, "sig_Multiplicity_ratio", "Multiplicity Signal Region #frac{pPb}{pp} ratio; Multiplicity ratio;");
    TH1D* Multiplicity_bkg_quotient = divide_histograms1D(Multiplicity_bkg_dividend, Multiplicity_bkg_divisor, "bkg_Multiplicity_ratio", "Multiplicity Background Region #frac{pPb}{pp} ratio; Multiplicity ratio;");
    TH1D* Multiplicity_adj_quotient = divide_histograms1D(Multiplicity_adj_dividend, Multiplicity_adj_divisor, "hadj_Multiplicity_ratio", "Multiplicity #frac{pPb}{pp} ratio; Multiplicity ratio;");
    
    TFile *fout = new TFile(Form("%s", outfilename.c_str()), "RECREATE");
    
    XobsPb_sig_quotient->Write("sig_XobsPb_ratio");
    XobsPb_bkg_quotient->Write("bkg_XobsPb_ratio");
    XobsPb_adj_quotient->Write("hadj_XobsPb_ratio");
    dPhi_sig_quotient->Write("sig_dPhi_ratio");
    dPhi_bkg_quotient->Write("bkg_dPhi_ratio");
    dPhi_adj_quotient->Write("hadj_dPhi_ratio");
    Xj_sig_quotient->Write("sig_Xj_ratio");
    Xj_bkg_quotient->Write("bkg_Xj_ratio");
    Xj_adj_quotient->Write("hadj_Xj_ratio");
    pTD_sig_quotient->Write("sig_pTD_ratio");
    pTD_bkg_quotient->Write("bkg_pTD_ratio");
    pTD_adj_quotient->Write("hadj_pTD_ratio");
    Multiplicity_sig_quotient->Write("sig_Multiplicity_ratio");
    Multiplicity_bkg_quotient->Write("bkg_Multiplicity_ratio");
    Multiplicity_adj_quotient->Write("hadj_Multiplicity_ratio");
    
    dividenddata->Close();
    divisordata->Close();
    fout->Close();
}
