// This macro separates out the truth variable(s) from a ROOT histogram so that CompareHistograms can put them on the same histogram as the measured one from the same run

// fdataname: the data file from where the truth histogram(s) is(are) to be drawn; outfilename: the new ROOT file with the histograms
void truthseparator(std::string fdataname, std::string outfilename) {
    TFile *fdata = new TFile(fdataname.c_str(), "READ");
    
    TH1D* XobsPb_truth_sig = 0;
    TH1D* XobsPb_truth_bkg = 0;
    
    TH1D* Xj_truth_sig = 0;
    TH1D* Xj_truth_bkg = 0;

    fdata->GetObject("sig_XobsPb_truth", XobsPb_truth_sig);
    fdata->GetObject("bkg_XobsPb_truth", XobsPb_truth_bkg);
    
    fdata->GetObject("sig_Xj_truth", Xj_truth_sig);
    fdata->GetObject("bkg_Xj_truth", Xj_truth_bkg);
    
    TFile *fout = new TFile(Form("%s", outfilename.c_str()), "RECREATE");
    
    XobsPb_truth_sig->Write("sig_XobsPb");
    XobsPb_truth_bkg->Write("bkg_XobsPb");
    Xj_truth_sig->Write("sig_Xj");
    Xj_truth_bkg->Write("bkg_Xj");
    
    fdata->Close();
    fout->Close();
}
