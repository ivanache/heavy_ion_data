// This program runs the CompareHistograms.py program on the dPhi, XobsPb, Xj, pTD and Multiplicity variables
#include <iostream>
#include <cstdlib>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <THStack.h>
#include <TProfile.h>

const int num_of_types = 3;
const std::string hist_types[num_of_types] = {"hadj", "sig", "bkg"};

// Two required arguments: the name of the input file, and the name of the tag to be applied to the file name
int main(int argc, char *argv[]) {
    if (argc < 3) {
        exit(EXIT_FAILURE);
    }
    
    // Get the name of the inputfile and filetag
    std::string filename = (std::string) argv[1];
    std::string tagname = (std::string) argv[2];
    
    for (int i=0; i < num_of_types; i++){
    // loop over hadj histograms, hSR histograms, hBR histograms
        system(Form("python CompareHistograms.py --hists %s_dPhi_ratio -m %s -t %s -x \"#Delta #phi (rads)\" -y \"#frac{pPb}{pp} correlation ratio\" -a 3.5 -z 0 -i %s", hist_types[i].c_str(), tagname.c_str(), tagname.c_str(), filename.c_str()));
        system(Form("python CompareHistograms.py --hists %s_XobsPb_ratio -m %s -t %s -x \"x_{obs}^{Pb}:= #frac{p_{T}^{#gamma}e^{-#eta^{#gamma}} + p_{T}^{jet}e^{-#eta^{jet}}}{2E_{Pb}}\" -y \"#frac{pPb}{pp} correlation ratio\" -a 3.5 -z 0 -i %s", hist_types[i].c_str(), tagname.c_str(), tagname.c_str(), filename.c_str()));
        system(Form("python CompareHistograms.py --hists %s_Xj_ratio -m %s -t %s -x \"#frac{p_{T}^{jet}}{p_{T}^{#gamma}}\" -y \"#frac{pPb}{pp} correlation ratio\" -a 3.5 -z 0 -i %s", hist_types[i].c_str(), tagname.c_str(), tagname.c_str(), filename.c_str()));
        system(Form("python CompareHistograms.py --hists %s_pTD_ratio -m %s -t %s -x \"p_{T}D\" -y \"#frac{pPb}{pp} correlation ratio\" -a 3.5 -z 0 -i %s", hist_types[i].c_str(), tagname.c_str(), tagname.c_str(), filename.c_str()));
        system(Form("python CompareHistograms.py --hists %s_Multiplicity_ratio -m %s -t %s -x \"Multiplicity\" -y \"#frac{pPb}{pp} correlation ratio\" -a 3.5 -z 0 -i %s", hist_types[i].c_str(), tagname.c_str(), tagname.c_str(), filename.c_str()));
    }
    
}
