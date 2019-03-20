 // This C++ file takes a 1D histogram from a ROOT file and integrates it over an inputted x-range
// Author: Ivan Chernyshev; Date: 9/12/2018
#include <TH1.h>
#include <TFile.h>

#include <TCanvas.h>

#include <iostream>

int main(int argc, char *argv[]) {
    if(argc < 5) {
        std::cout << "Error: syntax is ./Integrate [ROOT filename] [histogram name] [x-minimum] [x-maximum]" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Get filename
    std::string filename = (std::string)argv[1];
    TFile* inputfile = new TFile(filename.c_str(), "READ");
    if (inputfile == NULL){
        std::cout << "Error: file " << filename << " not found" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Get histogram
    std::string histname = (std::string)argv[2];
    TH1* histogram = 0;
    inputfile->GetObject(histname.c_str(), histogram);
    if (histogram == NULL){
        std::cout << "Error: histogram " << histname << " not found" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //Get the minimum and maximum indices and integrate
    double minval = atof(argv[3]);
    double maxval = atof(argv[4]);
    std::cout << "Integral of Histogram " << histname << " from file " << filename << " from " << minval << " to " << maxval << ": " << histogram->Integral(histogram->FindBin(minval), histogram->FindBin(maxval)) << std::endl;
    
    
    std::cout << "ending" << std::endl;
    return(EXIT_SUCCESS);
}
