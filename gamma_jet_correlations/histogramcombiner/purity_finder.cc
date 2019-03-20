// This program outputs the purity of a given data set, in the form of a ROOT output file
// Requires that the data posesses a TOT_clusterpt histogram, which contains a distribution of all of the purities of the data sample, and that the data's cluster pT is between 12.5 and 40 GeV
// Author: Ivan Chernyshev; Date: 10/16/2018

#include <TH1.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>

const int num_of_datatypes = 2;
const int num_of_id_vars = 2;
//enum photon_IDVARS {LAMBDA_0, DNN, EMAX_OVER_ECLUSTER};
//photon_IDVARS ID_possibilities[num_of_id_vars] = {LAMBDA_0, DNN, EMAX_OVER_ECLUSTER}
const std::string datatype_identifiers[num_of_datatypes] = {"pp", "pPb"};
const std::string Idvar_identifiers[num_of_id_vars] = {"lambda_0", "DNN"};

// Cluster pT intervals looked at
const int num_of_pt_intervals = 8;
const double cluspT_intervals[num_of_pt_intervals + 1] = {12.5, 13.5, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0};
const double purities[num_of_datatypes][num_of_id_vars][num_of_pt_intervals] = {{{22.7, 30.7, 27.6, 39.6, 48.6, 47.1, 49.1, 51.1}, {24.8, 31.8, 33.1, 41.1, 44.1, 45.1, 46.1, 37.1}}, {{23.7, 26.7, 28.6, 38.6, 42.6, 44.1, 50.1, 53.1}, {25.8, 29.8, 32.1, 37.1, 41.1, 47.1, 49.1, 46.1}}};
const double purityerrors[num_of_datatypes][num_of_id_vars][num_of_pt_intervals] = {{{2, 3, 2, 2, 2, 2, 4, 3}, {2, 2, 2, 2, 2, 2, 1, 4}}, {{2, 2, 1, 2, 2, 1, 2, 2}, {1, 2, 2, 1, 2, 1, 2, 2}}};

/*
const double purities_pp[num_of_id_vars][num_of_pt_intervals + 1] = {{27, 35, 33, 45, 54, 55, 57, 59}, {34, 41, 42, 50, 53, 56, 57, 48}, {30, 36, 36, 44, 46, 46, 53, 45}}
const double purities_pPb[num_of_id_vars][num_of_pt_intervals + 1] = {{28, 31, 34, 44, 48, 52, 58, 61}, {35, 39, 41, 46, 50, 58, 60, 57}, {30, 32, 33, 41, 41, 49, 55, 61}}
const double purityerrors_pp = {{2, 3, 2, 2, 2, 2, 4, 3}, {2, 2, 2, 2, 2, 2, 1, 4}, {2, 3, 2, 2, 2, 3, 3, 6}};
const double purityerrors_pPb = {{2, 2, 1, 2, 2, 1, 2, 2}, {1, 2, 2, 1, 2, 1, 2, 2}, {1 ,1, 1, 2, 2, 2, 3, 3}};
*/
 
/*
const double purities_DNN_pp[num_of_pt_intervals] = {34, 41, 42, 50, 53, 56, 57, 48};
const double purities_Lambda_pp[num_of_pt_intervals] = {27, 35, 33, 45, 54, 55, 57, 59};
const double purities_Emax_Ecluster_pp[num_of_pt_intervals] = {30, 36, 36, 44, 46, 46, 53, 45};
const double purities_DNN_pPb[num_of_pt_intervals] = {35, 39, 41, 46, 50, 58, 60, 57};
const double purities_Lambda_pPb[num_of_pt_intervals] = {28, 31, 34, 44, 48, 52, 58, 61};
const double purities_Emax_Ecluster_pPb[num_of_pt_intervals] = {30, 32, 33, 41, 41, 49, 55, 61};

const double purityerrors_DNN_pp[num_of_pt_intervals] = {2, 2, 2, 2, 2, 2, 1, 4};
const double purityerrors_Lambda_pp[num_of_pt_intervals] = {2, 3, 2, 2, 2, 2, 4, 3};
const double purityerrors_Emax_Ecluster_pp[num_of_pt_intervals] = {2, 3, 2, 2, 2, 3, 3, 6};
const double purityerrors_DNN_pPb[num_of_pt_intervals] = {1, 2, 2, 1, 2, 1, 2, 2};
const double purityerrors_Lambda_pPb[num_of_pt_intervals] = {2, 2, 1, 2, 2, 1, 2, 2};
const double purityerrors_Emax_Ecluster_pPb[num_of_pt_intervals] = {1 ,1, 1, 2, 2, 2, 3, 3};
*/

// The cutting function
void SetCut(TH1* h, double min, double max){
    //make a selection on the chosen variable
    double width = h->GetXaxis()->GetBinWidth(1);
    int binmin = h->GetXaxis()->FindBin(min);
    int binmax = h->GetXaxis()->FindBin(max);
    h->GetXaxis()->SetRange(binmin, binmax - 1);
    return;
}

// Integration function, based on the Integrate.cc function
double integrate(std::string filename, std::string histname, double minval, double maxval) {
    TFile* inputfile = new TFile(filename.c_str(), "READ");
    if (inputfile == NULL){
        std::cout << "Error: file " << filename << " not found" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Get histogram
    TH1* histogram = 0;
    inputfile->GetObject(histname.c_str(), histogram);
    if (histogram == NULL){
        std::cout << "Error: histogram " << histname << " not found" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //Get the minimum and maximum indices and integrate
    double integral_result = histogram->Integral(histogram->FindBin(minval), histogram->FindBin(maxval));
    std::cout << "Integral of Histogram " << histname << " from file " << filename << " from " << minval << " to " << maxval << ": " << integral_result << std::endl;
    
    
    return(integral_result);
}

int main(int argc, char *argv[]) {
    if(argc < 6) {
        std::cout << "Error: syntax is ./Integrate [ROOT filename] [Cluster pT min] [Cluster pT max] [Type of data (pp or pPb)] [Shower-shape selection variable]" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Get filename
    std::string filename = (std::string)argv[1];
    TFile* inputfile = new TFile(filename.c_str(), "READ");
    if (inputfile == NULL){
        std::cout << "Error: file " << filename << " not found" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Get the histogram for the cluster pT
    TH1D* cluspT_hist = 0;
    inputfile->GetObject("TOT_clusterpt", cluspT_hist);
    if (cluspT_hist == NULL){
        std::cout << "Error: histogram TOT_clusterpt not found" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Set the data type used
    std::string name_of_datatype = (std::string)argv[4];
    int datatype_index = -1;
    for (int i = 0; i < num_of_datatypes; i++) {
        if(name_of_datatype == datatype_identifiers[i]) {
            datatype_index = i;
            break;
        }
    }
    if (datatype_index == -1) {
        std::cout << "ERROR: Data type selection variable must be: \""  << datatype_identifiers[0] << "\", or \"" << datatype_identifiers[1] << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Set the variable specifying the shower-shape variable used
    std::string name_of_showershapevar = (std::string)argv[5];
    int showershape_index = -1;
    for (int i = 0; i < num_of_id_vars; i++) {
        if(name_of_showershapevar == Idvar_identifiers[i]) {
            showershape_index = i;
            break;
        }
    }
    if (showershape_index == -1) {
        std::cout << "ERROR: Shower shape selection variable must be: \"" << Idvar_identifiers[0] << "\" or \"" << Idvar_identifiers[2] << std::endl;
        exit(EXIT_FAILURE);
    }
        
    // Caluculate the purity and its uncertainty
    double total_purity = 0;
    double total_delta_purity = 0;
    double minval = atof(argv[2]);
    double maxval = atof(argv[3]);
    
    int i = 1;
    
    // Quit the program if the minimum cluster pT submitted is smaller than the minimum of the pT intervals or larger than the maximum of the pT intervals
    if (minval < cluspT_intervals[0]) {
        std::cout << "ERROR: minimum cluster p_{T} value must be at least " << cluspT_intervals[0] << " GeV" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (maxval > cluspT_intervals[num_of_pt_intervals]) {
        std::cout << "ERROR: maximum cluster p_{T} value must be at most " << cluspT_intervals[num_of_pt_intervals] << " GeV" << std::endl;
        exit(EXIT_FAILURE);
    }
    // A failsafe, in case submitted min pt is not smaller than submitted max pt
    if(minval >= maxval){
        std::cout << "ERROR: minimum cluster pT must be smaller than maximum cluster pT" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Iterate over all pT intervals that at least partially overlap with the specified datarange, adding up purities and purity uncertainties (latter is added in quadrature)
    while(minval > cluspT_intervals[i])
        i++;
    
    if (cluspT_intervals[i] < maxval) {
        double num_in_interval = integrate(filename, "TOT_clusterpt", minval, cluspT_intervals[i]);
        total_purity += (num_in_interval*purities[datatype_index][showershape_index][i-1]);
        total_delta_purity += (num_in_interval*purityerrors[datatype_index][showershape_index][i-1])*(num_in_interval*purityerrors[datatype_index][showershape_index][i-1]);
        i++;
        
        while(cluspT_intervals[i] < maxval) {
            num_in_interval = integrate(filename, "TOT_clusterpt", cluspT_intervals[i - 1], cluspT_intervals[i]);
            total_purity += (num_in_interval*purities[datatype_index][showershape_index][i-1]);
            total_delta_purity += (num_in_interval*purityerrors[datatype_index][showershape_index][i-1])*(num_in_interval*purityerrors[datatype_index][showershape_index][i-1]);
            i++;
        }
        num_in_interval = integrate(filename, "TOT_clusterpt", cluspT_intervals[i - 1], maxval);
        total_purity += (num_in_interval*purities[datatype_index][showershape_index][i-1]);
        total_delta_purity += (num_in_interval*purityerrors[datatype_index][showershape_index][i-1])*(num_in_interval*purityerrors[datatype_index][showershape_index][i-1]);
        i++;
    }
    else {
        double num_in_interval = integrate(filename, "TOT_clusterpt", minval, maxval);
        total_purity += (num_in_interval*purities[datatype_index][showershape_index][i-1]);
        total_delta_purity += (num_in_interval*purityerrors[datatype_index][showershape_index][i-1])*(num_in_interval*purityerrors[datatype_index][showershape_index][i-1]);
        i++;
    }
    
    // Now, normalize the total purity and the error (and take a square root of the latter before doing so, because we're adding in quadrature)
    total_delta_purity = TMath::Sqrt(total_delta_purity);
    double num_total = integrate(filename, "TOT_clusterpt", minval, maxval);
    total_purity = total_purity/num_total;
    total_delta_purity = total_delta_purity/num_total;
    
    // Now print out the results
    std::cout << "Resultant purity for dataset " << filename << " with " << name_of_datatype << " data processed by the " << name_of_showershapevar << " shower-shape cut over cluster pT range " << minval << " to " << maxval << " is: " << total_purity << " +/- " << total_delta_purity << std::endl;
    
    std::cout << "ending" << std::endl;
    return(EXIT_SUCCESS);
}

