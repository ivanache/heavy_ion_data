/**
   This program produces energy response plots from Monte-Carlo simulations
*/

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
#include <iostream>
#include <fstream>
#include <TGraphAsymmErrors.h>

#define NTRACK_MAX (1U << 15)

#include <vector>
#include <math.h>
#include <set>

const int MAX_INPUT_LENGTH = 200;

// Energy of lead, in GeV
const double EPb = 1560;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};
enum photon_IDVARS {LAMBDA_0, DNN, EMAX_OVER_ECLUSTER};

// Function to calculate the bin width of a TH1 graph
double calculatebinwidth(int numofbins, double binmin, double binmax){
    return (binmax - binmin)/numofbins;
}

int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }

    // Read configuration file for the cut variables
    FILE* config = fopen("GammaJet_config.yaml", "r");
    if (config == NULL)  std::cout<<"no config"<<std::endl;
    // Default values of various variables used in the file (actual values are to be determined by the configuration file)
    // Cut variables
    double primary_vertex_max = 10.0;
    double SIG_DNN_min = 0.55;
    double SIG_DNN_max = 0.85;
    double BKG_DNN_min = 0.0;
    double BKG_DNN_max = 0.3;
    double SIG_lambda_min = 0.0;
    double SIG_lambda_max = 0.4;
    double BKG_lambda_min = 0.5;
    double BKG_lambda_max = 2.0;
    double SIG_Emax_over_Ecluster_min = 0.0;
    double SIG_Emax_over_Ecluster_max = 0.5;
    double BKG_Emax_over_Ecluster_min = 0.5;
    double BKG_Emax_over_Ecluster_max = 1.0;
    double clus_pT_min = 10;
    double clus_pT_max = 16;
    double track_pT_max = 1;
    double jet_pT_min = 10.0;
    double Eta_max = 0.5;
    double Cluster_ncell_min = 2;
    double Cluster_locmaxima_max = 2.0;
    double Cluster_distobadchannel = 2.0;
    double EcrossoverE_min = 0.05;
    
    // The bounds for the events to fal into the isolation and nonisolation areas
    double iso_max = 1.0;
    double noniso_min = 2.0;
    double noniso_max = 10.0;
    
    // Delta eta
    double deta_max = 0.5;
    
    // Number of bins in correlation functions
    int xjbins = 10;
    int phibins = 5;
    int etabins = 20;
    
    // Which branch should be used to determine whether a cluster should fall into iso, noniso, or neither
    isolationDet determiner = CLUSTER_ISO_ITS_04;
    photon_IDVARS photon_identifier = LAMBDA_0;
    
    // Truth cuts
    int rightpdgcode = 22;
    int rightparentpdgcode = 22;
    
    // Number of events
    int nevents = 0;
    // Loop through config file
    char line[MAX_INPUT_LENGTH];
    while (fgets(line, MAX_INPUT_LENGTH, config) != NULL) {
        if (line[0] == '#') {
            continue;
        }
        
        // Declare char arrays needed to read the line
        char key[MAX_INPUT_LENGTH];
        char dummy[MAX_INPUT_LENGTH];
        char value[MAX_INPUT_LENGTH];
        
        // Cap off key[0] and value[0] with null characters and load the key, dummy-characters, and value of the line into their respective arrays
        key[0] = '\0';
        value[0] = '\0';
        sscanf(line, "%[^:]:%[ \t]%100[^\n]", key, dummy, value);
        
        // Use if statements to detect, based on key, which variable the line's content should be used to fill and fill that variable
        if (strcmp(key, "primary_vertex_max") == 0) {
            // Assign primary_vertex_max to the double-converted version of value
            primary_vertex_max = atof(value);
            std::cout << "primary_vertex_max is " << primary_vertex_max << std::endl;
        }
        else if (strcmp(key, "SIG_DNN_min") == 0) {
            // Assign SIG_DNN_min to the double-converted version of value
            SIG_DNN_min = atof(value);
            std::cout << "SIG_DNN_min is " << SIG_DNN_min << std::endl;
        }
        else if (strcmp(key, "SIG_DNN_max") == 0) {
            // Assign SIG_DNN_max to the double-converted version of value
            SIG_DNN_max = atof(value);
            std::cout << "SIG_DNN_max is " << SIG_DNN_max << std::endl;
        }
        else if (strcmp(key, "BKG_DNN_min") == 0) {
            // Assign BKG_DNN_min to the double-converted version of value
            BKG_DNN_min = atof(value);
            std::cout << "BKG_DNN_min is " << BKG_DNN_min << std::endl;
        }
        else if (strcmp(key, "BKG_DNN_max") == 0) {
            // Assign BKG_DNN_max to the double-converted version of value
            BKG_DNN_max = atof(value);
            std::cout << "BKG_DNN_max is " << BKG_DNN_max << std::endl;
        }
        else if (strcmp(key, "SIG_lambda_min") == 0) {
            // Assign SIG_lambda_min to the double-converted version of value
            SIG_lambda_min = atof(value);
            std::cout << "SIG_lambda_min is " << SIG_lambda_min << std::endl;
        }
        else if (strcmp(key, "SIG_lambda_max") == 0) {
            // Assign SIG_lambda_max to the double-converted version of value
            SIG_lambda_max = atof(value);
            std::cout << "SIG_lambda_max is " << SIG_lambda_max << std::endl;
        }
        else if (strcmp(key, "BKG_lambda_min") == 0) {
            // Assign BKG_lambda_min to the double-converted version of value
            BKG_lambda_min = atof(value);
            std::cout << "BKG_lambda_min is " << BKG_lambda_min << std::endl;
        }
        else if (strcmp(key, "BKG_lambda_max") == 0) {
            // Assign BKG_lambda_max to the double-converted version of value
            BKG_lambda_max = atof(value);
            std::cout << "BKG_lambda_max is " << BKG_lambda_max << std::endl;
        }
        else if (strcmp(key, "SIG_Emax_over_Ecluster_min") == 0) {
            // Assign SIG_lambda_min to the double-converted version of value
            SIG_Emax_over_Ecluster_min = atof(value);
            std::cout << "SIG_Emax_over_Ecluster_min is " << SIG_Emax_over_Ecluster_min << std::endl;
        }
        else if (strcmp(key, "SIG_Emax_over_Ecluster_max") == 0) {
            // Assign SIG_lambda_max to the double-converted version of value
            SIG_Emax_over_Ecluster_max = atof(value);
            std::cout << "SIG_Emax_over_Ecluster_max is " << SIG_Emax_over_Ecluster_max << std::endl;
        }
        else if (strcmp(key, "BKG_Emax_over_Ecluster_min") == 0) {
            // Assign BKG_lambda_min to the double-converted version of value
            BKG_Emax_over_Ecluster_min = atof(value);
            std::cout << "BKG_Emax_over_Ecluster_min is " << BKG_Emax_over_Ecluster_min << std::endl;
        }
        else if (strcmp(key, "BKG_Emax_over_Ecluster_max") == 0) {
            // Assign BKG_lambda_max to the double-converted version of value
            BKG_Emax_over_Ecluster_max = atof(value);
            std::cout << "BKG_Emax_over_Ecluster_max is " << BKG_Emax_over_Ecluster_max << std::endl;
        }
        else if (strcmp(key, "clus_pT_min") == 0) {
            clus_pT_min = atof(value);
            std::cout << "clus_pT_min is " << clus_pT_min << std::endl;
        }
        else if (strcmp(key, "clus_pT_max") == 0) {
            clus_pT_max = atof(value);
            std::cout << "clus_pT_max is " << clus_pT_max << std::endl;
        }
        else if (strcmp(key, "track_pT_max") == 0) {
            track_pT_max = atof(value);
            std::cout << "track_pT_max is " << track_pT_max << std::endl;
        }
        else if (strcmp(key, "jet_pT_min") == 0) {
            jet_pT_min = atof(value);
            std::cout << "jet_pT_min is " << jet_pT_min << std::endl;
        }
        else if (strcmp(key, "Eta_max") == 0) {
            Eta_max = atof(value);
            std::cout << "Eta_max is " << Eta_max << std::endl;
        }
        else if (strcmp(key, "Cluster_ncell_min") == 0) {
            Cluster_ncell_min = atof(value);
            std::cout << "Cluster_ncell_min is " << Cluster_ncell_min << std::endl;
        }
        else if (strcmp(key, "Cluster_locmaxima_max") == 0) {
            Cluster_locmaxima_max = atof(value);
            std::cout << "Cluster_locmaxima_max is " << Cluster_locmaxima_max << std::endl;
        }
        else if (strcmp(key, "Cluster_distobadchannel") == 0) {
            Cluster_distobadchannel = atof(value);
            std::cout << "Cluster_distobadchannel is " << Cluster_distobadchannel << std::endl;
        }
        else if (strcmp(key, "EcrossoverE_min") == 0) {
            EcrossoverE_min = atof(value);
            std::cout << "EcrossoverE_min is " << EcrossoverE_min << std::endl;
        }
        else if (strcmp(key, "iso_max") == 0) {
            iso_max = atof(value);
            std::cout << "iso_max is " << iso_max << std::endl;
        }
        else if (strcmp(key, "noniso_min") == 0) {
            noniso_min = atof(value);
            std::cout << "noniso_min is " << noniso_min << std::endl;
        }
        else if (strcmp(key, "noniso_max") == 0) {
            noniso_max = atof(value);
            std::cout << "noniso_max is " << noniso_max << std::endl;
        }
        else if (strcmp(key, "deta_max") == 0) {
            deta_max = atof(value);
            std::cout << "deta_max is " << deta_max << std::endl;
        }
        else if (strcmp(key, "phi_func_bins") == 0) {
            phibins = atoi(value);
            std::cout << "Bins in a phi function: " << phibins << std::endl;
        }
        else if (strcmp(key, "eta_func_bins") == 0) {
            etabins = atoi(value);
            std::cout << "Bins in an eta function: " << etabins << std::endl;
        }
        else if (strcmp(key, "xj_func_bins") == 0) {
            xjbins = atoi(value);
            std::cout << "Bins in an xj function: " << xjbins << std::endl;
        }
        else if (strcmp(key, "photon_idvar") == 0) {
            if (strcmp(value, "lambda_0") == 0){
                photon_identifier = LAMBDA_0;
                std::cout << "lambda_0 will determine photon selection" << std::endl;
            }
            else if (strcmp(value, "DNN") == 0){
                photon_identifier = DNN;
                std::cout << "Deep Neural Net will determine photon selection" << std::endl;
            }
            else if (strcmp(value, "Emax_over_Ecluster") == 0){
                photon_identifier = EMAX_OVER_ECLUSTER;
                std::cout << "#frac{E_{max}}{E_{cluster}} will determine photon selection" << std::endl;
            }
            else {
                std::cout << "ERROR: Photon selection determinant in configuration file must be \"lambda_0\", \"DNN\", or \"Emax_over_Ecluster\"" << std::endl << "Aborting the program" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(key, "Cluster_isolation_determinant") == 0) {
            if (strcmp(value, "cluster_iso_tpc_04") == 0){
                determiner = CLUSTER_ISO_TPC_04;
                std::cout << "cluster_iso_tpc_04 will determine the isolation and non-isolation placement" << std::endl;
            }
            else if (strcmp(value, "cluster_iso_its_04") == 0){
                determiner = CLUSTER_ISO_ITS_04;
                std::cout << "cluster_iso_its_04 will determine the isolation and non-isolation placement" << std::endl;
            }
            else if (strcmp(value, "cluster_frixione_tpc_04_02") == 0){
                determiner = CLUSTER_FRIXIONE_TPC_04_02;
                std::cout << "cluster_frixione_tpc_04_02 will determine the isolation and non-isolation placement" << std::endl;
            }
            else if (strcmp(value, "cluster_frixione_its_04_02") == 0){
                determiner = CLUSTER_FRIXIONE_ITS_04_02;
                std::cout << "cluster_frixione_its_04_02 will determine the isolation and non-isolation placement" << std::endl;
            }
            else {
                std::cout << "ERROR: Cluster_isolation_determinant in configuration file must be \"cluster_iso_tpc_04\", \"cluster_iso_its_04\", \"cluster_frixione_tpc_04_02\", or \"cluster_frixione_its_04_02\"" << std::endl << "Aborting the program" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(key, "pdg_code") == 0) {
            rightpdgcode = atoi(value);
            std::cout << "Right pdg_code: " << rightpdgcode << std::endl;
        }
        else if (strcmp(key, "parent_pdg_code") == 0) {
            rightparentpdgcode = atoi(value);
            std::cout << "Right parent_pdg_code: " << rightparentpdgcode << std::endl;
        }
        else if (strcmp(key, "Num_events") == 0) {
            nevents = atoi(value);
            std::cout << "Num_events: " << nevents << std::endl;
        }
        else {
            std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
        }
    }
    fclose(config);
  
  int dummyc = 1;
  char **dummyv = new char *[1];


  std::cout << " Number of events requested " << nevents << std::endl;
  std::cout << " ptmin " << clus_pT_min << " ptmax = " << clus_pT_max << std::endl;
  std::cout << "minimum pt jet " << jet_pT_min<< std::endl;

  dummyv[0] = strdup("main");
  TApplication application("", &dummyc, dummyv);    
  std::cout <<" Number of arguments " << argc << std::endl; 


  TH1D h_zvertex("h_zvertex","vertex z " , 100, -20.0, 20.0);
  double h_zvertex_binwidth = calculatebinwidth(100, -20.0, 20.0);

  TH1D h_cutflow("h_cutflow","cut flow for photons", 10, -0.5,9.5);
  TH1D h_evtcutflow("h_evtcutflow","Event cut flow", 7, -0.5,6.5);
  TH1D h_trkcutflow("h_trkcutflow","Track cut flow", 4, -0.5,3.5);
  TH1D h_jetcutflow("h_jetcutflow","Jet cut flow", 6, -0.5,5.5);
  TH1D h_evt_rho("h_evt_rho", "average UE density, rho", 100, 0, 10.0);
    double h_cutflow_binwidth = calculatebinwidth(10, -0.5,9.5);
    double h_evtcutflow_binwidth = calculatebinwidth(7, -0.5,6.5);
    double h_evt_rho_binwidth = calculatebinwidth(100, 0, 10.0);
   
  TH1D h_evt_rhoITS("h_evt_rhoITS", "average UE density, rho, from ITS", 100, 0, 10.0);
  TH1D h_evt_rhoTPC("h_evt_rhoTPC", "average UE density, rho, from TPC", 100, 0, 10.0);
    double h_evt_rhoITS_binwidth = calculatebinwidth(100, 0, 10.0);
    double h_evt_rhoTPC_binwidth = calculatebinwidth(100, 0, 10.0);
   
  TH1D h_reco("h_reco", "reco photons filled with pt reco", 50, 0, 50);    
  TH1D h_reco_truthpt("h_reco_truthpt", "reco photons filled with truthpt reco", 50, 0, 50); 
  TH1D h_truth("h_truth", "truth photons", 50, 0, 50);
    double h_reco_binwidth = calculatebinwidth(50, 0, 50);
    double h_reco_truthpt_binwidth = calculatebinwidth(50, 0, 50);
    double h_truth_binwidth = calculatebinwidth(50, 0, 50);

  TH1D hSR_njet("hSR_njet", "Number of associated jets, signal region" , 5, -0.5, 4.5);
  TH1D hBR_njet("hBR_njet", "Number of associated jets, bkg region" , 5, -0.5, 4.5);
  double hSR_njet_binwidth = calculatebinwidth(5, -0.5, 4.5);
  double hBR_njet_binwidth = calculatebinwidth(5, -0.5, 4.5);

  TH1D TOT_clusterpt("TOT_clusterpt", "Isolated cluster pt [GeV], all data", 80, 10.0, 30.0);
  TH1D hBR_clusterpt("hBR_clusterpt", "Isolated cluster pt [GeV], bkg region", 80, 10.0, 30.0);
  TH1D hSR_clusterpt("hSR_clusterpt", "Isolated cluster pt [GeV], signal region", 80, 10.0, 30.0);
    double TOT_clusterpt_binwidth = calculatebinwidth(80, 10.0, 30.0);
    double hBR_clusterpt_binwidth = calculatebinwidth(80, 10.0, 30.0);
    double hSR_clusterpt_binwidth = calculatebinwidth(80, 10.0, 30.0);

  TH1D hBR_clustereta("hBR_clustereta", "Isolated cluster eta, bkg region", 40, -1.0, 1.0);
  TH1D hSR_clustereta("hSR_clustereta", "Isolated cluster eta, signal region", 40, -1.0, 1.0);
  TH1D hBR_clusterphi("hBR_clusterphi", "Isolated cluster phi, bkg region", 40, -1.0*TMath::Pi(), TMath::Pi());
  TH1D hSR_clusterphi("hSR_clusterphi", "Isolated cluster phi, signal region", 40, -1.0*TMath::Pi(), TMath::Pi());
    double hSR_clustereta_binwidth = calculatebinwidth(40, -1.0, 1.0);
    double hBR_clustereta_binwidth = calculatebinwidth(40, -1.0, 1.0);
    double hSR_clusterphi_binwidth = calculatebinwidth(40, -1.0*TMath::Pi(), TMath::Pi());
    double hBR_clusterphi_binwidth = calculatebinwidth(40, -1.0*TMath::Pi(), TMath::Pi());

  TH1D h_clustereta("h_clustereta", "all photon cluster eta", 40, -1.0, 1.0);
  TH1D h_clustereta_iso("h_clustereta_iso", "all photon cluster eta, passing isolation", 40, -1.0, 1.0);
    double h_clustereta_binwidth = calculatebinwidth(40, -1.0, 1.0);
    double h_clustereta_iso_binwidth = calculatebinwidth(40, -1.0, 1.0);

  TH1D h_clusterphi("h_clusterphi", "all photon cluster phi", 100, -1.0*TMath::Pi(), TMath::Pi());
  TH1D h_clusterphi_iso("h_clusterphi_iso", "all photon cluster phi, passing isolation", 100, -1.0*TMath::Pi(), TMath::Pi());
  TH1D h_trackphi("h_trackphi", " track phi" , 100,  -1.0*TMath::Pi(), TMath::Pi());
  TH1D h_jetphi("h_jetphi", " jetphi phi" , 100,  -1.0*TMath::Pi(), TMath::Pi());
    double h_clusterphi_binwidth = calculatebinwidth(100, -1.0*TMath::Pi(), TMath::Pi());
    double h_clusterphi_iso_binwidth = calculatebinwidth(100, -1.0*TMath::Pi(), TMath::Pi());
    double h_trackphi_binwidth = calculatebinwidth(100, -1.0*TMath::Pi(), TMath::Pi());
    double h_jetphi_binwidth = calculatebinwidth(100, -1.0*TMath::Pi(), TMath::Pi());

    
  TH1D TOT_jetpt("TOT_jetpt",   "Associated jet pt spectrum (reco), all data", 30, 0, 30);
  TH1D hBR_jetpt("hBR_jetpt",   "Associated jet pt spectrum (reco), bkg region", 30, 0, 30);
  TH1D hSR_jetpt("hSR_jetpt",   "Associated jet pt spectrum (reco), signal region", 30, 0, 30);
  TH1D hBR_jeteta("hBR_jeteta", "Associated jet eta spectrum (reco), bkg region", 20, -1.0, 1.0);
  TH1D hSR_jeteta("hSR_jeteta", "Associated jet eta spectrum (reco), signal region", 20, -1.0, 1.0); 
  TH1D hBR_jetphi("hBR_jetphi", "Associated jet phi spectrum (reco), bkg region", 20, -1.0*TMath::Pi(), TMath::Pi());
  TH1D hSR_jetphi("hSR_jetphi", "Associated jet phi spectrum (reco), signal region", 20, -1.0*TMath::Pi(), TMath::Pi());
    double TOT_jetpt_binwidth = calculatebinwidth(30, 0, 30);
    double hBR_jetpt_binwidth = calculatebinwidth(30, 0, 30);
    double hSR_jetpt_binwidth = calculatebinwidth(30, 0, 30);
    double hBR_jeteta_binwidth = calculatebinwidth(20, -1.0, 1.0);
    double hSR_jeteta_binwidth = calculatebinwidth(20, -1.0, 1.0);
    double hBR_jetphi_binwidth = calculatebinwidth(20, -1.0*TMath::Pi(), TMath::Pi());
    double hSR_jetphi_binwidth = calculatebinwidth(20, -1.0*TMath::Pi(), TMath::Pi());


  TH1D hBR_jetpt_truth("hBR_jetpt_truth",   "Associated jet pt spectrum (truth), bkg region", 30, 0, 30);
  TH1D hSR_jetpt_truth("hSR_jetpt_truth",   "Associated jet pt spectrum (truth), signal region", 30, 0, 30);
  TH1D hBR_jeteta_truth("hBR_jeteta_truth", "Associated jet eta spectrum (truth), bkg region", 20, -1.0, 1.0);
  TH1D hSR_jeteta_truth("hSR_jeteta_truth", "Associated jet eta spectrum (truth), signal region", 20, -1.0, 1.0);
  TH1D hBR_jetphi_truth("hBR_jetphi_truth", "Associated jet phi spectrum (truth), bkg region", 20, -1.0*TMath::Pi(), TMath::Pi());
  TH1D hSR_jetphi_truth("hSR_jetphi_truth", "Associated jet phi spectrum (truth), signal region", 20, -1.0*TMath::Pi(), TMath::Pi());
    double hBR_jetpt_truth_binwidth = calculatebinwidth(30, 0, 30);
    double hSR_jetpt_truth_binwidth = calculatebinwidth(30, 0, 30);
    double hBR_jeteta_truth_binwidth = calculatebinwidth(20, -1.0, 1.0);
    double hSR_jeteta_truth_binwidth = calculatebinwidth(20, -1.0, 1.0);
    double hBR_jetphi_truth_binwidth = calculatebinwidth(20, -1.0*TMath::Pi(), TMath::Pi());
    double hSR_jetphi_truth_binwidth = calculatebinwidth(20, -1.0*TMath::Pi(), TMath::Pi());

  TH1D h_jetpt_truth("h_jetpt_truth", "truth jet pt", 30, 0, 30);
  TH1D h_jetpt_truthreco("h_jetpt_truthreco", "reco jet truth pt (numerator of efficiency)", 30, 0, 30);
  TH1D h_jetpt_reco("h_jetpt_reco", "reco jet reco pt", 30, 0, 30);
    double h_jetpt_truth_binwidth = calculatebinwidth(30, 0, 30);
    double h_jetpt_truthreco_binwidth = calculatebinwidth(30, 0, 30);
    double h_jetpt_reco_binwidth = calculatebinwidth(30, 0, 30);
 

  TH1D hSR_Xj("hSR_Xj", "Xj distribution, Signal region", xjbins, 0.0,2.0);
  TH1D hBR_Xj("hBR_Xj", "Xj distribution, BKG region", xjbins, 0.0,2.0);
  TH1D hSR_pTD("hSR_pTD", "pTD distribution, Signal region", 5, 0.0,1.0);
  TH1D hBR_pTD("hBR_pTD", "pTD distribution, BKG region", 5, 0.0,1.0);
  TH1D hSR_Multiplicity("hSR_Multiplicity", "Jet Multiplicity distribution, Signal region", 10, 0.0 , 20.0);
  TH1D hBR_Multiplicity("hBR_Multiplicity", "Jet Multiplicity distribution, BKG region", 10, 0.0, 20.0);
  TH1D hSR_jetwidth("hSR_jetwidth", "jet width distribution, Signal region", 20, -10, 0);
  TH1D hBR_jetwidth("hBR_jetwidth", "jet width distribution, BKG region", 20, -10, 0);
    double hSR_Xj_binwidth = calculatebinwidth(xjbins, 0.0,2.0);
    double hBR_Xj_binwidth = calculatebinwidth(xjbins, 0.0,2.0);
    double hSR_pTD_binwidth = calculatebinwidth(5, 0.0,1.0);
    double hBR_pTD_binwidth = calculatebinwidth(5, 0.0,1.0);
    double hSR_Multiplicity_binwidth = calculatebinwidth(20, 0.0, 20.0);
    double hBR_Multiplicity_binwidth = calculatebinwidth(20, 0.0, 20.0);
    double hSR_jetwidth_binwidth = calculatebinwidth(20, -10, 0);
    double hBR_jetwidth_binwidth = calculatebinwidth(20, -10, 0);


  TH1D hSR_Xj_truth("hSR_Xj_truth", "True Xj distribution, Signal region", xjbins, 0.0,2.0);
  TH1D hBR_Xj_truth("hBR_Xj_truth", "True Xj distribution, BKG region",xjbins, 0.0,2.0);
  TH2D h_Xj_Matrix("h_Xj_Matrix", "Truth/Reco matrix", xjbins, 0.0,2.0, xjbins, 0.0,2.0);
  TH1D h_Xj_truth("h_Xj_truth", "Xj truth distribution", xjbins, 0.0,2.0);
    double hSR_Xj_truth_binwidth = calculatebinwidth(xjbins, 0.0,2.0);
    double hBR_Xj_truth_binwidth = calculatebinwidth(xjbins, 0.0,2.0);
    double h_Xj_truth_binwidth = calculatebinwidth(xjbins, 0.0,2.0);

  TH1D hSR_dPhi("hSR_dPhi", "delta phi gamma-jet signal region", phibins, 0, TMath::Pi());
  TH1D hBR_dPhi("hBR_dPhi", "delta phi gamma-jet background region", phibins, 0, TMath::Pi());
  TH1D hSR_dPhi_truth("hSR_dPhi_truth", "delta phi gamma-jet signal region, truth", phibins, 0, TMath::Pi());
  TH1D hBR_dPhi_truth("hBR_dPhi_truth", "delta phi gamma-jet background region, truth", phibins, 0, TMath::Pi());
    double hSR_dPhi_binwidth = calculatebinwidth(phibins, 0, TMath::Pi());
    double hBR_dPhi_binwidth = calculatebinwidth(phibins, 0, TMath::Pi());
    double hSR_dPhi_truth_binwidth = calculatebinwidth(phibins, 0, TMath::Pi());
    double hBR_dPhi_truth_binwidth = calculatebinwidth(phibins, 0, TMath::Pi());


  TH1D hSR_dEta("hSR_dEta", "delta eta gamma-jet signal region", etabins, -1.2, 1.2);
  TH1D hBR_dEta("hBR_dEta", "delta eta gamma-jet background region", etabins, -1.2, 1.2);
  TH1D hSR_dEta_truth("hSR_dEta_truth", "delta eta gamma-jet signal region, truth", etabins, -1.2, 1.2);
  TH1D hBR_dEta_truth("hBR_dEta_truth", "delta eta gamma-jet background region, truth", etabins, -1.2, 1.2);
    double hSR_dEta_binwidth = calculatebinwidth(etabins, -1.2, 1.2);
    double hBR_dEta_binwidth = calculatebinwidth(etabins, -1.2, 1.2);
    double hSR_dEta_truth_binwidth = calculatebinwidth(etabins, -1.2, 1.2);
    double hBR_dEta_truth_binwidth = calculatebinwidth(etabins, -1.2, 1.2);

  TH1D hSR_AvgEta("hSR_AvgEta", "Average eta gamma-jet signal region", 2*etabins, -1.2, 1.2);
  TH1D hBR_AvgEta("hBR_AvgEta", "Average eta gamma-jet background region", 2*etabins, -1.2, 1.2);
  TH1D hSR_AvgEta_truth("hSR_AvgEta_truth", "Average eta gamma-jet signal region, truth", 2*etabins, -1.2, 1.2);
  TH1D hBR_AvgEta_truth("hBR_AvgEta_truth", "Average eta gamma-jet background region, truth", 2*etabins, -1.2, 1.2);
    double hSR_AvgEta_binwidth = calculatebinwidth(2*etabins, -1.2, 1.2);
    double hBR_AvgEta_binwidth = calculatebinwidth(2*etabins, -1.2, 1.2);
    double hSR_AvgEta_truth_binwidth = calculatebinwidth(2*etabins, -1.2, 1.2);
    double hBR_AvgEta_truth_binwidth = calculatebinwidth(2*etabins, -1.2, 1.2);

  TH1D hSR_XobsPb("hSR_XobsPb", "x_{pPb}^{obs} distribution: signal region; x_{pPb}^{obs}; #frac{d #sigma}{dx^{obs}_{pPb}}", 7, 0.004, 0.024);
  TH1D hBR_XobsPb("hBR_XobsPb", "x_{pPb}^{obs} distribution: background region; x_{pPb}^{obs}; #frac{d #sigma}{dx^{obs}_{pPb}}", 7, 0.004, 0.024);
  TH1D hSR_XobsPb_truth("hSR_XobsPb_truth", "x_{pPb}^{obs} distribution: signal region, truth; x_{pPb}^{obs}; #frac{d #sigma}{dx^{obs}_{pPb}}", 7, 0.004, 0.024);
  TH1D hBR_XobsPb_truth("hBR_XobsPb_truth", "x_{pPb}^{obs} distribution: background region, truth; x_{pPb}^{obs}; #frac{d #sigma}{dx^{obs}_{pPb}}", 7, 0.004, 0.024);
    double hSR_XobsPb_binwidth = calculatebinwidth(7, 0.004, 0.024);
    double hBR_XobsPb_binwidth = calculatebinwidth(7, 0.004, 0.024);
    double hSR_XobsPb_truth_binwidth = calculatebinwidth(7, 0.004, 0.024);
    double hBR_XobsPb_truth_binwidth = calculatebinwidth(7, 0.004, 0.024);

  TH1D h_dPhi_truth("h_dPhi_truth", "delta phi gamma-jet truth MC", phibins, 0, TMath::Pi());
    double h_dPhi_truth_binwidth = calculatebinwidth(phibins, 0, TMath::Pi());


  TH1D h_weights("h_weights", "weights; bin1 is for SR and bin2 is for BR", 2, -0.5, 1.5);
    double h_weights_binwidth = calculatebinwidth(2, -0.5, 1.5);
    
    TH1D h_lambda_0("h_lambda_0", "#lambda_{0}^2 distribution; #lambda_{0}^2; # of photons", 70, 0, 1.4);
    TH1D h_DNN("h_DNN", "Deep Neural Net distribution; DNN; # of photons", 50, 0, 1.0);
    TH1D h_EmaxOverEcluster("h_EmaxOverEcluster", "#frac{E_{max}}{E_{cluster}} distribution; #frac{E_{max}}{E_{cluster}} ; # of photons", 50, 0, 1.0);
    

  h_clusterphi.Sumw2();
  h_clusterphi_iso.Sumw2();

  h_clustereta.Sumw2();
  h_clustereta_iso.Sumw2();

  h_trackphi.Sumw2();
  h_jetphi.Sumw2();

  hSR_pTD.Sumw2();
  hBR_pTD.Sumw2();
  hSR_Multiplicity.Sumw2();
  hBR_Multiplicity.Sumw2();
  hSR_jetwidth.Sumw2();
  hBR_jetwidth.Sumw2();

  h_evt_rho.Sumw2();
  hSR_Xj.Sumw2();
  hBR_Xj.Sumw2();
  hSR_Xj_truth.Sumw2();
  hBR_Xj_truth.Sumw2();
  h_Xj_truth.Sumw2();     

  hSR_njet.Sumw2();
  hBR_njet.Sumw2();

  hSR_dPhi.Sumw2();
  hBR_dPhi.Sumw2();
  hSR_dPhi_truth.Sumw2();
  hBR_dPhi_truth.Sumw2();

  hSR_dEta.Sumw2();
  hBR_dEta.Sumw2();
  hSR_dEta_truth.Sumw2();
  hBR_dEta_truth.Sumw2();

  hSR_AvgEta.Sumw2();
  hBR_AvgEta.Sumw2();
  hSR_AvgEta_truth.Sumw2();
  hBR_AvgEta_truth.Sumw2();


  h_jetpt_truth.Sumw2();
  h_jetpt_truthreco.Sumw2();
  h_jetpt_reco.Sumw2();

  h_dPhi_truth.Sumw2();
    
  TOT_clusterpt.Sumw2();
  hSR_clusterpt.Sumw2();
  hBR_clusterpt.Sumw2();
  hSR_clustereta.Sumw2();
  hBR_clustereta.Sumw2();
  hSR_clusterphi.Sumw2();
  hBR_clusterphi.Sumw2();

  TOT_jetpt.Sumw2();
  hBR_jetpt.Sumw2();
  hSR_jetpt.Sumw2();
  hBR_jeteta.Sumw2();
  hSR_jeteta.Sumw2();
  hBR_jetphi.Sumw2();
  hSR_jetphi.Sumw2();
    
  hBR_jetpt_truth.Sumw2();
  hSR_jetpt_truth.Sumw2();
  hBR_jeteta_truth.Sumw2();
  hSR_jeteta_truth.Sumw2();
  hBR_jetphi_truth.Sumw2();
  hSR_jetphi_truth.Sumw2();

  hSR_XobsPb.Sumw2();
  hBR_XobsPb.Sumw2();
  hSR_XobsPb_truth.Sumw2();
  hBR_XobsPb_truth.Sumw2();

  h_evt_rhoITS.Sumw2();
  h_evt_rhoTPC.Sumw2();

  h_Xj_Matrix.Sumw2();

  hSR_Xj.SetTitle("; X_{j} ; 1/N_{#gamma} dN_{J#gamma}/dX_{j}");
  hBR_Xj.SetTitle("; X_{j} ; 1/N_{#gamma} dN_{J#gamma}/dX_{j}");   
  h_Xj_truth.SetTitle("; X_{j}^{true} ; counts");

    
  TCanvas* canvas = new TCanvas();
  Float_t N_SR = 0; //float because it might be weighted in MC
  Float_t N_BR = 0;
  Float_t N_eventpassed = 0;
  Float_t N_truth = 0;
    
  TH1D hBR("hBR", "Isolated cluster, bkg region", 40, 10.0, 50.0);
  TH1D hweight("hweight", "Isolated cluster, signal region", 40, 10.0, 50.0);
  
    int num_sig_dPhi = 0;
    int num_bkg_dPhi = 0;
    int num_sig_XobsPb = 0;
    int num_bkg_XobsPb = 0;
 
  Bool_t isRealData = true;
  for (int iarg = 1; iarg < argc; iarg++) {
      std::string filestring = (std::string)argv[iarg];
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
    TFile *file = TFile::Open((TString)argv[iarg]);
    
    if (file == NULL) {
        std::cout << " fail; could not open file" << std::endl;
        exit(EXIT_FAILURE);
    }
    file->Print();
    
    // Get all the TTree variables from the file to open, I guess
    TTree *_tree_event = NULL;
    std::cout << " About to try getting the ttree" << std::endl;
    _tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
    if (_tree_event == NULL) {
        std::cout << "First try did not got trying again" << std::endl;
        _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
        if (_tree_event == NULL) {
            std::cout << " fail; could not find _tree_event " << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (_tree_event == NULL) {
        std::cout << " fail; the _tree_event is NULL " << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //
    
    std::cout <<"_tree_event->GetEntries() " << _tree_event->GetEntries() << std::endl;
    
    //you define variables
    Double_t primary_vertex[3];
    Bool_t is_pileup_from_spd_5_08;
    Bool_t is_pileup_from_spd_3_08;
    Float_t ue_estimate_its_const;
    Float_t ue_estimate_tpc_const;
   
    UInt_t ntrack;
    Float_t track_e[NTRACK_MAX];
    Float_t track_pt[NTRACK_MAX];
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    UChar_t track_quality[NTRACK_MAX];
        
    UInt_t ncluster;
    Float_t cluster_e[NTRACK_MAX];
    Float_t cluster_e_cross[NTRACK_MAX];
    Float_t cluster_e_max[NTRACK_MAX];
    Float_t cluster_pt[NTRACK_MAX];
    Float_t cluster_eta[NTRACK_MAX];
    Float_t cluster_phi[NTRACK_MAX];
    
    Float_t cluster_iso_tpc_04[NTRACK_MAX];
    Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
    
    Float_t cluster_iso_its_04[NTRACK_MAX];
    
    Float_t cluster_iso_its_04_ue[NTRACK_MAX];
  

    Float_t cluster_frixione_its_04_02[NTRACK_MAX];
    Float_t cluster_s_nphoton[NTRACK_MAX][4];
    UChar_t cluster_nlocal_maxima[NTRACK_MAX];
    Float_t cluster_distance_to_bad_channel[NTRACK_MAX];   
 
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    Int_t cluster_ncell[NTRACK_MAX];
    UShort_t  cluster_cell_id_max[NTRACK_MAX];
    Float_t cluster_lambda_square[NTRACK_MAX][2];
    Float_t cell_e[17664];

    //Jets reco 
    UInt_t njet_ak04its;
    Float_t jet_ak04its_pt_raw[NTRACK_MAX];
    Float_t jet_ak04its_eta_raw[NTRACK_MAX];
    Float_t jet_ak04its_phi[NTRACK_MAX];
    
    Float_t jet_ak04its_pt_truth[NTRACK_MAX];
    Float_t jet_ak04its_eta_truth[NTRACK_MAX];
    Float_t jet_ak04its_phi_truth[NTRACK_MAX];
  
    //The z_reco is defined as the fraction of the true jet that ended up in this reco jet 
    //There are two entries and indices, the first is the best. 
    Int_t   jet_ak04its_truth_index_z_reco[NTRACK_MAX][2];
    Float_t jet_ak04its_truth_z_reco[NTRACK_MAX][2];
    Float_t jet_ak04its_ptd_raw[NTRACK_MAX];
    Float_t jet_ak04its_width_sigma[NTRACK_MAX][2];
    UShort_t jet_ak04its_multiplicity[NTRACK_MAX];

    //Truth Jets
    UInt_t njet_truth_ak04;
    Float_t jet_truth_ak04_pt[NTRACK_MAX];
    Float_t jet_truth_ak04_eta[NTRACK_MAX];
    Float_t jet_truth_ak04_phi[NTRACK_MAX];
   
    //Int_t eg_ntrial;

    Float_t eg_cross_section;
    Int_t   eg_ntrial;
    
    //MC
    unsigned int nmc_truth;
    Float_t mc_truth_pt[NTRACK_MAX];
    Float_t mc_truth_eta[NTRACK_MAX];
    Float_t mc_truth_phi[NTRACK_MAX];
    short mc_truth_pdg_code[NTRACK_MAX];
    short mc_truth_first_parent_pdg_code[NTRACK_MAX];
    char mc_truth_charge[NTRACK_MAX];
    UChar_t mc_truth_status[NTRACK_MAX];
        
    Float_t mc_truth_first_parent_e[NTRACK_MAX];
    Float_t mc_truth_first_parent_pt[NTRACK_MAX];
    Float_t mc_truth_first_parent_eta[NTRACK_MAX];
    Float_t mc_truth_first_parent_phi[NTRACK_MAX];
    
    ULong64_t trigger_mask[2];
        
    // Set the branch addresses of the branches in the TTrees
    //    _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("is_pileup_from_spd_5_08", &is_pileup_from_spd_5_08);
    _tree_event->SetBranchAddress("is_pileup_from_spd_3_08", &is_pileup_from_spd_3_08);
    _tree_event->SetBranchAddress("ue_estimate_its_const", &ue_estimate_its_const);
    _tree_event->SetBranchAddress("ue_estimate_tpc_const", &ue_estimate_tpc_const);

    _tree_event->SetBranchAddress("trigger_mask", &trigger_mask);


    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_e", track_e);
    _tree_event->SetBranchAddress("track_pt", track_pt);
    _tree_event->SetBranchAddress("track_eta", track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    _tree_event->SetBranchAddress("track_quality", track_quality);
        
    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_e_max", cluster_e_max);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt); // here
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton); // here
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
   
    _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
    _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
    
    _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
    _tree_event->SetBranchAddress("cluster_iso_its_04_ue",cluster_iso_its_04_ue);

    _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);
    _tree_event->SetBranchAddress("cluster_nlocal_maxima", cluster_nlocal_maxima);        
    _tree_event->SetBranchAddress("cluster_distance_to_bad_channel", cluster_distance_to_bad_channel);

    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);
        
    _tree_event->SetBranchAddress("nmc_truth", &nmc_truth);
    _tree_event->SetBranchAddress("mc_truth_pdg_code", mc_truth_pdg_code);
    _tree_event->SetBranchAddress("mc_truth_pt", mc_truth_pt);
    _tree_event->SetBranchAddress("mc_truth_phi", mc_truth_phi);
    _tree_event->SetBranchAddress("mc_truth_eta", mc_truth_eta);
    _tree_event->SetBranchAddress("mc_truth_status", mc_truth_status);        
    _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code",mc_truth_first_parent_pdg_code);

    _tree_event->SetBranchAddress("eg_cross_section",&eg_cross_section);
    _tree_event->SetBranchAddress("eg_ntrial",&eg_ntrial);
    

    //jets
    _tree_event->SetBranchAddress("njet_ak04its", &njet_ak04its);
    _tree_event->SetBranchAddress("jet_ak04its_pt_raw", jet_ak04its_pt_raw);
    _tree_event->SetBranchAddress("jet_ak04its_eta_raw", jet_ak04its_eta_raw);
    _tree_event->SetBranchAddress("jet_ak04its_phi", jet_ak04its_phi);
    _tree_event->SetBranchAddress("jet_ak04its_pt_truth", jet_ak04its_pt_truth);
    _tree_event->SetBranchAddress("jet_ak04its_eta_truth", jet_ak04its_eta_truth);
    _tree_event->SetBranchAddress("jet_ak04its_phi_truth", jet_ak04its_phi_truth);

    //quark-gluon discriminator variables
    _tree_event->SetBranchAddress("jet_ak04its_ptd_raw", jet_ak04its_ptd_raw);
    _tree_event->SetBranchAddress("jet_ak04its_width_sigma", jet_ak04its_width_sigma);
    _tree_event->SetBranchAddress("jet_ak04its_multiplicity_raw", jet_ak04its_multiplicity);



    _tree_event->SetBranchAddress("jet_ak04its_truth_index_z_reco",     jet_ak04its_truth_index_z_reco);
    _tree_event->SetBranchAddress("jet_ak04its_truth_z_reco", jet_ak04its_truth_z_reco);    

    //truth jets
    _tree_event->SetBranchAddress("njet_truth_ak04", &njet_truth_ak04);
    _tree_event->SetBranchAddress("jet_truth_ak04_pt", jet_truth_ak04_pt);
    _tree_event->SetBranchAddress("jet_truth_ak04_phi", jet_truth_ak04_phi);
    _tree_event->SetBranchAddress("jet_truth_ak04_eta", jet_truth_ak04_eta);    
    // Loop over events

    _tree_event->GetEntry(1);
    if(nmc_truth>0) isRealData= false;
    else isRealData = true;
     
 
    std::cout<<" About to start looping over events to get weights" << std::endl;

    if( not(nevents>0)){
      nevents = _tree_event->GetEntries();
    }

   

    if(isRealData){
      for(Long64_t ievent = 0; ievent < nevents ; ievent++){
	    if (ievent % 100000 == 0) std::cout << " event " << ievent << std::endl;

            _tree_event->GetEntry(ievent);
            if(not( TMath::Abs(primary_vertex[2])<primary_vertex_max)) continue; //vertex z position
            if(not (primary_vertex[2]!=0.00 )) continue; //removes default of vertex z = 0
            if(is_pileup_from_spd_5_08) continue; //removes pileup


            h_zvertex.Fill(primary_vertex[2]);
            //fill UE: 
	    h_evt_rhoITS.Fill(ue_estimate_its_const);
            h_evt_rhoTPC.Fill(ue_estimate_tpc_const);



            ULong64_t one1 = 1;
            ULong64_t triggerMask_13data = (one1 << 17) | (one1 << 18) | (one1 << 19) | (one1 << 20); //EG1 or EG2 or EJ1 or EJ2
            //if(triggerMask_13data & trigger_mask[0] == 0) continue; //trigger selection

            for (ULong64_t n = 0; n < ncluster; n++) {
            
            double isolation;
            //remove UE subtraction
            if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n] + cluster_iso_its_04_ue[n];
            else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n];
            else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n] + cluster_iso_its_04_ue[n];
            else isolation = cluster_frixione_its_04_02[n] + cluster_iso_its_04_ue[n];
                
	        isolation = isolation - ue_estimate_its_const*0.4*0.4*TMath::Pi(); //Use rhoxA subtraction

                if( not(cluster_pt[n]>clus_pT_min)) continue; //select pt of photons
                if( not(cluster_pt[n]<clus_pT_max)) continue;
                if( not(cluster_ncell[n]>Cluster_ncell_min)) continue;   //removes clusters with 1 or 2 cells
                if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters
                if( not(cluster_nlocal_maxima[n]<= Cluster_locmaxima_max)) continue; //require to have at most 2 local maxima.
                if( not(cluster_distance_to_bad_channel[n]>=Cluster_distobadchannel)) continue;
                if( not(isolation < iso_max)) continue;

                Bool_t inSignalRegion;
                Bool_t inBkgRegion;
                
                if (photon_identifier == DNN) {
                    inSignalRegion = ((cluster_s_nphoton[n][1] > SIG_DNN_min) and (cluster_s_nphoton[n][1]<SIG_DNN_max));
                    inBkgRegion    = ((cluster_s_nphoton[n][1]>BKG_DNN_min) and (cluster_s_nphoton[n][1]<BKG_DNN_max));
                }
                else if (photon_identifier == LAMBDA_0) {
                    inSignalRegion = ((cluster_lambda_square[n][0]>SIG_lambda_min) and (cluster_lambda_square[n][0]<SIG_lambda_max));
                    inBkgRegion    = ((cluster_lambda_square[n][0]>BKG_lambda_min) and (cluster_lambda_square[n][0]<BKG_lambda_max));
                }
                else {
                    float eratio = cluster_e_max[n]/cluster_e[n];
                    inSignalRegion = (( eratio > SIG_Emax_over_Ecluster_min) and (eratio < SIG_Emax_over_Ecluster_max));
                    inBkgRegion    = ((eratio > BKG_Emax_over_Ecluster_min) and (eratio < BKG_Emax_over_Ecluster_max));
                }
           
                if(inSignalRegion){
		    hweight.Fill(cluster_pt[n]);
		}
                else if(inBkgRegion){
		    hBR.Fill(cluster_pt[n]);
		}
	    } //cluster
	}//events
	hweight.Divide(&hBR);
	std::cout << " Weights " << std::endl;
        for(int i=0 ; i< hweight.GetNbinsX() ; i++) std::cout <<" i" << i << " weight= " << hweight.GetBinContent(i) << std::endl;
    }//end loop over events to get weights for background region
    
    std::cout<<" About to start looping over events" << std::endl;
    h_cutflow.GetXaxis()->SetBinLabel(1, "All");
    h_cutflow.GetXaxis()->SetBinLabel(2, Form("%2.2f<pt<%2.2f GeV", clus_pT_min, clus_pT_max));
    h_cutflow.GetXaxis()->SetBinLabel(3, Form("N_{cell}>%2.2f", Cluster_ncell_min));
    h_cutflow.GetXaxis()->SetBinLabel(4, "E_{cross}/E_{cell}");
    h_cutflow.GetXaxis()->SetBinLabel(5, Form("NLM <%2.2f", Cluster_locmaxima_max));
    h_cutflow.GetXaxis()->SetBinLabel(6, Form("Distance-to-bad chanel>=%2.2f", Cluster_distobadchannel));
    h_cutflow.GetXaxis()->SetBinLabel(7, "Isolation < 2 GeV");
    h_cutflow.GetXaxis()->SetBinLabel(8, Form("Isolation < %2.2f GeV", iso_max));
    h_cutflow.GetXaxis()->SetBinLabel(9, "In Signal Region");
    h_cutflow.GetXaxis()->SetBinLabel(10, "In Background Region");

    h_evtcutflow.GetXaxis()->SetBinLabel(1, "All");
    h_evtcutflow.GetXaxis()->SetBinLabel(2, "Vertex |z| < 10 cm");
    h_evtcutflow.GetXaxis()->SetBinLabel(3, "z!=0.0 cm");

    h_evtcutflow.GetXaxis()->SetBinLabel(4, "Pileup rejection");
    h_evtcutflow.GetXaxis()->SetBinLabel(5, "Trigger Selection");
    
    h_trkcutflow.GetXaxis()->SetBinLabel(1, "All");
    h_trkcutflow.GetXaxis()->SetBinLabel(2, Form("Track pt<%2.2f", track_pT_max));
    h_trkcutflow.GetXaxis()->SetBinLabel(3, "Selection 3");
      
    h_jetcutflow.GetXaxis()->SetBinLabel(1, "All");
    h_jetcutflow.GetXaxis()->SetBinLabel(2, Form("Jet p_{T}>%2.2f", jet_pT_min));
    h_jetcutflow.GetXaxis()->SetBinLabel(3, Form("Jet eta<%2.2f", Eta_max));
    h_jetcutflow.GetXaxis()->SetBinLabel(4, "In Signal Region");
    h_jetcutflow.GetXaxis()->SetBinLabel(5, "In Background Region");
      
    for(Long64_t ievent = 0; ievent < nevents ; ievent++){
      if(ievent%2) continue;
      _tree_event->GetEntry(ievent);
      h_evtcutflow.Fill(0);
      //Eevent Selection: 
      if(not( TMath::Abs(primary_vertex[2])<primary_vertex_max)) continue; //vertex z position
      h_evtcutflow.Fill(1);
      if(not (primary_vertex[2]!=0.00 )) continue; //removes default of vertex z = 0
      h_evtcutflow.Fill(2);
      if(is_pileup_from_spd_5_08) continue; //removes pileup
      h_evtcutflow.Fill(3);     
      ULong64_t one1 = 1;
      ULong64_t triggerMask_13data = (one1 << 17) | (one1 << 18) | (one1 << 19) | (one1 << 20); //EG1 or EG2 or EJ1 or EJ2
      //if(isRealData and (triggerMask_13data & trigger_mask[0]) == 0) continue; //trigger selection
      h_evtcutflow.Fill(4);

      N_eventpassed +=1;


      double weight = 1.0;
      if(not isRealData){
	if(eg_cross_section>0 and eg_ntrial>0){
          weight = eg_cross_section/(double)eg_ntrial;
        }
          //17g6a1 weights
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat1_ptmin12.0_Nevent_500000.root"){
              if(ievent == 0)
                  std::cout << "PtHat 1 of 17g6a1 series opened" << std::endl;
              weight = 1.60e-11;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat2_ptmin12.0_Nevent_500000.root"){
              if(ievent == 0)
                  std::cout << "PtHat 2 of 17g6a1 series opened" << std::endl;
              weight = 2.72e-12;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat3_ptmin12.0_Nevent_500000.root"){
              if(ievent == 0)
                  std::cout << "PtHat 3 of 17g6a1 series opened" << std::endl;
              weight = 3.69e-13;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat4_ptmin12.0_Nevent_500000.root"){
              if(ievent == 0)
                  std::cout << "PtHat 4 of 17g6a1 series opened" << std::endl;
              weight = 6.14e-14;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat5_ptmin12.0_Nevent_500000.root"){
              if(ievent == 0)
                  std::cout << "PtHat 5 of 17g6a1 series opened" << std::endl;
              weight = 1.27e-14;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat1_ptmin12.0_Nevent_300000.root"){
              if(ievent == 0)
                  std::cout << "PtHat 1 of 17g6a1 series opened" << std::endl;
              weight = 1.60e-11;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat2_ptmin12.0_Nevent_300000.root"){
              if(ievent == 0)
                  std::cout << "PtHat 2 of 17g6a1 series opened" << std::endl;
              weight = 2.72e-12;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat3_ptmin12.0_Nevent_300000.root"){
              if(ievent == 0)
                  std::cout << "PtHat 3 of 17g6a1 series opened" << std::endl;
              weight = 3.69e-13;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat4_ptmin12.0_Nevent_300000.root"){
              if(ievent == 0)
                  std::cout << "PtHat 4 of 17g6a1 series opened" << std::endl;
              weight = 6.14e-14;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat5_ptmin12.0_Nevent_300000.root"){
              if(ievent == 0)
                  std::cout << "PtHat 5 of 17g6a1 series opened" << std::endl;
              weight = 1.27e-14;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/17g6a1_pthat1.root"){
              if(ievent == 0)
                  std::cout << "PtHat 1 of 17g6a1 series opened" << std::endl;
              weight = 1.60e-11;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/17g6a1_pthat2.root"){
              if(ievent == 0)
                  std::cout << "PtHat 2 of 17g6a1 series opened" << std::endl;
              weight = 2.72e-12;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/17g6a1_pthat3.root"){
              if(ievent == 0)
                  std::cout << "PtHat 3 of 17g6a1 series opened" << std::endl;
              weight = 3.69e-13;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/17g6a1_pthat4.root"){
              if(ievent == 0)
                  std::cout << "PtHat 4 of 17g6a1 series opened" << std::endl;
              weight = 6.14e-14;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/17g6a1_pthat5.root"){
              if(ievent == 0)
                  std::cout << "PtHat 5 of 17g6a1 series opened" << std::endl;
              weight = 1.27e-14;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat2_4L_allruns_ptmin15.0.root"){
              if(ievent == 0)
                  std::cout << "PtHat 2 of 17g6a1 series opened" << std::endl;
              weight = 2.72e-12;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat3_4L_allruns_ptmin15.0.root"){
              if(ievent == 0)
                  std::cout << "PtHat 3 of 17g6a1 series opened" << std::endl;
              weight = 3.69e-13;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat4_4L_allruns_ptmin15.0.root"){
              if(ievent == 0)
                  std::cout << "PtHat 4 of 17g6a1 series opened" << std::endl;
              weight = 6.14e-14;
          }
          if(filestring == "/project/projectdirs/alice/NTuples/MC/17g6a1/Skimmed_17g6a1_pthat5_4L_allruns_ptmin15.0.root"){
              if(ievent == 0)
                  std::cout << "PtHat 5 of 17g6a1 series opened" << std::endl;
              weight = 1.27e-14;
          }
      }
      //std::cout << " weight " << weight << std::endl;        

      h_evt_rho.Fill(ue_estimate_its_const, weight);

      //loop over tracks
      const int TrackCutBit =16;
      for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {
          h_trkcutflow.Fill(0);
	if(track_pt[itrack] < track_pT_max) continue; //1GeV Tracks
          h_trkcutflow.Fill(1);
	if((track_quality[itrack]&TrackCutBit)==0) continue; //select only tracks that pass selection 3
	h_trackphi.Fill(track_phi[itrack],weight);
          h_trkcutflow.Fill(2);
      }

      //loop over jets
      for (ULong64_t ijet = 0; ijet < njet_ak04its; ijet++) { //start loop over jets
	if(not (jet_ak04its_pt_raw[ijet]>jet_pT_min)) continue;
	if(not (TMath::Abs(jet_ak04its_eta_raw[ijet]) <Eta_max)) continue;
        h_jetphi.Fill(jet_ak04its_phi[ijet], weight);
      }


      //loop over clusters
      for (ULong64_t n = 0; n < ncluster; n++) {
        //Photon Selection
       
          double isolation;
          //remove UE subtraction
          if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n] + cluster_iso_its_04_ue[n];
          else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n] + cluster_iso_its_04_ue[n];
          else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n] + cluster_iso_its_04_ue[n];
          else isolation = cluster_frixione_its_04_02[n] + cluster_iso_its_04_ue[n];
          
          isolation = isolation - ue_estimate_its_const*0.4*0.4*TMath::Pi(); //Use rhoxA subtraction
          
          // Photon Identification
          Bool_t inSignalRegion;
          Bool_t inBkgRegion;
          float eratio = cluster_e_max[n]/cluster_e[n];
          if (photon_identifier == DNN) {
              inSignalRegion = ((cluster_s_nphoton[n][1] > SIG_DNN_min) and (cluster_s_nphoton[n][1]<SIG_DNN_max));
              inBkgRegion    = ((cluster_s_nphoton[n][1]>BKG_DNN_min) and (cluster_s_nphoton[n][1]<BKG_DNN_max));
          }
          else if (photon_identifier == LAMBDA_0) {
              inSignalRegion = ((cluster_lambda_square[n][0]>SIG_lambda_min) and (cluster_lambda_square[n][0]<SIG_lambda_max));
              inBkgRegion    = ((cluster_lambda_square[n][0]>BKG_lambda_min) and (cluster_lambda_square[n][0]<BKG_lambda_max));
          }
          else {
              inSignalRegion = ((eratio > SIG_Emax_over_Ecluster_min) and (eratio < SIG_Emax_over_Ecluster_max));
              inBkgRegion    = ((eratio > BKG_Emax_over_Ecluster_min) and (eratio < BKG_Emax_over_Ecluster_max));
          }

	h_cutflow.Fill(0);
        if( not(cluster_pt[n]>clus_pT_min)) continue; //select pt of photons
        if( not(cluster_pt[n]<clus_pT_max)) continue;
        h_cutflow.Fill(1);
        if( not(cluster_ncell[n]>Cluster_ncell_min)) continue;   //removes clusters with 1 or 2 cells
        h_cutflow.Fill(2);
	if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters
        h_cutflow.Fill(3);
        if( not(cluster_nlocal_maxima[n]< Cluster_locmaxima_max)) continue; //require to have at most 2 local maxima.
        h_cutflow.Fill(4);
	if( not(cluster_distance_to_bad_channel[n]>=Cluster_distobadchannel)) continue;
        h_cutflow.Fill(5);
	h_clusterphi.Fill(cluster_phi[n], weight);
        h_clustereta.Fill(cluster_eta[n], weight);
        if( not(isolation < 2.0)) continue;
        h_cutflow.Fill(6);
	if( not(isolation < iso_max)) continue;
        h_clusterphi_iso.Fill(cluster_phi[n], weight);
        h_clustereta_iso.Fill(cluster_eta[n], weight);
        h_cutflow.Fill(7);
    h_lambda_0.Fill(cluster_lambda_square[n][0]);
    h_DNN.Fill(cluster_s_nphoton[n][1]);
          h_EmaxOverEcluster.Fill(eratio);
    if (inSignalRegion)
        h_cutflow.Fill(8);
    if(inBkgRegion)
        h_cutflow.Fill(9);

	Bool_t isTruePhoton = false;
        Float_t truth_pt = -999.0;
        Float_t truth_eta = -999.0;
        Float_t truth_phi = -999.0;

	for(int counter = 0 ; counter<32; counter++){
	  unsigned short index = cluster_mc_truth_index[n][counter];                   
          if(isTruePhoton) break;
          if(index==65535) continue;
          if(mc_truth_pdg_code[index]!=rightpdgcode) continue;
          if(mc_truth_first_parent_pdg_code[index]!=rightparentpdgcode) continue;
          if( not (mc_truth_status[index] >0)) continue;        
          isTruePhoton = true;
          truth_pt     = mc_truth_pt[index];
          truth_phi    =  mc_truth_phi[index];
          truth_eta    =  mc_truth_eta[index];
	}//end loop over indices

        //if (not isTruePhoton){ std::cout << " photon is not true " << std::endl;} 
	// if((not isRealData) and (not isTruePhoton) and weight!=1.0 ) continue; //17g samples don't have weight FIX ME: of course this cut only works for GJ and not JJ
	
        if( not isRealData){
	  if(cluster_phi[n]<0 and cluster_phi[n]>-2.0){
	    weight = weight*0.45; //adhoc weighting for DCAL acceptance
	  }
	}
	
	//start jet loop


        if(inSignalRegion){
           N_SR +=weight;
	}
        else if(inBkgRegion){
            double bkg_weight = hweight.GetBinContent(hweight.FindBin(cluster_pt[n]));
	    if( isRealData) weight = weight*bkg_weight; //pt-dependent weight for background;					      
            N_BR +=weight;
	}

        Int_t njets_SR = 0; 
        Int_t njets_BR = 0;
	for (ULong64_t ijet = 0; ijet < njet_ak04its; ijet++) { //start loop over jets
          h_jetcutflow.Fill(0);
          if(not (jet_ak04its_pt_raw[ijet]>jet_pT_min)) continue;
          h_jetcutflow.Fill(1);
          if(not (TMath::Abs(jet_ak04its_eta_raw[ijet]) <Eta_max)) continue;
          h_jetcutflow.Fill(2);
          if(inSignalRegion)
              h_jetcutflow.Fill(3);
          if(inBkgRegion)
              h_jetcutflow.Fill(4);
          Float_t dphi = TMath::Abs(TVector2::Phi_mpi_pi(jet_ak04its_phi[ijet] - cluster_phi[n]));
	  Float_t deta = jet_ak04its_eta_raw[ijet] - cluster_eta[n];
	  Float_t dphi_truth = 0;
          Float_t deta_truth = 0;

          if(eg_ntrial>0){
              dphi_truth = TMath::Abs(TVector2::Phi_mpi_pi(jet_ak04its_phi_truth[ijet] - truth_phi));
	      deta_truth = jet_ak04its_eta_truth[ijet] -  truth_eta;
	  }
          if(inSignalRegion){
	    hSR_dPhi.Fill(dphi,weight);
            hSR_dPhi_truth.Fill(dphi_truth,weight);
              num_sig_dPhi++;
	 
	  }
	  else if(inBkgRegion){
	    hBR_dPhi.Fill(dphi,weight);
            hBR_dPhi_truth.Fill(dphi_truth,weight);
              num_bkg_dPhi++;
	  }
 	  //if(not (dphi>0.4)) continue; 
          if( not (dphi>TMath::Pi()/2.0)) continue;

          //counts jets associated with clusters   
          if(inSignalRegion){
	    njets_SR =+1 ;
	  }
          else if(inBkgRegion){
            njets_BR =+1; 
	  }
          

	  Float_t xj = jet_ak04its_pt_raw[ijet]/cluster_pt[n];
	  //std::cout <<"truthptjet: " << jet_ak04its_pt_truth[ijet] << "reco pt jet" << jet_ak04its_pt_raw[ijet] << std::endl;           
          Float_t xj_truth = jet_ak04its_pt_truth[ijet]/truth_pt; 
        TOT_jetpt.Fill(jet_ak04its_pt_raw[ijet], weight);
	  if( inSignalRegion){
            hSR_Xj.Fill(xj, weight);
	    hSR_dEta.Fill(deta, weight);
            hSR_AvgEta.Fill(0.5*(jet_ak04its_eta_raw[ijet] + cluster_eta[n]), weight);

            hSR_pTD.Fill(jet_ak04its_ptd_raw[ijet],weight);
	    //	    std::cout << jet_ak04its_multiplicity[ijet] << " " << jet_ak04its_width_sigma[ijet] << " " << jet_ak04its_ptd_raw[ijet] << std::endl;
            hSR_Multiplicity.Fill(jet_ak04its_multiplicity[ijet],weight);
            hSR_jetwidth.Fill(TMath::Log(jet_ak04its_width_sigma[ijet][0]), weight);

            //associated jet rate
	    hSR_jetpt.Fill(jet_ak04its_pt_raw[ijet], weight);
	    hSR_jeteta.Fill(jet_ak04its_eta_raw[ijet], weight);
	    hSR_jetphi.Fill(jet_ak04its_phi[ijet],weight);

	    hSR_jetpt_truth.Fill(jet_ak04its_pt_truth[ijet], weight);
            hSR_jeteta_truth.Fill(jet_ak04its_eta_truth[ijet], weight);
            hSR_jetphi_truth.Fill(jet_ak04its_phi_truth[ijet],weight);

            hSR_Xj_truth.Fill(xj_truth, weight);
            hSR_dEta_truth.Fill(deta_truth,weight);
            hSR_AvgEta_truth.Fill(0.5*(jet_ak04its_eta_truth[ijet] +  truth_eta), weight);
          
          //std::cout << "Signal region Xj truth: " << xj_truth << std::endl;
          
        hSR_XobsPb.Fill(((cluster_pt[n]*TMath::Exp(-cluster_eta[n]))+(jet_ak04its_pt_raw[ijet]*TMath::Exp(-jet_ak04its_eta_raw[ijet])))/(2*EPb), weight);
          
            hSR_XobsPb_truth.Fill(((truth_pt*TMath::Exp(-truth_eta))+(jet_ak04its_pt_truth[ijet]*TMath::Exp(-jet_ak04its_eta_truth[ijet])))/(2*EPb), weight);
          
          //std::cout << "Signal region XobsPb truth: " << ((truth_pt*TMath::Exp(-truth_eta))+(jet_ak04its_pt_truth[ijet]*TMath::Exp(-jet_ak04its_eta_truth[ijet])))/(2*EPb) << std::endl;
          
          num_sig_XobsPb++;
          
            h_Xj_Matrix.Fill(xj_truth, xj, weight);
	    //std::cout<<" xj " << xj << " " << " xj_truth "<< xj_truth << std::endl;

	  }
          else if(inBkgRegion){
	    hBR_Xj.Fill(xj,weight);
	    hBR_dEta.Fill(deta,weight);
	    hBR_AvgEta.Fill(0.5*(jet_ak04its_eta_raw[ijet] + cluster_eta[n]), weight);

	    hBR_pTD.Fill(jet_ak04its_ptd_raw[ijet],weight);
	    hBR_Multiplicity.Fill(jet_ak04its_multiplicity[ijet],weight);
            hBR_jetwidth.Fill(TMath::Log(jet_ak04its_width_sigma[ijet][0]), weight);

	    //Associated jet rate
	    hBR_jetpt.Fill(jet_ak04its_pt_raw[ijet], weight);
	    hBR_jeteta.Fill(jet_ak04its_eta_raw[ijet], weight);
            hBR_jetphi.Fill(jet_ak04its_phi[ijet],weight);

	    hBR_jetpt_truth.Fill(jet_ak04its_pt_truth[ijet], weight);
            hBR_jeteta_truth.Fill(jet_ak04its_eta_truth[ijet], weight);
            hBR_jetphi_truth.Fill(jet_ak04its_phi_truth[ijet],weight);

            hBR_Xj_truth.Fill(xj_truth,weight);
	    hBR_dEta_truth.Fill(deta_truth,weight);
	    hSR_AvgEta_truth.Fill(0.5*(jet_ak04its_eta_truth[ijet] +  truth_eta), weight);
              
              //std::cout << "Background region Xj truth: " << xj_truth << std::endl;
              
        hBR_XobsPb.Fill(((cluster_pt[n]*TMath::Exp(-cluster_eta[n]))+(jet_ak04its_pt_raw[ijet]*TMath::Exp(-jet_ak04its_eta_raw[ijet])))/(2*EPb), weight);
              
        hBR_XobsPb_truth.Fill(((truth_pt*TMath::Exp(-truth_eta))+(jet_ak04its_pt_truth[ijet]*TMath::Exp(-jet_ak04its_eta_truth[ijet])))/(2*EPb), weight);
        
        //std::cout << "Background region XobsPb truth: " << ((truth_pt*TMath::Exp(-truth_eta))+(jet_ak04its_pt_truth[ijet]*TMath::Exp(-jet_ak04its_eta_truth[ijet])))/(2*EPb) << std::endl;
              
              num_bkg_XobsPb++;
	  }



	}//end loop over jets
      TOT_clusterpt.Fill(cluster_pt[n], weight);
        if(inSignalRegion){
	  hSR_clusterpt.Fill(cluster_pt[n], weight);
	  hSR_clustereta.Fill(cluster_eta[n], weight);
	  hSR_clusterphi.Fill(cluster_phi[n], weight);
          hSR_njet.Fill(njets_SR, weight); 
	}
        else if(inBkgRegion){
	  hBR_clusterpt.Fill(cluster_pt[n], weight);
          hBR_clustereta.Fill(cluster_eta[n], weight);
          hBR_clusterphi.Fill(cluster_phi[n], weight);
	  hBR_njet.Fill(njets_BR, weight);
	}
	//fill in this histogram only photons that can be traced to a generated non-decay photon.	
        h_reco_truthpt.Fill(truth_pt,weight);
        h_reco.Fill(cluster_pt[n],weight); 

      }//end loop on clusters


      //** Study of jet reconstruction efficiency:

      //loop over truth jets
      for (ULong64_t ijet = 0; ijet < njet_truth_ak04; ijet++) {
	if(not(TMath::Abs(jet_truth_ak04_eta[ijet])<Eta_max)) continue;
         h_jetpt_truth.Fill(jet_truth_ak04_pt[ijet], weight);
      }
    
      std::set<int> temp; //to store truth indices associated with reco jets
      for (ULong64_t ijet = 0; ijet < njet_ak04its; ijet++) { 
	if(not (jet_ak04its_pt_raw[ijet]>jet_pT_min)) continue;
	if(not (TMath::Abs(jet_ak04its_eta_raw[ijet])  <Eta_max ) ) continue;
	h_jetpt_reco.Fill(jet_ak04its_pt_raw[ijet], weight);
	temp.insert(jet_ak04its_truth_index_z_reco[ijet][0]);
      } //end loop over reco jets
      for(auto& index: temp){
	if(index>0){
	  if(not(TMath::Abs(jet_truth_ak04_eta[index])<Eta_max)) continue;
          h_jetpt_truthreco.Fill(jet_truth_ak04_pt[index],weight);
	}
      }//end loop over indices of reco jets

      for (ULong64_t nmc = 0; nmc < nmc_truth; nmc++) {
        if(not(mc_truth_pt[nmc]>clus_pT_min)) continue;
	if(not(mc_truth_pt[nmc]<clus_pT_max)) continue;
	if(mc_truth_pdg_code[nmc]==rightpdgcode && int(mc_truth_status[nmc])>0 &&  mc_truth_first_parent_pdg_code[nmc]==rightparentpdgcode){
	  //std::cout << mc_truth_pt[nmc] << "phi " << mc_truth_phi[nmc] << " eta " << mc_truth_eta[nmc] << 
	  //  " code: " << mc_truth_pdg_code[nmc] << " status " << int(mc_truth_status[nmc]) << " parentpdg " << mc_truth_first_parent_pdg_code[nmc] << std::endl;    
	  h_truth.Fill(mc_truth_pt[nmc],weight);
          N_truth +=1;
	  //std::cout << " number of truth jets " << njet_truth_ak04 << std::endl;
	  for (ULong64_t ijet = 0; ijet < njet_truth_ak04; ijet++) {
            if(not(jet_truth_ak04_pt[ijet]>jet_pT_min)) continue;
	    if(not(TMath::Abs(jet_truth_ak04_eta[ijet])<Eta_max)) continue;
	    Float_t dphi_truth = TMath::Abs(TVector2::Phi_mpi_pi(jet_truth_ak04_phi[ijet] - mc_truth_phi[nmc]));
	    //std::cout<< dphi_truth << std::endl;
            
	    h_dPhi_truth.Fill(dphi_truth,weight);
            //if( not(dphi_truth>0.4)) continue;
            if( not(dphi_truth>TMath::Pi()/2.0)) continue;
	    Float_t xj_truth = jet_truth_ak04_pt[ijet]/mc_truth_pt[nmc];
            h_Xj_truth.Fill(xj_truth,weight);
	  }//end loop over truth jets
	}
      }//end loop over mc particles

     // Create the file label, to be used within the filenames, to represent the source file
     

      h_truth.SetLineColor(2);
      THStack* hs = new THStack("hs","stack histo for plotting");
      hs->Add(&h_reco);
      hs->Add(&h_truth);

      if (ievent % 10000 == 0) {
	//SR_Xj.Draw("e1x0nostack");
        hSR_dPhi.Draw("e1x0nostack");
	std::cout << ievent << " " << _tree_event->GetEntries() << std::endl;
        canvas->Update();
      } 

     }//end over events
  } // end over files
    std::cout << " Numbers of events passing selection " << N_eventpassed << std::endl;
    std::cout << " Number of clusters in signal region " << N_SR << std::endl;
    std::cout << " Number of clusters in background region " << N_BR << std::endl;
    std::cout << " Number of truth photons " << N_truth << std::endl;

    std::string opened_files = "";
    for (int iarg = 1; iarg < argc; iarg++) {
        std::string filepath = argv[iarg];
        opened_files += "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
    }
    
    std::string photonselectionvar;
    if (photon_identifier == DNN) {
        photonselectionvar = "DNN";
    }
    else if (photon_identifier == LAMBDA_0) {
        photonselectionvar = "Lambda0";
    }
    else {
        photonselectionvar = "EmaxOverEcluster";
    }
    //TFile* fout = new TFile(Form("GammaJet_config_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_%s_PHOTONSELECT_%s.root", clus_pT_min, clus_pT_max, jet_pT_min, opened_files.c_str(), photonselectionvar.c_str()),"RECREATE");
    TFile* fout = new TFile(Form("GammaJet_config_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_MC17g6a1_PHOTONSELECT_%s.root", clus_pT_min, clus_pT_max, jet_pT_min, photonselectionvar.c_str()),"RECREATE");
    //TFile* fout = new TFile(Form("GammaJet_config_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_MCdijet_PHOTONSELECT_%s.root", clus_pT_min, clus_pT_max, jet_pT_min, photonselectionvar.c_str()),"RECREATE");
    //TFile* fout = new TFile(Form("GammaJet_config_clusptmin%2.1f_clusptmax%2.1f_JETPTMIN_%2.1f_DATANAME_MCgammajet_PHOTONSELECT_%s.root", clus_pT_min, clus_pT_max, jet_pT_min, photonselectionvar.c_str()),"RECREATE");
    fout->Print();

    //Save the sum of weights 
    h_weights.SetBinContent(1, N_SR);
    h_weights.SetBinContent(2, N_BR);
    
    h_zvertex.Write("zvertex");
    h_evt_rhoITS.Write("h_evt_rhoITS");
    h_evt_rhoTPC.Write("h_evt_rhoTPC");

    h_weights.Write("h_weights");


    ///weighting
    hSR_njet.Scale(1.0/(N_SR*hSR_njet_binwidth));
    hBR_njet.Scale(1.0/(N_BR*hBR_njet_binwidth));

    hSR_Xj.Scale(1.0/(N_SR*hSR_Xj_binwidth));
    hBR_Xj.Scale(1.0/(N_BR*hBR_Xj_binwidth));
    hSR_dPhi.Scale(1.0/(N_SR*hSR_dPhi_binwidth));
    hBR_dPhi.Scale(1.0/(N_BR*hBR_dPhi_binwidth));
    hSR_dEta.Scale(1.0/(N_SR*hSR_dEta_binwidth));
    hBR_dEta.Scale(1.0/(N_BR*hBR_dEta_binwidth));

    hSR_AvgEta.Scale(1.0/(N_SR*hSR_AvgEta_binwidth));
    hBR_AvgEta.Scale(1.0/(N_BR*hBR_AvgEta_binwidth));

    hSR_jetpt.Scale(1.0/(N_SR*hSR_jetpt_binwidth));
    hBR_jetpt.Scale(1.0/(N_BR*hBR_jetpt_binwidth));
    hSR_jeteta.Scale(1.0/(N_SR*hSR_jeteta_binwidth));
    hBR_jeteta.Scale(1.0/(N_BR*hBR_jeteta_binwidth));
    hSR_jetphi.Scale(1.0/(N_SR*hSR_jetphi_binwidth));
    hBR_jetphi.Scale(1.0/(N_BR*hBR_jetphi_binwidth));


    hSR_jetpt_truth.Scale(1.0/(N_SR*hSR_jetpt_truth_binwidth));
    hBR_jetpt_truth.Scale(1.0/(N_BR*hBR_jetpt_truth_binwidth));
    hSR_jeteta_truth.Scale(1.0/(N_SR*hSR_jeteta_truth_binwidth));
    hBR_jeteta_truth.Scale(1.0/(N_BR*hBR_jeteta_truth_binwidth));
    hSR_jetphi_truth.Scale(1.0/(N_SR*hSR_jetphi_truth_binwidth));
    hBR_jetphi_truth.Scale(1.0/(N_BR*hBR_jetphi_truth_binwidth));


    hSR_Xj_truth.Scale(1.0/(N_SR*hSR_Xj_truth_binwidth));
    hBR_Xj_truth.Scale(1.0/(N_BR*hBR_Xj_truth_binwidth));
    hSR_dPhi_truth.Scale(1.0/(N_SR*hSR_dPhi_truth_binwidth));
    hBR_dPhi_truth.Scale(1.0/(N_BR*hBR_dPhi_truth_binwidth));
    hSR_dEta_truth.Scale(1.0/(N_SR*hSR_dEta_truth_binwidth));
    hBR_dEta_truth.Scale(1.0/(N_BR*hBR_dEta_truth_binwidth));
    hSR_AvgEta_truth.Scale(1.0/(N_SR*hSR_AvgEta_truth_binwidth));
    hBR_AvgEta_truth.Scale(1.0/(N_BR*hBR_AvgEta_truth_binwidth));


    hSR_pTD.Scale(1.0/(N_SR*hSR_pTD_binwidth));
    hBR_pTD.Scale(1.0/(N_BR*hBR_pTD_binwidth));
    hSR_Multiplicity.Scale(1.0/(N_SR*hSR_Multiplicity_binwidth));
    hBR_Multiplicity.Scale(1.0/(N_BR*hBR_Multiplicity_binwidth));
    hSR_jetwidth.Scale(1.0/(N_SR*hSR_jetwidth_binwidth));
    hBR_jetwidth.Scale(1.0/(N_BR*hBR_jetwidth_binwidth));

    
    hSR_clustereta.Scale(1.0/(N_SR*hSR_clustereta_binwidth));
    hBR_clustereta.Scale(1.0/(N_BR*hBR_clustereta_binwidth));
    hSR_clusterphi.Scale(1.0/(N_SR*hSR_clusterphi_binwidth));
    hBR_clusterphi.Scale(1.0/(N_BR*hBR_clusterphi_binwidth));

    hSR_XobsPb.Scale(1.0/(N_SR*hSR_XobsPb_binwidth));
    hBR_XobsPb.Scale(1.0/(N_BR*hBR_XobsPb_binwidth));
    hSR_XobsPb_truth.Scale(1.0/(N_SR*hSR_XobsPb_truth_binwidth));
    hBR_XobsPb_truth.Scale(1.0/(N_BR*hBR_XobsPb_truth_binwidth));














    h_dPhi_truth.Scale(1.0/(N_truth*h_dPhi_truth_binwidth));
    h_Xj_truth.Scale(1.0/(N_truth*h_Xj_truth_binwidth));


    h_evtcutflow.Write("EventCutFlow");
    h_cutflow.Write("ClusterCutFlow");
    h_trkcutflow.Write("TrackCutFlow");
    h_jetcutflow.Write("JetCutFlow");
    h_evt_rho.Write("h_evt_rho");

    //number of jets
    hSR_njet.Write("sig_njet");
    hBR_njet.Write("bkg_njet");

    //Xj
    hSR_Xj.Write("sig_Xj");
    hBR_Xj.Write("bkg_Xj");
    if (not isRealData) {
        hSR_Xj_truth.Write("sig_Xj_truth");
        hBR_Xj_truth.Write("bkg_Xj_truth");
    }
    //dPhi
    std::cout << "Number of dPhi entries: singal: " << num_sig_dPhi << " background: " << num_bkg_dPhi << std::endl;
    hSR_dPhi.Write("sig_dPhi");
    hBR_dPhi.Write("bkg_dPhi");
    if (not isRealData) {
        hSR_dPhi_truth.Write("sig_dPhi_truth");
        hBR_dPhi_truth.Write("bkg_dPhi_truth");
    }
    //dEta
    hBR_dEta.Write("bkg_dEta");
    hSR_dEta.Write("sig_dEta");
    if (not isRealData) {
        hBR_dEta_truth.Write("bkg_dEta_truth");
        hSR_dEta_truth.Write("sig_dEta_truth");
    }
    //Average Eta  
    hBR_AvgEta.Write("bkg_AvgEta");
    hSR_AvgEta.Write("sig_AvgEta");
    if (not isRealData) {
        hBR_AvgEta_truth.Write("bkg_AvgEta_truth");
        hSR_AvgEta_truth.Write("sig_AvgEta_truth");
    }

    //Flavor variables
    hSR_pTD.Write("sig_pTD");
    hBR_pTD.Write("bkg_pTD");
    hSR_Multiplicity.Write("sig_Multiplicity");
    hBR_Multiplicity.Write("bkg_Multiplicity");
    hSR_jetwidth.Write("sig_jetwidth");
    hBR_jetwidth.Write("bkg_jetwidth");
    //Associated jet histograms
    TOT_jetpt.Write("TOT_jetpt");
    hSR_jetpt.Write("sig_jetpt");
    hBR_jetpt.Write("bkg_jetpt");
    if (not isRealData) {
        hSR_jetpt_truth.Write("sig_jetpt_truth");
        hBR_jetpt_truth.Write("bkg_jetpt_truth");
    }

    hSR_jeteta.Write("sig_jeteta");
    hBR_jeteta.Write("bkg_jeteta");
    if (not isRealData) {
        hSR_jeteta_truth.Write("sig_jeteta_truth");
        hBR_jeteta_truth.Write("bkg_jeteta_truth");
    }

    hSR_jetphi.Write("sig_jetphi");
    hBR_jetphi.Write("bkg_jetphi");
    if (not isRealData) {
        hSR_jetphi_truth.Write("sig_jetphi_truth");
        hBR_jetphi_truth.Write("bkg_jetphi_truth");
    }
    
    
    std::cout << "Number of XobsPb entries: singal: " << num_sig_XobsPb << " background: " << num_bkg_XobsPb << std::endl;
    hSR_XobsPb.Write("sig_XobsPb");
    hBR_XobsPb.Write("bkg_XobsPb");
    if (not isRealData) {
        hSR_XobsPb_truth.Write("sig_XobsPb_truth");
        hBR_XobsPb_truth.Write("bkg_XobsPb_truth");
    }


    //Cluster pt
    TOT_clusterpt.Write("TOT_clusterpt");
    hSR_clusterpt.Write("sig_clusterpt");
    hBR_clusterpt.Write("bkg_clusterpt");
    hSR_clustereta.Write("sig_clustereta");
    hBR_clustereta.Write("bkg_clustereta");
    hSR_clusterphi.Write("sig_clusterphi");
    hBR_clusterphi.Write("bkg_clusterphi");

    h_clusterphi.Write("h_clusterphi");
    h_clusterphi_iso.Write("h_clusterphi_iso");
    h_clusterphi_iso.Divide(&h_clusterphi);
    h_clusterphi_iso.Write("h_clusterphi_isoratio");

    h_clustereta.Write("h_clustereta");
    h_clustereta_iso.Write("h_clustereta_iso");
    h_clustereta_iso.Divide(&h_clustereta);
    h_clustereta_iso.Write("h_clustereta_isoratio");

    //MC truth
    if (not isRealData) {
        h_dPhi_truth.Write("h_dPhi_truth");
    }
    h_Xj_Matrix.Write("xj_matrix");    
    h_Xj_truth.Write("h_Xj_truth");
    if (not isRealData) {
        h_jetpt_truth.Write("h_jetpt_truth");
        h_jetpt_truthreco.Write("h_jetpt_truthreco");
    }
    h_jetpt_reco.Write("h_jetpt");


    h_trackphi.Write("h_trackphi");
    h_jetphi.Write("h_jetphi");
   
    TH1D* jet_eff = (TH1D*)h_jetpt_truthreco.Clone();
    jet_eff->Divide(&h_jetpt_truth);
    jet_eff->Write("jetefficiency");
 
    std::cout << " ending " << std::endl;
    fout->Close();
  //end of arguments
  return EXIT_SUCCESS;
 
}
