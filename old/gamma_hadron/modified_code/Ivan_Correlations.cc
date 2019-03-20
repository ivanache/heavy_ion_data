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

#define NTRACK_MAX (1U << 14)

#include <vector>
#include <math.h>

const int MAX_INPUT_LENGTH = 200;

enum isolationDet {CLUSTER_ISO_TPC_04, CLUSTER_ISO_ITS_04, CLUSTER_FRIXIONE_TPC_04_02, CLUSTER_FRIXIONE_ITS_04_02};

int main(int argc, char *argv[])
{
  if (argc < 2) {
    exit(EXIT_FAILURE);
  }
  int dummyc = 1;
  char **dummyv = new char *[1];

  dummyv[0] = strdup("main");
  
  // Read configuration file for the DNN_min and DNN_max variables
  FILE* config = fopen("Corr_config.yaml", "r");
  
  // Default values of various variables used in the file (actual values are to be determined by the configuration file)
  // Cut variables
  double DNN_min = 0.75;
  double DNN_max = 0.85;
  double pT_min = 12;
  double pT_max = 14;
  double Eta_max = 0.6;
  double Cluster_min = 2;
  double EcrossoverE_min = 0.03;
  int selection_number = 3;
    
  // The bounds for the events to fal into the isolation and nonisolation areas
  double iso_max = 2.0;
  double noniso_min = 5.0;
  double noniso_max = 15.0;
  
  // Delta eta
  double deta_max = 0.6;
  
  // Zt bins
  int nztbins = 7;
  float* ztbins;
  ztbins = new float[nztbins+1];
  ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
   
  // Number of bins in correlation functions
  int n_correlationbins = 18;
    
  // Which branch should be used to determine whether a cluster should fall into iso, noniso, or neither
  isolationDet determiner = CLUSTER_ISO_TPC_04;
  
  // Loop over the config file's lines
  // '#' means comment
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
      if (strcmp(key, "DNN_min") == 0) {
          // Assign DNN_min to the double-converted version of value
          DNN_min = atof(value);
          std::cout << "DNN_min is " << DNN_min << std::endl;
      }
      else if (strcmp(key, "DNN_max") == 0) {
          // Assign DNN_max to the double-converted version of value
          DNN_max = atof(value);
          std::cout << "DNN_max is " << DNN_max << std::endl;
      }
      else if (strcmp(key, "pT_min") == 0) {
          pT_min = atof(value);
          std::cout << "pT_min is " << pT_min << std::endl;
      }
      else if (strcmp(key, "pT_max") == 0) {
          pT_max = atof(value);
          std::cout << "pT_max is " << pT_max << std::endl;
      }
      else if (strcmp(key, "Eta_max") == 0) {
          Eta_max = atof(value);
          std::cout << "Eta_max is " << Eta_max << std::endl;
      }
      else if (strcmp(key, "Cluster_min") == 0) {
          Cluster_min = atof(value);
          std::cout << "Cluster_min is " << Cluster_min << std::endl;
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
      else if (strcmp(key, "Correlation_func_bins") == 0) {
          n_correlationbins = atoi(value);
          std::cout << "Bins in a correlation function: " << n_correlationbins << std::endl;
      }
      else if (strcmp(key, "track_selection_num") == 0) {
          selection_number = atoi(value);
          std::cout << "Number of Selection that tracks must pass: " << selection_number << std::endl;
      }
      else if (strcmp(key, "Zt_bins") == 0) {
          nztbins = -1;
          for (const char *v = value; *v != ']';) {
              while (*v != ']' && !isdigit(*v)) {
                  v++;
              }
              
              nztbins++;
              
              while (*v != ']' && (isdigit(*v) || *v == '.')) {
                  v++;
              }
          }
          ztbins = new float[nztbins + 1];
          int i = 0;
          for (const char *v = value; *v != ']' ;) {
              while (*v != ']' && !isdigit(*v)) {
                  v++;
              }
              
              ztbins[i] = atof(v);
              i++;
              
              while (*v != ']' && (isdigit(*v) || *v == '.')) {
                  v++;
              }
          }
          std::cout << "Number of Zt bins: " << nztbins << std::endl << "Zt bins: {";
          for (int i = 0; i <= nztbins; i++)
              std::cout << ztbins[i] << ", ";
          std::cout << "}\n";
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
      // Warning message if the key contents are unrecognized
      else {
          std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
      }
  }
  fclose(config);

  // Loop over files to open, I guess. I may be wrong though.
  for (int iarg = 1; iarg < argc; iarg++) {
    std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
    TFile *file = TFile::Open((TString)argv[iarg]);

    if (file == NULL) {
      std::cout << " fail" << std::endl;
      exit(EXIT_FAILURE);
    }
    file->Print();
    
    // Get all the TTree variables from the file to open, I guess
    TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

    if (_tree_event == NULL) {
      std::cout << " fail " << std::endl;
      exit(EXIT_FAILURE);
    }  
    //_tree_event->Print();

    TApplication application("", &dummyc, dummyv);
    
   // Create the TCanvas and the histograms
    TCanvas canvas("canvas", "");
    TH1D histogram0("histogram0", "", 16, 8.0, 16.0);
    //TH2D histogram1("histogram1", "", 30, -1.5, 1.5, 18, -0.5, 1.5);
    //TH1D histogram2("histogram2", "", 18, -0.5,1.5);
    TH1D histogram3("histogram3", "", 18, -0.5,1.5);
    TH1D h_ntrig("h_ntrig", "", 2, -0.5,1.0);
      
    // Create the histogram for the 2D plots
    // TH2D histogram2D0("histogram2D0", "", );
      
    // Zt bins
    //const int nztbins = 7;
    //const float ztbins[nztbins+1] = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
    
    // Function declarations of h_dPhi_iso and h_dPhi_noniso
    TH1F* h_dPhi_iso[nztbins];
    TH1F* h_dPhi_noniso[nztbins];

    // Function initializations of h_dPhi_iso and h_dPhi_noniso for each and every zt bin
    for (int izt = 0; izt<nztbins; izt++){
      h_dPhi_iso[izt] = new TH1F(Form("dPhi_iso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]) ,"", n_correlationbins,-0.5,1.5);
      h_dPhi_noniso[izt] = new TH1F(Form("dPhi_noniso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt], 10*ztbins[izt+1]), "", n_correlationbins,-0.5,1.5);
      h_dPhi_iso[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
      h_dPhi_noniso[izt]->SetTitle("; #Delta#phi/#pi [rad]; entries");
      h_dPhi_iso[izt]->Sumw2();
      h_dPhi_noniso[izt]->Sumw2();
    }

    //histogram2.Sumw2();
    //histogram3.Sumw2();

// TTree variables of each track
    //variables
    Double_t primary_vertex[3];
    UInt_t ntrack;
    Float_t track_e[NTRACK_MAX];
    Float_t track_pt[NTRACK_MAX];
    Float_t track_eta[NTRACK_MAX];
    Float_t track_phi[NTRACK_MAX];
    UChar_t track_quality[NTRACK_MAX];

    UInt_t ncluster;
    Float_t cluster_e[NTRACK_MAX];
    Float_t cluster_e_cross[NTRACK_MAX];
    Float_t cluster_pt[NTRACK_MAX];
    Float_t cluster_eta[NTRACK_MAX];
    Float_t cluster_phi[NTRACK_MAX]; 
    Float_t cluster_iso_tpc_04[NTRACK_MAX];
    Float_t cluster_iso_its_04[NTRACK_MAX];
    Float_t cluster_frixione_tpc_04_02[NTRACK_MAX];
    Float_t cluster_frixione_its_04_02[NTRACK_MAX];
    Float_t cluster_s_nphoton[NTRACK_MAX][4];
    unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
    Int_t cluster_ncell[NTRACK_MAX];
    UShort_t  cluster_cell_id_max[NTRACK_MAX];
    Float_t cluster_lambda_square[NTRACK_MAX][2];   
    Float_t cell_e[17664];
      
    Long64_t Mix_Events[10];

    //MC
    unsigned int nmc_truth;
    Float_t mc_truth_pt[NTRACK_MAX];
    Float_t mc_truth_eta[NTRACK_MAX];  
    Float_t mc_truth_phi[NTRACK_MAX];
    short mc_truth_pdg_code[NTRACK_MAX];
    short mc_truth_first_parent_pdg_code[NTRACK_MAX];
    char mc_truth_charge[NTRACK_MAX];

    Float_t mc_truth_first_parent_e[NTRACK_MAX];
    Float_t mc_truth_first_parent_pt[NTRACK_MAX];
    Float_t mc_truth_first_parent_eta[NTRACK_MAX];
    Float_t mc_truth_first_parent_phi[NTRACK_MAX];
    UChar_t mc_truth_status[NTRACK_MAX];
 
    
    // Set the branch addresses of the branches in the TTrees
    _tree_event->SetBranchStatus("*mc*", 0);
  
    _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
    _tree_event->SetBranchAddress("ntrack", &ntrack);
    _tree_event->SetBranchAddress("track_e", track_e);
    _tree_event->SetBranchAddress("track_pt", track_pt);
    _tree_event->SetBranchAddress("track_eta", track_eta);
    _tree_event->SetBranchAddress("track_phi", track_phi);
    _tree_event->SetBranchAddress("track_quality", track_quality);

    _tree_event->SetBranchAddress("ncluster", &ncluster);
    _tree_event->SetBranchAddress("cluster_e", cluster_e);
    _tree_event->SetBranchAddress("cluster_e_cross", cluster_e_cross);
    _tree_event->SetBranchAddress("cluster_pt", cluster_pt); // here
    _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
    _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
    _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton); // here
    _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
    _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
    _tree_event->SetBranchAddress("cluster_iso_tpc_04",cluster_iso_tpc_04);
    _tree_event->SetBranchAddress("cluster_iso_its_04",cluster_iso_its_04);
    _tree_event->SetBranchAddress("cluster_frixione_tpc_04_02",cluster_frixione_tpc_04_02);
    _tree_event->SetBranchAddress("cluster_frixione_its_04_02",cluster_frixione_its_04_02);

    _tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
    _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
    _tree_event->SetBranchAddress("cell_e", cell_e);
      
    _tree_event->SetBranchAddress("Mix_Events", Mix_Events);
      // Other two variables: mc_truth_pdg_code and mc_truth_pt
     
    
    std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;

   // Loop over events
    for(Long64_t ievent = 0; ievent < _tree_event->GetEntries() ; ievent++){     
      //for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
      _tree_event->GetEntry(ievent);
        
      //loop over clusters
      for (ULong64_t n = 0; n < ncluster; n++) {
    // Apply cuts
	if( not(cluster_pt[n]>pT_min and cluster_pt[n]<pT_max)) continue; //select pt of photons
	if( not(cluster_s_nphoton[n][1]>DNN_min and cluster_s_nphoton[n][1]<DNN_max)) continue; //select deep-photons
	if( not(TMath::Abs(cluster_eta[n])<Eta_max)) continue; //cut edges of detector
	if( not(cluster_ncell[n]>Cluster_min)) continue;   //removes clusters with 1 or 2 cells
	if( not(cluster_e_cross[n]/cluster_e[n]>EcrossoverE_min)) continue; //removes "spiky" clusters
    
    // This part uses the determiner varable, acquired from the configuration file's Cluster_isolation_determinant value, to determine isolation
    // where isolation is the variable to be used to fill the isolated deep-photon pt spectra and ntriggers
    float isolation;
    if (determiner == CLUSTER_ISO_TPC_04) isolation = cluster_iso_tpc_04[n];
    else if (determiner == CLUSTER_ISO_ITS_04) isolation = cluster_iso_its_04[n];
    else if (determiner == CLUSTER_FRIXIONE_TPC_04_02) isolation = cluster_frixione_tpc_04_02[n];
    else isolation = cluster_frixione_its_04_02[n];
          
    // Use the variable "isolation" to fill the isolated deep-photon pt spectra and the ntriggers
	if(isolation<iso_max){
	  histogram0.Fill(cluster_pt[n]); //isolated deep-photon pt spectra
	  h_ntrig.Fill(0);
	}
	if(isolation>noniso_min && isolation<noniso_max){
	  h_ntrig.Fill(0.5);
	}
    // Loop over tracks
	for (ULong64_t itrack = 0; itrack < ntrack; itrack++) {            
	  if((track_quality[itrack]&selection_number)==0) continue; //select only tracks that pass selection 3
      // Find the track's zt, delta eta, and delta phi with the cluster
	  Double_t zt = track_pt[itrack]/cluster_pt[n];
	  Float_t deta =  cluster_eta[n]-track_eta[itrack];
	  Float_t dphi =  TVector2::Phi_mpi_pi(cluster_phi[n]-track_phi[itrack])/TMath::Pi();
	  if(!(TMath::Abs(deta)<deta_max)) continue; // delta eta cut
	  if(dphi<-0.5) dphi +=2;
        
      // Loop over zt bins
	  for(int izt = 0; izt<nztbins ; izt++){
	    if(zt>ztbins[izt] and  zt<ztbins[izt+1])
	      {
              // Where the  h_dPhi_iso and h_dPhi_noniso bins are filled
		if(isolation< iso_max)    h_dPhi_iso[izt]->Fill(dphi);
		if(isolation> noniso_min && isolation<noniso_max) h_dPhi_noniso[izt]->Fill(dphi);
	      }
	  } // end loop over bins
	}//end loop over tracks

      }//end loop on clusters. 
      if (ievent % 25000 == 0) {
	histogram0.Draw("e1x0");
	canvas.Update();
	std::cout << "Event # " << ievent << " / " << _tree_event->GetEntries() << std::endl;
      }
    } //end loop over events

    // Write to fout
    std::string opened_files = "";
    for (int iarg = 1; iarg < argc; iarg++) {
        std::string filepath = argv[iarg];
        
        opened_files += "_" + filepath.substr(filepath.find_last_of("/")+1, filepath.find_last_of(".")-filepath.find_last_of("/")-1);
    }
    TFile* fout = new TFile(Form("fout_Corr_config%s.root", opened_files.c_str()),"RECREATE");
    histogram0.Write("DeepPhotonSpectra");
    h_ntrig.Write("ntriggers");
   
    for (int izt = 0; izt<nztbins; izt++){
      h_dPhi_iso[izt]->SetMinimum(0.0);
      h_dPhi_iso[izt]->Write();
      h_dPhi_noniso[izt]->SetMinimum(0.0);
      h_dPhi_noniso[izt]->Write();  
    }
    fout->Close();     
      
  }//end loop over samples

  std::cout << " ending " << std::endl;
  return EXIT_SUCCESS;
}
