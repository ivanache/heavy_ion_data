/**
   This program clones an NTuple, then uses data contained in text files to addmixed events to the clone
*/
// Author: Ivan Chernyshev; Date: 6/18/2018

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
#include <sstream>

#define NTRACK_MAX (1U << 15)

#include <vector>
#include <math.h>

float c_eta_array[17664];
float c_phi_array[17664];


unsigned int GetSuperModule(const unsigned int n){
    unsigned int sm = n < 11520 ? n / 1152 :
    n < 12288 ? 10 + (n - 11520) / 384 :
    n < 16896 ? 12 + (n - 12288) / 768 :
    18 + (n - 16896) / 384;
    
    
    return sm;
}

void cell_5_5(unsigned int n_5_5[], const unsigned int n,
              const unsigned int ld = 5)
{
    const unsigned int sm = n < 11520 ? n / 1152 :
    n < 12288 ? 10 + (n - 11520) / 384 :
    n < 16896 ? 12 + (n - 12288) / 768 :
    18 + (n - 16896) / 384;
    const unsigned int nphi =
    sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;
    
    n_5_5[0 * ld + 0] = n - 2 * nphi - 4;
    n_5_5[0 * ld + 1] = n - 2 * nphi - 2;
    n_5_5[0 * ld + 2] = n - 2 * nphi;
    n_5_5[0 * ld + 3] = n - 2 * nphi + 2;
    n_5_5[0 * ld + 4] = n - 2 * nphi + 4;
    if (n % 2 == 0) {
        n_5_5[1 * ld + 0] = n - 3;
        n_5_5[1 * ld + 1] = n - 1;
        n_5_5[1 * ld + 2] = n + 1;
        n_5_5[1 * ld + 3] = n + 3;
        n_5_5[1 * ld + 4] = n + 5;
    }
    else {
        n_5_5[1 * ld + 0] = n - 2 * nphi - 5;
        n_5_5[1 * ld + 1] = n - 2 * nphi - 3;
        n_5_5[1 * ld + 2] = n - 2 * nphi - 1;
        n_5_5[1 * ld + 3] = n - 2 * nphi + 1;
        n_5_5[1 * ld + 4] = n - 2 * nphi + 3;
    }
    n_5_5[2 * ld + 0] = n - 4;
    n_5_5[2 * ld + 1] = n - 2;
    n_5_5[2 * ld + 2] = n;
    n_5_5[2 * ld + 3] = n + 2;
    n_5_5[2 * ld + 4] = n + 4;
    if (n % 2 == 0) {
        n_5_5[3 * ld + 0] = n + 2 * nphi - 3;
        n_5_5[3 * ld + 1] = n + 2 * nphi - 1;
        n_5_5[3 * ld + 2] = n + 2 * nphi + 1;
        n_5_5[3 * ld + 3] = n + 2 * nphi + 3;
        n_5_5[3 * ld + 4] = n + 2 * nphi + 5;
    }
    else {
        n_5_5[3 * ld + 0] = n - 5;
        n_5_5[3 * ld + 1] = n - 3;
        n_5_5[3 * ld + 2] = n - 1;
        n_5_5[3 * ld + 3] = n + 1;
        n_5_5[3 * ld + 4] = n + 3;
    }
    n_5_5[4 * ld + 0] = n + 2 * nphi - 4;
    n_5_5[4 * ld + 1] = n + 2 * nphi - 2;
    n_5_5[4 * ld + 2] = n + 2 * nphi;
    n_5_5[4 * ld + 3] = n + 2 * nphi + 2;
    n_5_5[4 * ld + 4] = n + 2 * nphi + 4;
}


int main(int argc, char *argv[])
{
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }
    int dummyc = 1;
    char **dummyv = new char *[1];
    
    dummyv[0] = strdup("main");
    
    for (int iarg = 1; iarg < argc; iarg++) {
        std::cout << "Opening: " << (TString)argv[iarg] << std::endl;
        TFile *file = TFile::Open((TString)argv[iarg]);
        
        if (file == NULL) {
            std::cout << " fail" << std::endl;
            exit(EXIT_FAILURE);
        }
        file->Print();
        
        TTree *_tree_event = NULL;
        _tree_event = dynamic_cast<TTree *> (dynamic_cast<TDirectoryFile *>   (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
        if (_tree_event == NULL) {
            std::cout << "First try did not got (AliAnalysisTaskNTGJ does not exist, trying again" << std::endl;
            _tree_event = dynamic_cast<TTree *> (file->Get("_tree_event"));
            if (_tree_event == NULL) {
                std::cout << " fail " << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        //_tree_event->Print();
        std::cout<<"TTree successfully acquired" << std::endl;
        
        // Existing branch items
        UInt_t ncluster;
        UInt_t cluster_nmc_truth[NTRACK_MAX];
        Float_t cluster_e[NTRACK_MAX];
        Float_t cluster_pt[NTRACK_MAX];
        Float_t cluster_eta[NTRACK_MAX];
        Float_t cluster_phi[NTRACK_MAX];
        
        UShort_t  cluster_cell_id_max[NTRACK_MAX];
        
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
        
        Double_t primary_vertex[3];
        
        Float_t cluster_s_nphoton[NTRACK_MAX][4];
        unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
        Int_t cluster_ncell[NTRACK_MAX];
        
        Float_t cluster_lambda_square[NTRACK_MAX][2];
        Float_t cluster_lambda_square_angle[NTRACK_MAX][2];
        
        UChar_t track_quality[NTRACK_MAX];
        UInt_t ntrack;
        Float_t track_e[NTRACK_MAX];
        Float_t track_pt[NTRACK_MAX];
        Float_t track_eta[NTRACK_MAX];
        Float_t track_phi[NTRACK_MAX];
        UChar_t track_its_ncluster[NTRACK_MAX];
        Float_t track_dca_xy[NTRACK_MAX];
        Float_t track_dca_z[NTRACK_MAX];
        Float_t track_its_chi_square[NTRACK_MAX];
        
        _tree_event->SetBranchAddress("ntrack", &ntrack);
        _tree_event->SetBranchAddress("track_e", track_e);
        _tree_event->SetBranchAddress("track_pt", track_pt);
        _tree_event->SetBranchAddress("track_eta", track_eta);
        _tree_event->SetBranchAddress("track_phi", track_phi);
        _tree_event->SetBranchAddress("track_quality", track_quality);
        _tree_event->SetBranchAddress("track_dca_xy", track_dca_xy);
        _tree_event->SetBranchAddress("track_dca_z", track_dca_z);
        _tree_event->SetBranchAddress("track_its_chi_square", track_its_chi_square);
        _tree_event->SetBranchAddress("track_its_ncluster", track_its_ncluster);
        
        _tree_event->SetBranchAddress("primary_vertex", primary_vertex);
        _tree_event->SetBranchAddress("ncluster", &ncluster);
        _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
        _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
        _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
        _tree_event->SetBranchAddress("cluster_e", cluster_e);
        _tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
        _tree_event->SetBranchAddress("cluster_lambda_square",  cluster_lambda_square);
        _tree_event->SetBranchAddress("cluster_lambda_square_angle",  cluster_lambda_square_angle);
        _tree_event->SetBranchAddress("cluster_nmc_truth", cluster_nmc_truth);
        _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
        
        Float_t cell_e[17664];
        Float_t cell_eta[17664];
        Float_t cell_phi[17664];
        
        _tree_event->SetBranchAddress("cell_e", cell_e);
        _tree_event->SetBranchAddress("cell_eta", cell_eta);
        _tree_event->SetBranchAddress("cell_phi", cell_phi);
        
        _tree_event->SetBranchAddress("nmc_truth",&nmc_truth);
        _tree_event->SetBranchAddress("mc_truth_pt",mc_truth_pt);
        _tree_event->SetBranchAddress("mc_truth_eta",mc_truth_eta);
        _tree_event->SetBranchAddress("mc_truth_phi",mc_truth_phi);
        _tree_event->SetBranchAddress("mc_truth_charge",mc_truth_charge);
        _tree_event->SetBranchAddress("mc_truth_pdg_code",mc_truth_pdg_code);
        _tree_event->SetBranchAddress("mc_truth_first_parent_pdg_code", mc_truth_first_parent_pdg_code);
        _tree_event->SetBranchAddress("mc_truth_first_parent_e", mc_truth_first_parent_e);
        _tree_event->SetBranchAddress("mc_truth_first_parent_pt", mc_truth_first_parent_pt);
        _tree_event->SetBranchAddress("mc_truth_first_parent_eta", mc_truth_first_parent_eta);
        _tree_event->SetBranchAddress("mc_truth_first_parent_phi", mc_truth_first_parent_phi);
        _tree_event->SetBranchAddress("mc_truth_status",  mc_truth_status);
        
        
        std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;
        
        // New file
        TFile *newfile = new TFile("13def_mixedadded.root", "RECREATE");
        TTree *newtree = _tree_event->CloneTree(0);
        
        //new branch: mixed_events
        Float_t mixed_events[NTRACK_MAX];
        newtree->Branch("mixed_events", mixed_events, "mixed_events[300]/L"); // One more entry needed for this to work
        
        std::cout<< "New branch successfully created " <<std::endl;
        
        // Get the mixed event textfiles
        std::ifstream mixed_0_19("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_0_19.txt");
        std::cout<< "First textfile successfully accessed" <<std::endl;
        std::ifstream mixed_20_39("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_20_39.txt");
        std::ifstream mixed_40_59("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_40_59.txt");
        std::ifstream mixed_60_79("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_60_79.txt");
        std::ifstream mixed_80_99("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_80_99.txt");
        std::ifstream mixed_100_119("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_100_119.txt");
        std::ifstream mixed_120_139("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_120_139.txt");
        std::ifstream mixed_140_159("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_140_159.txt");
        std::ifstream mixed_160_179("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_160_179.txt");
        std::ifstream mixed_180_199("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_180_199.txt");
        std::ifstream mixed_200_219("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_200_219.txt");
        std::ifstream mixed_220_239("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_220_239.txt");
        std::ifstream mixed_240_259("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_220_239.txt");
        std::ifstream mixed_260_279("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_220_239.txt");
        std::ifstream mixed_280_299("/project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/InputData/4GeVpairs_v1_280_299.txt");
        
        
        const Long64_t nevents = _tree_event->GetEntries();
        // Loop over events
        for(Long64_t ievent = 0; ievent < nevents ; ievent++){
            _tree_event->GetEntry(ievent);
            // Get the appropriate line from each file
            std::string eventline_0_19;
            getline(mixed_0_19, eventline_0_19);
            std::string eventline_20_39;
            getline(mixed_20_39, eventline_20_39);
            std::string eventline_40_59;
            getline(mixed_40_59, eventline_40_59);
            std::string eventline_60_79;
            getline(mixed_60_79, eventline_60_79);
            std::string eventline_80_99;
            getline(mixed_80_99, eventline_80_99);
            std::string eventline_100_119;
            getline(mixed_100_119, eventline_100_119);
            std::string eventline_120_139;
            getline(mixed_120_139, eventline_120_139);
            std::string eventline_140_159;
            getline(mixed_140_159, eventline_140_159);
            std::string eventline_160_179;
            getline(mixed_160_179, eventline_160_179);
            std::string eventline_180_199;
            getline(mixed_180_199, eventline_180_199);
            std::string eventline_200_219;
            getline(mixed_200_219, eventline_200_219);
            std::string eventline_220_239;
            getline(mixed_220_239, eventline_220_239);
            std::string eventline_240_259;
            getline(mixed_240_259, eventline_240_259);
            std::string eventline_260_279;
            getline(mixed_260_279, eventline_260_279);
            std::string eventline_280_299;
            getline(mixed_280_299, eventline_280_299);
            
            
            std::string mixednum_string;
            long mixednum;
            std::istringstream parser1(eventline_0_19);
            // Loop over mixed events, fill the mixed_events histogram while at it
            for(int m = 0; m <20; m++) {
                getline(parser1, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser2(eventline_20_39);
            for(int m = 20; m <40; m++) {
                getline(parser2, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser3(eventline_40_59);
            for(int m = 40; m <60; m++) {
                getline(parser3, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser4(eventline_60_79);
            for(int m = 60; m <80; m++) {
                getline(parser4, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser5(eventline_80_99);
            for(int m = 80; m <100; m++) {
                getline(parser5, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser6(eventline_100_119);
            for(int m = 100; m <120; m++) {
                getline(parser6, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser7(eventline_120_139);
            for(int m = 120; m <140; m++) {
                getline(parser7, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser8(eventline_140_159);
            for(int m = 140; m <160; m++) {
                getline(parser8, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser9(eventline_160_179);
            for(int m = 160; m <180; m++) {
                getline(parser9, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser10(eventline_180_199);
            for(int m = 180; m <200; m++) {
                getline(parser10, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser11(eventline_200_219);
            for(int m = 200; m <220; m++) {
                getline(parser11, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser12(eventline_220_239);
            for(int m = 220; m <240; m++) {
                getline(parser12, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser13(eventline_240_259);
            for(int m = 240; m <260; m++) {
                getline(parser13, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser14(eventline_260_279);
            for(int m = 260; m <280; m++) {
                getline(parser14, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            std::istringstream parser15(eventline_280_299);
            for(int m = 280; m <300; m++) {
                getline(parser15, mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            
            newtree->Fill();
            int numentries = newtree->GetEntries();
            if (ievent % 10000 == 0) {
                std::cout << "Event number: " << ievent << " GetEntries entry:" << numentries << std::endl;
            }
        }
        newtree->AutoSave();
        delete newfile;
        delete newtree;
    }
    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
