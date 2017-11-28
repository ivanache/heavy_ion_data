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


void to_sm_nphi(unsigned int &sm, unsigned int &nphi,
		unsigned int n)
{
  sm = n < 11520 ? n / 1152 :
    n < 12288 ? 10 + (n - 11520) / 384 :
    n < 16896 ? 12 + (n - 12288) / 768 :
    18 + (n - 16896) / 384;
  nphi = sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;
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

void cell_3_3(unsigned int n_3_3[], const unsigned int n,
	      const unsigned int ld = 3)
{
  unsigned int sm;
  unsigned int nphi;

  to_sm_nphi(sm, nphi, n);

  if (n % 2 == 0) {
    n_3_3[0 * ld + 0] = n - 1;
    n_3_3[0 * ld + 1] = n + 1;
    n_3_3[0 * ld + 2] = n + 3;
  }
  else {
    n_3_3[0 * ld + 0] = n - 2 * nphi - 3;
    n_3_3[0 * ld + 1] = n - 2 * nphi - 1;
    n_3_3[0 * ld + 2] = n - 2 * nphi + 1;
  }
  n_3_3[1 * ld + 0] = n - 2;
  n_3_3[1 * ld + 1] = n;
  n_3_3[1 * ld + 2] = n + 2;
  if (n % 2 == 0) {
    n_3_3[2 * ld + 0] = n + 2 * nphi - 1;
    n_3_3[2 * ld + 1] = n + 2 * nphi + 1;
    n_3_3[2 * ld + 2] = n + 2 * nphi + 3;
  }
  else {
    n_3_3[2 * ld + 0] = n - 3;
    n_3_3[2 * ld + 1] = n - 1;
    n_3_3[2 * ld + 2] = n + 1;
  }
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

        TTree *_tree_event = dynamic_cast<TTree *>(file->Get("_tree_event"));

        //TTree *_tree_event = dynamic_cast<TTree *>
        //    (dynamic_cast<TDirectoryFile *>
        //     (file->Get("AliAnalysisTaskNTGJ"))->Get("_tree_event"));
	if (_tree_event == NULL) {
	  std::cout << " fail " << std::endl;
	  exit(EXIT_FAILURE);
        }  
        //_tree_event->Print();

	TApplication application("", &dummyc, dummyv);
	TCanvas canvas("canvas", "");//, 960 + 4, 720 + 28);
	//canvas.SetLogy();

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

        Float_t cluster_s_nphoton[NTRACK_MAX][4];
        unsigned short cluster_mc_truth_index[NTRACK_MAX][32];
        Int_t cluster_ncell[NTRACK_MAX];
        UShort_t  cluster_cell_id_max[NTRACK_MAX];
        Float_t cluster_lambda_square[NTRACK_MAX][2];   


        Float_t cell_e[17664];

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
        _tree_event->SetBranchAddress("cluster_pt", cluster_pt);
        _tree_event->SetBranchAddress("cluster_eta", cluster_eta);
        _tree_event->SetBranchAddress("cluster_phi", cluster_phi);
	_tree_event->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton);
        _tree_event->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index);
        _tree_event->SetBranchAddress("cluster_lambda_square", cluster_lambda_square);
 	_tree_event->SetBranchAddress("cluster_ncell", cluster_ncell);
        _tree_event->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max);
        _tree_event->SetBranchAddress("cell_e", cell_e);

	 
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
    

	//for(Long64_t ievent = 0; ievent < _tree_event->GetEntries() ; ievent++){     
        for(Long64_t ievent = 0; ievent < 1000 ; ievent++){
             _tree_event->GetEntry(ievent);

	      for (ULong64_t n = 0; n < ncluster; n++) {
	          if(TMath::Abs(primary_vertex[2])>10) continue;
	          if(!(cluster_pt[n]>10 and cluster_pt[n]<14)) continue;
                  if(TMath::Abs(cluster_eta[n])>0.6) continue;
                  if(cluster_ncell[n]<3) continue; 
                  if(cluster_e_cross[n]/cluster_e[n]<0.03) continue;
	          unsigned int cell_id_5_5[25];
	          cell_5_5(cell_id_5_5, cluster_cell_id_max[n]);
		  std::cout << Form("%2.1f, %2.2f", cluster_pt[n], cluster_s_nphoton[n][1]);
	          
                  for (size_t i = 0; i < 25; i++) {
		    std::cout <<"," <<Form("%2.2f",100*cell_e[cell_id_5_5[i]]/cluster_e[n]);
		  }
		  std::cout<<std::endl;
	    }//end loop on clusters. 
	   

	} //end loop over events

	TFile* fout = new TFile("fout.root","RECREATE");
	

	fout->Close();     
	  
    }//end loop over samples

    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
