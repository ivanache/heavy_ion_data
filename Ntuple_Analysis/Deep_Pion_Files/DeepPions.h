//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 27 13:33:22 2017 by ROOT version 6.08/00
// from TTree _tree_event/
// found on file: 16c3b_small.root
//////////////////////////////////////////////////////////

#ifndef DeepPions_h
#define DeepPions_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DeepPions {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  Double_t        primary_vertex[3];
  Double_t        primary_vertex_sigma[3];
  Int_t           primary_vertex_ncontributor;
  Double_t        primary_vertex_spd[3];
  Double_t        primary_vertex_spd_sigma[3];
  Int_t           primary_vertex_spd_ncontributor;
  Int_t           npileup_vertex_spd;
  Int_t           pileup_vertex_spd_ncontributor;
  Double_t        pileup_vertex_spd_min_z_distance;
  Float_t         eg_primary_vertex[3];
  UInt_t          ncluster;
  Float_t         cluster_e[62];   //[ncluster]
  Float_t         cluster_pt[62];   //[ncluster]
  Float_t         cluster_eta[62];   //[ncluster]
  Float_t         cluster_phi[62];   //[ncluster]
  Float_t         cluster_lambda_square[62][2];   //[ncluster]
  Float_t         cluster_tof[62];   //[ncluster]
  Int_t           cluster_ncell[62];   //[ncluster]
  UShort_t        cluster_cell_id_max[62];   //[ncluster]
  Float_t         cluster_e_max[62];   //[ncluster]
  Float_t         cluster_e_cross[62];   //[ncluster]
  UInt_t          cluster_nmc_truth[62];   //[ncluster]
  UShort_t        cluster_mc_truth_index[62][32];   //[ncluster]
  Float_t         cluster_iso_tpc_01[62];   //[ncluster]
  Float_t         cluster_iso_tpc_02[62];   //[ncluster]
  Float_t         cluster_iso_tpc_03[62];   //[ncluster]
  Float_t         cluster_iso_tpc_04[62];   //[ncluster]
  Float_t         cluster_iso_its_01[62];   //[ncluster]
  Float_t         cluster_iso_its_02[62];   //[ncluster]
  Float_t         cluster_iso_its_03[62];   //[ncluster]
  Float_t         cluster_iso_its_04[62];   //[ncluster]
  Float_t         cluster_frixione_tpc_04_02[62];   //[ncluster]
  Float_t         cluster_frixione_tpc_04_05[62];   //[ncluster]
  Float_t         cluster_frixione_tpc_04_10[62];   //[ncluster]
  Float_t         cluster_frixione_its_04_02[62];   //[ncluster]
  Float_t         cluster_frixione_its_04_05[62];   //[ncluster]
  Float_t         cluster_frixione_its_04_10[62];   //[ncluster]
  Float_t         cluster_iso_01_truth[62];   //[ncluster]
  Float_t         cluster_iso_02_truth[62];   //[ncluster]
  Float_t         cluster_iso_03_truth[62];   //[ncluster]
  Float_t         cluster_iso_04_truth[62];   //[ncluster]
  Float_t         cluster_frixione_04_02_truth[62];   //[ncluster]
  Float_t         cluster_frixione_04_05_truth[62];   //[ncluster]
  Float_t         cluster_frixione_04_10_truth[62];   //[ncluster]
  Float_t         cluster_s_nphoton[62][4];   //[ncluster]
  Float_t         cluster_s_ncharged_hadron[62][4];   //[ncluster]
  UShort_t        cell_mc_truth_index[17664];
  UInt_t          ntrack;
  Float_t         track_e[182];   //[ntrack]
  Float_t         track_pt[182];   //[ntrack]
  Float_t         track_eta[182];   //[ntrack]
  Float_t         track_phi[182];   //[ntrack]
  Char_t          track_charge[182];   //[ntrack]
  UChar_t         track_quality[182];   //[ntrack]
  Float_t         track_tpc_dedx[182];   //[ntrack]
  Float_t         track_tpc_length_active_zone[182];   //[ntrack]
  UChar_t         track_tpc_xrow[182];   //[ntrack]
  UChar_t         track_tpc_ncluster[182];   //[ntrack]
  UChar_t         track_tpc_ncluster_dedx[182];   //[ntrack]
  UChar_t         track_tpc_ncluster_findable[182];   //[ntrack]
  UChar_t         track_its_ncluster[182];   //[ntrack]
  Float_t         track_its_chi_square[182];   //[ntrack]
  Float_t         track_dca_xy[182];   //[ntrack]
  Float_t         track_dca_z[182];   //[ntrack]
  UShort_t        track_mc_truth_index[182];   //[ntrack]
  Float_t         track_voronoi_area[182];   //[ntrack]
  UInt_t          nmuon_track;
  UShort_t        muon_track_mc_truth_index[10];   //[nmuon_track]
  UInt_t          nmc_truth;
  Float_t         mc_truth_e[747];   //[nmc_truth]
  Float_t         mc_truth_pt[747];   //[nmc_truth]
  Float_t         mc_truth_eta[747];   //[nmc_truth]
  Float_t         mc_truth_phi[747];   //[nmc_truth]
  Char_t          mc_truth_charge[747];   //[nmc_truth]
  Short_t         mc_truth_pdg_code[747];   //[nmc_truth]
  UChar_t         mc_truth_status[747];   //[nmc_truth]
  UChar_t         mc_truth_generator_index[747];   //[nmc_truth]
  Short_t         mc_truth_first_parent_pdg_code[747];   //[nmc_truth]
  Float_t         mc_truth_first_parent_e[747];   //[nmc_truth]
  Float_t         mc_truth_first_parent_pt[747];   //[nmc_truth]
  Float_t         mc_truth_first_parent_eta[747];   //[nmc_truth]
  Float_t         mc_truth_first_parent_phi[747];   //[nmc_truth]
  UShort_t        mc_truth_sibling_index[747];   //[nmc_truth]
  Float_t         cluster_NN1[62];   //[ncluster]
  Float_t         cluster_NN2[62];   //[ncluster]
  Float_t         cluster_Lambda[62];   //[ncluster]

  // List of branches
  TBranch        *b_primary_vertex;   //!
  TBranch        *b_primary_vertex_sigma;   //!
  TBranch        *b_primary_vertex_ncontributor;   //!
  TBranch        *b_primary_vertex_spd;   //!
  TBranch        *b_primary_vertex_spd_sigma;   //!
  TBranch        *b_primary_vertex_spd_ncontributor;   //!
  TBranch        *b_npileup_vertex_spd;   //!
  TBranch        *b_pileup_vertex_spd_ncontributor;   //!
  TBranch        *b_pileup_vertex_spd_min_z_distance;   //!
  TBranch        *b_eg_primary_vertex;   //!
  TBranch        *b_ncluster;   //!
  TBranch        *b_cluster_e;   //!
  TBranch        *b_cluster_pt;   //!
  TBranch        *b_cluster_eta;   //!
  TBranch        *b_cluster_phi;   //!
  TBranch        *b_cluster_lambda_square;   //!
  TBranch        *b_cluster_tof;   //!
  TBranch        *b_cluster_ncell;   //!
  TBranch        *b_cluster_cell_id_max;   //!
  TBranch        *b_cluster_e_max;   //!
  TBranch        *b_cluster_e_cross;   //!
  TBranch        *b_cluster_nmc_truth;   //!
  TBranch        *b_cluster_mc_truth_index;   //!
  TBranch        *b_cluster_iso_tpc_01;   //!
  TBranch        *b_cluster_iso_tpc_02;   //!
  TBranch        *b_cluster_iso_tpc_03;   //!
  TBranch        *b_cluster_iso_tpc_04;   //!
  TBranch        *b_cluster_iso_its_01;   //!
  TBranch        *b_cluster_iso_its_02;   //!
  TBranch        *b_cluster_iso_its_03;   //!
  TBranch        *b_cluster_iso_its_04;   //!
  TBranch        *b_cluster_frixione_tpc_04_02;   //!
  TBranch        *b_cluster_frixione_tpc_04_05;   //!
  TBranch        *b_cluster_frixione_tpc_04_10;   //!
  TBranch        *b_cluster_frixione_its_04_02;   //!
  TBranch        *b_cluster_frixione_its_04_05;   //!
  TBranch        *b_cluster_frixione_its_04_10;   //!
  TBranch        *b_cluster_iso_01_truth;   //!
  TBranch        *b_cluster_iso_02_truth;   //!
  TBranch        *b_cluster_iso_03_truth;   //!
  TBranch        *b_cluster_iso_04_truth;   //!
  TBranch        *b_cluster_frixione_04_02_truth;   //!
  TBranch        *b_cluster_frixione_04_05_truth;   //!
  TBranch        *b_cluster_frixione_04_10_truth;   //!
  TBranch        *b_cluster_s_nphoton;   //!
  TBranch        *b_cluster_s_ncharged_hadron;   //!
  TBranch        *b_cell_mc_truth_index;   //!
  TBranch        *b_ntrack;   //!
  TBranch        *b_track_e;   //!
  TBranch        *b_track_pt;   //!
  TBranch        *b_track_eta;   //!
  TBranch        *b_track_phi;   //!
  TBranch        *b_track_charge;   //!
  TBranch        *b_track_quality;   //!
  TBranch        *b_track_tpc_dedx;   //!
  TBranch        *b_track_tpc_length_active_zone;   //!
  TBranch        *b_track_tpc_xrow;   //!
  TBranch        *b_track_tpc_ncluster;   //!
  TBranch        *b_track_tpc_ncluster_dedx;   //!
  TBranch        *b_track_tpc_ncluster_findable;   //!
  TBranch        *b_track_its_ncluster;   //!
  TBranch        *b_track_its_chi_square;   //!
  TBranch        *b_track_dca_xy;   //!
  TBranch        *b_track_dca_z;   //!
  TBranch        *b_track_mc_truth_index;   //!
  TBranch        *b_track_voronoi_area;   //!
  TBranch        *b_nmuon_track;   //!
  TBranch        *b_muon_track_mc_truth_index;   //!
  TBranch        *b_nmc_truth;   //!
  TBranch        *b_mc_truth_e;   //!
  TBranch        *b_mc_truth_pt;   //!
  TBranch        *b_mc_truth_eta;   //!
  TBranch        *b_mc_truth_phi;   //!
  TBranch        *b_mc_truth_charge;   //!
  TBranch        *b_mc_truth_pdg_code;   //!
  TBranch        *b_mc_truth_status;   //!
  TBranch        *b_mc_truth_generator_index;   //!
  TBranch        *b_mc_truth_first_parent_pdg_code;   //!
  TBranch        *b_mc_truth_first_parent_e;   //!
  TBranch        *b_mc_truth_first_parent_pt;   //!
  TBranch        *b_mc_truth_first_parent_eta;   //!
  TBranch        *b_mc_truth_first_parent_phi;   //!
  TBranch        *b_mc_truth_sibling_index;   //!
  TBranch        *b_cluster_NN1;   //!
  TBranch        *b_cluster_NN2;   //!
  TBranch        *b_cluster_Lambda;   //!

  DeepPions(TTree *tree=0);
  virtual ~DeepPions();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  virtual void     Efficiency();
};

#endif

#ifdef DeepPions_cxx
DeepPions::DeepPions(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("13d_small.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("13d_small.root");
    }
    f->GetObject("_tree_event",tree);

  }
  Init(tree);
}

DeepPions::~DeepPions()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t DeepPions::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t DeepPions::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void DeepPions::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("primary_vertex", primary_vertex, &b_primary_vertex);
  fChain->SetBranchAddress("primary_vertex_sigma", primary_vertex_sigma, &b_primary_vertex_sigma);
  fChain->SetBranchAddress("primary_vertex_ncontributor", &primary_vertex_ncontributor, &b_primary_vertex_ncontributor);
  fChain->SetBranchAddress("primary_vertex_spd", primary_vertex_spd, &b_primary_vertex_spd);
  fChain->SetBranchAddress("primary_vertex_spd_sigma", primary_vertex_spd_sigma, &b_primary_vertex_spd_sigma);
  fChain->SetBranchAddress("primary_vertex_spd_ncontributor", &primary_vertex_spd_ncontributor, &b_primary_vertex_spd_ncontributor);
  fChain->SetBranchAddress("npileup_vertex_spd", &npileup_vertex_spd, &b_npileup_vertex_spd);
  fChain->SetBranchAddress("pileup_vertex_spd_ncontributor", &pileup_vertex_spd_ncontributor, &b_pileup_vertex_spd_ncontributor);
  fChain->SetBranchAddress("pileup_vertex_spd_min_z_distance", &pileup_vertex_spd_min_z_distance, &b_pileup_vertex_spd_min_z_distance);
  fChain->SetBranchAddress("eg_primary_vertex", eg_primary_vertex, &b_eg_primary_vertex);
  fChain->SetBranchAddress("ncluster", &ncluster, &b_ncluster);
  fChain->SetBranchAddress("cluster_e", cluster_e, &b_cluster_e);
  fChain->SetBranchAddress("cluster_pt", cluster_pt, &b_cluster_pt);
  fChain->SetBranchAddress("cluster_eta", cluster_eta, &b_cluster_eta);
  fChain->SetBranchAddress("cluster_phi", cluster_phi, &b_cluster_phi);
  fChain->SetBranchAddress("cluster_lambda_square", cluster_lambda_square, &b_cluster_lambda_square);
  fChain->SetBranchAddress("cluster_tof", cluster_tof, &b_cluster_tof);
  fChain->SetBranchAddress("cluster_ncell", cluster_ncell, &b_cluster_ncell);
  fChain->SetBranchAddress("cluster_cell_id_max", cluster_cell_id_max, &b_cluster_cell_id_max);
  fChain->SetBranchAddress("cluster_e_max", cluster_e_max, &b_cluster_e_max);
  fChain->SetBranchAddress("cluster_e_cross", cluster_e_cross, &b_cluster_e_cross);
  fChain->SetBranchAddress("cluster_nmc_truth", cluster_nmc_truth, &b_cluster_nmc_truth);
  fChain->SetBranchAddress("cluster_mc_truth_index", cluster_mc_truth_index, &b_cluster_mc_truth_index);
  fChain->SetBranchAddress("cluster_iso_tpc_01", cluster_iso_tpc_01, &b_cluster_iso_tpc_01);
  fChain->SetBranchAddress("cluster_iso_tpc_02", cluster_iso_tpc_02, &b_cluster_iso_tpc_02);
  fChain->SetBranchAddress("cluster_iso_tpc_03", cluster_iso_tpc_03, &b_cluster_iso_tpc_03);
  fChain->SetBranchAddress("cluster_iso_tpc_04", cluster_iso_tpc_04, &b_cluster_iso_tpc_04);
  fChain->SetBranchAddress("cluster_iso_its_01", cluster_iso_its_01, &b_cluster_iso_its_01);
  fChain->SetBranchAddress("cluster_iso_its_02", cluster_iso_its_02, &b_cluster_iso_its_02);
  fChain->SetBranchAddress("cluster_iso_its_03", cluster_iso_its_03, &b_cluster_iso_its_03);
  fChain->SetBranchAddress("cluster_iso_its_04", cluster_iso_its_04, &b_cluster_iso_its_04);
  fChain->SetBranchAddress("cluster_frixione_tpc_04_02", cluster_frixione_tpc_04_02, &b_cluster_frixione_tpc_04_02);
  fChain->SetBranchAddress("cluster_frixione_tpc_04_05", cluster_frixione_tpc_04_05, &b_cluster_frixione_tpc_04_05);
  fChain->SetBranchAddress("cluster_frixione_tpc_04_10", cluster_frixione_tpc_04_10, &b_cluster_frixione_tpc_04_10);
  fChain->SetBranchAddress("cluster_frixione_its_04_02", cluster_frixione_its_04_02, &b_cluster_frixione_its_04_02);
  fChain->SetBranchAddress("cluster_frixione_its_04_05", cluster_frixione_its_04_05, &b_cluster_frixione_its_04_05);
  fChain->SetBranchAddress("cluster_frixione_its_04_10", cluster_frixione_its_04_10, &b_cluster_frixione_its_04_10);
  fChain->SetBranchAddress("cluster_iso_01_truth", cluster_iso_01_truth, &b_cluster_iso_01_truth);
  fChain->SetBranchAddress("cluster_iso_02_truth", cluster_iso_02_truth, &b_cluster_iso_02_truth);
  fChain->SetBranchAddress("cluster_iso_03_truth", cluster_iso_03_truth, &b_cluster_iso_03_truth);
  fChain->SetBranchAddress("cluster_iso_04_truth", cluster_iso_04_truth, &b_cluster_iso_04_truth);
  fChain->SetBranchAddress("cluster_frixione_04_02_truth", cluster_frixione_04_02_truth, &b_cluster_frixione_04_02_truth);
  fChain->SetBranchAddress("cluster_frixione_04_05_truth", cluster_frixione_04_05_truth, &b_cluster_frixione_04_05_truth);
  fChain->SetBranchAddress("cluster_frixione_04_10_truth", cluster_frixione_04_10_truth, &b_cluster_frixione_04_10_truth);
  fChain->SetBranchAddress("cluster_s_nphoton", cluster_s_nphoton, &b_cluster_s_nphoton);
  fChain->SetBranchAddress("cluster_s_ncharged_hadron", cluster_s_ncharged_hadron, &b_cluster_s_ncharged_hadron);
  fChain->SetBranchAddress("cell_mc_truth_index", cell_mc_truth_index, &b_cell_mc_truth_index);
  fChain->SetBranchAddress("ntrack", &ntrack, &b_ntrack);
  fChain->SetBranchAddress("track_e", track_e, &b_track_e);
  fChain->SetBranchAddress("track_pt", track_pt, &b_track_pt);
  fChain->SetBranchAddress("track_eta", track_eta, &b_track_eta);
  fChain->SetBranchAddress("track_phi", track_phi, &b_track_phi);
  fChain->SetBranchAddress("track_charge", track_charge, &b_track_charge);
  fChain->SetBranchAddress("track_quality", track_quality, &b_track_quality);
  fChain->SetBranchAddress("track_tpc_dedx", track_tpc_dedx, &b_track_tpc_dedx);
  fChain->SetBranchAddress("track_tpc_length_active_zone", track_tpc_length_active_zone, &b_track_tpc_length_active_zone);
  fChain->SetBranchAddress("track_tpc_xrow", track_tpc_xrow, &b_track_tpc_xrow);
  fChain->SetBranchAddress("track_tpc_ncluster", track_tpc_ncluster, &b_track_tpc_ncluster);
  fChain->SetBranchAddress("track_tpc_ncluster_dedx", track_tpc_ncluster_dedx, &b_track_tpc_ncluster_dedx);
  fChain->SetBranchAddress("track_tpc_ncluster_findable", track_tpc_ncluster_findable, &b_track_tpc_ncluster_findable);
  fChain->SetBranchAddress("track_its_ncluster", track_its_ncluster, &b_track_its_ncluster);
  fChain->SetBranchAddress("track_its_chi_square", track_its_chi_square, &b_track_its_chi_square);
  fChain->SetBranchAddress("track_dca_xy", track_dca_xy, &b_track_dca_xy);
  fChain->SetBranchAddress("track_dca_z", track_dca_z, &b_track_dca_z);
  fChain->SetBranchAddress("track_mc_truth_index", track_mc_truth_index, &b_track_mc_truth_index);
  fChain->SetBranchAddress("track_voronoi_area", track_voronoi_area, &b_track_voronoi_area);
  fChain->SetBranchAddress("nmuon_track", &nmuon_track, &b_nmuon_track);
  fChain->SetBranchAddress("muon_track_mc_truth_index", muon_track_mc_truth_index, &b_muon_track_mc_truth_index);
  fChain->SetBranchAddress("nmc_truth", &nmc_truth, &b_nmc_truth);
  fChain->SetBranchAddress("mc_truth_e", mc_truth_e, &b_mc_truth_e);
  fChain->SetBranchAddress("mc_truth_pt", mc_truth_pt, &b_mc_truth_pt);
  fChain->SetBranchAddress("mc_truth_eta", mc_truth_eta, &b_mc_truth_eta);
  fChain->SetBranchAddress("mc_truth_phi", mc_truth_phi, &b_mc_truth_phi);
  fChain->SetBranchAddress("mc_truth_charge", mc_truth_charge, &b_mc_truth_charge);
  fChain->SetBranchAddress("mc_truth_pdg_code", mc_truth_pdg_code, &b_mc_truth_pdg_code);
  fChain->SetBranchAddress("mc_truth_status", mc_truth_status, &b_mc_truth_status);
  fChain->SetBranchAddress("mc_truth_generator_index", mc_truth_generator_index, &b_mc_truth_generator_index);
  fChain->SetBranchAddress("mc_truth_first_parent_pdg_code", mc_truth_first_parent_pdg_code, &b_mc_truth_first_parent_pdg_code);
  fChain->SetBranchAddress("mc_truth_first_parent_e", mc_truth_first_parent_e, &b_mc_truth_first_parent_e);
  fChain->SetBranchAddress("mc_truth_first_parent_pt", mc_truth_first_parent_pt, &b_mc_truth_first_parent_pt);
  fChain->SetBranchAddress("mc_truth_first_parent_eta", mc_truth_first_parent_eta, &b_mc_truth_first_parent_eta);
  fChain->SetBranchAddress("mc_truth_first_parent_phi", mc_truth_first_parent_phi, &b_mc_truth_first_parent_phi);
  fChain->SetBranchAddress("mc_truth_sibling_index", mc_truth_sibling_index, &b_mc_truth_sibling_index);
  fChain->SetBranchAddress("cluster_NN1", cluster_NN1, &b_cluster_NN1);
  fChain->SetBranchAddress("cluster_NN2", cluster_NN2, &b_cluster_NN2);
  fChain->SetBranchAddress("cluster_Lambda", cluster_Lambda, &b_cluster_Lambda);
  Notify();
}

Bool_t DeepPions::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void DeepPions::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t DeepPions::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef DeepPions_cxx

