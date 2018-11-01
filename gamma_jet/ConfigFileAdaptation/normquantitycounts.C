// This ROOT macro takes a set of output gamma-jet ROOT files and counts the number of events, clusters, and gamma-jet pairs
// Author: Ivan Chernyshev; Date: 11/1/2018

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <string>

const int num_of_files = 6;
const std::string files_to_analyze[num_of_files] = {"GammaJet_config_clusptmin15.0_clusptmax30.0_JETPTMIN_10.0_DATANAME__13d_0GeVTrack_paired_13e_0GeVTrack_paired_13f_0GeVTrack_paired_PHOTONSELECT_DNN.root", "GammaJet_config_clusptmin15.0_clusptmax30.0_JETPTMIN_10.0_DATANAME__13d_0GeVTrack_paired_13e_0GeVTrack_paired_13f_0GeVTrack_paired_PHOTONSELECT_EmaxOverEcluster.root", "GammaJet_config_clusptmin15.0_clusptmax30.0_JETPTMIN_10.0_DATANAME__13d_0GeVTrack_paired_13e_0GeVTrack_paired_13f_0GeVTrack_paired_PHOTONSELECT_Lambda0.root", "GammaJet_config_clusptmin15.0_clusptmax30.0_JETPTMIN_10.0_DATANAME__17q_0GeVTrack_paired_PHOTONSELECT_DNN.root", "GammaJet_config_clusptmin15.0_clusptmax30.0_JETPTMIN_10.0_DATANAME__17q_0GeVTrack_paired_PHOTONSELECT_EmaxOverEcluster.root", "GammaJet_config_clusptmin15.0_clusptmax30.0_JETPTMIN_10.0_DATANAME__17q_0GeVTrack_paired_PHOTONSELECT_Lambda0.root"};

void normquantitycounts() {
    for (std::string filename : files_to_analyze) {
        
        // Get file name, and set event, cluster, and jet cutflows
        TFile* file = new TFile(filename.c_str(), "READ");
        TH1D* event_cutflow = 0;
        TH1D* cluster_cutflow = 0;
        TH1D* jet_cutflow = 0;
        file->GetObject("EventCutFlow", event_cutflow);
        file->GetObject("ClusterCutFlow", cluster_cutflow);
        file->GetObject("JetCutFlow", jet_cutflow);
        
        // Get num of events, signal region clusters, background region clusters, signal region pairs (or jets, because there's really no difference)
        double num_of_events = event_cutflow->GetBinContent(5);
        double num_of_sig_clusters = cluster_cutflow->GetBinContent(9);
        double num_of_bkg_clusters = cluster_cutflow->GetBinContent(10);
        double num_of_sig_pairs = jet_cutflow->GetBinContent(4);
        double num_of_bkg_pairs = jet_cutflow->GetBinContent(5);
        
        std::cout << "Statistics for " << filename << ":" << std::endl << "Events: " << num_of_events << std::endl << "Signal clusters: " << num_of_sig_clusters << " Background clusters: " << num_of_bkg_clusters << std::endl << "Signal pairs: " << num_of_sig_pairs << " Background pairs: " << num_of_bkg_pairs << std::endl;
    }
}
