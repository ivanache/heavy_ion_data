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

int num_of_files = 15;


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

        
        std::cout << " Total Number of entries in TTree: " << _tree_event->GetEntries() << std::endl;
        
        // New file
        TFile *newfile = new TFile("13def_mixedadded.root", "RECREATE");
        TTree *newtree = _tree_event->CloneTree(0);
        
        //new branch: mixed_events
        Float_t mixed_events[NTRACK_MAX];
        newtree->Branch("mixed_events", mixed_events, "mixed_events[300]/L"); // One more entry needed for this to work
        
        std::cout<< "New branch successfully created " <<std::endl;
        
        // Get the mixed event textfiles
        std::ifstream mixed_textfiles[num_of_files];
        for(int i = 0; i < num_of_files; i++) {
            std::ostringstream filename;
            filename << Form("4GeVpairs_v1_%i_%i.txt", i*20, ((i+1)*20)-1);
            mixed_textfiles[i].open(filename.str());
        }
        
        
        const Long64_t nevents = _tree_event->GetEntries();
        // Loop over events
        for(Long64_t ievent = 0; ievent < nevents ; ievent++){
            _tree_event->GetEntry(ievent);
            // Get the appropriate line from each file, break out of the loop if you hit an empty file
            std::string eventlines[num_of_files];
            bool event_end = false;
            for(int i = 0; i < num_of_files; i++) {
                getline(mixed_textfiles[i], eventlines[i]);
                if (eventlines[i] == "") {
                    event_end = true;
                    break;
                }
            }
            if(event_end)
                break;
            
            //try {
            std::string mixednum_string;
            long mixednum;
            std::istringstream parsers[num_of_files];
            for(int i = 0; i < num_of_files; i++) {
                parsers[i].str(eventlines[i]);
            }
            int currentindex;
            // Loop over mixed events, fill the mixed_events histogram while at it
            for(int m = 0; m <300; m++) {
                currentindex = m/20;
                getline(parsers[currentindex], mixednum_string, '\t');
                mixed_events[m] = stol(mixednum_string);
            }
            //}
            //catch(std::invalid_argument) {
                //std::cout << "std_invalid_argument thrown at Event " << ievent << std::endl;
            //}
            newtree->Fill();
            int numentries = newtree->GetEntries();
            if (ievent % 10000 == 0) {
                std::cout << "Event number: " << ievent << " GetEntries entry:" << numentries << std::endl;
            }
        }
        std::cout << "Successfully exited the eventloop" << std::endl;
        newtree->AutoSave();
        std::cout << "Successful autosave" <<std::endl;
        delete newfile;
        delete file;
        std::cout << "Deleted newfile" << std::endl;
    }
    std::cout << " ending " << std::endl;
    return EXIT_SUCCESS;
}
