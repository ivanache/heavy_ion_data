// The code that outputs PNG plots of the Fernando code correlations
// Author: Ivan Chernyshev

#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

#include "TH1F.h"

const int MAX_INPUT_LENGTH = 200;

void fernando_corr_converter() {
    TFile* frixione_corr = new TFile("fout_mixed_frixione.root", "READ");
    
    // Get Atlas Style and define the Canvas
    TCanvas* c = new TCanvas();
    gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();
    
    /**
     Set up Zt bins
     */
    
    // Zt bins
    // Read configuration file for the DNN_min and DNN_max variables
    FILE* config = fopen("Corr_config.yaml", "r");
    
    int nztbins = 7;
    float* ztbins;
    ztbins = new float[nztbins+1];
    ztbins[0] = 0.0; ztbins[1] = 0.1; ztbins[2] = 0.2; ztbins[3] = 0.4; ztbins[4] = 0.6; ztbins[5] = 0.8; ztbins[6] = 1.0; ztbins[7] = 1.2;
    
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
        
        if (strcmp(key, "Zt_bins") == 0) {
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
        else {
            std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
        }
    }
    fclose(config);
    
    
    // For all Zt bins, get the data histograms and graph them
    for (int izt = 0; izt<nztbins; izt++) {
        TH1F* iso_phimap = 0;
        std::cout << Form("dPhi_iso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]) << std::endl;
        frixione_corr->GetObject(Form("dPhi_iso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]), iso_phimap);
        TH1F* noniso_phimap = 0;
        frixione_corr->GetObject(Form("dPhi_iso_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]), noniso_phimap);
        TH2D* corrmap = 0;
        frixione_corr->GetObject(Form("Correlation_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]), corrmap);
        TH2D* corrmap_iso = 0;
        frixione_corr->GetObject(Form("IsoCorrelation_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]), corrmap_iso);
        TH2D* corrmap_noniso = 0;
        frixione_corr->GetObject(Form("AntiIsoCorrelation_ztmin%1.0f_ztmax%1.0f",10*ztbins[izt],10*ztbins[izt+1]), corrmap_noniso);
        
        iso_phimap->SetTitle(Form("Isolated Spectrum for Zt %1.1f - %1.1f; #Delta#phi/#pi [rad]; entries", ztbins[izt], ztbins[izt+1]));
        iso_phimap->Draw();
        myText(0.2, 0.92, 1, Form("Isolated Spectrum for Zt %1.1f - %1.1f", ztbins[izt], ztbins[izt+1]));
        c->SaveAs(Form("dPhi_iso_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        c->Clear();
        
        noniso_phimap->SetTitle(Form("Nonisolated Spectrum for Zt %1.1f - %1.1f; #Delta#phi/#pi [rad]; entries", ztbins[izt], ztbins[izt+1]));
        noniso_phimap->Draw();
        myText(0.2, 0.92, 1, Form("NonIsolated Spectrum for Zt %1.1f - %1.1f", ztbins[izt], ztbins[izt+1]));
        c->SaveAs(Form("dPhi_iso_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        c->Clear();
        
        corrmap->SetTitle(Form("Correlation Function for Zt %1.1f - %1.1f; #phi (#frac{rad}{#pi}); #eta",ztbins[izt], ztbins[izt+1]));
        corrmap->Draw("COLZ");
        myText(0.2, 0.92, 1, Form("Correlation Function for Zt %1.1f - %1.1f",ztbins[izt], ztbins[izt+1]));
        c->SaveAs(Form("Correlation_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        c->Clear();
        
        corrmap_iso->SetTitle(Form("Isolated Correlation Function for Zt %1.1f - %1.1f; #phi (#frac{rad}{#pi}); #eta",ztbins[izt], ztbins[izt+1]));
        corrmap_iso->Draw("COLZ");
        myText(0.2, 0.92, 1, Form("Isolated Correlation Function for Zt %1.1f - %1.1f",ztbins[izt], ztbins[izt+1]));
        c->SaveAs(Form("IsoCorrelation_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        c->Clear();
        
        corrmap_noniso->SetTitle(Form("Nonisolated Correlation Function for Zt %1.1f - %1.1f; #phi (#frac{rad}{#pi}); #eta",ztbins[izt], ztbins[izt+1]));
        corrmap_noniso->Draw("COLZ");
        myText(0.2, 0.92, 1, Form("Nonisolated Correlation Function for Zt %1.1f - %1.1f",ztbins[izt], ztbins[izt+1]));
        c->SaveAs(Form("AntiIsoCorrelation_ztmin%1.0f_ztmax%1.0f.png",10*ztbins[izt],10*ztbins[izt+1]));
        c->Clear();
    }
    c->Close();
}
