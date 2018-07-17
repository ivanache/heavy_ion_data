// This program takes histograms from a root output file and makes them into a graph
// Author: Ivan Chernyshev
// Date Created: 5/11/2018

#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"

// Usage: Insert file name as the first argument, number of rows and columns in the canvas are next, followed by all of the histograms. Syntax for histogram insertion: there is one argument per canvas to be created, with a limit of 4 canvases. For graphs within a canvas: ; delimits graphs and , delimits histograms
void grapher(std::string filename, int canvas_rows, int canvas_cols, std::string canvas1, std::string canvas2 = NULL, std::string canvas3 = NULL, std::string canvas4 = NULL) {
    TFile* fIn = new TFile(filename.c_str(), "READ");
    
    // The vectors of histogram names and graphs to be graphed
    std::vector<char*> canvas1graphs;
    std::vector<char*> canvas2graphs;
    std::vector<char*> canvas3graphs;
    std::vector<char*> canvas4graphs;
    
    std::vector<std::vector<char*>> canvas1hists;
    std::vector<std::vector<char*>> canvas2hists;
    std::vector<std::vector<char*>> canvas3hists;
    std::vector<std::vector<char*>> canvas4hists;
    
    // Read the contents of the canvas arguments into the corresponding vectors. Note that it does require at least one element in the canvas argument, so a little re-working will be needed for implementation outside canvas1
    char* canvasgraphs;
    canvasgraphs = strtok((char*) canvas1.c_str(), ";");
    canvas1graphs.push_back(canvasgraphs);
    
    canvasgraphs = strtok(NULL, ";");
    while(canvasgraphs != NULL) {
        canvas1graphs.push_back(canvasgraphs);
        canvasgraphs = strtok(NULL, ";");
    }
    
    for(char* graph : canvas1graphs) {
        std::vector<char*> temp;
        char* canvashists;
        canvashists = strtok(graph, ",");
        temp.push_back(canvashists);
        
        canvashists = strtok(NULL, ",");
        while(canvashists != NULL) {
            temp.push_back(canvashists);
            canvashists = strtok(NULL, ",");
        }
        
        canvas1hists.push_back(temp);
    }
    
    std::cout << "Canvas 1: graphs and histograms contained within the canvas: ";
    for(std::vector<char*> graph : canvas1hists) {
        for (char* histogram : graph) {
            std::cout << histogram << ",";
        }
        std::cout<< ";";
    }
    std::cout << std::endl;
}
