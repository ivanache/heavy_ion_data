// This program takes in 1-D correlation functions and graphs them on one graph, usually signal and background
// It obtains the functions from output root files from a same-event correlator, a mixed-event correlator, or a combined-event correlator
// Author: Ivan Chernyshev; Date: 5/3/2018
#include <TFile.h>

const int MAX_INPUT_LENGTH = 1000;

void correlationgrapher() {
    // Use config file to get the files to graph
    FILE* config = fopen("Graph_config.yaml", "r");
    std::string*** graphs_to_graph = 0;
    std::string** graph_titles = 0;
    int* colors = 0;
    
    int hists_per_graph = 0;
    int graphs_per_canvas = 0;
    int canvases_to_graph = 0;
    int tiles_per_canvas = 0;
    int canvases_to_title = 0;
    int num_of_colors = 0;
    
    int canvas_rows = 0;
    int canvas_cols = 0;
    
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
        
        if (strcmp(key, "Graphs_to_open") == 0) {
            canvases_to_graph = -1;
            for(const char *v = value; *v != '}';) {
                graphs_per_canvas = -1;
                while (*v != '}' && *v != '{') {
                    v++;
                }
                canvases_to_graph++;
                for(; *v != '}';) {
                    hists_per_graph = -1;
                    while (*v != '}' && *v != '{') {
                        v++;
                    }
                    graphs_per_canvas++;
                    for(; *v != '}';) {
                        while (*v != '}' && *v != '"') {
                            v++;
                        }
                        if (*v != '}')
                            break;
                        hists_per_graph++;
                        while (*v != '}' && *v != '"') {
                            v++;
                        }
                        if (*v != '}')
                            break;
                        v++;
                    }
                    v++;
                }
                v++;
            }
            
            graphs_to_graph = std::string[canvases_to_graph + 1][graphs_per_canvas + 1][hists_per_graph +];
            int i = 0;
            int j = 0;
            int k = 0;
            int index = 0;
            for(const char *v = value; *v != '}';) {
                graphs_per_canvas = -1;
                while (*v != '}' && *v != '{') {
                    v++;
                    index++;
                }
                for(; *v != '}';) {
                    hists_per_graph = -1;
                    while (*v != '}' && *v != '{') {
                        v++;
                        index++;
                    }
                    for(; *v != '}';) {
                        while (*v != '}' && *v != '"') {
                            v++;
                            index++;
                        }
                        if (*v != '}')
                            break;
                        int pos = index;
                        v++;
                        index++;
                        while (*v != '}' && *v != '"') {
                            v++;
                            index++;
                        }
                        if (*v != '}')
                            break;
                        int length = index - pos;
                        graphs_to_graph[i][j][k] = std::string(value, pos, length);
                        k++;
                        v++;
                        index++;
                    }
                    j++;
                    v++;
                    index++;
                }
                i++;
                v++;
                index++;
            }
            std::cout << "Number of canvases to graph: " << canvases_to_graph << "\nNumber of graphs per canvas: " << graphs_per_canvas << "\Number of histograms per canvas: " << hists_per_graph << std::endl << "Titles: {";
            for(int i = 0; i <= canvases_to_title; i++) {
                for(int j = 0; j <= titles_per_canvas; i++) {
                    for(int k = 0; j <= titles_per_canvas; i++) {
                        std::cout << graphs_to_graph[i][j][k] << ", ";
                    }
                    std::cout << "}, ";
                }
                std::cout << "}, ";
            }
            std::cout << "}\n";
        }
        
        else if (strcmp(key, "Graph_titles") == 0) {
            canvases_to_title = -1;
            for(const char *v = value; *v != '}';) {
                tiles_per_canvas = -1;
                while (*v != '}' && *v != '{') {
                    v++;
                }
                canvases_to_title++;
                for(; *v != '}';) {
                    while (*v != '}' && *v != '"') {
                        v++;
                    }
                    if (*v != '}')
                        break;
                    titles_per_canvas++;
                    while (*v != '}' && *v != '"') {
                        v++;
                    }
                    if (*v != '}')
                        break;
                    v++;
                }
                v++;
            }
            graph_titles = new std::string[canvases_to_title+1][tiles_per_canvas+1]
            int i = 0;
            int j = 0;
            int index = 0;
            for(const char *v = value; *v != '}';) {
                while (*v != '}' && *v != '{') {
                    v++;
                    index++;
                }
                if ((*v != '}')
                    break;

                    for(; *v != '}';) {
                        while (*v != '}' && *v != '"') {
                            v++;
                            index++;
                        }
                        if (*v != '}')
                            break;
                        int pos = index;
                        v++;
                        index++;
                        while (*v != '}' && *v != '"') {
                            v++;
                            index++;
                        }
                        if (*v != '}')
                            break;
                        int length = index - pos;
                        graph_titles[i][j] = std::string(value, pos, length);
                        j++;
                        v++;
                        index++;
                    }
                    i++;
                    v++;
                    index++;
            }
            std::cout << "Number of canvases to title: " << canvases_to_title << "\nNumber of titles per canvas: " << titles_per_canvas << std::endl << "Titles: {";
            for(int i = 0; i <= canvases_to_title; i++) {
                for(int j = 0; j <= titles_per_canvas; i++) {
                    std::cout << graph_titles[i][j] << ", ";
                }
                std::cout << "}, ";
            }
            std::cout << "}\n";
        }
        
        else if (strcmp(key, "Graph_Colors") == 0) {
            num_of_colors = -1;
            for (const char *v = value; *v != '}';) {
                while (*v != '}' && !isdigit(*v)) {
                    v++;
                }
                
                num_of_colors++;
                
                while (*v != '}' && isdigit(*v)) {
                    v++;
                }
            }
            
            colors = new int[num_of_colors + 1];
            int i = 0;
            
            for (const char *v = value; *v != '}';) {
                while (*v != '}' && !isdigit(*v)) {
                    v++;
                }
                
                colors[i] = atoi(v);
                i++;
                
                while (*v != '}' && isdigit(*v)) {
                    v++;
                }
            }
            std::cout << "Number of colors in graphs: " << num_of_colors << std::endl << "Colors: {";
            for (int i = 0; i <= num_of_colors; i++) {
                std::cout << colors[i] << ", ";
            }
            std::cout << "}\n";
        }
                    
        else if (strcmp(key, "Canvas_rows") == 0) {
            canvas_rows = atoi(value);
            std::cout << "Number of graphing rows in the canvas is: " << canvas_rows << std::endl;
        }
        
        else if (strcmp(key, "Canvas_columns") == 0) {
            canvas_cols = atoi(value);
            std::cout << "Number of graphing columns in the canvas is: " << canvas_cols << std::endl;
        }
                    
        // Warning message if the key contents are unrecognized
        else {
            std::cout << "WARNING: Unrecognized keyvariable " << key << std::endl;
        }
    }
    fclose(config);
}
