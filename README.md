# heavy_ion_data
Heavy Ion Data Version 1.1
Creator: Ivan Chernyshev, June 20 2017

## General Info
- Heavy Ion Data is a basic data analysis tool designed to analyze collision data from heavy ion collisions to make it easier to draw conclusions from the data
- Depends heavily on CERN's data analysis framework ROOT; data is in the form of THnSparse ROOT files
- Uses LaTeX to generate report files and to put data on an easy-to-read surface
- Outputs data as graphs in ATLAS style
- Currently, the program is incomplete, so errors may occur and all users are assumed to understand this fact 

## Techical Specifications
- gamma_jet: Everything involving gamma-jet study
- jet_jet: everything involving jet-jet study
- Photon_Correlations: The segment of the data analysis that uses Ntuples to construct photon-hadron correlations
- Ntuple_Analysis: The segment of the data analysis responsible for deep-pion efficiency analysis, which compares Dr. Yue Shi Lai's neural net to the existing lambda-based method of discerning signal from background
- Pion_Hadron Correlations: The segment of the data analysis that deals with pion-hadron correlations
- Pion Data: the folder with all of the data that is directly from pions
- Derived_Photon_data: the folder with all of the data that is about the photons with which the pions are measured and that at least partly comes out of the pion data analysis
- Direct_Photon_data: the folder with all of the data about photons that comes directly from the THnSparses that I start with and is not associated with pions
- Presentations: the folder with previously-presented presentations
- old: a storage compartment for data created by older versions of the code, as well as old versions of the code
- atlasstyle-00-03-05: the folder with the libraries for ATLAS's style. This is the master copy.
- THnSparses: the data format from where the data is drawn. Versions include (in order of release): Ntuple.root, THnSparses.root (not in use), THnSparses_060717.root, THnSparses_062017.root (first one to be used outside of Pion_data), THnSparses_062717.root, THnSparses_062817.root, THnSparses_070517 (first one to be used in Direct_Pion_data), THnSparses_071217.root, THnSparses_071417.root (only used in Direct_Photon_Data), THnSparses_071617.root (first to include more than one run, has Ncells>1, distance to border > 0, and distance to bad cell > 1 cuts built-in to its pion data), THnSparses_072417.root (An improvement over THnSparses_071617.root. Includes smaller bin widths and a fix of the double-counting problem in THnSparses_071617.root, as well as a built-in R < 0.4 and R < 0.2 cuts on some variables in hCluster), and THnSparses_LHC13d.root ( It relaxes some cuts, should include the entire data sample, and is the first one do direct photon analysis (specifically the distance to border and bad cell cuts) with a momentum cut of pT > 8 GeV), THnSparses_LHC13d_090717.root (an updated version of the previous one), THnSparses_LHC13e_090717.root THnSparses_LHC13f_090717.root (the previous 2 are just like THnSparses_LHC13d_090717.root, except they represent different runs), and THnSparses_LHC13d_101517.root (Used in Pion_Hadron_Correlations)
- monte_carlo: versions of the rest of the data analysis that uses a THnSparses file with Monte Carlo simulation data. These THnSparses include: THnSparses_080117_MC.root (which unfortunately is too big for GitHub to hold) and THnSparses_JetJetMC.root. Changes: Gaussian peak, quadric background is the default mass-pion data fit. Is divided into different folders for various types of simulations, such as background, signal, and signal+background
- combination: plots ALICE data superimposed on simulation plots. pion_combiner is the macro in use, but is currently being replaced by the Compare_Histograms data
- AtlasCommonUtils.py, CompareHistograms.py, Legend.py, and input.txt: Files needed to run the Compare_Histograms algorithm, which combines data histograms for comparison. input.txt specifies the file coordinates, CompareHistograms.py is the file to be called.

## Guide to working with heavy_ion_data
- To start, start up ROOT by entering "root" into Terminal after navigating to the desired directory (make sure that you have ROOT 6 installed and hooked up to Terminal). 
- To view the data inside a THnSparses, first enter the constructor TFile* <file object (your choice) name here> = new TFile(<insert the name of the THnSparses here (make sure you are in the THnSparses's directory)>);. Alternatively, you may enter the name of the THnSparses right after the "root" command when you launch ROOT. Then enter the command TBrowser* <insert the name of the browser object (your choice) here> = new TBrowser();. A  Browser window will pop up, and your THnSparses will be one of the files on the menu to the left. Click on it to access the data inside.
- To run a macro (one of the C++ programs which processes the data), type in .x <enter macro name here, followed by any arguments that the macro requires>.

**Flowchart of the data processing:**<br />
THnSparses<br />
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;||<br />
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                   ||<br />
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; \\/<br />
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Data processing macro<br />
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;||&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;||&nbsp;&nbsp;/\ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;||<br />
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ||&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;||&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;||<br />
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; \\/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\/&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\/<br />
.png graph ROOT output Terminal print and<br />
files &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; files &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; to-file output<br />
    ||&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;        ||<br />
    \\/ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; \\/ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; \\/<br />
Presentations and other communication<br />
