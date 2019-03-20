# General Information
- This repository contatins the code necessary to replicate Ivan Chernyshev's research results taken between June 2017 and February 2019
- gamma_jet_correlations contains the gamma-jet correlations which Ivan has been working on since July 2019
- general_tools contains several useful tools: Integrating certain histograms over their variables, plotting histograms over various variables (both ROOT and HDF5) into a pdf plot, conversion from ROOT to HDF5 files, taking the ratio of the data in one ROOT file to one in another root file, injecting a mixed event list into a ROOT file, setting the plot style of an output plot, and merging the outputs of 3 different ROOT files
- Root files contain TH1D, TH1F, TH2D, THNSparses, etc. type plots, that can be accessed from Terminal via:
&nbsp;&nbsp; root
&nbsp;&nbsp; TBrowser <insert a random name for the TBrowser object that will display the plots>
or via plotting software like what there is in our files
- HDF5 files are a certain data type that organizes data efficiently for puropses of event-mixing

