# Histogram Combiner
## This directory contains tools for taking histograms generated in the main gamma jet folder and combine and format them into the appropriate histograms
- CompareHistograms.py: takes in a set of ROOT files with histograms and creates a graph of superimposed plots with the same name, one from each ROOT file. The ROOT file names come from a test file (default name is input.txt, although alternative names can be inputted as arguments to the program). Name of the plots must be inputted as an argument. Arguments to CompareHistograms include the text file with the names of the ROOT files, the maxima and minima of the x- and y- axes, etc. In the terminal, arguments are inputted after the command python CompareHistograms.py, and come in the form of -<letter index> <input>. Letter indices for each purpose can be found at the top of CompareHistograms.py.
- loopovervariables.cc : given an input text file for CompareHistograms, calls CompareHistograms to create graphs for all of the following variables: \Delta \phi, x^{obs}_{Pb}, X_{j}, p_{T}D, and multiplicity, for shower shape signal (sig), shower shape background (bkg), and signal with background subtracted out (hadj)

- background_purity_subtraction: the macro which subtracts the shower shape background from the signal in order to reduce pion background
- autorun.C: scripts which repeatedly run background_purity_subtraction through all relevant histograms for gamma-jet correlation comparison (at least the ones I've been doing). If you want to know how to use background_purity_subtraction, go here
- comparehistsforallvars.sh: a shell script that runs  loopovervariables.cc for each relevant input text file, (i.e. for each data-set comparison group)
- purity_finder.cc: finds the purity of a given ROOT file's data, with the help of purity values from the analysis note
-truthseparator.C: separates the truth data sets from the rest of the data in order to enable the architecture of the HistogramCombiner programs

## Format of an input text file:
On each row: <ROOT file name (file that will be access)> <label for legend>

Example:
GammaJet_config_clusptmin15.0_clusptmax30.0_JETPTMIN_10.0_DATANAME__13f_PHOTONSELECT_Lambda0.root, 13f(Pb-p)
GammaJet_config_clusptmin15.0_clusptmax30.0_JETPTMIN_10.0_DATANAME__13e_PHOTONSELECT_Lambda0.root, 13e(p-Pb)
GammaJet_config_clusptmin15.0_clusptmax30.0_JETPTMIN_10.0_DATANAME__13d_PHOTONSELECT_Lambda0.root, 13d(p-Pb)
