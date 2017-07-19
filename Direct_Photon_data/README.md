# Direct_Photon_data
## Technical specifications
- THnSparses_070517.root: the first direct source of the data and the only current source of data for Cut Comparisons
- THnSparses_071217.root: an old direct source of the data for Modeler
- THnSparses_071617.root: the current source of the data for Modeler
- auto_run.C: the macro that auto-runs all possible combinations of the other macros
- atlasstyle-00-03-05: the folder with all of the atlas style libraries
- old: files that are probably not going to be currently used

### Cut Comparisons, the part of the program that is responsible for comparing data with different distance to bad cuts and distance to botder cuts
- photon_grapher.C: the macro that graphs the data for both the photon spectra and the ratios between data of various distance to border and distance to bad cells cuts
- ratio_combiner: the LaTeX beamer that combines spectra and ratio plots from all data folders
- 0.1<lambda<0.4: the folder with data that has been cut to 0.1<lambda<0.4
- 0.4<lambda: the folder with data that has been cut to 0.4<lambda
- 4<pT<20, 10<pT<20: folders inside 0.1<lambda<0.4 and 0.4<lambda which contain data whose transverse momentum range is 4-20 GeV/c and 10-20 GeV/c, respectively

### Modeler, the part of the program that is responsible for examining the photon data in order to find out how many photons are in the sample, as well as generating basic plots
- photon_modeler.C: the macro that graphs and examines the pion data spectrum. Initially based on a specific iteration of photon_grapher.C from Cut Comparisons
- photon_Elambda.C: the macro that graphs a correlation of energy readings and lambda
