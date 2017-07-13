# Pion Data Folder
## Technical Specifications
- Different_models: the data of the comparisons of various data models
- Gaussian Model Cuts: the first tests of what effects data cuts have on the data, done before the switch from the Gaussian model to the Crystal Ball model. Includes no source code. Includes: lambda02 cut, asymmetry cut, angle cut, Ncells cut
- THnSparses 060717 Crystal Ball Cuts: the tests of what effects data cuts have on the data, done after the switch from the Gaussian Model to the Crystal Ball Model but before the switch from THnSparses_060717 to THnSparses_062017. In addition to the Gaussian Model Cuts's cuts, this includes a Matched Tracks cut. 
- THnSparses 062017 Crystal Ball Cuts: Essentially the same as THnSparses 060717 Crystal Ball Cuts, except that it uses THnSparses 062017 as its THnSparses and includes a few minor upgrades done after the switch to THnSparses 062017
- my_code: the main generator of the data produced by this software (my_code_old_comparison for Different_models)
- pt_combiner: a comparator which compares graohs created by my_code (pt_combiner_old_comparison for Different_models")
- my_code_old: a copy of the oldest version of my_code, and the only one that uses Ntuple.root
- THnSparses_062717.root output: just like THnSparses 062017 Crystal Ball Cuts, but with THnSparses_062717.root. Used as a source of data for parts of Photon_data.
- THnSparses_062817_output: just like THnSparses 062017 Crystal Ball Cuts, but with THnSparses_062817.root as its source of data. First to replace matched track cut with distance to charged particles cut. Also used to analyze the effects of various cuts. Is the first to output data to .tex files for conversion into tables
- THnSparses_070517_output: has the same job description as THnSparses_062817_output, but it uses THnSparses_070517.root as its THnSparses and has perfected the C++ side of outputting data to .tex files
- THnSparses_071217_output: has the same job description as THnSparses_070517_output, but it uses THnSparses_071217.root as its THnSparses. Also includes the exoticity and time cuts
- table_file<other name parts go here>.tex - LaTeX files that contain table content data from my_code to be used by files in Presentations in tables
