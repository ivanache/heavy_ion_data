# Technical Specifications

-13e_clusv1: the folder containing correlation analysis of data from the data NTuple 13e_clusv1_small.root
-NTuples: the folder with raw data from several NTuples
-Energy_Response: the folder with the C++ file in charge of the Energy response data
-Early_SameEvent_Correlations: the folder with the C++ file run on PDSF to find the signal and background spectra there, and associated files, run back in January to create 13e clusv1 clusterizer same-event correlations
-Correlations_With_Mixing: Correlations that include mixed events. Right now these inclued 13d from the v2 clusterizer mixed with 13d, 13d with the v2 clusterizer mixed with minimum bias 13c
-SameEvent_Correlations: Same-event Correlations started in late April-early May of 2018. So far these include 17p_passv1_wSDD pp correlations, simulated same event correlations with the lambda cut, and a mixed sample of 13d and 13e.
-PhotonEfficiency: the folder with the C++ file used to find efficiency of cluster photon processes, and associated files
- Fernando_code: the correlation code from Fernando. It uses event mixing for background.
-Corr_config.yaml: the configuration file used to change the values of various constants in Correlations.cc
- All .root files here are output, and their names indicate the settings of the configuration file constants and/or the Ntuples used to find them


