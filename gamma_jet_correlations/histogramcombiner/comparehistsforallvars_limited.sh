#!/bin/bash

./loopovervariables input_showershapecompare_pp.txt pp_data
./loopovervariables input_showershapecompare_pPb.txt pPb_data

./loopovervariables input_pPb-pp_comparison_DNN.txt DNN_pPb_pp_comparison
./loopovervariables input_pPb-pp_comparison_EmaxOverEcluster.txt Emax_over_Ecluster_pPb_pp_comparison
./loopovervariables input_pPb-pp_comparison_Lambda0.txt Lambda0_pPb_pp_comparison
