// This macro is my first attempt at pion-hadron correlations, inspired by the sample that Miguel sent me, Plotting.C in the Pion_hadron_correlations/Miguel_sample_code
//Author: Ivan Chernyhsev
//Date: 10/19/17

#include "TList.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TVirtualFitter.h"
#include <iostream>

#include "atlasstyle-00-03-05/AtlasStyle.h"
#include "atlasstyle-00-03-05/AtlasStyle.C"
#include "atlasstyle-00-03-05/AtlasUtils.h"
#include "atlasstyle-00-03-05/AtlasUtils.C"
#include "atlasstyle-00-03-05/AtlasLabels.h"
#include "atlasstyle-00-03-05/AtlasLabels.C"

//variables of hPion
const int axis_pion_Cen         = 0;
const int axis_pion_Zvtx        = 1;
const int axis_pionMass         = 2;
const int axis_pionPt           = 3;
const int axis_pionRapidity     = 4;
const int axis_pion_asymmetry   = 5;
const int axis_pion_Ph1_Pt      = 6;
const int axis_pion_Ph2_Pt      = 7;
const int axis_pionOpeningAngle = 8;
const int axis_pion_Ph1_lambda02= 9;
const int axis_pion_Ph2_lambda02= 10;
const int axis_pion_Ph1_dR      = 11;
const int axis_pion_Ph2_dR      = 12;
