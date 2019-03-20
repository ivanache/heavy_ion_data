// This macro auto-loops through all possible permutations of photon_grapher
// Author: Ivan Chernyshev; Date: 7/11/17

void auto_run() {
    gROOT->ProcessLine(".x photon_grapher.C(\"0.1-0.4\", \"4-20\")");
    gROOT->ProcessLine(".x photon_grapher.C(\"0.1-0.4\", \"10-20\")");
    gROOT->ProcessLine(".x photon_grapher.C(\">0.4\", \"4-20\")");
    gROOT->ProcessLine(".x photon_grapher.C(\">0.4\", \"10-20\")");
}
