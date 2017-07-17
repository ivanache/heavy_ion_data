// This program auto-runs every permutation of my_code, followed by pt_combiner
// Prog rammer: Ivan Chernyshev; Date: 7/6/2017

void auto_run() {
    gROOT->ProcessLine(".x set_atlas_style.C()");
    gROOT->ProcessLine(".x my_code.C(0, \"quadric\")");
    gROOT->ProcessLine(".x my_code.C(1, \"quadric\")");
    gROOT->ProcessLine(".x my_code.C(2, \"quadric\")");
    gROOT->ProcessLine(".x my_code.C(3, \"quadric\")");
    gROOT->ProcessLine(".x my_code.C(4, \"quadric\")");
    gROOT->ProcessLine(".x pt_combiner.C()");
}
