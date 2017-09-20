// This program auto-runs every permutation of my_code, followed by pt_combiner
// Prog rammer: Ivan Chernyshev; Date: 7/6/2017

void auto_run() {
    gROOT->ProcessLine(".x set_atlas_style.C()");
    gROOT->ProcessLine(".x my_code.C(0)");
    gROOT->ProcessLine(".x my_code.C(1)");
    gROOT->ProcessLine(".x my_code.C(2)");
    gROOT->ProcessLine(".x my_code.C(3)");
    gROOT->ProcessLine(".x my_code.C(4)");
    gROOT->ProcessLine(".x my_code.C(5)");
    gROOT->ProcessLine(".x my_code.C(6)");
    gROOT->ProcessLine(".x my_code.C(7)");
    gROOT->ProcessLine(".x pt_combiner.C()");
}
