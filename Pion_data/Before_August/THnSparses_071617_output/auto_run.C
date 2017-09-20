// This program auto-runs every non-deprecated permutation of my_code, followed by the default pt_combiner setting (which should be equal to the total number of my_code permutation runs done
// Programmer: Ivan Chernyshev; Date: 7/6/2017

void auto_run() {
    gROOT->ProcessLine(".x set_atlas_style.C()");
    gROOT->ProcessLine(".x my_code.C(0)");
    gROOT->ProcessLine(".x my_code.C(1)");
    gROOT->ProcessLine(".x my_code.C(2)");
    gROOT->ProcessLine(".x my_code.C(3)");
    //gROOT->ProcessLine(".x my_code.C(4)");
    gROOT->ProcessLine(".x pt_combiner.C(4)");
}
