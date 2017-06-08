// Runs every possible permutation of my_code.C with one macro call
void run_all_permutations(){
    gROOT->ProcessLine(".x my_code.C(true, true, true)");
    gROOT->ProcessLine(".x my_code.C(true, true, false)");
    gROOT->ProcessLine(".x my_code.C(true, false, true)");
    gROOT->ProcessLine(".x my_code.C(true, false, false)");
    gROOT->ProcessLine(".x my_code.C(false, true, true)");
    gROOT->ProcessLine(".x my_code.C(false, true, false)");
    gROOT->ProcessLine(".x my_code.C(false, false, true)");
    gROOT->ProcessLine(".x my_code.C(false, false, false)");
}
