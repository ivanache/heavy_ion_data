// Runs every possible permutation of my_code.C with one macro call
void run_all_permutations(){
    gROOT->ProcessLine(".x my_code.C(0)");
    gROOT->ProcessLine(".x my_code.C(1)");
    gROOT->ProcessLine(".x my_code.C(2)");
    gROOT->ProcessLine(".x my_code.C(3)");
    gROOT->ProcessLine(".x my_code.C(4)");
}
