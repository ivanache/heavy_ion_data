// Runs every possible permutation of my_code.C with one macro call
#include "my_code.C"
void run_all_permutations(){
    gROOT->ProcessLine(".x my_code.C(0), 'quadric'");
    gROOT->ProcessLine(".x my_code.C(1), 'quadric'");
    gROOT->ProcessLine(".x my_code.C(2), 'quadric'");
    gROOT->ProcessLine(".x my_code.C(3), 'quadric'");
    gROOT->ProcessLine(".x my_code.C(4), 'quadric'");
}
