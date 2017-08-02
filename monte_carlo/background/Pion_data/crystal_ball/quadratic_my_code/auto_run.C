// This program auto-runs every non-deprecated permutation of my_code
// Programmer: Ivan Chernyshev; Date: 7/24/2017

void auto_run() {
  gROOT->ProcessLine(".x set_atlas_style.C()");
  gROOT->ProcessLine(".x my_code.C()");
  gROOT->ProcessLine(".x my_code.C(\"TESTSENSITIVITY\")");
  gROOT->ProcessLine(".x my_code.C(\"DOUBLETESTSENSITIVITY\")");
}
