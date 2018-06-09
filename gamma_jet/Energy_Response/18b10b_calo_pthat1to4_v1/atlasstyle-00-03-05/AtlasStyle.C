//
// ATLAS Style, based on a style file from BaBar
//

#include <iostream>

#include "AtlasStyle.h"

#include "TROOT.h"

void SetAtlasStyle ()
{
  static TStyle* atlasStyle = 0;
  std::cout << "\nApplying ATLAS style settings...\n" << std::endl ;
  if ( atlasStyle==0 ) atlasStyle = AtlasStyle();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
}

TStyle* AtlasStyle() 
{
  TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetFrameFillColor(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetStatColor(icol);
  //atlasStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  atlasStyle->SetPaperSize(20,26);

  // set margin sizes
  atlasStyle->SetPadTopMargin(0.10);
  atlasStyle->SetPadRightMargin(0.10);
  atlasStyle->SetPadBottomMargin(0.150);
  atlasStyle->SetPadLeftMargin(0.150);

  // set title offsets (for axis label)
  atlasStyle->SetTitleXOffset(1.50);
  atlasStyle->SetTitleYOffset(1.50);
  atlasStyle->SetLabelOffset(0.018,"X");
  atlasStyle->SetLabelOffset(0.018,"Y");
  atlasStyle->SetTitleAlign();
  //atlasStyle->SetLabelOffset(10.4);
  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica42
  Double_t tsize=0.04;
  atlasStyle->SetTextFont(font);

  atlasStyle->SetTextSize(tsize);
  atlasStyle->SetLabelFont(font,"x");
  atlasStyle->SetTitleFont(font,"x");
  atlasStyle->SetLabelFont(font,"y");
  atlasStyle->SetTitleFont(font,"y");
  atlasStyle->SetLabelFont(font,"z");
  atlasStyle->SetTitleFont(font,"z");
  
  atlasStyle->SetLabelSize(tsize,"x");
  atlasStyle->SetTitleSize(tsize,"x");
  atlasStyle->SetLabelSize(tsize,"y");
  atlasStyle->SetTitleSize(tsize,"y");
  atlasStyle->SetLabelSize(tsize,"z");
  atlasStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(1);
  atlasStyle->SetHistLineWidth(2);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  atlasStyle->SetErrorX(0.001);
  // get rid of error bar caps
 

  // atlasStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  atlasStyle->SetOptTitle(0);
  //atlasStyle->SetOptStat(1111);
  atlasStyle->SetOptStat(0);
  //atlasStyle->SetOptFit(1111);
  atlasStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(0);
  atlasStyle->SetPadTickY(0);
  //atlasStyle->SetTickLength(0.00,"x");
  //atlasStyle->SetTickLength(0.00,"y");
  //atlasStyle->SetTickLength(0.00,"z");
  atlasStyle->SetPalette(kViridis); //I want 112
  atlasStyle->SetOptStat("");
  atlasStyle->SetTitleAlign(4);

  //atlasStyle->SetNdivisions(220,"x");
  return atlasStyle;

}

