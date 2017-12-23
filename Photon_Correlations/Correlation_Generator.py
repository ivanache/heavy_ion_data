import ROOT 
import numpy

import sys
sys.path.append('/mnt/c/Users/marratia/Linux/PionHadronCorr/CompareHistograms')
from ROOT import TLatex

c = ROOT.TCanvas("c","c")

# Get the isolation signal and nonisolation background from fout.root, the output of the selection and dPhi spectrum generating file Correlations.cc, and generate files that contain the graphs for the components of
# The signal is to be stored in the file that will eventually have the correlation
fRaw = ROOT.TFile("fout.root")
fRaw.Print()

corr_ztmin0_ztmax1   = fRaw.Get("dPhi_iso_ztmin0_ztmax1")
corr_ztmin0_ztmax1.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin0_ztmax1.Draw()
c.SaveAs("Deep_Photon_Isolation_ztmin_00_ztmax_01.pdf")
c.Clear()

corr_ztmin1_ztmax2   = fRaw.Get("dPhi_iso_ztmin1_ztmax2")
corr_ztmin1_ztmax2.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin1_ztmax2.Draw()
c.SaveAs("Deep_Photon_Isolation_ztmin_01_ztmax_02.pdf")
c.Clear()

corr_ztmin2_ztmax4   = fRaw.Get("dPhi_iso_ztmin2_ztmax4")
corr_ztmin2_ztmax4.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin2_ztmax4.Draw()
c.SaveAs("Deep_Photon_Isolation_ztmin_02_ztmax_04.pdf")
c.Clear()

corr_ztmin4_ztmax6   = fRaw.Get("dPhi_iso_ztmin4_ztmax6")
corr_ztmin4_ztmax6.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin4_ztmax6.Draw()
c.SaveAs("Deep_Photon_Isolation_ztmin_04_ztmax_06.pdf")
c.Clear()

corr_ztmin6_ztmax8   = fRaw.Get("dPhi_iso_ztmin6_ztmax8")
corr_ztmin6_ztmax8.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin6_ztmax8.Draw()
c.SaveAs("Deep_Photon_Isolation_ztmin_06_ztmax_08.pdf")
c.Clear()

corr_ztmin8_ztmax10  = fRaw.Get("dPhi_iso_ztmin8_ztmax10")
corr_ztmin8_ztmax10.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin8_ztmax10.Draw()
c.SaveAs("Deep_Photon_Isolation_ztmin_08_ztmax_10.pdf")
c.Clear()

corr_ztmin10_ztmax12 = fRaw.Get("dPhi_iso_ztmin10_ztmax12")
corr_ztmin10_ztmax12.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin10_ztmax12.Draw()
c.SaveAs("Deep_Photon_Isolation_ztmin_10_ztmax_12.pdf")
c.Clear()


dPhi_noniso_ztmin0_ztmax1   = fRaw.Get("dPhi_noniso_ztmin0_ztmax1")
dPhi_noniso_ztmin0_ztmax1.SetTitle("; \Delta \phi / \pi; entries")
dPhi_noniso_ztmin0_ztmax1.Draw()
c.SaveAs("Deep_Photon_Nonisolation_ztmin_00_ztmax_01.pdf")
c.Clear()

dPhi_noniso_ztmin1_ztmax2   = fRaw.Get("dPhi_noniso_ztmin1_ztmax2")
dPhi_noniso_ztmin1_ztmax2.SetTitle("; \Delta \phi / \pi; entries")
dPhi_noniso_ztmin1_ztmax2.Draw()
c.SaveAs("Deep_Photon_Nonisolation_ztmin_01_ztmax_02.pdf")
c.Clear()

dPhi_noniso_ztmin2_ztmax4   = fRaw.Get("dPhi_noniso_ztmin2_ztmax4")
dPhi_noniso_ztmin2_ztmax4.SetTitle("; \Delta \phi / \pi; entries")
dPhi_noniso_ztmin2_ztmax4.Draw()
c.SaveAs("Deep_Photon_Nonisolation_ztmin_02_ztmax_04.pdf")
c.Clear()

dPhi_noniso_ztmin4_ztmax6   = fRaw.Get("dPhi_noniso_ztmin4_ztmax6")
dPhi_noniso_ztmin4_ztmax6.SetTitle("; \Delta \phi / \pi; entries")
dPhi_noniso_ztmin4_ztmax6.Draw()
c.SaveAs("Deep_Photon_Nonisolation_ztmin_04_ztmax_06.pdf")
c.Clear()

dPhi_noniso_ztmin6_ztmax8   = fRaw.Get("dPhi_noniso_ztmin6_ztmax8")
dPhi_noniso_ztmin6_ztmax8.SetTitle("; \Delta \phi / \pi; entries")
dPhi_noniso_ztmin6_ztmax8.Draw()
c.SaveAs("Deep_Photon_Nonisolation_ztmin_06_ztmax_08.pdf")
c.Clear()

dPhi_noniso_ztmin8_ztmax10  = fRaw.Get("dPhi_noniso_ztmin8_ztmax10")
dPhi_noniso_ztmin8_ztmax10.SetTitle("; \Delta \phi / \pi; entries")
dPhi_noniso_ztmin8_ztmax10.Draw()
c.SaveAs("Deep_Photon_Nonisolation_ztmin_08_ztmax_10.pdf")
c.Clear()

dPhi_noniso_ztmin10_ztmax12 = fRaw.Get("dPhi_noniso_ztmin10_ztmax12")
dPhi_noniso_ztmin10_ztmax12.SetTitle("; \Delta \phi / \pi; entries")
dPhi_noniso_ztmin10_ztmax12.Draw()
c.SaveAs("Deep_Photon_Nonisolation_ztmin_10_ztmax_12.pdf")
c.Clear()


#One by one, generate the correlation functions, save them in pT bin-labeled files
corr_ztmin0_ztmax1.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin0_ztmax1.Divide(dPhi_noniso_ztmin0_ztmax1)
corr_ztmin0_ztmax1.Draw()
c.SaveAs("Deep_Photon_Correlation_ztmin_00_ztmax_01.pdf")

c.Clear()
corr_ztmin1_ztmax2.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin1_ztmax2.Divide(dPhi_noniso_ztmin1_ztmax2)
corr_ztmin1_ztmax2.Draw()
c.SaveAs("Deep_Photon_Correlation_ztmin_01_ztmax_02.pdf")

c.Clear()
corr_ztmin2_ztmax4.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin2_ztmax4.Divide(dPhi_noniso_ztmin2_ztmax4)
corr_ztmin2_ztmax4.Draw()
c.SaveAs("Deep_Photon_Correlation_ztmin_02_ztmax_04.pdf")

c.Clear()
corr_ztmin4_ztmax6.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin4_ztmax6.Divide(dPhi_noniso_ztmin4_ztmax6)
corr_ztmin4_ztmax6.Draw()
c.SaveAs("Deep_Photon_Correlation_ztmin_04_ztmax_06.pdf")

c.Clear()
corr_ztmin6_ztmax8.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin6_ztmax8.Divide(dPhi_noniso_ztmin6_ztmax8)
corr_ztmin6_ztmax8.Draw()
c.SaveAs("Deep_Photon_Correlation_ztmin_06_ztmax_08.pdf")

c.Clear()
corr_ztmin8_ztmax10.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin8_ztmax10.Divide(dPhi_noniso_ztmin8_ztmax10)
corr_ztmin8_ztmax10.Draw()
c.SaveAs("Deep_Photon_Correlation_ztmin_08_ztmax_10.pdf")

c.Clear()
corr_ztmin10_ztmax12.SetTitle("; \Delta \phi / \pi; entries")
corr_ztmin10_ztmax12.Divide(dPhi_noniso_ztmin10_ztmax12)
corr_ztmin10_ztmax12.Draw()
c.SaveAs("Deep_Photon_Correlation_ztmin_10_ztmax_12.pdf")
