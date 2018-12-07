#!/usr/bin/env python
from AtlasCommonUtils import *
import time


def my_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

import optparse,os
op = optparse.OptionParser(usage="%prog [opts] ")

op.add_option("","--outputDir",dest="outdir",
              default=".",
              help="Output Directory")

op.add_option("-i","--inputFileList",dest="inlist",
              default='input.txt',
              help="Input File list with xsections")

op.add_option("-n","--Normalized",dest="Normalized",
              default=False,
              help="Normalized")

op.add_option("-l","--RatioLimit",dest="RatioLimit",
              default=4,
              help="The limit for the ratio plot (in %)")

op.add_option("-t","--Tag",dest="Tag",
              default="",
              help="tag to add to plots")

op.add_option("-m","--legendname",dest="legendname",
              default="",
              help="legend name")

op.add_option("-x","--xaxislabel",dest="xaxislabel",
              default="",
              help="x axis label")

op.add_option("-y","--yaxislabel",dest="yaxislabel",
              default="",
              help="y axis label")

op.add_option("-a","--ymaximum",dest="ymaximum",
              default="",
              help="y axis maximum")

op.add_option("-z","--yminimum",dest="yminimum",
              default="",
              help="y axis minimum")

op.add_option("","--hists",type='string', dest = "hists",
              action='callback', callback=my_callback,
              help="comma separated list of histograms you would like to Compare, e.g. hists,hists_i,hists_o,hists_ss" )


(p_ops, p_args) = op.parse_args()


output_dir = p_ops.outdir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
else:
    print "Using output dir ",output_dir,", overwriting existing contents"


if not p_ops.inlist:
    raise RuntimeError("Must specify infile list")


## Array to store a list of names of graphs to be compared
hists = []
if not p_ops.hists:
    hists = ["nhits"]
else:
    hists = p_ops.hists


RatioLimit = float(p_ops.RatioLimit)
Tag = p_ops.Tag

print 'The ratio plot will be plotted at maximum' , RatioLimit
import ROOT 
from Legend import Legend
from ROOT import TLatex
SetAtlasStyle()


color = [1,2,4,7,9,8]
style = [20,21,24,26,27,32] 

def DivideGraphs(den,num): 
    x = ROOT.Double() 
    y1= ROOT.Double()
    y2= ROOT.Double()
    ratio = ROOT.TGraphErrors()
    for i in range(den.GetN()):
        num.GetPoint(i,x,y2)
        den.GetPoint(i,x,y1)
        if(y1!=0):
            #print "SETTING POINT ", y2/y1
            ratio.SetPoint(i,x,(y2/y1-1)*100)
            ratio.SetPointError(i,0,0)
    #num.Print()
    #den.Print()
    #ratio.Print()
    return ratio

def TransformTH1toTGraph(histo):
    graph = ROOT.TGraphErrors(histo)
    return graph
            
def NormalizeGraph(graph):

    x = ROOT.Double()
    y= ROOT.Double()
    integral = 0
    for i in range(graph.GetN()):
        graph.GetPoint(i, x, y)
        integral = integral + y
  
    if integral==0: return

    for i in range(graph.GetN()):
        graph.GetPoint(i, x, y)
        
        graph.SetPoint(i, x, y/integral)
        graph.SetPointError(i, graph.GetErrorX(i), graph.GetErrorY(i)/integral)

    return graph


def plotHists(hists,files, xtitle, ytitle, title, pdfname,ly=False):
    c = ROOT.TCanvas()
    start = 0
    label = Legend(p_ops.legendname)
    multi = ROOT.TMultiGraph()

    for i in range(len(hists)): #This loops over the files. nhists(f1,f2,f3) is the format of hists.
        try:
            temp = files[i].GetName()
        
            if(p_ops.Normalized): hists[i] = NormalizeGraph(hists[i])
            hists[i].SetLineColor(color[i])
            hists[i].SetMarkerColor(color[i])
            hists[i].SetMarkerStyle(style[i])
            hists[i].SetMarkerSize(1)
            hists[i].SetLineWidth(2)
            # print 'i ' , i , ' MINIMUM:  ' , hists[i].GetHistogram().GetMinimum() , ' MAXIMUM : ' , hists[i].GetHistogram().GetMaximum(),
            label.Add(hists[i],tag[i],"L")
            multi.Add(hists[i])
        except:
            print "Skipping over missing file or histogram"

    if (multi.GetListOfGraphs().GetSize() == 0):
        print "Multigraph is empty; aborting"
        return
    cluslabel = Legend("15 GeV < p_{T_{trigger}} < 30 GeV")
    jetlabel = Legend("p_{T}^{jet} > 10 GeV ; R=0.4")

    #multi.SetMaximum(15)
    #multi.SetMinimum(5)
    multi.Draw("ALP")
    ROOT.gStyle.SetErrorX(0.0001);
    multi.SetMaximum(1.02*multi.GetHistogram().GetMaximum())
    multi.SetTitle(title)
    multi.GetXaxis().SetTitle(xtitle)
    multi.GetYaxis().SetTitle(ytitle)
    #multi.GetXaxis().SetTitle('x_{obs}^{Pb}:= #frac{p_{T}^{#gamma}e^{-#eta^{#gamma}} + p_{T}^{jet}e^{-#eta^{jet}}}{2E_{Pb}}')
    multi.GetXaxis().SetTitleSize(0.04)
    multi.GetXaxis().SetTitleOffset(1.6)
    #multi.GetYaxis().SetTitle('#frac{1}{N_{trig}}#frac{dN}{dx_{obs}^{Pb}}')
    #multi.GetXaxis().SetTitle('#Delta #phi (rads)')
    multi.GetXaxis().SetTitle(p_ops.xaxislabel)
    multi.GetYaxis().SetTitle(p_ops.yaxislabel)
    multi.GetYaxis().SetRangeUser(float(p_ops.yminimum), float(p_ops.ymaximum))
#multi.GetYaxis().SetTitle('#frac{1}{N_{trig}}#frac{dN}{d#Delta #phi}}')
#multi.GetXaxis().SetRangeUser(0.005, 0.020)
    multi.GetXaxis().SetNdivisions(10)

    label.Draw(0.6,.96)
    cluslabel.Draw(0.2, 0.9)
    jetlabel.Draw(0.2, 0.80)

    if ly==True:
        ROOT.gPad.SetLogy(1)
    latex = TLatex()
    latex.SetNDC()

    
    c.SaveAs(pdfname)
    c.Clear()
    print "END OF PLOTTING FUNCTION " 

def get_hists(files,hists):
    allhists = []
    for h in hists:
        hlist = []
        for f in files:
            f.Print()
            f.cd()
            histo = f.Get(h)
            if not histo:
                continue
            #histo.Print()
            if 'TH1' in str(type(histo)): 
                print 'Found TH1, transforming it into TGraphErrors:'
                histo = TransformTH1toTGraph(histo)
                 
            hlist.append(histo)
        allhists.append(hlist)
    return allhists


files = []
tag = []
xsecs = []
legs = []
weight = []



########### START OF MAIN FUNCTION #########################################
print "Opening file file " 
inlist = open(p_ops.inlist)

for line in inlist:
    if line[0] == "#": ##skipping comments
       continue
    print line
    tag.append(line.split(",")[1])
    weight.append(1)
    line = line.split(",")[0] #the line after the ROOT file name is the label 
    print 'Opening file', line
    files.append(ROOT.TFile(line.strip("\n")))
    line = line.split("hists")[0].split("_")[-1]
    tagOutputPlots = line
    
print "Getting all hists"
c = ROOT.TCanvas()

allhists = get_hists(files, hists) #gets all histograms 
temp = get_hists(files,hists)
print "Waiting..."
for i in range(len(hists)):
    print "Waiting..."
    try:
        plotHists(temp[i], files, allhists[i][0].GetXaxis().GetTitle(),allhists[i][0].GetYaxis().GetTitle() ,  allhists[i][0].GetTitle(),  output_dir + "/"+hists[i]+Tag+'_Comparison.pdf',ly=False)
    except IndexError:
        print "Skipping Comparison step due to data not being there."

print 'there are a total of ', len(tag)

cRatio = ROOT.TCanvas()

for i in range(len(hists)):
    try:
        ref = allhists[i][0]
        ratio =  DivideGraphs(ref, allhists[i][1])
        ratio.Draw()
        cRatio.SaveAs(output_dir + '/'+hists[i]+Tag+'_Ratio.pdf')
    except IndexError:
        print "Skipping Ratio step due to data not being there."






