from ROOT import *

class Legend:
    def __init__(self,title):
        self.title = title
        self.latex = TLatex()
        self.latex.SetNDC(True)
        self.entry = []
        self.entryText = []
        self.entryOption = []

    def Add(self,graph,text,option="P"):
        marker = TMarker()
        line = TLine()
        if (option=="P"):
            marker.SetMarkerColor(graph.GetMarkerColor())
            marker.SetMarkerStyle(graph.GetMarkerStyle())
            marker.SetMarkerSize(graph.GetMarkerSize()+0.3)
            self.entry.append(marker)
            self.entryText.append(text)
            self.entryOption.append(option)
        elif (option=="L"):
            line.SetLineColor(graph.GetLineColor())
            line.SetLineWidth(graph.GetLineWidth()+1)
            self.entry.append(line)
            self.entryText.append(text)
            self.entryOption.append(option)
        else:
            print "Error: Unknown legend option"

    def Draw(self,x,y):
        #self.latex.DrawLatex(x,y,"#font[72]{#scale[0.9]{"+self.title+"}}")
        self.latex.DrawLatex(x,y,"#font[62]{#scale[0.9]{"+self.title+"}}")
        stepY = 0.05
        for i in range(len(self.entry)):
            if (self.entryOption[i] == "P"):
                m = self.entry[i]
                m.SetX(x+0.02)
                m.SetY(y-0.05-i*stepY)
                m.SetNDC(True)
                m.Draw()
                text = self.entryText[i]
            elif (self.entryOption[i] == "L"):
                l = self.entry[i]
                lineLength = 0.015
                l.DrawLineNDC(x+0.02-lineLength,y-0.05-i*stepY,x+0.02+lineLength,y-0.05-i*stepY)
                text = self.entryText[i]
            self.latex.DrawLatex(x+0.05,y-i*stepY-0.065,"#scale[0.7]{"+text+"}")

