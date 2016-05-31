from ROOT import *
import math

def getPerformanceCurve(fFileS1, fFileB1, pTmin, pTmax, fXMin=0, fXMax=0):
    #get files and histograms
    file_S1 = TFile(fFileS1)
    file_B1 = TFile(fFileB1)
    
    # h2_S_new = file_S1.Get(fPlotS1)
    # h2_B_new = file_B1.Get(fPlotB1)
    
    h2_S_new = TH1F("hBDTGDisc_S","",10000,-1.,1.)
    h2_B_new = TH1F("hBDTGDisc_B","",10000,-1.,1.)
    h2_S_subjet = TH1F("hSubjetDisc_S","",10000,-1.,1.)
    h2_B_subjet = TH1F("hSubjetDisc_B","",10000,-1.,1.)
    h2_S_fatjet = TH1F("hFatjetDisc_S","",10000,-1.,1.)
    h2_B_fatjet = TH1F("hFatjetDisc_B","",10000,-1.,1.)
    
    S_bin1 = TH1F("S_bin1","",10000,-5.,5.)
    S_bin2 = TH1F("S_bin2","",10000,-5.,5.)
    S_bin3 = TH1F("S_bin3","",10000,-5.,5.)
    B_bin1 = TH1F("B_bin1","",10000,-5.,5.)
    B_bin2 = TH1F("B_bin2","",10000,-5.,5.)
    B_bin3 = TH1F("B_bin3","",10000,-5.,5.)
    
    S_bin1_subjet = TH1F("S_bin1_subjet","",10000,-5.,5.)
    S_bin2_subjet = TH1F("S_bin2_subjet","",10000,-5.,5.)
    S_bin3_subjet = TH1F("S_bin3_subjet","",10000,-5.,5.)
    B_bin1_subjet = TH1F("B_bin1_subjet","",10000,-5.,5.)
    B_bin2_subjet = TH1F("B_bin2_subjet","",10000,-5.,5.)
    B_bin3_subjet = TH1F("B_bin3_subjet","",10000,-5.,5.)
    
    
    
    
    treeS=file_S1.Get('tree') 
    treeB=file_B1.Get('tree') 
    Sentry = treeS.GetEntries() 
    Bentry = treeB.GetEntries() 
    
   

    for event in treeS:
      subjetcsv =event.SubJet_csv
      if(event.SubJet_csv < -1):
         subjetcsv = -1
      elif(event.SubJet_csv >1):
         subjetcsv = -1 # if (abs(event.flavour)==5 and event.nbHadrons>1 and event.ptPruned > pTmin and event.ptPruned <= pTmax):
      fatjetcsv =event.FatJet_csv
      if(event.FatJet_csv < -1):
        fatjetcsv = -1
      elif(event.FatJet_csv >1):
         fatjetcsv = -1
      if (event.ptPruned > pTmin and event.ptPruned <= pTmax and abs(event.flavour)==5 and event.nbHadrons>1):
        h2_S_subjet.Fill( subjetcsv )
        h2_S_fatjet.Fill( fatjetcsv )
        h2_S_new.Fill( event.BDTG )
        if (event.ptPruned > 200. and event.ptPruned <= 350.):
          S_bin1.Fill( event.BDTG )
          #S_bin1_subjet.Fill( event.SubJet_csv )
        if (event.ptPruned > 350. and event.ptPruned <= 500.):
          S_bin2.Fill( event.BDTG )
         # S_bin2_subjet.Fill( event.SubJet_csv )
        if (event.ptPruned > 500. and event.ptPruned <= 2000.):
          S_bin3.Fill( event.BDTG )    
          #S_bin3_subjet.Fill( event.SubJet_csv )
    
    for event in treeB:
      subjetcsv =event.SubJet_csv
      if(event.SubJet_csv < -1):
         subjetcsv = -1
      elif(event.SubJet_csv >1):
         subjetcsv = -1 # if (abs(event.flavour)==5 and event.nbHadrons>1 and event.ptPruned > pTmin and event.ptPruned <= pTmax):
      fatjetcsv =event.FatJet_csv
      if(event.FatJet_csv < -1):
         fatjetcsv = -1
      elif(event.FatJet_csv >1):
        fatjetcsv = -1
      if (event.ptPruned > pTmin and event.ptPruned <= pTmax):
        h2_B_subjet.Fill( subjetcsv )
        h2_B_fatjet.Fill( fatjetcsv )
        h2_B_new.Fill( event.BDTG )
        if (event.ptPruned > 200. and event.ptPruned <= 350.):
          B_bin1.Fill( event.BDTG )
         # B_bin1_subjet.Fill( event.SubJet_csv )
        if (event.ptPruned > 350. and event.ptPruned <= 500.):
          B_bin2.Fill( event.BDTG )
        #  B_bin2_subjet.Fill( event.SubJet_csv )
        if (event.ptPruned > 500. and event.ptPruned <= 2000.):
          B_bin3.Fill( event.BDTG )
       #   B_bin3_subjet.Fill( event.SubJet_csv )
    
    
    
    
    mg = []    

    #total jet count for denominator of efficiency calculation
    denom_S_1 = float( h2_S_new.Integral(0,10000) )
    denom_B_1 = float( h2_B_new.Integral(0,10000) )
    #denomB_ = denom_B_1
    #denomS_ = denom_S_1 	
    g_new = TGraph(10000)

    for i in range(10000):
        num_S = float( h2_S_new.Integral(10000-i*1,10000) )
        num_B = float( h2_B_new.Integral(10000-i*1,10000) )
        g_new.SetPoint( i,(num_S/denom_S_1),(num_B/denom_B_1) )
        # g_new.SetPoint( i,(num_B/denom_B_1),(num_S/denom_S_1) )

    
    
    g_subjet = TGraph(10000)

    for i in range(10000):
        num_S = float( h2_S_subjet.Integral(10000-i*1,10000) )
        num_B = float( h2_B_subjet.Integral(10000-i*1,10000) )
        g_subjet.SetPoint( i,(num_S/denom_S_1),(num_B/denom_B_1) )
    g_fatjet = TGraph(10000)

    for i in range(10000):
        num_S = float( h2_S_fatjet.Integral(10000-i*1,10000) )
        num_B = float( h2_B_fatjet.Integral(10000-i*1,10000) )
        g_fatjet.SetPoint( i,(num_S/denom_S_1),(num_B/denom_B_1) )
        # g_fatjet.SetPoint( i,(num_B/denom_B_1),(num_S/denom_S_1) )

                
    mg.append(g_new)     
    #mg.append(g1)   
    #mg.append(g2)
    #mg.append(g3)         
    #mg.append(g_baseline)       
    mg.append(g_subjet) 
    mg.append(g_fatjet) 

    return mg

def formatGraph(graph, graphNum):
    #colors = [ kBlue+1, kAzure+1, kAzure+2, kBlack, kRed+2, kGreen+3, kCyan ]
    #style= [1,2,3,4,1,1,1]	
    colors = [ kBlue+1, kRed+2, kGreen+3, kCyan ]
    style= [1,1,1,1]
	
    graphColor = colors[graphNum % 7]
    lineStyle = style[graphNum % 7]
    graph.SetLineColor(graphColor)
    graph.SetLineStyle(lineStyle)
    graph.SetLineWidth(2)

def plotPerformanceCurves(graphs, ordering, fTitle, fXAxisTitle, fYAxisTitle, fExtraInfo, fOutputFile, fXmin, fXmax, fYmin, fYmax, fLogy=0):

    # gROOT.SetBatch(kTRUE)
    gStyle.SetGridColor(kGray)
    gStyle.SetOptStat(kFALSE)
    gStyle.SetPadTopMargin(0.07)
    gStyle.SetPadBottomMargin(0.13)
    gStyle.SetPadLeftMargin(0.14)
    gStyle.SetPadRightMargin(0.06)
    gROOT.ForceStyle()

    c = TCanvas("c", "",800,800)
    c.cd()
    bkg = TH2D("bkg","",100,fXmin,fXmax,100,fYmin,fYmax)
    bkg.GetXaxis().SetTitle(fXAxisTitle)
    bkg.GetYaxis().SetTitle(fYAxisTitle)
    bkg.SetTitleOffset(1.2,"X")
    bkg.SetTitleOffset(1.5,"Y")
    bkg.Draw()
    c.SetGridx()
    c.SetGridy()

    legend = TLegend(.16,.79,.36,.93)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.033)

    graphCounter = 0
    for g in range(0,len(ordering)):
        graph = graphs[g]
        legend.AddEntry(graph, ordering[g],"l")
        formatGraph(graph,graphCounter)
        graph.Draw("L")
        graphCounter += 1

    if (fLogy):
        c.SetLogy()
    legend.Draw()
    l1 = TLatex()
    l1.SetTextAlign(13)
    l1.SetTextFont(42)
    l1.SetNDC()
    l1.SetTextSize(0.04)
    l1.DrawLatex(0.14+0.03,0.25, fTitle)

    l1.SetTextAlign(12)
    l1.SetTextSize(0.035)
    l1.SetTextFont(42)
    l1.DrawLatex(0.82,0.96, "(13 TeV)")
    
    l1.SetTextAlign(12)
    l1.SetTextSize(0.035)
    l1.SetTextFont(61)
    l1.DrawLatex(0.15,0.96, "CMS")
    l1.SetTextSize(0.03)
    l1.SetTextFont(52)
    l1.DrawLatex(0.23,0.955, "Simulation Preliminary")

    l1.SetTextFont(42)
    l1.SetTextSize(0.03)
    l1.DrawLatex(0.17,0.75, "AK8")
    l1.DrawLatex(0.17,0.7, fExtraInfo)

    afOutputFile=fOutputFile+".png"
    c.SaveAs(afOutputFile)
    bfOutputFile=fOutputFile+".pdf"
    c.SaveAs(bfOutputFile)
    name = fOutputFile+".root"
    f = TFile(name,"UPDATE");
    c.Write()
    f.Close()


def makePlots():


   ordering = [] # vectors storing the order of legend entries
   mg = [] # maps to hold legend entries and TGraph*s
   mg = getPerformanceCurve("validation_None_bulkGrav800-4000_forTraining.root","validation_None_qcd_forTraining.root",500.,800.)
   
   # ordering.append("Inclusive QCD")
#    ordering.append("Baseline")
#    ordering.append("Subjet CSV")
   ordering.append("double-b-tag ")
   #ordering.append("double-b-tag, 200 GeV < p_{T} < 350 GeV")
   #ordering.append("double-b-tag, 350 GeV < p_{T} < 500 GeV")
   #ordering.append("double-b-tag, p_{T} > 500 GeV")
   ordering.append("Subjet CSVv2")
   ordering.append("Fatjet CSVv2")

   plotPerformanceCurves(mg,ordering,"","Tagging efficiency (H#rightarrowb#bar{b})","Mistagging efficiency ","70 < m < 200 GeV , 500 < p_{T} < 800 GeV","subjet_500to800_all",0, 1, 1E-3, 1, 1)
   
   


if __name__ == "__main__":
    makePlots()
