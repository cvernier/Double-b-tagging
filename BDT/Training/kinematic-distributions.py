from ROOT import *
import time

files = []

# files.append('tmp_reweighted_grav/TMVA_reweighted_grav.root')
# files.append('tmp_unweighted_grav/TMVA_unweighted_grav.root')
files.append('validation_None_bulkGrav800-4000_forTraining.root')
files.append('validation_None_qcd_forTraining.root')





gStyle.SetGridColor(kGray)
gStyle.SetOptStat(kFALSE)
gStyle.SetPadTopMargin(0.07)
gStyle.SetPadBottomMargin(0.13)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadRightMargin(0.06)
gROOT.ForceStyle()

histos = []
cat = False


for j in range(0,16):
  hname = "histos_%d" % (j) # Each histogram must have a unique name
  htitle = ""
  histos.append( TH1F(hname, htitle, 800, 300., 2000.) )
          
    
for t in files:
  f = TFile(t,"read") 
  # mytree=f.Get('TrainTree')
  mytree=f.Get('tree') 
  nentry = mytree.GetEntries()


# #########for sig eff == 70%
  if(f.GetName().find("") != -1):  
      cutBDTCat8 = 0.3 #83
      cutBDTG =  0.6
      cutG = 0.9	

      
  print f.GetName()  
  print cutBDTG
  print cutBDTCat8
  
  for event in mytree:
    # ptEtaWeight = event.weight
    
    if(f.GetName().find("") != -1):
      # if event.classID==0:
      if(f.GetName().find("Grav") != -1):
        histos[0].Fill( event.ptPruned )
        histos[6].Fill( event.massPruned )
        if(event.BDTG>= cutBDTG):
          histos[1].Fill(event.ptPruned)
	  histos[7].Fill(event.massPruned)
	if(event.BDTG>= cutBDTCat8):
          histos[8].Fill(event.massPruned)
          histos[2].Fill(event.ptPruned)
	if(event.BDTG>= cutG):
          histos[12].Fill(event.massPruned)
          histos[13].Fill(event.ptPruned)
      # if event.classID==1:
      if(f.GetName().find("qcd") != -1):
        histos[3].Fill(event.ptPruned)
        histos[9].Fill( event.massPruned )
        if(event.BDTG>= cutBDTG):
          histos[4].Fill(event.ptPruned)
	  histos[10].Fill(event.massPruned)
        if(event.BDTG>= cutBDTCat8):
          histos[11].Fill(event.massPruned)
          histos[5].Fill(event.ptPruned)
	if(event.BDTG>= cutG):
          histos[14].Fill(event.massPruned)
          histos[15].Fill(event.ptPruned)


#histos[8].Draw()          
#histos[9].Draw("same")
histos[1].Rebin(40)
histos[0].Rebin(40)
histos[4].Rebin(40)
histos[3].Rebin(40)
histos[2].Rebin(40)
histos[5].Rebin(40)
histos[13].Rebin(40)
histos[15].Rebin(40)

EffVsPt_Rew_grav =  TGraphAsymmErrors() 
EffVsPt_Rew_grav.Divide(histos[1],histos[0],"cl=0.683 b(1,1) mode") 
EffVsPt_Rew_grav.SetTitle("") 
EffVsPt_Rew_grav.GetXaxis().SetTitle("p_{T} (GeV)") 
EffVsPt_Rew_grav.GetXaxis().SetLimits(300.,2000.) 
EffVsPt_Rew_grav.GetYaxis().SetTitle("Tagging Efficiency (H#rightarrowb#bar{b})") 
EffVsPt_Rew_grav.GetYaxis().SetTitleOffset(1.26) 
EffVsPt_Rew_grav.GetHistogram().SetMaximum(1.1)       
EffVsPt_Rew_grav.GetHistogram().SetMinimum(0.)   
EffVsPt_Rew_grav.SetMarkerStyle(20) 
EffVsPt_Rew_grav.SetMarkerColor(kBlue+1) 


EffVsPt_Rew_grav2 =  TGraphAsymmErrors()
EffVsPt_Rew_grav2.Divide(histos[2],histos[0],"cl=0.683 b(1,1) mode")
EffVsPt_Rew_grav2.SetTitle("")
EffVsPt_Rew_grav2.GetXaxis().SetTitle("p_{T} [GeV]")
EffVsPt_Rew_grav2.GetXaxis().SetLimits(300.,2500.)
EffVsPt_Rew_grav2.GetYaxis().SetTitle("Tagging Efficiency (H#rightarrowb#bar{b})")
EffVsPt_Rew_grav2.GetYaxis().SetTitleOffset(1.16)
EffVsPt_Rew_grav2.GetHistogram().SetMaximum(1.1)
EffVsPt_Rew_grav2.GetHistogram().SetMinimum(0.)
EffVsPt_Rew_grav2.SetMarkerStyle(20)
EffVsPt_Rew_grav2.SetMarkerColor(2)

EffVsPt_Rew_grav3 =  TGraphAsymmErrors()
EffVsPt_Rew_grav3.Divide(histos[13],histos[0],"cl=0.683 b(1,1) mode")
EffVsPt_Rew_grav3.SetTitle("")
EffVsPt_Rew_grav3.GetXaxis().SetTitle("p_{T} (GeV)")
EffVsPt_Rew_grav3.GetXaxis().SetLimits(300.,2000.)
EffVsPt_Rew_grav3.GetYaxis().SetTitle("Tagging Efficiency (H#rightarrowb#bar{b})")
EffVsPt_Rew_grav3.GetYaxis().SetTitleOffset(1.16)
EffVsPt_Rew_grav3.GetHistogram().SetMaximum(1.1)
EffVsPt_Rew_grav3.GetHistogram().SetMinimum(0.)
EffVsPt_Rew_grav3.SetMarkerStyle(20)
EffVsPt_Rew_grav3.SetMarkerColor(kGreen+2)

c = TCanvas("c", "",800,800)
c.cd()
EffVsPt_Rew_grav.Draw("AP")
EffVsPt_Rew_grav2.Draw("PSAME")
EffVsPt_Rew_grav3.Draw("PSAME")
#c.SetGridx()
#c.SetGridy()

legend = TLegend(.6,.85,.81,.75)
# legend = TLegend(.16,.17,.3,.25)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.AddEntry(EffVsPt_Rew_grav2,"Double-b Loose",'p')
legend.AddEntry(EffVsPt_Rew_grav,"Double-b Medium",'p')
legend.AddEntry(EffVsPt_Rew_grav3,"Double-b Tight",'p')

legend.Draw()

l1 = TLatex()
l1.SetTextAlign(13)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.04)

#l1.DrawLatex(0.14+0.03,0.85, "R(0.6 to 2. TeV)")

l1.SetTextAlign(12)
l1.SetTextSize(0.035)
l1.SetTextFont(42)
l1.DrawLatex(0.82,0.96, "(13 TeV)")

l1.SetTextAlign(12)
l1.SetTextSize(0.035)
l1.SetTextFont(61)
l1.DrawLatex(0.15,0.96, "CMS")
l1.SetTextSize(0.033)
l1.SetTextFont(52)
l1.DrawLatex(0.23,0.955, "Simulation Preliminary")

l1.SetTextFont(42)
l1.SetTextSize(0.035)
#l1.DrawLatex(0.17,0.8, "20 % mistag efficiency")
l1.DrawLatex(0.17,0.85, "AK8, 70 < m < 200 GeV")# , p_{T} > 300 GeV")

c.Print("sigeff_cat.pdf")
EffVsPt_Rew_qcd =  TGraphAsymmErrors() 
EffVsPt_Rew_qcd.Divide(histos[4],histos[3],"cl=0.683 b(1,1) mode") 
EffVsPt_Rew_qcd.SetTitle("") 
EffVsPt_Rew_qcd.GetXaxis().SetTitle("p_{T} (GeV)") 
EffVsPt_Rew_qcd.GetXaxis().SetLimits(300.,2000.) 
EffVsPt_Rew_qcd.GetYaxis().SetTitle("Mistagging Efficiency") 
EffVsPt_Rew_qcd.GetYaxis().SetTitleOffset(1.6) 
EffVsPt_Rew_qcd.GetHistogram().SetMaximum(0.2)       
EffVsPt_Rew_qcd.GetHistogram().SetMinimum(0.)   
EffVsPt_Rew_qcd.SetMarkerStyle(20) 
EffVsPt_Rew_qcd.SetMarkerColor(kBlue+1) 

EffVsPt_Rew_qcd2 =  TGraphAsymmErrors()
EffVsPt_Rew_qcd2.Divide(histos[5],histos[3],"cl=0.683 b(1,1) mode")
EffVsPt_Rew_qcd2.SetTitle("")
EffVsPt_Rew_qcd2.GetXaxis().SetTitle("p_{T} [GeV]")
EffVsPt_Rew_qcd2.GetXaxis().SetLimits(300.,1500.)
EffVsPt_Rew_qcd2.GetYaxis().SetTitle("Tagging Efficiency (udscg)")
EffVsPt_Rew_qcd2.GetYaxis().SetTitleOffset(1.16)
EffVsPt_Rew_qcd2.GetHistogram().SetMaximum(0.35)
EffVsPt_Rew_qcd2.GetHistogram().SetMinimum(0.)
EffVsPt_Rew_qcd2.SetMarkerStyle(20)
EffVsPt_Rew_qcd2.SetMarkerColor(2)

EffVsPt_Rew_qcd3 =  TGraphAsymmErrors()
EffVsPt_Rew_qcd3.Divide(histos[15],histos[3],"cl=0.683 b(1,1) mode")
EffVsPt_Rew_qcd3.SetTitle("")
EffVsPt_Rew_qcd3.GetXaxis().SetTitle("p_{T} [GeV]")
EffVsPt_Rew_qcd3.GetXaxis().SetLimits(300.,1500.)
EffVsPt_Rew_qcd3.GetYaxis().SetTitle("Tagging Efficiency (udscg)")
EffVsPt_Rew_qcd3.GetYaxis().SetTitleOffset(1.16)
EffVsPt_Rew_qcd3.GetHistogram().SetMaximum(0.15)
EffVsPt_Rew_qcd3.GetHistogram().SetMinimum(0.)
EffVsPt_Rew_qcd3.SetMarkerStyle(20)
EffVsPt_Rew_qcd3.SetMarkerColor(kGreen+2)


c2 = TCanvas("c2", "",800,800)
c2.cd()
EffVsPt_Rew_qcd.Draw("AP")
EffVsPt_Rew_qcd2.Draw("PSAME")
EffVsPt_Rew_qcd3.Draw("PSAME")
#c2.SetGridx()
#c2.SetGridy()
legend.Draw()
l1.SetTextAlign(13)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.04)

l1.SetTextAlign(13)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.04)

#l1.DrawLatex(0.47+0.03,0.55, "R(0.6 to 2. TeV)")

l1.SetTextAlign(12)
l1.SetTextSize(0.035)
l1.SetTextFont(42)
l1.DrawLatex(0.82,0.96, "(13 TeV)")

l1.SetTextAlign(12)
l1.SetTextSize(0.035)
l1.SetTextFont(61)
l1.DrawLatex(0.15,0.96, "CMS")
l1.SetTextSize(0.033)
l1.SetTextFont(52)
l1.DrawLatex(0.23,0.955, "Simulation Preliminary")

l1.SetTextFont(42)
l1.SetTextSize(0.035)
#l1.DrawLatex(0.17,0.8, "20 % mistag efficiency")
l1.DrawLatex(0.17,0.85, "AK8, 70 < m < 200 GeV")# , p_{T} > 300 GeV")

c2.Print("eff_nocat.pdf")



time.sleep(200)
f.Close()
