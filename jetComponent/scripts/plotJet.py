#!/usr/bin/env python
import os
import ROOT
import rootTools

CMSSWBASE = os.getenv('CMSSW_BASE')
PLUGINBASE = os.path.join(CMSSWBASE,'src','sidm-bs','jetComponent')
#ROOT.gROOT.Macro(PLUGINBASE+"/scripts/rootlogon.C")
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)


f = ROOT.TFile(os.path.join(PLUGINBASE, 'data', 'output.root'))

c0 = ROOT.TCanvas('c0','',500,500)
c1 = ROOT.TCanvas('c1','',500,500)

Gen_slimmedGenJets = f.Get('jetComponent/Gen_slimmedGenJets')
Pat_slimmedJets    = f.Get('jetComponent/PAT_slimmedJets')

h_genPt_  = ROOT.TH1F('h_genPt_','',100,0,100)
h_genEta_ = ROOT.TH1F('h_genEta_','',100,-5,5)
for _ in Gen_slimmedGenJets:
    h_genPt_.Fill(Gen_slimmedGenJets.pt)
    h_genEta_.Fill(Gen_slimmedGenJets.eta)


h_patPt_  = ROOT.TH1F('h_patPt_','',100,0,100)
h_patEta_ = ROOT.TH1F('h_patEta_','',100,-5,5)
for _ in Pat_slimmedJets:
    h_patPt_.Fill(Pat_slimmedJets.pt)
    h_patEta_.Fill(Pat_slimmedJets.eta)


h_genPt = rootTools.addOverflow(h_genPt_)
h_patPt = rootTools.addOverflow(h_patPt_)
h_genPt.SetFillColor(ROOT.kAzure+3)
h_genPt.SetLineColor(ROOT.kAzure+3)
h_genPt.SetFillStyle(3004)
h_patPt.SetFillColor(ROOT.kOrange+7)
h_patPt.SetLineColor(ROOT.kOrange+7)
h_patPt.SetFillStyle(3005)
rootTools.styleAxisLabel([h_genPt, h_patPt], 132, 0.02)

h_genEta = addOverflow(h_genEta_)
h_patEta = addOverflow(h_patEta_)
#h_genEta.SetMarkerColor(ROOT.kAzure+3)
#h_genEta.SetMarkerStyle(27) # empty diamond
#h_patEta.SetMarkerColor(ROOT.kOrange+7)
#h_patEta.SetMarkerStyle(24) # empty circle
h_genEta.SetFillColor(ROOT.kAzure+3)
h_genEta.SetLineColor(ROOT.kAzure+3)
h_genEta.SetFillStyle(3004)
h_patEta.SetFillColor(ROOT.kOrange+7)
h_patEta.SetLineColor(ROOT.kOrange+7)
h_patEta.SetFillStyle(3005)
rootTools.styleAxisLabel([h_genEta, h_patEta], 132, 0.02)


c0.cd()
ROOT.gPad.SetLogy()
h_genPt.Draw()
h_patPt.Draw(h_patPt.GetOption()+'sames')
rootTools.styleStatsBox([h_genPt, h_patPt],'fill')
c0.Update()
xt_0, yt_0 = rootTools.addAxisTitle("Pt [GeV]", "No. Of Entries / 1GeV")
xt_0.Draw()
yt_0.Draw()
c0.Update()
ht_0 = rootTools.addHistTitle("pt of genJet and reconstructed jet")
ht_0.Draw()
c0.Update()
leg_0 = ROOT.TLegend(0.62, 0.8, 0.76, 0.9)
leg_0.AddEntry("h_genPt", "gen jets", 'f')
leg_0.AddEntry("h_patPt", "pat jets", 'f')
leg_0.SetTextSize(0.02)
leg_0.SetTextAlign(22)
leg_0.Draw()
c0.Update()

c1.cd()
#ROOT.gPad.SetLogy()
h_genEta.Draw(h_genEta.GetOption())
h_patEta.Draw(h_patEta.GetOption()+'sames')
rootTools.styleStatsBox([h_genEta, h_patEta],'fill')
c1.Update()
xt_1, yt_1 = rootTools.addAxisTitle("Eta", "No. Of Entries")
xt_1.Draw()
yt_1.Draw()
c1.Update()
ht_1 = rootTools.addHistTitle("Eta of genJet and reconstructed jet")
ht_1.Draw()
c1.Update()
leg_1 = ROOT.TLegend(0.1, 0.8, 0.24, 0.9,"", "NDC")
leg_1.AddEntry("h_genEta", "gen jets", 'f')
leg_1.AddEntry("h_patEta", "pat jets", 'f')
leg_1.SetTextSize(0.02)
leg_1.SetTextAlign(22)
leg_1.Draw()
c1.Update()


c0.Print(PLUGINBASE+'/plots/output.pdf(')
c1.Print(PLUGINBASE+'/plots/output.pdf)')

cmd = "cp {0}/plots/output.pdf /publicweb/w/wsi/public".format(PLUGINBASE)
os.system(cmd)
