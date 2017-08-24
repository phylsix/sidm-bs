#!/usr/bin/env python
import os
import math
import ROOT
import rootTools

CMSSWBASE = os.getenv('CMSSW_BASE')
PLUGINBASE = os.path.join(CMSSWBASE,'src','sidm-bs','jetComponent')
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPalette(ROOT.kRedBlue)


f = ROOT.TFile(os.path.join(PLUGINBASE, 'mydata', 'mcvalidation0822.root'))

c0 = ROOT.TCanvas('c0','',500,500)

event = f.Get('mcValidation/eventTree')

h_num_jetE = ROOT.TH2F('h_num_jetE','',10,0,10,15,0,15)

for _ in event:
    h_num_jetE.Fill(event.numberOfElectrons, event.numberOfPatJets)

rootTools.styleAxisLabel([h_num_jetE], 132, 0.02)
h_num_jetE.SetStats(0)

c0.cd()
h_num_jetE.Draw('colz')
xt_0, yt_0 = rootTools.addAxisTitle("Num.Of.Electrons", "Num.Of.Jets")
xt_0.Draw()
yt_0.Draw()
ht_0 = rootTools.addHistTitle("Number of Reconstructed Jets vs Electrons")
ht_0.Draw()
c0.Update()


c0.Print(PLUGINBASE+'/plots/numOfJetE0822.pdf')
cmd = "cp {0}/plots/numOfJetE0822.pdf /publicweb/w/wsi/public".format(PLUGINBASE)
os.system(cmd)
