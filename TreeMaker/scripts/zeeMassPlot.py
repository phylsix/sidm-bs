# -*- coding: utf-8 -*-
#!/usr/bin/env python
import sys
sys.argv.append('-b')
import ROOT; ROOT.TCanvas
sys.argv.remove('-b')

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

def styleStatsBox(hlist):
    for ih, h in enumerate(hlist):
        ROOT.gPad.Update()
        st = h.FindObject("stats")
        st.SetX1NDC(0.76)
        st.SetY1NDC(0.8-0.1*ih)
        st.SetX2NDC(0.9)
        st.SetY2NDC(0.9-0.1*ih)
        st.SetTextSize(0.02)
        st.SetTextColor(h.GetFillColor())

def addAxisTitle(xtitle, ytitle):
    yt = ROOT.TPaveText(0.021,0.7,0.048,0.88,"NDC")
    ytt = yt.AddText(ytitle)
    yt.SetTextAlign(22)
    ytt.SetTextAngle(90)
    yt.SetFillColor(0)
    yt.SetTextSize(0.03)
    yt.SetTextFont(42)
    yt.SetBorderSize(0)

    xt = ROOT.TPaveText(0.76,0.031,0.9,0.058,"NDC")
    xt.AddText(xtitle)
    xt.SetTextAlign(22)
    xt.SetFillColor(0)
    xt.SetBorderSize(0)
    xt.SetTextSize(0.03)
    xt.SetTextFont(42)
    return (xt,yt)

def addHistTitle(htitle):
    ht = ROOT.TPaveText(0.1,0.901,0.6,0.95,"NDC")
    ht.AddText(htitle)
    ht.SetTextAlign(12)
    ht.SetTextFont(62)
    ht.SetFillColor(0)
    ht.SetTextSize(0.034)
    ht.SetBorderSize(0)
    return ht

def addOverflow(h):
    nx = h.GetNbinsX()+1
    x1 = h.GetBinLowEdge(1)
    bw = h.GetBinWidth(nx)
    x2 = h.GetBinLowEdge(nx)+bw
    r = ROOT.TH1F(h.GetName()[0:-1],'',nx,x1,x2)
    for i in range(1, nx+1):
        r.Fill(r.GetBinCenter(i), h.GetBinContent(i))
    r.Fill(x1-1, h.GetBinContent(0))
    r.SetEntries(h.GetEntries())
    r.SetOption('HIST')
    return r

inf = ROOT.TFile('../data/output.root')

c0 = ROOT.TCanvas('c0','',500,500)
c1 = ROOT.TCanvas('c1','',500,500)
c2 = ROOT.TCanvas('c2','',500,500)
c3 = ROOT.TCanvas('c3','',500,500)
c4 = ROOT.TCanvas('c4','',500,500)
c5 = ROOT.TCanvas('c5','',500,500)

################################################################
genZeeTree_ = inf.Get('treeMaker/genZeeTree')
################################################################
c0.cd()
genZ1PtHist_ = ROOT.TH1F('genZ1PtHist_','',60,0,300)
genZ2PtHist_ = ROOT.TH1F('genZ2PtHist_','',60,0,300)

for _ in genZeeTree_:
    genZ1PtHist_.Fill(genZeeTree_.Z1_Pt)
    genZ2PtHist_.Fill(genZeeTree_.Z2_Pt)

genZ1PtHist = addOverflow(genZ1PtHist_)
genZ2PtHist = addOverflow(genZ2PtHist_)
genZ1PtHist.SetFillColorAlpha(ROOT.kAzure+3, 0.5)
genZ2PtHist.SetFillColorAlpha(ROOT.kOrange+7, 0.5)
genZ2PtHist.Draw()
genZ1PtHist.Draw(genZ1PtHist.GetOption()+'sames')

styleStatsBox([genZ1PtHist, genZ2PtHist])
c0.Update()
xt_0, yt_0 = addAxisTitle("Pt [GeV]","No. of Entries")
xt_0.Draw()
yt_0.Draw()
c0.Update()
ht_0 = addHistTitle("reco::GenParticle Z pt distribution")
ht_0.Draw()
c0.Update()

################################################################
c1.cd()
genM1eeHist_ = ROOT.TH1F('genM1eeHist_','',40,80,102)
genM2eeHist_ = ROOT.TH1F('genM2eeHist_','',40,80,102)

for _ in genZeeTree_:
    genM1eeHist_.Fill(genZeeTree_.m1ee)
    genM2eeHist_.Fill(genZeeTree_.m2ee)

genM1eeHist = addOverflow(genM1eeHist_)
genM2eeHist = addOverflow(genM2eeHist_)
genM1eeHist.SetFillColorAlpha(ROOT.kAzure+3, 0.5)
genM2eeHist.SetFillColorAlpha(ROOT.kOrange+7, 0.5)
genM1eeHist.Draw()
genM2eeHist.Draw(genM1eeHist.GetOption()+'sames')

styleStatsBox([genM1eeHist, genM2eeHist])
c1.Update()
xt_1, yt_1 = addAxisTitle("M [GeV]","No. of Entries")
xt_1.Draw()
yt_1.Draw()
c1.Update()
ht_1 = addHistTitle("Invariant mass reconstructed from two gen e distribution")
ht_1.Draw()
c1.Update()

################################################################
c2.cd()
genE1PtHist_ = ROOT.TH1F('genE1PtHist_','',40,0,200)
genE2PtHist_ = ROOT.TH1F('genE2PtHist_','',40,0,200)
genE3PtHist_ = ROOT.TH1F('genE3PtHist_','',40,0,200)
genE4PtHist_ = ROOT.TH1F('genE4PtHist_','',40,0,200)

for _ in genZeeTree_:
    genE1PtHist_.Fill(genZeeTree_.e1_Pt)
    genE2PtHist_.Fill(genZeeTree_.e2_Pt)
    genE3PtHist_.Fill(genZeeTree_.e3_Pt)
    genE4PtHist_.Fill(genZeeTree_.e4_Pt)

genE2PtHist = addOverflow(genE2PtHist_)
genE1PtHist = addOverflow(genE1PtHist_)
genE3PtHist = addOverflow(genE3PtHist_)
genE4PtHist = addOverflow(genE4PtHist_)
genE1PtHist.SetFillColorAlpha(ROOT.kAzure+3,0.5)
genE2PtHist.SetFillColorAlpha(ROOT.kOrange+7,0.5)
genE3PtHist.SetFillColorAlpha(ROOT.kPink-2,0.5)
genE4PtHist.SetFillColorAlpha(ROOT.kGreen+3,0.5)
genE2PtHist.Draw()
genE1PtHist.Draw(genE1PtHist.GetOption()+'sames')
genE3PtHist.Draw(genE3PtHist.GetOption()+'sames')
genE4PtHist.Draw(genE4PtHist.GetOption()+'sames')

styleStatsBox([genE1PtHist, genE2PtHist, genE3PtHist, genE4PtHist])
c2.Update()
xt_0.Draw()
yt_0.Draw()
c1.Update()
ht_2 = addHistTitle("gen e1, e2, e3, e4 pt distribution")
ht_2.Draw()
c2.Update()

################################################################
patZeeTree_ = inf.Get('treeMaker/patZeeTree')
################################################################
c3.cd()
patZ1PtHist_ = ROOT.TH1F('patZ1PtHist_','',60,0,300)
patZ2PtHist_ = ROOT.TH1F('patZ2PtHist_','',60,0,300)

for _ in patZeeTree_:
    patZ1PtHist_.Fill(patZeeTree_.Z1_Pt)
    patZ2PtHist_.Fill(patZeeTree_.Z2_Pt)

patZ1PtHist = addOverflow(patZ1PtHist_)
patZ2PtHist = addOverflow(patZ2PtHist_)
patZ1PtHist.SetFillColorAlpha(ROOT.kAzure+3, 0.5)
patZ2PtHist.SetFillColorAlpha(ROOT.kOrange+7, 0.5)
patZ2PtHist.Draw()
patZ1PtHist.Draw(patZ1PtHist.GetOption()+'sames')

styleStatsBox([patZ1PtHist, patZ2PtHist])
c3.Update()
xt_0.Draw()
yt_0.Draw()
c3.Update()
ht_3 = addHistTitle("pat Z pt distribution")
ht_3.Draw()
c3.Update()

################################################################
c4.cd()
patM1eeHist_ = ROOT.TH1F('patM1eeHist_','',40,80,102)
patM2eeHist_ = ROOT.TH1F('patM2eeHist_','',40,80,102)

for _ in patZeeTree_:
    patM1eeHist_.Fill(patZeeTree_.m1ee)
    patM2eeHist_.Fill(patZeeTree_.m2ee)

patM1eeHist = addOverflow(patM1eeHist_)
patM2eeHist = addOverflow(patM2eeHist_)
patM1eeHist.SetFillColorAlpha(ROOT.kAzure+3, 0.5)
patM2eeHist.SetFillColorAlpha(ROOT.kOrange+7, 0.5)
patM2eeHist.Draw()
patM1eeHist.Draw(patM1eeHist.GetOption()+'sames')

styleStatsBox([patM1eeHist, patM2eeHist])
c4.Update()
xt_1.Draw()
yt_1.Draw()
c4.Update()
ht_4 = addHistTitle("Invariant mass reconstructed from two pat e distribution")
ht_4.Draw()
c4.Update()

################################################################
c5.cd()
patE1PtHist_ = ROOT.TH1F('patE1PtHist_','',40,0,200)
patE2PtHist_ = ROOT.TH1F('patE2PtHist_','',40,0,200)
patE3PtHist_ = ROOT.TH1F('patE3PtHist_','',40,0,200)
patE4PtHist_ = ROOT.TH1F('patE4PtHist_','',40,0,200)

for _ in patZeeTree_:
    patE1PtHist_.Fill(patZeeTree_.e1_Pt)
    patE2PtHist_.Fill(patZeeTree_.e2_Pt)
    patE3PtHist_.Fill(patZeeTree_.e3_Pt)
    patE4PtHist_.Fill(patZeeTree_.e4_Pt)

patE1PtHist = addOverflow(patE1PtHist_)
patE2PtHist = addOverflow(patE2PtHist_)
patE3PtHist = addOverflow(patE3PtHist_)
patE4PtHist = addOverflow(patE4PtHist_)

patE1PtHist.SetFillColorAlpha(ROOT.kAzure+3,0.5)
patE2PtHist.SetFillColorAlpha(ROOT.kOrange+7,0.5)
patE3PtHist.SetFillColorAlpha(ROOT.kPink-2,0.5)
patE4PtHist.SetFillColorAlpha(ROOT.kGreen+3,0.5)
patE4PtHist.Draw()
patE1PtHist.Draw(patE1PtHist.GetOption()+'sames')
patE2PtHist.Draw(patE2PtHist.GetOption()+'sames')
patE3PtHist.Draw(patE3PtHist.GetOption()+'sames')

styleStatsBox([patE1PtHist, patE2PtHist, patE3PtHist, patE4PtHist])
c5.Update()
xt_0.Draw()
yt_0.Draw()
c5.Update()
ht_5 = addHistTitle("pat e1, e2, e3, e4 pt distribution")
ht_5.Draw()
c5.Update()

################################################################


c0.Print('../plots/output.pdf(')
c1.Print('../plots/output.pdf')
c2.Print('../plots/output.pdf')
c3.Print('../plots/output.pdf')
c4.Print('../plots/output.pdf')
c5.Print('../plots/output.pdf)')
