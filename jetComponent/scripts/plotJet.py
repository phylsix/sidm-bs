#!/usr/bin/env python
import os
import ROOT

CMSSWBASE = os.getenv('CMSSW_BASE')
PLUGINBASE = os.path.join(CMSSWBASE,'src','sidm-bs','jetComponent')
#ROOT.gROOT.Macro(PLUGINBASE+"/scripts/rootlogon.C")
ROOT.gROOT.SetBatch()

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


f = ROOT.TFile(os.path.join(PLUGINBASE, 'data', 'output.root'))

c0 = ROOT.TCanvas('c0','',500,500)
c1 = ROOT.TCanvas('c1','',500,500)

Gen_slimmedGenJets = f.Get('jetComponent/Gen_slimmedGenJets')
Pat_slimmedJets = f.Get('jetComponent/PAT_slimmedJets')

h_genPt_ = ROOT.TH1F('h_genPt_','',100,0,200)
h_genEta_ = ROOT.TH1F('h_genEta_','',100,-5,5)
for _ in Gen_slimmedGenJets:
    h_genPt_.Fill(Gen_slimmedGenJets.pt)
    h_genEta_.Fill(Gen_slimmedGenJets.eta)


h_patPt_ = ROOT.TH1F('h_patPt_','',100,0,200)
h_patEta_ = ROOT.TH1F('h_patEta_','',100,-5,5)
for _ in Pat_slimmedJets:
    h_patPt_.Fill(Pat_slimmedJets.pt)
    h_patEta_.Fill(Pat_slimmedJets.eta)


h_genPt = addOverflow(h_genPt_)
h_patPt = addOverflow(h_patPt_)

h_genEta = addOverflow(h_genEta_)
h_patEta = addOverflow(h_patEta_)

c0.cd()
h_genPt.Draw()
h_patPt.Draw(h_patPt.GetOption()+'sames')
c0.Update()

c1.cd()
h_genEta.Draw()
h_patEta.Draw(h_patEta.GetOption()+'sames')
c1.Update()

c0.Print(PLUGINBASE+'/plots/output.pdf(')
c1.Print(PLUGINBASE+'/plots/output.pdf)')
