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

INPUTTAG = 'electronFinder_100MeV_10E-13GeV'

f = ROOT.TFile(os.path.join(PLUGINBASE, 'mydata', '{0}.root'.format(INPUTTAG)))

c0 = ROOT.TCanvas('c0','',500,500)
c1 = ROOT.TCanvas('c1','',500,500)
c2 = ROOT.TCanvas('c2','',500,500)

darkPhotons = f.Get('electronFinder/darkPhotonFromGenElectrons')

h_ep_vxvy = ROOT.TH2F('h_ep_vxvy_','',250,-500,500,250,-500,500)
h_ep_vz   = ROOT.TH1F('h_ep_vz_','',100,-500,500)
h_ep_r    = ROOT.TH1F('h_ep_r_','',100,0,500) # range to be tuned


for _ in darkPhotons:
    h_ep_vxvy.Fill(darkPhotons.dv_x, darkPhotons.dv_y)
    h_ep_vz.Fill(darkPhotons.dv_z)
    h_ep_r.Fill(math.sqrt(
        darkPhotons.dv_x*darkPhotons.dv_x \
        + darkPhotons.dv_y*darkPhotons.dv_y \
        + darkPhotons.dv_z*darkPhotons.dv_z))


hs_1D = [h_ep_vz, h_ep_r]
hs_2D = [h_ep_vxvy]

for h in hs_1D:
    h.SetFillColor(ROOT.kAzure+3)
    h.SetLineColor(ROOT.kAzure+3)
    h.SetFillStyle(3004)

rootTools.styleAxisLabel(hs_1D, 132, 0.02)
rootTools.styleAxisLabel(hs_2D, 132, 0.02)

for h in hs_2D:
    h.SetStats(0)



c0.cd()
h_ep_vxvy.Draw('colz')
xt_0, yt_0 = rootTools.addAxisTitle("x", "y")
xt_0.Draw()
yt_0.Draw()
ht_0 = rootTools.addHistTitle("Electron Pair Vertex in x-y")
ht_0.Draw()
c0.Update()

c1.cd()
h_ep_vz.Draw()
xt_1, yt_1 = rootTools.addAxisTitle("z", "Num.Of.Entries")
xt_1.Draw()
yt_1.Draw()
ht_1 = rootTools.addHistTitle("Electron Pair Vertex in z")
ht_1.Draw()
rootTools.styleStatsBox([h_ep_vz],'fill')
c1.Update()

c2.cd()
h_ep_r.Draw()
xt_2, yt_2 = rootTools.addAxisTitle("r", "Num.Of.Entries")
xt_2.Draw()
yt_2.Draw()
ht_2 = rootTools.addHistTitle("Electron Pair Vertex in r")
ht_2.Draw()
rootTools.styleStatsBox([h_ep_r],'fill')
c2.Update()



c0.Print(PLUGINBASE+'/plots/{0}.pdf('.format(INPUTTAG))
c1.Print(PLUGINBASE+'/plots/{0}.pdf'.format(INPUTTAG))
c2.Print(PLUGINBASE+'/plots/{0}.pdf)'.format(INPUTTAG))

cmd = "cp {0}/plots/{1}.pdf /publicweb/w/wsi/public".format(PLUGINBASE, INPUTTAG)
os.system(cmd)
