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
c1 = ROOT.TCanvas('c1','',500,500)
c2 = ROOT.TCanvas('c2','',500,500)
c3 = ROOT.TCanvas('c3','',500,500)
c4 = ROOT.TCanvas('c4','',500,500)
c5 = ROOT.TCanvas('c5','',500,500)
c6 = ROOT.TCanvas('c6','',500,500)
c7 = ROOT.TCanvas('c7','',500,500)
c8 = ROOT.TCanvas('c8','',500,500)

darkPhotons = f.Get('mcValidation/darkPhoton_reco')
pscalars    = f.Get('mcValidation/pscalar_reco')

h_zp_m = ROOT.TH1F('h_zp_m_','',20,0.,0.2) # 0.1GeV
h_ps_m = ROOT.TH1F('h_ps_m_','',30,98.5,101.5) # 100GeV
h_ep_dEtadPhi = ROOT.TH2F('h_ep_dEtadPhi_','',100,0,0.1,100,0,0.1)
h_ep_dR = ROOT.TH1F('h_ep_dR_','',100,0,0.15)
h_ep_vxvy = ROOT.TH2F('h_ep_vxvy_','',260,-130,130,260,-130,130)
h_ep_vz = ROOT.TH1F('h_ep_vz_','',100,-314,314)
h_ep_r = ROOT.TH1F('h_ep_r_','',100,0,300) # range to be tuned
h_ps_dEtadPhi = ROOT.TH2F('h_ps_dEtadPhi_','',100,0,7,100,0,7)
h_ps_dR = ROOT.TH1F('h_ps_dR_','',100,0,7)


for _ in darkPhotons:
    h_zp_m.Fill(darkPhotons.invm)
    h_ep_dEtadPhi.Fill(darkPhotons.dEta_ep, darkPhotons.dPhi_ep)
    h_ep_dR.Fill(darkPhotons.dR_ep)
    h_ep_vxvy.Fill(darkPhotons.dv_x, darkPhotons.dv_y)
    h_ep_vz.Fill(darkPhotons.dv_z)
    h_ep_r.Fill(math.sqrt(
        darkPhotons.dv_x*darkPhotons.dv_x \
        + darkPhotons.dv_y*darkPhotons.dv_y \
        + darkPhotons.dv_z*darkPhotons.dv_z))

for _ in pscalars:
    h_ps_m.Fill(pscalars.invm)
    h_ps_dEtadPhi.Fill(pscalars.dEta, pscalars.dPhi)
    h_ps_dR.Fill(pscalars.dR)


hs_1D = [h_zp_m, h_ps_m, h_ep_dR, h_ep_vz, h_ep_r, h_ps_dR]
hs_2D = [h_ep_dEtadPhi, h_ep_vxvy, h_ps_dEtadPhi]

for h in hs_1D:
    h.SetFillColor(ROOT.kAzure+3)
    h.SetLineColor(ROOT.kAzure+3)
    h.SetFillStyle(3004)

rootTools.styleAxisLabel(hs_1D, 132, 0.02)
rootTools.styleAxisLabel(hs_2D, 132, 0.02)

for h in hs_2D:
    h.SetStats(0)


c0.cd()
h_zp_m.Draw()
xt_0, yt_0 = rootTools.addAxisTitle("M [GeV]", "Num.Of.Entries / 0.01GeV")
xt_0.Draw()
yt_0.Draw()
ht_0 = rootTools.addHistTitle("Inv Mass of Dark Photon")
ht_0.Draw()
rootTools.styleStatsBox([h_zp_m],'fill')
c0.Update()

c1.cd()
h_ps_m.Draw()
xt_1, yt_1 = rootTools.addAxisTitle("M [GeV]", "Num.Of.Entries / 0.1GeV")
xt_1.Draw()
yt_1.Draw()
ht_1 = rootTools.addHistTitle("Inv Mass of Pseudo-scalar")
ht_1.Draw()
rootTools.styleStatsBox([h_ps_m],'fill')
c1.Update()

c2.cd()
h_ep_dEtadPhi.Draw('colz')
xt_2, yt_2 = rootTools.addAxisTitle("dEta", "dPhi")
xt_2.Draw()
yt_2.Draw()
ht_2 = rootTools.addHistTitle("Electron Pair Separation in dEta-dPhi")
ht_2.Draw()
c2.Update()

c3.cd()
h_ep_dR.Draw()
xt_3, yt_3 = rootTools.addAxisTitle("dR", "Num.Of.Entries")
xt_3.Draw()
yt_3.Draw()
ht_3 = rootTools.addHistTitle("Electron Pair Separation in dR")
ht_3.Draw()
rootTools.styleStatsBox([h_ep_dR],'fill')
c3.Update()

c4.cd()
h_ep_vxvy.Draw('colz')
xt_4, yt_4 = rootTools.addAxisTitle("x", "y")
xt_4.Draw()
yt_4.Draw()
ht_4 = rootTools.addHistTitle("Electron Pair Vertex in x-y")
ht_4.Draw()
c4.Update()

c5.cd()
h_ep_vz.Draw()
xt_5, yt_5 = rootTools.addAxisTitle("z", "Num.Of.Entries")
xt_5.Draw()
yt_5.Draw()
ht_5 = rootTools.addHistTitle("Electron Pair Vertex in z")
ht_5.Draw()
rootTools.styleStatsBox([h_ep_vz],'fill')
c5.Update()

c6.cd()
h_ep_r.Draw()
xt_6, yt_6 = rootTools.addAxisTitle("r", "Num.Of.Entries")
xt_6.Draw()
yt_6.Draw()
ht_6 = rootTools.addHistTitle("Electron Pair Vertex in r")
ht_6.Draw()
rootTools.styleStatsBox([h_ep_r],'fill')
c6.Update()

c7.cd()
h_ps_dEtadPhi.Draw('colz')
xt_7, yt_7 = rootTools.addAxisTitle("dEta", "dPhi")
xt_7.Draw()
yt_7.Draw()
ht_7 = rootTools.addHistTitle("Dark Photons Separation in dEta-dPhi")
ht_7.Draw()
c7.Update()

c8.cd()
h_ps_dR.Draw()
xt_8, yt_8 = rootTools.addAxisTitle("dR", "Num.Of.Entries")
xt_8.Draw()
yt_8.Draw()
ht_8 = rootTools.addHistTitle("Dark Photons Separation in dR")
ht_8.Draw()
rootTools.styleStatsBox([h_ps_dR],'fill')
c8.Update()


c0.Print(PLUGINBASE+'/plots/mcValidation0822.pdf(')
c1.Print(PLUGINBASE+'/plots/mcValidation0822.pdf')
c2.Print(PLUGINBASE+'/plots/mcValidation0822.pdf')
c3.Print(PLUGINBASE+'/plots/mcValidation0822.pdf')
c4.Print(PLUGINBASE+'/plots/mcValidation0822.pdf')
c5.Print(PLUGINBASE+'/plots/mcValidation0822.pdf')
c6.Print(PLUGINBASE+'/plots/mcValidation0822.pdf')
c7.Print(PLUGINBASE+'/plots/mcValidation0822.pdf')
c8.Print(PLUGINBASE+'/plots/mcValidation0822.pdf)')

cmd = "cp {0}/plots/mcValidation0822.pdf /publicweb/w/wsi/public".format(PLUGINBASE)
os.system(cmd)
