#!/usr/bin/env python
import ROOT

def styleStatsBox(hs, c):
    for ih, h in enumerate(hs):
        ROOT.gPad.Update()
        st = h.FindObject("stats")
        st.SetX1NDC(0.76)
        st.SetY1NDC(0.8-0.1*ih)
        st.SetX2NDC(0.9)
        st.SetY2NDC(0.9-0.1*ih)
        st.SetTextSize(0.02)
        if c == 'fill':
            st.SetTextColor(h.GetFillColor())
        if c == 'marker':
            st.SetTextColor(h.GetMarkerColor())

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

def styleAxisLabel(hs, f, s):
    for h in hs:
        h.GetXaxis().SetLabelFont(f)
        h.GetXaxis().SetLabelSize(s)
        h.GetYaxis().SetLabelFont(f)
        h.GetYaxis().SetLabelSize(s)


def addHistTitle(htitle):
    ht = ROOT.TPaveText(0.1,0.901,0.6,0.95,"NDC")
    ht.AddText(htitle)
    ht.SetTextAlign(12)
    ht.SetTextFont(132)
    ht.SetFillColor(0)
    ht.SetTextSize(0.034)
    ht.SetBorderSize(1)
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
