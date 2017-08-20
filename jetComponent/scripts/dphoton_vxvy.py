#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ROOT
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

ROOT.gROOT.SetBatch()
matplotlib.style.use('bmh')

f = ROOT.TFile('/Users/WNSI/Code/Work/sidm-bs/mcvalidation0819.root')

darkPhotons = f.Get('mcValidation/darkPhoton_reco')
#pscalars    = f.Get('mcValidation/pscalar_reco')

nullfmt = NullFormatter()  # no labels

# definitions for the axes
left, width = 0.1, 0.7
bottom, height = 0.1, 0.7
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx   = [left, bottom_h, width, 0.12]
rect_histy   = [left_h, bottom, 0.12, height]

# start with a rectangular Figure
plt.figure(1, figsize=(7.5, 7.5))

axScatter = plt.axes(rect_scatter)
axHistx   = plt.axes(rect_histx)
axHisty   = plt.axes(rect_histy)

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

vx = np.array([p.dv_x for p in darkPhotons])
vy = np.array([p.dv_y for p in darkPhotons])

# scatter plot
axScatter.scatter(vx, vy, alpha=0.8, edgecolors='w')
axScatter.set_xlabel('X')
axScatter.set_ylabel('Y')
#plt.grid(False)

# projection hists plot
axHistx.hist(vx, bins=50)
axHisty.hist(vy, bins=50, orientation='horizontal')

# stats box
df = pd.DataFrame({'x':vx, 'y':vy}, columns=['x', 'y'])
stats = pd.DataFrame( np.round([df.count(), df.mean(), df.std()],3),
                      index=['Entry','Mean','Std'])
axScatter.table(cellText=stats.values, loc='upper right',
                colWidths=[0.12]*len(df.columns),
                rowLabels=stats.index,
                colLabels=df.columns)


fig = plt.gcf()
fig.suptitle("Electron Pair Vertex Position", fontsize=14)
plt.show()