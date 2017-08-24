#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ROOT
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch()
matplotlib.style.use('bmh')

f = ROOT.TFile('/Users/WNSI/Code/Work/sidm-bs/mcvalidation0819.root')

#darkPhotons = f.Get('mcValidation/darkPhoton_reco')
pscalars    = f.Get('mcValidation/pscalar_reco')

a = np.array([ x.invm for x in pscalars ])
df = pd.DataFrame({'a':a}, columns=['a'])
plt.hist(a, bins='auto',label=None)

stats = pd.DataFrame( np.round([df.count(), df.mean(), df.std()],3),
                      index=['Entry','Mean','Std'])
plt.table(cellText=stats.values, loc='upper right', colWidths=[0.15],
          rowLabels=stats.index)

plt.title("Pseudo-scalar's invariant mass")
plt.xlabel('M [GeV]')
plt.ylabel('Num.Of.Entries')
plt.grid(True)

fig = plt.gcf()
plt.show()