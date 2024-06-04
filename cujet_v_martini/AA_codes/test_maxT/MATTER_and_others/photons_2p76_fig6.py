#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.interpolate import InterpolatedUnivariateSpline
from dictionaries import multiplicity, channel_colors, channel_linestyles

#mpl.use('Qt5Agg')
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 25,
    "lines.linewidth": 2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size" : 8,
    "ytick.major.size" : 8,
    "xtick.minor.size" : 4,
    "ytick.minor.size" : 4,
    "axes.spines.right": False,
    "axes.spines.top" : False,
    "legend.frameon":False
}
plt.rcParams.update(my_rcParams)


inel_Xsec = 62.03948 
expdata = {}
pT_low, pT_high = 3, 20
fnames = {}
fnames[200] = '../../../jetscape_data/sqrt_s_2760/test_MaxT_200/PbPb_2760/PbPb2760_00-05_photon_spec_0.80.csv'
fnames[20] = '../../../jetscape_data/sqrt_s_2760/martini/PbPb_2760/PbPb2760_00-05_photon_spec_0.80.csv'
jetscape = {}
oversampling_factors = {20:1000000, 200:10000}
for maxT in fnames:
    jetscape[maxT] = {}
    fname = fnames[maxT]
    factor = oversampling_factors[maxT]
    tmp = pd.read_csv(fname,comment='#')
    xmin = tmp['ptmin']
    xmax = tmp['ptmax']
    x = 0.5*(xmin+xmax)
    dx = xmax - xmin
    tmp['conv']  *= 1.0/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
    tmp['dconv'] *= 1.0/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
    tmp['brem']  *= 1.0/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
    tmp['dbrem'] *= 1.0/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
    tmp['pT'] = x
    tmp['dpT'] = dx
    jetscape[maxT] = tmp

fig, ax = plt.subplots(1,figsize=(15,8),
                        gridspec_kw={
                                     #'height_ratios':(3,1),
                                     'top':0.945,'bottom':0.105,
                                     'left':0.12,'right':0.965,
                                     'hspace':0.05,'wspace':0.01},
                        sharey='row', sharex=True)
colors = {20:'#7cae7a', 200:'#d64045'}
linestyles = {'brem':'solid','conv':'dashed'}
for item in jetscape:
    print(f"Dealing with {item}")
    tmp = jetscape[item]
    color = colors[item]
    for ch in ['conv','brem']:
        y, dy = tmp[ch], tmp[f'd{ch}']
        ax.plot(tmp['pT'],y, color=color, linestyle=linestyles[ch])
        ax.fill_between(tmp['pT'], y-dy, y+dy, color=color, linestyle=linestyles[ch])

labels = [Line2D([],[],label=r'MATTER: $\tau_{\mathrm{max}}=$'+ f'{t}' + '$\mathrm{fm}/c$ + MARTINI', color=c) for t, c in colors.items()]
for item in linestyles:
    lstyle = linestyles[item]
    labels.append(Line2D([],[],label=f'{item}.'.capitalize(),color='black',linestyle=lstyle))

ax.legend(loc='best', handles=labels)
ax.set_ylabel(r'$E_{\gamma} \frac{d N^{\gamma}}{d^{3}p}$, (GeV${}^{-2}$)')
ax.set_xlabel(r'$p^{\gamma}_{T}$, (GeV)')
ax.set_yscale('log')
plt.show()