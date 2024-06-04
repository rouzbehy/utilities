#!/usr/bin/env python3
import numpy as np
from pandas import read_csv, DataFrame
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d
from matplotlib import ticker
from COLORS import module_colors, rate_set_colors
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
linestyles = {'200':'dashed', '2p76':'solid'}

loc = '~/Documents/research/jetscape_project/v2/jetscape_data/photon_sp/'
template = 'photons_00-20_{s}_{m}_jetscape.csv'

data_jetscape = {s: {m: read_csv(loc + template.format(s=s,m=m)) for m in ['cujet','martini']} for s in ['200','2p76']}

loc = './photon_sp/'
template = 'photons_00-20_rset_{r}_{s}_martini_alone.csv'
data_alone = {s: read_csv(loc+template.format(r=1,s=s),comment='#') for s in ['200','2p76']}

fig, ax = plt.subplots(1, 1, figsize=(16,9))


for s in data_jetscape:
    lstyle = linestyles[s]
    vals = data_jetscape[s]
    for r in vals:
        color = module_colors['MATTER+'+r.upper()]
        result = vals[r]['jmed']/vals[r]['total']
        ax.plot(vals[r]['pT'], result, color=color, linestyle=lstyle)
for s in data_alone:
    lstyle = linestyles[s]
    vals   = data_alone[s]
    color  = module_colors['PYTHIA+MARTINI']
    result = vals['jmed']/vals['total']
    ax.plot(vals['pT'], result, color=color, linestyle=lstyle)

labels = [Line2D([],[],color=c, label=l) for l, c in module_colors.items() if l != 'MATTER']

systs = {'200': r'Au-Au, $\sqrt{s}=200$ AGeV, $0$-$20\%$, $|\eta|<0.35$',
         '2p76':r'Pb-Pb, $\sqrt{s}=2.76$ ATeV, $0$-$20\%$, $|\eta|<0.8$'}

for s, l in linestyles.items():
    labels.append(Line2D([],[],color='black',linestyle=l,label=systs[s]))
ax.legend(loc='best', handles=labels, fontsize=20)

ax.set_ylabel('Jet-Medium/Total', fontsize=25)
ax.set_xlabel(r'$p^{\gamma}_T$ (GeV)', fontsize=25)
plt.show()