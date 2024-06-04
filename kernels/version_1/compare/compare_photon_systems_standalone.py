## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle
from scipy.interpolate import interp1d as interpolate
## my custom modules
import util
import dictionaries as my_dicts
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
import COLORS 
rate_colours = COLORS.rate_set_colors
rate_names = {1:'LO',2:'NLO',3:'NP'}

loc = './photon_sp/'
template = 'photons_00-20_rset_{r}_{s}_martini_alone.csv'

data = {s: {r:pd.read_csv(loc+template.format(r=r,s=s),comment='#') for r in rate_names} for s in ['2p76']}

fig , ax = plt.subplots(1,1,figsize=(16,9))

linestyles = {'200':'dashed', '2p76':'solid'}

for s in data:
    lstyle = linestyles[s]
    vals = data[s]
    for r in vals:
        color = rate_colours[rate_names[r]]
        result = vals[r]['jmed']/vals[r]['total']
        ax.plot(vals[r]['pT'], result, color=color, linestyle=lstyle)
labels = [Line2D([],[],color=c, label=l) for l, c in rate_colours.items()]

systs = {'200': r'Au-Au, $\sqrt{s}=200$ AGeV, $0$-$20\%$, $|\eta|<0.35$',
         '2p76':r'Pb-Pb, $\sqrt{s}=2.76$ ATeV, $0$-$20\%$, $|\eta|<0.8$',}
#for s, l in linestyles.items():
#    labels.append(Line2D([],[],color='black',linestyle=l,label=systs[s]))
ax.legend(loc='best', handles=labels)
ax.text(0.05,0.05,'Pb-Pb, $\sqrt{s}=2.76$ ATeV, $0$-$20\%$, $|\eta|<0.8$', transform=ax.transAxes)
ax.set_ylabel('Jet-Medium/Total', fontsize=25)
ax.set_xlabel(r'$p^{\gamma}_T$ (GeV)', fontsize=25)
plt.show()
