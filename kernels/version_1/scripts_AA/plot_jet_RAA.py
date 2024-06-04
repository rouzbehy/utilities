## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
#from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.interpolate import interp1d as interpolate
## my custom modules
import util
import dictionaries as my_hdicts
import jetDicts as my_jdicts
# my_rcParams={
#     "text.usetex": True,
#     "font.family": "Georgia",
#     "font.size": 30,
#     "lines.linewidth": 2,
#     "xtick.direction": "in",
#     "ytick.direction": "in",
#     "xtick.minor.visible": True,
#     "ytick.minor.visible": True,
#     "xtick.major.size" : 8,
#     "ytick.major.size" : 8,
#     "xtick.minor.size" : 4,
#     "ytick.minor.size" : 4,
#     "axes.spines.right": False,
#     "axes.spines.top" : False,
#     "legend.frameon":False
# }
plt.rcParams.update(util.my_rcParams)
plt.rcParams.update({'lines.linewidth':3})

pTlow  = 60
pThigh = 300

def get_spectrum(fname):
    tmp = pd.read_csv(fname, comment='#')
    tmp['pT'] = 0.5 * (tmp['pTmin'] + tmp['pTmax'])
    tmp['dpT'] = tmp['pTmax'] - tmp['pTmin']
    result = tmp[tmp['pT'].between(pTlow, pThigh)]
    return result

## read in the experimental RAA:
cms = {}
for R in ['0p2','0p3','0p4']:
    cms[R] = pd.read_csv(f'../../../exp_data/jets/CMS/Table_JETRAA_{R}.csv',comment='#').rename(columns=my_jdicts.colnames_CMS_RAA)
from COLORS import rate_set_colors as rate_colours

rate_names = {1:'LO', 2:'NLO',3:'NP'}

#### read in the calculated pp charged hadron spec
pp_loc  = '../martini_results/pp/'
AA_loc_hadronic = '../martini_results/final_PbPb_2p76/rset_{r}/cent_0_5/jet_spectra.csv'
AA_loc_partonic = '../martini_results/final_PbPb_2p76/rset_{r}/cent_0_5/partonic_jet_spectra.csv'
pp_dat_hadronic = get_spectrum(pp_loc+'jet_spectra.csv')
pp_dat_partonic = get_spectrum(pp_loc+'partonic_jet_spectra.csv')

AA_dat_hadronic = {r:get_spectrum(AA_loc_hadronic.format(r=r)) for r in rate_names}
AA_dat_partonic = {r:get_spectrum(AA_loc_partonic.format(r=r)) for r in rate_names}


fig, axes = plt.subplots(3,1,gridspec_kw={'top':0.985,'bottom':0.12,
                                      'left':0.08,'right':0.78,
                                      'hspace':0.1,'wspace':0.09},
                              sharex=True, figsize=(16,9), sharey=True)

axes = axes.flatten()
for ir, r in enumerate(cms):
    ax = axes[ir]
    util.plot_expr_data_on_axis(ax,cms[r],'*')
    rtxt = r'$R=$' + r.replace('p','.')
    ax.text(0.45, 0.1, rtxt, transform=ax.transAxes)

x = pp_dat_partonic['pT']
for rset in rate_names:
    color = rate_colours[rate_names[rset]]
    AA_had = AA_dat_hadronic[rset]
    AA_prt = AA_dat_partonic[rset]
    for ir, r in enumerate(cms):
        ax = axes[ir]
        Rtag = f"N{r}"
        DRtag = f"dN{r}"
        y1, dy1 = AA_had[Rtag], AA_had[DRtag]
        y2, dy2 = AA_prt[Rtag], AA_prt[DRtag]
        b1, db1 = pp_dat_hadronic[Rtag], pp_dat_hadronic[DRtag]
        b2, db2 = pp_dat_partonic[Rtag], pp_dat_partonic[DRtag]
        raa1 = (y1/b1)
        raa2 = (y2/b2)
        draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + (db1*db1)/(b1*b1))
        draa2 = raa2*np.sqrt(dy2*dy2/(y2*y2) + (db2*db2)/(b2*b2))
        ax.plot(x, raa1, color=color, linestyle='solid') ## Hadronic Jet RAA 
        ax.fill_between(x, raa1+draa1, raa1-draa1, color=color, alpha=0.2)

        # ax.plot(x, raa2, color=color, linestyle='dashed') ## Partonic Jet RAA
        # ax.fill_between(x, raa2+draa2, raa2-draa2, color=color, alpha=0.2)

expt_hands = []
hands = [Line2D([],[],label=l,color=c) for (l,c) in rate_colours.items()]
hands.append(Line2D([],[],color=css['black'], label='CMS, $|\eta^{\mathrm{jet}}|<2.0$', marker='*', markersize=10, linestyle="None"))
#'Partonic Jet','dashed',
for item in zip([ 'Hadronic Jet'],['solid' ]):
    hands.append(Line2D([],[],color=css['black'],label=item[0],linestyle=item[1]))
    
axes[0].legend(loc='upper left', handles=hands, bbox_to_anchor=(1.0, 0.9), prop={'size':25})


ax = axes[2]
ax.set_ylim(bottom=-0.01)
ax.set_xlabel(r'$p^{\mathrm{jet}}_T$ (GeV)', fontsize=30)
for ax in axes:
    ax.set_ylabel(r'$R^{\mathrm{jet}}_{\mathrm{AA}}$', fontsize=30)

axes[2].legend(loc='lower right', handles=expt_hands)
plt.show()
