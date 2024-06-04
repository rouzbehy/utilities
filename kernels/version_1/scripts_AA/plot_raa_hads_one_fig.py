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
import dictionaries as my_dicts
from COLORS import rate_set_colors as rate_colours
# my_rcParams={
#     "text.usetex": True,
#     "font.family": "Georgia",
#     "font.size": 25,
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


## read in the experimental RAA:
cms     = pd.read_csv('../../../exp_data/charged/CMS/Table5_1.csv',comment='#').rename(columns=my_dicts.colnames_raa_cms)
alice   = pd.read_csv('../../../exp_data/charged/ALICE/PbPb_RAA_cent_00_05_eta_0p8.csv',comment='#').rename(columns=my_dicts.colnames_raa_alice)
atlas_1 = pd.read_csv("../../../exp_data/charged/ATLAS/Table49.csv",comment='#').rename(columns=my_dicts.colnames_raa_atlas)
atlas_2 = pd.read_csv("../../../exp_data/charged/ATLAS/Table33.csv",comment='#').rename(columns=my_dicts.colnames_raa_atlas) 
atlas_1 = atlas_1[atlas_1['x'].between(8, 104)]
atlas_2 = atlas_2[atlas_2['x'].between(8, 104)]

rate_names = {1:'LO', 2:'NLO',3:'NP'}
fig, ax = plt.subplots(1,1,gridspec_kw={'top':0.985,'bottom':0.11,
                                      'left':0.08,'right':0.98,
                                      'hspace':0.11,'wspace':0.09},
                              sharex=True, figsize=(16,9), sharey=True)

pT_low_cut, pT_high_cut = 10, 100
cms = cms[cms['x'].between(pT_low_cut, pT_high_cut)]
alice = alice[alice['x'].between(pT_low_cut, pT_high_cut)]

#for ax in [ax1, ax2]:
util.plot_expr_data_on_axis(ax,cms,'*')
util.plot_expr_data_on_axis(ax,alice,'o')
util.plot_expr_data_on_axis(ax,atlas_1,'P')
util.plot_expr_data_on_axis(ax,atlas_2,'s')

pp_loc  = '../martini_results/pp/hadron_spectra.csv'
### read in the calculated pp charged hadron spec
pp_data = pd.read_csv(pp_loc, comment='#')
pp_data['pT'] = 0.5 * (pp_data['pTmin'] + pp_data['pTmax'])
pp_data['dpT'] = pp_data['pTmax'] - pp_data['pTmin']
pp_data = pp_data[pp_data['pT'].between(pT_low_cut,pT_high_cut)]

x = pp_data['pT']
for rate_set in rate_names:
        color = rate_colours[rate_names[rate_set]]
        AA = pd.read_csv(f"../martini_results/final_PbPb_2p76/rset_{rate_set}/cent_0_5/hadron_spectra.csv", comment='#')
        AA['pT'] = 0.5 * (AA['pTmin'] + AA['pTmax'])
        AA['dpT'] = AA['pTmax'] - AA['pTmin']
        AA = AA[AA['pT'].between(pT_low_cut,pT_high_cut)]

        y1, dy1 = AA['N'], AA['dN']
        raa1 = y1/pp_data['N'] 
        draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + pp_data['dN']*pp_data['dN']/(pp_data['N']*pp_data['N']))
        ax.plot(x, raa1, color=color) 
        ax.fill_between(x, raa1+draa1, raa1-draa1, color=color, alpha=0.2)


expt_hands = [Line2D([],[],color=css['black'], label='CMS $|\eta|<1.0$', marker='*', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ALICE $|\eta|<0.8$', marker='o', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ATLAS $|\eta|<1.0$', marker='P', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ATLAS $|\eta|<2.0$', marker='s', markersize=10, linestyle="None")]

theor_hands = [Line2D([],[],label=l,color=c) for (l,c) in rate_colours.items()]

artist = ax.legend(loc='lower right', handles=expt_hands, ncol=1, bbox_to_anchor=(0.95, 0.01))
ax.add_artist(artist)

ax.legend(loc='center left', handles=theor_hands)

ax.set_ylim(top=0.81,bottom=-0.01)
ax.set_xlabel(r'$p_T$ (GeV)')
ax.set_ylabel(r'$R^{\mathrm{h}^{\pm}}_{\mathrm{AA}}$')
ax.text(0.1, 0.9, r'Pb-Pb \@ $\sqrt{s}=2.76$ ATeV' +'\n'+r'$|\eta|<1$, $0$-$5\%$', transform=ax.transAxes) 
plt.show()
