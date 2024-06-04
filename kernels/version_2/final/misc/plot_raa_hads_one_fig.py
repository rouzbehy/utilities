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

plt.rcParams.update(util.my_rcParams)
pT_low_cut, pT_high_cut = 33, 100

## read in the experimental RAA:
exp_loc = '../../../exp_data/charged/'
cms_cents = {1:'0_5',2:'5_10',3:'10_30'}
cms     = {c:pd.read_csv(exp_loc+f'CMS/Table5_{i}.csv',comment='#').rename(columns=my_dicts.colnames_raa_cms)
           for i,c in cms_cents.items()}

alice_cents = {'00_05':'0_5','05_10':'5_10','10_20':'10_20'}
alice   = {c:pd.read_csv(exp_loc+f'ALICE/PbPb_RAA_cent_{i}_eta_0p8.csv',comment='#').rename(columns=my_dicts.colnames_raa_alice)
           for i,c in alice_cents.items()}

atlas_1 = pd.read_csv(exp_loc+"ATLAS/Table49.csv",comment='#').rename(columns=my_dicts.colnames_raa_atlas)
atlas_2_cents = {33:'0_5',34:'5_10',35:'10_20'}
atlas_2 = {c:pd.read_csv(exp_loc+f"ATLAS/Table{i}.csv",comment='#').rename(columns=my_dicts.colnames_raa_atlas)
           for i, c in atlas_2_cents.items()}

cms = {c: itm[itm['x'].between(pT_low_cut,pT_high_cut)] for c, itm in cms.items()}
alice = {c:itm[itm['x'].between(pT_low_cut,pT_high_cut)] for c, itm in alice.items()}

atlas_1 = atlas_1[atlas_1['x'].between(pT_low_cut,pT_high_cut)]
atlas_2 = {c:itm[itm['x'].between(pT_low_cut,pT_high_cut)] for c, itm in atlas_2.items()}

rate_names = {1:'LO', 2:'NLO',3:'NP'}
fig, axes = plt.subplots(1,3,gridspec_kw={'top':0.985,'bottom':0.11,
                                      'left':0.08,'right':0.98,
                                      'hspace':0.11,'wspace':0.02},
                              sharex=True, figsize=(16,9), sharey=True)



for cent, ax in zip(alice, axes):
    alice_spectrum = alice[cent]
    atlas_spectrum = atlas_2[cent]
    util.plot_expr_data_on_axis(ax, alice_spectrum, 'P', s=150)
    util.plot_expr_data_on_axis(ax, atlas_spectrum, 's', s=150)

for cent, ax in zip(cms, axes):
    cms_spectrum = cms[cent]
    util.plot_expr_data_on_axis(ax,cms_spectrum,'*', s=150)

util.plot_expr_data_on_axis(axes[0],atlas_1,'o')

pp_loc  = '../../calcs/pp_2p76/chgd_spec.csv'
### read in the calculated pp charged hadron spec
pp_data = pd.read_csv(pp_loc, comment='#')
pp_data['pT'] = 0.5 * (pp_data['pTmin'] + pp_data['pTmax'])
pp_data['dpT'] = pp_data['pTmax'] - pp_data['pTmin']
pp_data = pp_data[pp_data['pT'].between(pT_low_cut,pT_high_cut)]

x = pp_data['pT']
for rate_set in rate_names:
        color = rate_colours[rate_names[rate_set]]
        for cent, ax in zip(alice, axes):
            AA = pd.read_csv(f"../../calcs/final_runs/rset_{rate_set}/2p76/{cent}/hadron_spectra.csv", comment='#')
            AA['pT'] = 0.5 * (AA['pTmin'] + AA['pTmax'])
            AA['dpT'] = AA['pTmax'] - AA['pTmin']
            AA = AA[AA['pT'].between(pT_low_cut,pT_high_cut)]

            y1, dy1 = AA['N'], AA['dN']
            raa1 = y1/pp_data['N']
            draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + pp_data['dN']*pp_data['dN']/(pp_data['N']*pp_data['N']))
            ax.plot(x, raa1, color=color)
            ax.fill_between(x, raa1+draa1, raa1-draa1, color=color, alpha=0.2)

theo_hands = [Line2D([],[],label=l,color=c) for (l,c) in rate_colours.items()]
expt_hands = [Line2D([],[],color=css['black'], label='CMS $|\eta|<1.0$', marker='*', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ALICE $|\eta|<0.8$', marker='P', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ATLAS $|\eta|<1.0$', marker='o', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ATLAS $|\eta|<2.0$', marker='s', markersize=10, linestyle="None")]

artist = axes[2].legend(loc='lower right', handles=expt_hands,
                        ncol=1, bbox_to_anchor=(0.98, 0.01),handletextpad=0.02)
axes[2].add_artist(artist)

axes[0].legend(loc='lower left', handles=theo_hands)

axes[0].set_ylim(top=1.21,bottom=-0.19)
for cent, ax in zip(alice, axes):
    ax.set_xlabel(r'$p_T$ (GeV)')
    low, high = cent.split('_')
    centrality = f'${low}$-${high}$\%'
    ax.text(0.1,0.9,centrality, transform=ax.transAxes)
axes[1].text(0.05, 0.1, r'Pb-Pb,$\sqrt{s}=2.76$ ATeV' +'\n'+r'$|\eta|<1$', transform=axes[1].transAxes)
axes[2].text(0.05, 0.8, 'CMS: $10$-$30$\%', transform=axes[2].transAxes)
axes[0].set_ylabel(r'$R^{\mathrm{h}^{\pm}}_{\mathrm{AA}}$')
plt.savefig("../../plots/compare_coll_kernels/charged_hadron_raa.png", dpi=200)
plt.show()
