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
from COLORS import rate_set_colors as rate_colours

plt.rcParams.update(util.my_rcParams)
pTlow  = 60
pThigh = 300

## read in the experimental RAA:
cms = {}
cents = {'00-05':'0_5','05-10':'5_10','10-30':'10_30'}
radii = ['0p2','0p3','0p4']
for cent, label in cents.items():
    cms[label] = {}
    for R in radii:
        cms[label][R] = pd.read_csv(f'../../../exp_data/jets/CMS/jetscape/cent_{cent}_R_{R}.csv',comment='#').rename(columns=my_jdicts.colnames_CMS_RAA)

rate_names = {1:'LO', 2:'NLO',3:'NP'}

#### read in the calculated pp charged hadron spec
pp_loc  = '../../calcs/pp_2p76/jet_spec.csv'
pp_data = pd.read_csv(pp_loc, comment='#')
pp_data['pT'] = 0.5 * (pp_data['pTmin'] + pp_data['pTmax'])
pp_data['dpT'] = pp_data['pTmax'] - pp_data['pTmin']
pp_data = pp_data[pp_data['pT'].between(pTlow, pThigh)]


fig, axes = plt.subplots(3,3,
                         gridspec_kw={'top':0.985,'bottom':0.1,
                                      'left':0.08,'right':0.98,
                                      'hspace':0.05,'wspace':0.07},
                         sharex=True, figsize=(16,9), sharey=True)
#axes = axes.flatten()
axes[0][0].set_ylim(bottom=-0.09, top=1.2)
for i, r in enumerate(radii):
    rtxt = r'$R=$' + r.replace('p','.')
    axes[i][0].text(0.45, 0.1, rtxt, transform=axes[i][0].transAxes)
    for j, c in enumerate(cents.values()):
        ax = axes[i][j]
        util.plot_expr_data_on_axis(ax,cms[c][r],'*')

x = pp_data['pT']

## AA jet simulation results:
default_path = '../../calcs/final_runs/rset_{r}/2p76/{c}/jet_spectra.csv'
for rate in rate_names:
    for icol, c in enumerate(['0_5','5_10','10_20']):
        loc = default_path.format(r=rate, c=c)
        AA = pd.read_csv(loc,comment='#')
        color = rate_colours[rate_names[rate]]
        AA['pT'] = 0.5 * (AA['pTmin'] + AA['pTmax'])
        AA['dpT'] = AA['pTmax'] - AA['pTmin']
        AA = AA[AA['pT'].between(pTlow,pThigh)]
        for irow, r in enumerate(radii):
            ax = axes[irow][icol]
            Rtag = f"N{r}"
            DRtag = f"dN{r}"
            y, dy = AA[Rtag], AA[DRtag]
            raa = y/pp_data[Rtag]
            draa = raa*np.sqrt(dy*dy/(y*y) + pp_data[DRtag]*pp_data[DRtag]/(pp_data[Rtag]*pp_data[Rtag]))
            ax.plot(x, raa, color=color)
            ax.fill_between(x, raa+draa, raa-draa, color=color, alpha=0.2)

for ax, c in zip(axes[0],cents.values()):
    c = c if c!= '10_30' else '10_20'
    cent_low, cent_high = c.split('_')
    centrality = f'${cent_low}$-${cent_high}$\%'
    ax.text(0.45, 0.8, centrality, transform=ax.transAxes)

for ax in axes[2]:
    ax.set_xlabel(r'$p_T$ (GeV)')

for ax in axes[:,0]:
    ax.set_ylabel(r'$R^{\mathrm{jet}}_{\mathrm{AA}}$', fontsize=30)


axes[0][1].text(0.08, 0.1, 'PbPb, $\sqrt{s}=2.76$ATeV'+'\n'+'$|\eta|<2.0$', transform=axes[0][1].transAxes)
axes[0][2].text(0.08, 0.1, 'Data: $10$-$30$\%', transform=axes[0][2].transAxes)

expt_hands = [Line2D([],[],color=css['black'], label='CMS $|\eta|<1.0$', marker='*', markersize=10, linestyle="None")]
theor_hands = [Line2D([],[],label=l,color=c) for (l, c) in rate_colours.items()]
artist = axes[2][2].legend(loc='lower right', handles=expt_hands, ncol=2, bbox_to_anchor=(0.95, 0.2),prop={'size':20})
axes[2][2].add_artist(artist)
axes[1][1].legend(handles=theor_hands, loc='lower center', ncols=3, handletextpad=0.02, fontsize=20)

plt.savefig("../../plots/compare_coll_kernels/jet_raa.png", dpi=200)
plt.show()
