## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
## my custom modules
import util
import dictionaries as my_dicts
import getSpecs as spc 
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
pT_low_cut, pT_high_cut = 10, 100

## read in the experimental RAA:
cms     = pd.read_csv('../../../exp_data/charged/CMS/Table5_1.csv',comment='#').rename(columns=my_dicts.colnames_raa_cms)
alice   = pd.read_csv('../../../exp_data/charged/ALICE/PbPb_RAA_cent_00_05_eta_0p8.csv',comment='#').rename(columns=my_dicts.colnames_raa_alice)
atlas_1 = pd.read_csv("../../../exp_data/charged/ATLAS/Table49.csv",comment='#').rename(columns=my_dicts.colnames_raa_atlas)
atlas_2 = pd.read_csv("../../../exp_data/charged/ATLAS/Table33.csv",comment='#').rename(columns=my_dicts.colnames_raa_atlas) 
atlas_1 = atlas_1[atlas_1['x'].between(8, 104)]
atlas_2 = atlas_2[atlas_2['x'].between(8, 104)]
from COLORS import rate_set_colors
rate_colours = rate_set_colors 
rate_names = {1:'LO', 2:'NLO',3:'NP'}


pp_loc  = '../martini_results/pp/hadron_spectra.csv'
pp_data = spc.get_charged_spec(pp_loc, pT_low_cut, pT_high_cut )
## now read in fixed alphas where all rate sets have alphas=0.3
fname_tmp = '../martini_results/fixed_alphas_const/rset_{r}/hadron_spectra.csv'
#fname_tmp = '../martini_results/final_PbPb_2p76/rset_{r}/cent_0_5/hadron_spectra.csv'
AA_fxd_c = {r: spc.get_charged_spec(fname_tmp.format(r=r), pT_low_cut, pT_high_cut ) for r in rate_names}
fname_tmp = '../martini_results/fixed_alphas_fitted/rset_{r}/hadron_spectra.csv'
AA_fxd_f = {r: spc.get_charged_spec(fname_tmp.format(r=r), pT_low_cut, pT_high_cut ) for r in rate_names}
fname_tmp = '../martini_results/run_alphas_fit/rset_{r}/hadron_spectra.csv'
AA_run_f_v2 = {r: spc.get_charged_spec(fname_tmp.format(r=r), pT_low_cut, pT_high_cut ) for r in rate_names}
#fname_tmp = '../martini_results/rset_{r}/hadron_spectra.csv'
#AA_run_f_v1 = {r: spc.get_charged_spec(fname_tmp.format(r=r), pT_low_cut, pT_high_cut ) for r in rate_names}

fig, ax = plt.subplots(1,1,gridspec_kw={'top':0.985,'bottom':0.11,
                                      'left':0.08,'right':0.98,
                                      'hspace':0.11,'wspace':0.09},
                              sharex=True, figsize=(16,9), sharey=True)

cms = cms[cms['x'].between(pT_low_cut, pT_high_cut)]
alice = alice[alice['x'].between(pT_low_cut, pT_high_cut)]
util.plot_expr_data_on_axis(ax,cms,'*')
util.plot_expr_data_on_axis(ax,alice,'o')
util.plot_expr_data_on_axis(ax,atlas_1,'P')
util.plot_expr_data_on_axis(ax,atlas_2,'s')

x = pp_data['pT']
yp, dyp = pp_data['N'], pp_data['dN']
err_fac_pp = dyp*dyp/(yp*yp)
for rate_set in rate_names:
        color = rate_colours[rate_names[rate_set]]
        AA1 = AA_run_f_v2[rate_set]
        y1, dy1 = AA1['N'], AA1['dN']
        raa1 = y1/yp
        draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + err_fac_pp)
        ax.plot(x, raa1, color=color) 
        ax.fill_between(x, raa1+draa1, raa1-draa1, color=color, alpha=0.2)

        AA2 = AA_fxd_f[rate_set]
        y2, dy2 = AA2['N'], AA2['dN']
        raa2 = y2/yp
        draa2 = raa2*np.sqrt(dy2*dy2/(y2*y2) + err_fac_pp)
        ax.plot(x, raa2, color=color, linestyle='dashed') 
        ax.fill_between(x, raa2+draa2, raa2-draa2, color=color, alpha=0.2)

        AA3 = AA_fxd_c[rate_set]
        y3, dy3 = AA3['N'], AA3['dN']
        raa3 = y3/yp
        draa3 = raa3*np.sqrt(dy3*dy3/(y3*y3) + err_fac_pp)
        ax.plot(x, raa3, color=color, linestyle='dotted') 
        ax.fill_between(x, raa3+draa3, raa3-draa3, color=color, alpha=0.2)


expt_hands = [Line2D([],[],color=css['black'], label='CMS $|\eta|<1.0$', marker='*', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ALICE $|\eta|<0.8$', marker='o', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ATLAS $|\eta|<1.0$', marker='P', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ATLAS $|\eta|<2.0$', marker='s', markersize=10, linestyle="None")]


extra_labels = {r'$\alpha_s=\alpha_s(\mu)$':'solid',
                r'$\alpha_s=\alpha_{s,i}$':'dashed',
                r'$\alpha_s=0.3$':'dotted'}
theor_hands = [Line2D([],[],label=l,color=c) for (l,c) in rate_colours.items()]
for l, st in extra_labels.items():
    theor_hands.append(Line2D([],[],label=l,linestyle=st,color='black'))

artist = ax.legend(loc='lower right', handles=expt_hands, ncol=2, bbox_to_anchor=(0.95, 0.05),prop={'size':20})
ax.add_artist(artist)
ax.legend(loc='upper left', handles=theor_hands, ncol=1, fontsize=18)

ax.set_ylim(top=0.88,bottom=-0.01)
ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30)
ax.set_ylabel(r'$R^{\mathrm{h}^{\pm}}_{\mathrm{AA}}$', fontsize=30)
ax.text(0.75, 0.9, r'Pb-Pb \@ $\sqrt{s}=2.76$ ATeV' +'\n'+r'$|\eta|<1$, $0$-$5\%$', transform=ax.transAxes) 
plt.show()
