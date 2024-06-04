## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
import util
import dictionaries as my_hdicts
import jetDicts as my_jdicts
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

pTlow  = 60
pThigh = 300


## read in the experimental RAA:
cms = {}
for R in ['0p2','0p3','0p4']:
    cms[R] = pd.read_csv(f'../../../exp_data/jets/CMS/Table_JETRAA_{R}.csv',comment='#').rename(columns=my_jdicts.colnames_CMS_RAA)
from COLORS import rate_set_colors
rate_colours = rate_set_colors
rate_names = {1:'LO', 2:'NLO',3:'NP'}

#### read in the calculated pp jet spectra 
pp_loc  = '../martini_results/pp/jet_spectra.csv'
pp_data = spc.get_jet_spec(pp_loc, pTlow, pThigh )
## now read in fixed alphas where all rate sets have alphas=0.3
fname_tmp = '../martini_results/fixed_alphas_const/rset_{r}/jet_spectra.csv'
AA_fxd_c = {r: spc.get_charged_spec(fname_tmp.format(r=r), pTlow, pThigh ) for r in rate_names}
fname_tmp = '../martini_results/fixed_alphas_fitted/rset_{r}/jet_spectra.csv'
AA_fxd_f = {r: spc.get_charged_spec(fname_tmp.format(r=r), pTlow, pThigh ) for r in rate_names}
fname_tmp = '../martini_results/run_alphas_fit/rset_{r}/jet_spectra.csv'
AA_run_f_v2 = {r: spc.get_charged_spec(fname_tmp.format(r=r), pTlow, pThigh ) for r in rate_names}
#fname_tmp = '../martini_results/rset_{r}/jet_spectra.csv'
#AA_run_f_v1 = {r: spc.get_charged_spec(fname_tmp.format(r=r), pTlow, pThigh ) for r in rate_names}

fig, axes = plt.subplots(3,1,gridspec_kw={'top':0.985,'bottom':0.105,
                                      'left':0.08,'right':0.8,
                                      'hspace':0.1,'wspace':0.09},
                              sharex=True, figsize=(16,9), sharey=True)

axes = axes.flatten()
for ir, r in enumerate(cms):
    ax = axes[ir]
    util.plot_expr_data_on_axis(ax,cms[r],'*')
    rtxt = r'$R=$' + r.replace('p','.')
    ax.text(0.85, 0.95, rtxt, transform=ax.transAxes)

x = pp_data['pT']
for rset in rate_names:
    color = rate_colours[rate_names[rset]]
    #AA_1 = AA_run_f_v2[rset]
    AA_2 = AA_fxd_f[rset]
    AA_3 = AA_fxd_c[rset]
    for ir, r in enumerate(cms):
        ax = axes[ir]
        Rtag = f"N{r}"
        DRtag = f"dN{r}"
        #y1, dy1 = AA_1[Rtag], AA_1[DRtag]
        y2, dy2 = AA_2[Rtag], AA_2[DRtag]
        y3, dy3 = AA_3[Rtag], AA_3[DRtag] 
        p, dp   = pp_data[Rtag], pp_data[DRtag]
        err_fac = dp*dp/(p*p)
        #raa1 = y1/p
        raa2 = y2/p
        raa3 = y3/p
        #draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + err_fac)
        draa2 = raa2*np.sqrt(dy2*dy2/(y2*y2) + err_fac)
        draa3 = raa3*np.sqrt(dy3*dy3/(y3*y3) + err_fac)
        #ax.plot(x, raa1, color=color, linestyle='solid') ##Fitted, running 
        #ax.fill_between(x, raa1+draa1, raa1-draa1, color=color, alpha=0.2)
        ax.plot(x, raa2, color=color, linestyle='solid') ## fitted, constant
        ax.fill_between(x, raa2+draa2, raa2-draa2, color=color, alpha=0.2)
        ax.plot(x, raa3, color=color, linestyle='dashed') ## constant, same
        ax.fill_between(x, raa3+draa3, raa3-draa3, color=color, alpha=0.2)


handles = [Line2D([],[],color=css['black'], label='CMS, $|\eta^{\mathrm{jet}}|<2.0$, $0$-$5\%$', marker='*', markersize=10, linestyle="None")]
for (l,c) in rate_colours.items():
    handles.append(Line2D([],[],label=l,color=c))
extra_labels = {#r'$\alpha_s=\alpha_s(\mu)$':'solid',
                r'$\alpha_s=\alpha_{s,i}$':'solid',
                r'$\alpha_s=0.3$':'dashed'}
for l, st in extra_labels.items():
    handles.append(Line2D([],[],label=l,linestyle=st,color='black'))
axes[0].legend(loc='upper left', fontsize=18, handles=handles, bbox_to_anchor=(1.0,0.8))


ax = axes[2]
ax.set_ylim(bottom=-0.01, top=1.01)
ax.set_xlim(left=50, right=320)
ax.set_xlabel(r'$p^{\mathrm{jet}}_T$ (GeV)', fontsize=30)
for ax in axes:
    ax.set_ylabel(r'$R^{\mathrm{jet}}_{\mathrm{AA}}$', fontsize=30)


plt.show()