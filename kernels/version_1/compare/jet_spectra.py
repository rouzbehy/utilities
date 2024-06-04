## import scipy and matplotlib modules:
import numpy as np
from pandas import read_csv
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
from compareLabels import data as data_dirs
import sys
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
#
pTlow=30
pThigh=300

## type_of_jet = {momentum_cut : results}
parton_jet_success = {} 
parton_jet_failed  = {} 
hadron_jet         = {}

#setnum = sys.argv[1]
momentum_cut_and_colors = {#'0p2':'orange',#'0p5':'green',
                           #'1p0':'red','2p0':'blue',
                           '4T':'cyan','8T':'purple','10T':'black',
                           '4T_old':'green', '8T_old':'red', '10T_old':'blue'}
radii = ['0p2','0p3','0p4']

templ_loc_p_success = "./processed/data_set{id}/pcut_{pcut}/parton_jets_success.csv"
templ_loc_p_fail = "./processed/data_set{id}/pcut_{pcut}/parton_jets_failed.csv"
templ_loc_h = "./processed/data_set{id}/pcut_{pcut}/jet_spectra.csv"
for pcut in momentum_cut_and_colors:
    setnum = 9 if 'old' not in pcut else 8
    pcut_tmp = pcut.replace('_old','')
    #setnum = 4 if 'T' not in pcut else 6
    tmp = read_csv(templ_loc_p_success.format(id=setnum, pcut=pcut_tmp), comment='#')
    tmp['pT'] = 0.5*(tmp['pTmin']+tmp['pTmax'])
    tmp = tmp[tmp['pT'].between(pTlow,pThigh)]
    parton_jet_success[pcut] = tmp

    tmp  = read_csv(templ_loc_p_fail.format(id=setnum, pcut=pcut_tmp), comment='#')
    tmp['pT'] = 0.5*(tmp['pTmin']+tmp['pTmax'])
    tmp = tmp[tmp['pT'].between(pTlow,pThigh)]
    parton_jet_failed[pcut] = tmp

    tmp = read_csv(templ_loc_h.format(id=setnum, pcut=pcut_tmp), comment='#')
    tmp['pT'] = 0.5*(tmp['pTmin']+tmp['pTmax'])
    tmp = tmp[tmp['pT'].between(pTlow,pThigh)]
    hadron_jet[pcut] = tmp

    

fig, axes = plt.subplots(2,3, gridspec_kw={'height_ratios':(2,1),
                        'top':0.985, 'bottom':0.1, 'left':0.1,
                      'right':0.95, 'hspace':0.02,
                      'wspace':0.18},sharex=True,sharey='row')

for p in momentum_cut_and_colors:
    psuc = parton_jet_success[p]
    pfal = parton_jet_failed[p]
    had = hadron_jet[p]
    col = momentum_cut_and_colors[p]

    for ir, r in enumerate(radii):
        xp_s, yp_s, dyp_s = psuc['pT'],psuc[f'N{r}'],psuc[f'dN{r}']
        axes[0][ir].plot(xp_s, yp_s, color= col)
        axes[0][ir].fill_between(xp_s, yp_s-dyp_s,yp_s+dyp_s,color=col,alpha=0.2)

        xp_f, yp_f, dyp_f = pfal['pT'],pfal[f'N{r}'],pfal[f'dN{r}']
        axes[0][ir].plot(xp_f, yp_f, color= col, linestyle='dashed')
        axes[0][ir].fill_between(xp_f, yp_f-dyp_f,yp_f+dyp_f,color=col,alpha=0.2)

        xp_h, yp_h, dyp_h = had['pT'],had[f'N{r}'],had[f'dN{r}']
        axes[0][ir].plot(xp_h, yp_h, color= col, linestyle='dotted')
        axes[0][ir].fill_between(xp_h, yp_h-dyp_h,yp_h+dyp_h,color=col,alpha=0.2)

        ## Ratios:
        ratio_success_hadron = yp_h/yp_s
        ratio_fail_hadron = yp_f/yp_s

        axes[1][ir].plot(xp_s, ratio_success_hadron, color=col, linestyle='dotted')
        axes[1][ir].plot(xp_s, ratio_fail_hadron   , color=col, linestyle='dashed')
        R = r.replace('p','.')
        axes[0][ir].text(0.05, 0.06, f"R={R}",transform=axes[0][ir].transAxes)

axes[0][0].set_yscale('log')
axes[0][0].set_ylabel(r'$\sigma_{\mathrm{jet}}$, $|\eta|<2.0$')
axes[1][0].set_ylabel('(F.,H.)/S.')
fig.supxlabel(r'$p_T$ (GeV)')

## Legend work:
cut_legend = [Line2D([],[],label="pCut="+p.replace('p','.')+" GeV",color=l) for p,l in momentum_cut_and_colors.items()]
axes[0][2].legend(loc='upper right',handles=cut_legend)
line_legend = [Line2D([],[],label='Successful, partonic',color='black',linestyle='solid'),
               Line2D([],[],label='Failed, partonic', color='black', linestyle='dashed'),
               Line2D([],[],label='Hadronic', color='black', linestyle='dotted')]
axes[0][1].legend(loc='upper right', handles=line_legend)
axes[0][0].set_xscale('log')
plt.show()
