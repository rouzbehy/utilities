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
import dictionaries as my_dicts
from compareLabels import data as data_dict
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

pT_low_cut, pT_high_cut = 1, 900
momentum_cut_and_colors = {#'0p2':'orange',#'0p5':'green',
                           #'1p0':'red','2p0':'blue',
                           '4T':'cyan','8T':'purple','10T':'black',
                           '4T_old':'green', '8T_old':'red', '10T_old':'blue'}
line_style = {'f':'solid','g':'dashed'}

tmpl_loc_p = "./processed/data_set{id}/pcut_{pcut}/parton_spectra.csv"
templ_loc_h = "./processed/data_set{id}/pcut_{pcut}/hadron_spectra.csv"

hadronic_spec = {}
partonic_spec = {}
for pcut in momentum_cut_and_colors:
    setnum = 9 if 'old' not in pcut else 8
    pcut_tmp = pcut.replace('_old','')
    tmp = read_csv(tmpl_loc_p.format(id=setnum, pcut=pcut_tmp), comment='#')
    tmp['pT'] = 0.5*(tmp['pTmin']+tmp['pTmax'])
    tmp = tmp[tmp['pT'].between(pT_low_cut,pT_high_cut)]
    partonic_spec[pcut] = tmp

    tmp = read_csv(templ_loc_h.format(id=setnum, pcut=pcut_tmp), comment='#')
    tmp['pT'] = 0.5*(tmp['pTmin']+tmp['pTmax'])
    tmp = tmp[tmp['pT'].between(pT_low_cut,pT_high_cut)]
    hadronic_spec[pcut] = tmp



fig, axes = plt.subplots(2,2, gridspec_kw={'height_ratios':(2,1),
                        'top':0.985, 'bottom':0.1, 'left':0.1,
                      'right':0.95, 'hspace':0.02,
                      'wspace':0.18},sharex=True,sharey='row')
axes = axes.flatten()
axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[2].set_yscale('log')
for p in momentum_cut_and_colors:
    par = partonic_spec[p]
    had = hadronic_spec[p]
    col = momentum_cut_and_colors[p]

    xp_h, yp_h, dyp_h = had['pT'],had[f'N'],had[f'dN']

    for ax in axes[0:2]:
        ax.plot(xp_h, yp_h, color= col, linestyle='dotted')
        ax.fill_between(xp_h, yp_h-dyp_h,yp_h+dyp_h,color=col,alpha=0.2)

    for parton in ['f','g']:
        lstyle = line_style[parton]
        for indx, tag in enumerate(['Failed','Success']):
            ax = axes[indx]
            yname, dyname = f'{parton}{tag}', f'd{parton}{tag}'
            x, y, dy = par['pT'], par[yname], par[dyname]

            ax.plot(x, y, color= col, linestyle=lstyle)
            ax.fill_between(x, y-dy,y+dy,color=col,alpha=0.2)
            ax.text(0.2,0.8,tag+" Events",transform=ax.transAxes)

            ## for ratios, take parton/hadron
            ratio = y/yp_h
            ax = axes[indx+2]
            ax.plot(x, ratio, color=col, linestyle=lstyle)

fig.supxlabel(r'$p_T$ (GeV)')

axes[0].text(0.09,0.05,r'$|\eta|<1.0$', transform=axes[0].transAxes)
axes[0].set_ylabel(r'$\sigma_{h^{\pm}}$')
axes[2].set_ylabel(r'Ratio to $h^{\pm}$')

cut_legend = [Line2D([],[],label="pCut="+p.replace('p','.')+" GeV",color=l) for p,l in momentum_cut_and_colors.items()]
axes[1].legend(loc='upper right',handles=cut_legend, fontsize=16)
line_legend = [Line2D([],[],label=r'$q(\bar{q})$',color='black',linestyle='solid'),
               Line2D([],[],label=r'$g$', color='black', linestyle='dashed'),
               Line2D([],[],label=r'$h^{\pm}$', color='black', linestyle='dotted')]
axes[3].legend(loc='upper center', handles=line_legend,ncol=3)
plt.show() 
