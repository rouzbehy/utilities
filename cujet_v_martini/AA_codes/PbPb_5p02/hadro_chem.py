#!/usr/bin/env python3
import numpy as np
from pandas import read_csv, DataFrame
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
import util 
import dictionaries as ddicts

my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 18,
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

def get_shape_info(fname):
    tmp = read_csv(fname,comment='#',index_col=False)
    num_evt, num_jets = 1, 1
    with open(fname, 'r') as f:
        line = f.readline()
        #print(line)
        line = line.split()
        (num_evt, num_jets) = float(line[-3]), float(line[-1])
        #print(num_evt, num_jets)
    #tmp = tmp[tmp['rmax'] < 0.31]
    delta_r = tmp['rmax'] - tmp[' rmin']
    r = 0.5*(tmp['rmax'] + tmp[' rmin'])
    rho_meson   = tmp['meson']   /( num_jets )
    drho_meson  = tmp['dmeson']  /( num_jets )
    rho_baryon  = tmp['baryon']  /( num_jets )
    drho_baryon = tmp['dbaryon'] /( num_jets )
    norm_meson  = sum(rho_meson.to_list())
    norm_baryon = sum(rho_baryon.to_list())
    #print(norm_baryon, norm_meson)
    total_norm  = norm_baryon + norm_meson 
    rho_meson_norm   = rho_meson  /( delta_r* total_norm)
    drho_meson_norm  = drho_meson /( delta_r* total_norm)
    rho_baryon_norm  = rho_baryon /( delta_r* total_norm)
    drho_baryon_norm = drho_baryon/( delta_r* total_norm)
    result = DataFrame({'r':r, 'dr':delta_r, 
                        'meson':rho_meson_norm, 
                        'dmeson':drho_meson_norm, 
                        'baryon':rho_baryon_norm, 
                        'dbaryon':drho_baryon_norm})
    return result


data_loc = "../../jetscape_data/hadroChem/"
AA_prefix = '{eloss}_PbPb5020_{cent}_jet_'
pp_prefix = 'pp_5020_jet_'
fname_tail_1 = 'hadchem_LHC_cuts_jet_rad_{R}_-2.00_2.00.csv'
fname_tail_2 = 'hadchem_LHC_cuts_jet_rad_{R}_-2.00_2.00_pTWeighted_False.csv'

pp_names_1 = data_loc + pp_prefix + fname_tail_1
pp_names_2 = data_loc + pp_prefix + fname_tail_2
AA_names_pTweight_1 = data_loc + AA_prefix + fname_tail_1
AA_names_pTweight_2 = data_loc + AA_prefix + fname_tail_2

radii = ['0.2','0.3','0.4','0.6','0.8']
eloss = ['cujet','martini']
linestyles = {'baryon':'dashed', 'meson':'solid'}
xticks = {'0.2':[0., 0.05, 0.1, 0.15],
          '0.3':[0., 0.10, 0.15, 0.20, 0.25],
          '0.4':[0., 0.20, 0.4],
          '0.6':[0., 0.30, 0.6],
          '0.8':[0., 0.40, 0.8]}

pp_results_1 = {r:get_shape_info(pp_names_1.format(R=r)) for r in radii}
pp_results_2 = {r:get_shape_info(pp_names_2.format(R=r)) for r in radii}

fig, axes = plt.subplots(nrows=3, ncols=5, 
                        gridspec_kw={'left':0.08,'right':0.85,
                                     'bottom':0.1,'top':0.95},
                        figsize=(16,9), sharex='col',sharey='row')

fig2, axes2 = plt.subplots(nrows=3, ncols=5, 
                        gridspec_kw={'left':0.08,'right':0.85,
                                     'bottom':0.1,'top':0.95},
                        figsize=(16,9), sharex='col',sharey=True)

for icent, cent in enumerate(['00-10','10-20','30-50']):
    row = axes[icent]
    row[-1].text(1.05,0.8,f'Cent: {cent}\%', transform=row[-1].transAxes)
    row2 = axes2[icent]
    row2[-1].text(1.05,0.5,f'Cent: {cent}\%', transform=row2[-1].transAxes, fontsize=18)

    for ir, r in enumerate(radii):
        ax = row[ir]
        ax2 = row2[ir]
        ax.axhline(y=1, linestyle='dotted', color='black')
        ax2.axhline(y=1, linestyle='dotted', color='black')
        if icent == 0:
            ax.set_title(f'R={r}' )
            ax2.set_title(f'R={r}')
        if icent == 2:
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax2.xaxis.set_major_locator(MaxNLocator(4))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        pp_1 = pp_results_1[r]
        pp_2 = pp_results_2[r]
        for eloss in ['cujet','martini']:
            col = ddicts.eloss_colours[eloss]
            AA_1 = get_shape_info(AA_names_pTweight_1.format(eloss=eloss, R=r, cent=cent))
            AA_2 = get_shape_info(AA_names_pTweight_2.format(eloss=eloss, R=r, cent=cent))
            for part in ['meson','baryon']:
                lstyle = linestyles[part]
                raa_1 = AA_1[part]/pp_1[part]
                y1,dy1 = AA_1[part], AA_1[f'd{part}']
                y2,dy2 = pp_1[part], pp_1[f'd{part}']
                draa_1 = raa_1*np.sqrt(dy1*dy1/(y1*y1) + dy2*dy2/(y2*y2))
                ax.plot(AA_1['r'],raa_1,color=col, linestyle=lstyle)
                ax.fill_between(AA_1['r'], raa_1+draa_1, raa_1-draa_1, color=col, alpha=0.2)

                raa_2  = AA_2[part]/ pp_2[part]
                y1,dy1 = AA_2[part], AA_2[f'd{part}']
                y2,dy2 = pp_2[part], pp_2[f'd{part}']
                draa_2 = raa_2*np.sqrt(dy1*dy1/(y1*y1) + dy2*dy2/(y2*y2))
                ax2.plot(AA_2['r'],raa_2,color=col, linestyle=lstyle)
                ax2.fill_between(AA_2['r'], raa_2+draa_2, raa_2-draa_2, color=col, alpha=0.2)

theo_labels = [Line2D([],[],color=c, label=l.upper()) for l, c in ddicts.eloss_colours.items() if 'no' not in l]
for l, lstyle in linestyles.items():
    theo_labels.append(Line2D([],[],color='black',label=l,linestyle=lstyle))

axes[0][0].legend(handles=theo_labels, fontsize=16, loc='upper left')#bbox_to_anchor=(0.04,0.8), fontsize=16)
axes2[0][0].legend(handles=theo_labels, fontsize=16, loc='upper left')#bbox_to_anchor=(0.04,0.8), fontsize=16)
# line_labels = []
# axes[1][-1].legend(handles=line_labels,bbox_to_anchor=(1.0,0.8), fontsize=16)

fig.supxlabel('$r$')
fig.supylabel(r'$R_{\rho}= \rho_{\mathrm{AA}}(r)/\rho_{\mathrm{pp}}(r)$, momentum weighted')

fig2.supxlabel('$r$')
fig2.supylabel(r'$R_{\rho}= \zeta_{\mathrm{AA}}(r)/\zeta_{\mathrm{pp}}(r)$')
axes[0][2].text(0.05,0.9,'Pb-Pb @ 5.02 ATeV', transform=axes[0][2].transAxes)
plt.show()