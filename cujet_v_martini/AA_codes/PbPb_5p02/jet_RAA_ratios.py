#!/usr/bin/env python3
import numpy as np
from pandas import read_csv, DataFrame
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d
from matplotlib import ticker
import util 
import jetDicts as ddicts
from COLORS import module_colors
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
sigma = 62.049
experiment_colors = {'alice':"#7cae7a",
                     'atlas':"#d64045", 
                     'cms':"#eca400"}
#cut_min, cut_max = 7, 440
exp_loc="../../../expt/PbPb_5p02/Jets/"

pT_cut_low, pT_cut_high = 85, 810
pp_specs = {}
AA_specs = {}
for radius in [0.2, 0.3, 0.4, 0.6, 0.8]:
    fname_pp = "../../jetscape_data/max_time/maxT_200_highstat/pp_5020_jet_spec_jet_rad_{r}_2.00.csv"
    r = str(radius)#.replace('.','p')
    #eta_window = 0.9 - radius
    #eta = f'{eta_window:0.2f}'
    tmp = read_csv(fname_pp.format(r=r), comment='#')
    tmp['pT'] = 0.5*(tmp['ptmin'] + tmp['ptmax'])
    tmp = tmp[tmp['pT'].between(pT_cut_low, pT_cut_high)]
    pp_specs[r] = tmp

for eloss in ['cujet','martini']:
    AA_specs[eloss] = {}
    for cent in ['00-10','10-20','30-50']:
        AA_specs[eloss][cent] = {}
        for radius in [0.2, 0.3, 0.4, 0.6,0.8]:
            fname_AA = '../../jetscape_data/sqrt_s_5020/maxt_200/{eloss}/PbPb5020_{cent}_jet_spec_jet_rad_{r}_2.00.csv'
            r = str(radius)#.replace('.','p')
            #eta_window = 0.9 - radius
            #eta = f'{eta_window:0.2f}'
            tmp = read_csv(fname_AA.format(eloss=eloss, cent=cent, r=r),comment='#')
            tmp['pT'] = 0.5*(tmp['ptmin'] + tmp['ptmax'])
            tmp = tmp[tmp['pT'].between(pT_cut_low, pT_cut_high)]
            AA_specs[eloss][cent][r] = tmp

fig, axes = plt.subplots(2,2,sharex=True, sharey='row', figsize=(16,9))
axes = axes.flatten()
linestyles = {'00-10':'solid','10-20':'dashed','30-50':'dotted'}
for ir, r in enumerate(['0.3','0.4','0.6','0.8']):
    ax = axes[ir]
    ax.set_ylim(top=3.75,bottom=1.0)
    ax.set_xscale('log')

    for eloss in AA_specs:
        color =module_colors['MATTER+'+eloss.upper()] 
        for indx, cent in enumerate(AA_specs[eloss]):
            if cent != '00-10':
                continue
            line_style = linestyles[cent]
            pp = pp_specs[r]
            aa = AA_specs[eloss][cent][r]
            y_0p2, dy_0p2 = AA_specs[eloss][cent]['0.2']['Ncut'], AA_specs[eloss][cent]['0.2']['dNcut']
            raa_0p2 = y_0p2/pp['Ncut']
            draa_0p2 = raa_0p2*np.sqrt(dy_0p2*dy_0p2/(y_0p2*y_0p2)+ \
                               pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']))
            raa = aa['Ncut']/pp['Ncut']
            draa = raa*np.sqrt(aa['dNcut']*aa['dNcut']/(aa['Ncut']*aa['Ncut']) + \
                               pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']))
            raa_ratio = raa/raa_0p2
            draa_ratip = raa_ratio*np.sqrt(draa_0p2*draa_0p2/(raa_0p2*raa_0p2)+\
                                            draa*draa/(raa*raa))

            ax.plot(pp['pT'], raa_ratio, color=color, linestyle=line_style)
            ax.fill_between(pp['pT'], raa_ratio+draa_ratip, raa_ratio-draa_ratip, color=color, alpha=0.3)
            eta = 2.0 #0.9-float(r)
            ax.text(0.7,0.8,f'R={r}',transform=ax.transAxes)
for ax in axes:
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    # ax.xaxis.set_major_formatter(plt.FixedFormatter(['100','200', '400','800']))
    ax.set_xticks([100, 200, 400, 800])
    ax.set_xticklabels(['100','200','400','800'])

ax = axes[0]
label_set_1 = [Line2D([],[],label='MATTER+'+eloss.upper(), color=module_colors['MATTER+'+eloss.upper()]) \
                for eloss in ['cujet','martini']]

for cent, lstyle in linestyles.items():
    if cent != '00-10':
        continue
    label_set_1.append(Line2D([],[],label=cent+'$\%$',color='black',linestyle=lstyle))

artist = ax.legend(handles=label_set_1, fontsize=20, loc='upper left')#, bbox_to_anchor=(1.,0.7), fontsize=20)
ax.add_artist(artist)

for ax in axes[2:]:
    ax.set_xlabel(r'$p_T$ (GeV)', fontsize=35)
for ax in [axes[0], axes[2]]:
    ax.set_ylabel(r'$R^{\mathrm{jet}}_{\mathrm{AA}}(R)/R^{\mathrm{jet}}_{\mathrm{AA}}(0.2)$')

axes[2].text(0.1, 0.8, 
            r'Pb-Pb $\sqrt{s}=5.02$ ATeV,'+'\n'+r'$|\eta^{\mathrm{jet}}|<2.0$'+'\n'+'Inclusive Jets', 
            transform=axes[2].transAxes, fontsize=20)
plt.show()