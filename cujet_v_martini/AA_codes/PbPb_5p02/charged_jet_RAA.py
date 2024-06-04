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

pT_cut_low, pT_cut_high = 20, 200
pp_specs = {}
AA_specs = {}
for radius in [0.2, 0.3, 0.4, 0.6]:
    fname_pp = "../../jetscape_data/max_time/maxT_200_highstat/pp_5020_jet_spec_jet_rad_{r}_{eta}.csv"
    r = str(radius)#.replace('.','p')
    eta_window = 0.7 - radius
    eta = f'{eta_window:0.2f}'
    tmp = read_csv(fname_pp.format(r=r,eta=eta), comment='#')
    tmp['pT'] = 0.5*(tmp['ptmin'] + tmp['ptmax'])
    tmp = tmp[tmp['pT'].between(pT_cut_low, pT_cut_high)]
    pp_specs[r] = tmp

for eloss in ['cujet','martini']:
    AA_specs[eloss] = {}
    for cent in ['00-10','10-20','30-50']:
        AA_specs[eloss][cent] = {}
        for radius in [0.2, 0.3, 0.4, 0.6]:
            fname_AA = '../../jetscape_data/sqrt_s_5020/maxt_200/{eloss}/PbPb5020_{cent}_jet_spec_jet_rad_{r}_{eta}.csv'
            r = str(radius)#.replace('.','p')
            eta_window = 0.7 - radius
            eta = f'{eta_window:0.2f}'
            tmp = read_csv(fname_AA.format(eloss=eloss, cent=cent, r=r, eta=eta),comment='#')
            tmp['pT'] = 0.5*(tmp['ptmin'] + tmp['ptmax'])
            tmp = tmp[tmp['pT'].between(pT_cut_low, pT_cut_high)]
            AA_specs[eloss][cent][r] = tmp

fig, axes = plt.subplots(2,2,gridspec_kw={'top':0.98,'bottom':0.12,
                                          'left':0.10,'right':0.75,
                                          'hspace':0.1,'wspace':0.1},
                                          sharex=True, sharey=True, figsize=(16,9))
axes = axes.flatten()
linestyles = {'00-10':'solid','10-20':'dashed','30-50':'dotted'}
for ir, r in enumerate(['0.2','0.3','0.4','0.6']):
    ax = axes[ir]
    for eloss in AA_specs:
        color = ddicts.eloss_colours[eloss]
        for indx, cent in enumerate(AA_specs[eloss]):
            line_style = linestyles[cent]
            #ax = axes[indx]
            pp = pp_specs[r]
            aa = AA_specs[eloss][cent][r]
            raa = aa['Chg']/pp['Chg']
            draa = raa*np.sqrt(aa['dChg']*aa['dChg']/(aa['Chg']*aa['Chg']) + \
                               pp['dChg']*pp['dChg']/(pp['Chg']*pp['Chg']))

            ax.plot(pp['pT'], raa, color=color, linestyle=line_style)
            ax.fill_between(pp['pT'], raa+draa, raa-draa, color=color, alpha=0.3)
            eta = 0.7-float(r)
            ax.text(0.7,0.1,f'R={r}',transform=ax.transAxes)

ax = axes[1]
label_set_1 = [Line2D([],[],label='MATTER+'+eloss.upper(), color=col) for eloss, col in ddicts.eloss_colours.items()]
for cent, lstyle in linestyles.items():
    label_set_1.append(Line2D([],[],label=cent+'$\%$',color='black',linestyle=lstyle))

artist = ax.legend(handles=label_set_1, bbox_to_anchor=(1.,0.7), fontsize=20)
ax.add_artist(artist)

for ax in axes[2:4]:
    ax.set_xlabel(r'$p^{\mathrm{jet}^{\pm}}_T$', fontsize=35)
fig.supylabel(r'$R^{\mathrm{jet}^{\pm}}_{\mathrm{AA}}$')

axes[3].text(1.05, 0.5, 
            r'Pb-Pb $\sqrt{s}=5.02$ ATeV,'+'\n'+r'$|\eta^{\mathrm{jet^{\pm}}}|<0.7-$R', 
            transform=axes[3].transAxes, fontsize=20)
fig.savefig('../../../Thesis_Plots/sqrt_s_5p02/jets/charged_jet_RAA_0p7-r.png', dpi=200)

plt.show()
