#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.interpolate import InterpolatedUnivariateSpline
import util 
import dictionaries as ddicts
from COLORS import module_colors as colors
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 30,
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

elosses = ['matter','cujet','martini']
matter_options = [200]
pT_lower_lim = 8
pT_upper_lim = 175
inel_Xsec= 62.03948 

absol_path_expt='/Users/rmyazdi/Documents/research/jetscape_project/expt/PbPb_2p76/charged'
centralities = {0:'00-05'}

alice_experimental_data = {}
for icent, cent in centralities.items():
    indx = 16 + icent
    tmp = pd.read_csv(f"{absol_path_expt}/ALICE/Table{indx:d}.csv",comment='#').rename(columns=ddicts.colnames_raa_alice)
    tmp = tmp[tmp['xhigh'].between(pT_lower_lim,pT_upper_lim)]
    alice_experimental_data[cent] = tmp

atlas2_experimental_data = {}
for icent, cent in centralities.items():
    indx = 33 + icent
    tmp = pd.read_csv(f"{absol_path_expt}/ATLAS/Table{indx:d}.csv",comment='#').rename(columns=ddicts.colnames_raa_atlas)
    tmp = tmp[tmp['xhigh'].between(pT_lower_lim,pT_upper_lim)]
    atlas2_experimental_data[cent] = tmp

atlas1_experimental_data = {}
for icent, cent in enumerate(['00-05']):
    indx = 49
    tmp = pd.read_csv(f"{absol_path_expt}/ATLAS/Table{indx:d}.csv",comment='#').rename(columns=ddicts.colnames_raa_atlas)
    tmp = tmp[tmp['xhigh'].between(pT_lower_lim,pT_upper_lim)]
    atlas1_experimental_data[cent] = tmp

cms_experimental_data = {}
for icent, cent in enumerate(['00-05']):
    tmp = pd.read_csv(f"{absol_path_expt}/CMS/Table5_{cent}.csv",comment='#').rename(columns=ddicts.colnames_raa_cms)
    tmp = tmp[tmp['xhigh'].between(pT_lower_lim,pT_upper_lim)]
    cms_experimental_data[cent] = tmp

pp_spec = pd.read_csv("../../../jetscape_data/max_time/maxT_200_highstat/pp_2760_charged_hadrons_eta_cut_1.csv") ## eta < 1, ptmin, ptmax, N, dN^2
pp_spec['pT']  = 0.5*(pp_spec['pTmax']+pp_spec['pTmin'])
pp = pp_spec[pp_spec['pT'].between(pT_lower_lim,pT_upper_lim)]

AA_spec = {}
for eloss in elosses:
    AA_spec[eloss] = {}
    for opt in matter_options:
        AA_spec[eloss][opt] = {}
        for icent, cent in centralities.items():
            fname  = f'../../../jetscape_data/sqrt_s_2760/{eloss}'
            if eloss == 'matter':
                fname += f'/PbPb_2760_maxT_{opt}/PbPb2760_{cent}_charged_hadrons_eta_cut_1.csv'
            else:
                fname += f'/PbPb_2760/PbPb2760_{cent}_charged_hadrons_eta_cut_1.csv'
            tmp = pd.read_csv(fname,comment='#')
            tmp['pT']  = 0.5*(tmp['pTmax']+tmp['pTmin'])
            tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
            AA_spec[eloss][opt][cent] = tmp

fig, ax = plt.subplots(1,1,gridspec_kw={'top':0.985,'bottom':0.1,
                                          'left':0.10,'right':0.95,
                                          'hspace':0.1,'wspace':0.0},
                                          sharex=True, sharey=True, figsize=(16,9))

for icent, cent in centralities.items():
    alice = alice_experimental_data[cent]
    util.plot_expr_data_on_axis(ax, alice, marker=ddicts.markers['alice'], s=50)

    atlas2 = atlas2_experimental_data[cent]
    util.plot_expr_data_on_axis(ax, atlas2, marker=ddicts.markers['atlas2'], s=50)

    if cent in cms_experimental_data:
        cms = cms_experimental_data[cent]
        util.plot_expr_data_on_axis(ax, cms, marker=ddicts.markers['cms'], s=50) 
   
(ptmin, ptmax, proton, dproton) = util.combine_bins(pp['pTmin'].to_list(), pp['pTmax'].to_list(),
                                                    pp['N'].to_list(), pp['dN'].to_list())
pT = 0.5*(ptmax + ptmin)
#eloss = 'matter'
#for opt in matter_options:
#    for icent, cent in centralities.items():
#        aa = AA_spec[eloss][opt][cent]
#        (_, _, PbPb, dPbPb) = util.combine_bins(aa['pTmin'].to_list(), \
#                                                aa['pTmax'].to_list(),\
#                                                aa['N'].to_list(),\
#                                                aa['dN'].to_list())
#        raa = PbPb/proton
#        err = raa*np.sqrt(dPbPb*dPbPb/(PbPb*PbPb) + dproton*dproton/(proton*proton))
#        #coltag = f"{eloss}_{opt}"
#        color = ddicts.eloss_colours[coltag]
#        ax.plot(pT, raa, color=color)
#        ax.fill_between(pT, raa+err, raa-err, color=color, alpha=0.3)

for eloss in elosses:
    #if eloss == 'matter':
    #    continue
    for icent, cent in centralities.items():
        aa = AA_spec[eloss][opt][cent]
        (_, _, PbPb, dPbPb) = util.combine_bins(aa['pTmin'].to_list(), \
                                                aa['pTmax'].to_list(),\
                                                aa['N'].to_list(),\
                                                aa['dN'].to_list())
        raa = PbPb/proton
        err = raa*np.sqrt(dPbPb*dPbPb/(PbPb*PbPb) + dproton*dproton/(proton*proton))
        #coltag = f"{eloss}"
        coltag = 'MATTER'
        if eloss in ['cujet', 'martini']:
            coltag += f'+{eloss}'
            coltag=coltag.upper()
        color = colors[coltag]#ddicts.eloss_colours[coltag]
        ax.plot(pT, raa, color=color)
        ax.fill_between(pT, raa+err, raa-err, color=color, alpha=0.3)

ax.set_ylabel(r"$R^{h^{\pm}}_{\mathrm{AA}}$", fontsize=30)
ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30) 
handles_experiment = [Line2D([],[],marker=ddicts.markers['alice'] ,markersize=10,label=r'ALICE $|\eta|<1.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['atlas2'],markersize=10,label=r'ATLAS $|\eta|<2.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['atlas1'],markersize=10,label=r'ATLAS $|\eta|<1.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['cms'],markersize=10, label=r'CMS $|\eta|<1.0$',color='black')]
handles_theory = [Line2D([],[],color=ddicts.eloss_colours[eloss],label=label) for eloss,label in zip(['matter_20','matter_200'],['MATTER (maxT=20)', 'MATTER (maxT=200)'])]
handles_theory.append(Line2D([],[],color=ddicts.eloss_colours['martini'], label='MATTER+MARTINI'))
handles_theory.append(Line2D([],[],color=ddicts.eloss_colours['cujet'], label='MATTER+CUJET'))

handles_theory = [Line2D([],[],color=colors[l], label=l) for l in colors]

artist=ax.legend(handles=handles_experiment,loc='upper left')
ax.add_artist(artist)
ax.legend(handles=handles_theory,loc='upper right')
ax.set_ylim(bottom=-0.01, top=1.51)
ax.text(0.65,0.1,r'$|\eta|<1.0$', transform=ax.transAxes)
ax.text(0.65,0.2,r'Pb-Pb, $\sqrt{s}=2.76$ ATeV', fontsize=30, transform=ax.transAxes)
ax.text(0.65,0.15,r'$0$-$5\%$', transform=ax.transAxes)
plt.show()
