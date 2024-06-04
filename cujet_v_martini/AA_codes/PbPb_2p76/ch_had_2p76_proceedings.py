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

my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 25,
    "lines.linewidth": 4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size" : 12,
    "ytick.major.size" : 12,
    "xtick.minor.size" : 6,
    "ytick.minor.size" : 6,
    "axes.spines.right": False,
    "axes.spines.top" : False,
    "legend.frameon":False
}
plt.rcParams.update(my_rcParams)

elosses = ['cujet', 'martini']

pT_lower_lim = 19
pT_upper_lim = 175
inel_Xsec= 62.03948 
from COLORS import module_colors
absol_path_expt='/Users/rmyazdi/Documents/research/jetscape_project/expt/PbPb_2p76/charged'
centralities = {0:'00-05',1:'05-10',2:'10-20',
                3:'20-30',4:'30-40',5:'40-50'}

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
for icent, cent in enumerate(['00-05','05-10','10-30','30-50']):
    tmp = pd.read_csv(f"{absol_path_expt}/CMS/Table5_{cent}.csv",comment='#').rename(columns=ddicts.colnames_raa_cms)
    tmp = tmp[tmp['xhigh'].between(pT_lower_lim,pT_upper_lim)]
    cms_experimental_data[cent] = tmp

pp_spec = pd.read_csv("../../jetscape_data/max_time/maxT_200_highstat/pp_2760_charged_hadrons_eta_cut_1.csv") ## eta < 1, ptmin, ptmax, N, dN^2
pp_spec['pT']  = 0.5*(pp_spec['pTmax']+pp_spec['pTmin'])
pp = pp_spec[pp_spec['pT'].between(pT_lower_lim,pT_upper_lim)]

AA_spec = {}
for eloss in elosses:
    AA_spec[eloss] = {}
    for icent, cent in centralities.items():
        tmp = pd.read_csv(f"../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_charged_hadrons_eta_cut_1.csv",comment='#')
        tmp['pT']  = 0.5*(tmp['pTmax']+tmp['pTmin'])
        tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
        AA_spec[eloss][cent] = tmp



fig, axes = plt.subplots(1,1,sharex=True, sharey=True, figsize=(16,9))

for icent, cent in centralities.items():
    if cent != '00-05':
        continue
    alice = alice_experimental_data[cent]
    util.plot_expr_data_on_axis(axes, alice, marker=ddicts.markers['alice'], s=70)

    atlas2 = atlas2_experimental_data[cent]
    util.plot_expr_data_on_axis(axes, atlas2, marker=ddicts.markers['atlas2'], s=70)

    if cent in cms_experimental_data:
        cms = cms_experimental_data[cent]
        util.plot_expr_data_on_axis(axes, cms, marker=ddicts.markers['cms'], s=70) 
    if cent in atlas1_experimental_data:
        atlas = atlas1_experimental_data[cent]
        util.plot_expr_data_on_axis(axes, atlas, marker=ddicts.markers['atlas1'], s=70)
(ptmin, ptmax, proton, dproton) = util.combine_bins(pp['pTmin'].to_list(), pp['pTmax'].to_list(),
                                                    pp['N'].to_list(), pp['dN'].to_list())
pT = 0.5*(ptmax + ptmin)
for eloss in elosses:
    for icent, cent in centralities.items():
        if cent != '00-05':
            continue
        aa = AA_spec[eloss][cent]
        (_, _, PbPb, dPbPb) = util.combine_bins(aa['pTmin'].to_list(), \
                                                aa['pTmax'].to_list(),\
                                                aa['N'].to_list(),\
                                                aa['dN'].to_list())
        raa = PbPb/proton
        err = raa*np.sqrt(dPbPb*dPbPb/(PbPb*PbPb) + dproton*dproton/(proton*proton))
        axes.plot(pT, raa, color=module_colors['MATTER+'+eloss.upper()])
        axes.fill_between(pT, raa+err, raa-err, color=module_colors['MATTER+'+eloss.upper()], alpha=0.3)

axes.set_ylabel(r"$R^{h^{\pm}}_{\mathrm{AA}}$", fontsize=30)
axes.set_xlabel(r'$p_T$ (GeV)', fontsize=30) 
handles_experiment = [Line2D([],[],marker=ddicts.markers['alice'] ,markersize=10,label=r'ALICE $|\eta|<1.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['atlas2'],markersize=10,label=r'ATLAS $|\eta|<2.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['atlas1'],markersize=10,label=r'ATLAS $|\eta|<1.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['cms'],markersize=10, label=r'CMS $|\eta|<1.0$',color='black')]

handles_theory = [Line2D([],[],color=module_colors['MATTER+'+eloss.upper()],label='MATTER+'+eloss.upper()) for eloss in elosses]

artist = axes.legend(handles=handles_experiment,loc='lower right')
axes.add_artist(artist)
axes.legend(handles=handles_theory,loc='upper left')
axes.set_ylim(bottom=-0.01, top=1.71)
axes.text(0.1,0.7,r'Pb-Pb, $\sqrt{s}=2.76$ ATeV', fontsize=30, transform=axes.transAxes)
axes.text(0.1,0.6,r'$|\eta|<1.0$', transform=axes.transAxes)
axes.text(0.1,0.5,r'$0$-$5\%$', transform=axes.transAxes)
plt.show()
