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

plt.rcParams.update(util.my_rcParams)

elosses = ['cujet', 'martini']

pT_lower_lim = 8
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



fig, axes = plt.subplots(2,3,sharex=True, sharey=True, figsize=(16,9))
axes = axes.flatten()

for icent, cent in centralities.items():
    alice = alice_experimental_data[cent]
    util.plot_expr_data_on_axis(axes[icent], alice, marker=ddicts.markers['alice'], s=50)

    atlas2 = atlas2_experimental_data[cent]
    util.plot_expr_data_on_axis(axes[icent], atlas2, marker=ddicts.markers['atlas2'], s=50)

    if cent in cms_experimental_data:
        cms = cms_experimental_data[cent]
        util.plot_expr_data_on_axis(axes[icent], cms, marker=ddicts.markers['cms'], s=50) 
    if cent == '40-50':
        cms = cms_experimental_data['30-50']
        util.plot_expr_data_on_axis(axes[icent], cms, marker=ddicts.markers['cms'], s=50) 
    if cent in atlas1_experimental_data:
        atlas = atlas1_experimental_data[cent]
        util.plot_expr_data_on_axis(axes[icent], atlas, marker=ddicts.markers['atlas1'], s=50) 
(ptmin, ptmax, proton, dproton) = util.combine_bins(pp['pTmin'].to_list(), pp['pTmax'].to_list(),
                                                    pp['N'].to_list(), pp['dN'].to_list())
pT = 0.5*(ptmax + ptmin)
#pT = pp['pT']#0.5*(pp['pTmin'] + pp['pTmax'])
for eloss in elosses:
    for icent, cent in centralities.items():
        aa = AA_spec[eloss][cent]
        (_, _, PbPb, dPbPb) = util.combine_bins(aa['pTmin'].to_list(), \
                                                aa['pTmax'].to_list(),\
                                                aa['N'].to_list(),\
                                                aa['dN'].to_list())

        #raa = aa['N']/pp['N']
        #p, a = pp['N'], aa['N']
        #dp, da = pp['dN'], aa['dN']
        #err = raa*np.sqrt(dp*dp/(p*p) + da*da/(a*a))

        raa = PbPb/proton
        err = raa*np.sqrt(dPbPb*dPbPb/(PbPb*PbPb) + dproton*dproton/(proton*proton))
        axes[icent].plot(pT, raa, color=module_colors['MATTER+'+eloss.upper()])
        axes[icent].fill_between(pT, raa+err, raa-err, color=module_colors['MATTER+'+eloss.upper()], alpha=0.3)
        axes[icent].text(0.65, 0.1, f'{ddicts.nice_cent_label[cent]}', transform=axes[icent].transAxes)
        if cent == '40-50':
            axes[icent].text(0.1,0.1,'CMS: $30-50\%$', transform=axes[icent].transAxes, fontsize=25)

for (iax, ax) in enumerate(axes):
    if iax%3 == 0:
        ax.set_ylabel(r"$R^{h^{\pm}}_{\mathrm{AA}}$", fontsize=30)
    if iax >= 3:
        ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30) 
handles_experiment = [Line2D([],[],marker=ddicts.markers['alice'] ,markersize=10,label=r'ALICE $|\eta|<1.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['atlas2'],markersize=10,label=r'ATLAS $|\eta|<2.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['atlas1'],markersize=10,label=r'ATLAS $|\eta|<1.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['cms'],markersize=10, label=r'CMS $|\eta|<1.0$',color='black')]

handles_theory = [Line2D([],[],color=module_colors['MATTER+'+eloss.upper()],label='MATTER+'+eloss.upper()) for eloss in elosses]

axes[2].legend(handles=handles_experiment,loc='best')
axes[0].legend(handles=handles_theory,loc='upper left')
axes[0].set_ylim(bottom=-0.01, top=2.01)
axes[0].text(0.1,0.6,r'$|\eta|<1.0$', transform=axes[0].transAxes)
axes[1].text(0.07,0.8,r'Pb-Pb, $\sqrt{s}=2.76$ ATeV', fontsize=30, transform=axes[1].transAxes)
plt.show()
