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
import util 
import dictionaries as ddicts
from COLORS import module_colors
plt.rcParams.update(util.my_rcParams)
sigma = 62.049
experiment_colors = {'alice':"#7cae7a",
                     'atlas':"#d64045", 
                     'cms':"#eca400"}
cut_min, cut_max = 7, 440
exp_loc="../../../expt/PbPb_5p02/Charged/"

cms_data = {}
cms_centralities = {'00-05':8,'05-10':9,'10-30':10,'30-50':11}
for centrality in cms_centralities:
    table_num = cms_centralities[centrality]
    tmp = read_csv(exp_loc + f"CMS/Table{table_num}.csv", comment='#').rename(columns=ddicts.cols_CMS_chRAA)
    tmp = tmp[tmp['x'].between(cut_min, cut_max)]
    cms_data[centrality] = tmp

jetscapeLoc = "../../jetscape_data/sqrt_s_5020/maxt_200/{eloss}/PbPb5020_{centrality}_charged_hadrons_eta_cut_1.csv"
eloss_data = {}
for eloss in ['cujet','martini']:
    eloss_data[eloss] = {} 
    for cent in ['00-10','10-20','30-50']:
        tmp = read_csv(jetscapeLoc.format(eloss=eloss, centrality=cent), comment='#')
        tmp['pT'] = 0.5*(tmp['pTmin'] + tmp['pTmax'])
        tmp = tmp[tmp['pT'].between(cut_min, cut_max)]
        eloss_data[eloss][cent] = tmp

## pp base line:
#pp_spec = read_csv("../../jetscape_data/sqrt_s_5020/pp/pp_5020_charged_hadrons_eta_cut_1.csv", comment='#')
pp_spec = read_csv("../../jetscape_data/max_time/maxT_200_highstat/pp_5020_charged_hadrons_eta_cut_1.csv", comment='#')
pp_spec['pT'] = 0.5*(pp_spec['pTmin']+pp_spec['pTmax'])
pp_spec = pp_spec[pp_spec['pT'].between(cut_min, cut_max)]
#ppx,ppdx, ppy, ppdy = util.combine_bins(pp_spec['pTmin'], pp_spec['pTmax'], pp_spec['N'], pp_spec['dN'])

## ACTUAL PLOTTING
choice = {'00-05':0, '05-10':0,'10-30':1,'30-50':2}
centmarkers = markers = {'00-05':'P','05-10':'o','10-30':'*','30-50':'s'}

fig, axes = plt.subplots(3,1,sharex=True, sharey=True, figsize=(16,9))
exp_labels = {}

for cent in cms_data:
    axnum =  choice[cent]
    if axnum not in exp_labels:
        exp_labels[axnum] = []
    ax = axes[axnum]
    util.plot_expr_data_on_axis(ax, cms_data[cent], marker=centmarkers[cent],color='black',face='gray',s=90)
    exp_labels[axnum].append(Line2D([],[],marker=centmarkers[cent],color='black',markersize=10, label='CMS:'+ cent + '$\%$'))


## plot the theory
for eloss in eloss_data:
    color = module_colors['MATTER+'+eloss.upper()]
    for indx, cent in enumerate(['00-10','10-20','30-50']):
        ax = axes[indx]
        aa_spec = eloss_data[eloss][cent]
        #aax,aadx, aay, aady = util.combine_bins(aa_spec['pTmin'], aa_spec['pTmax'], aa_spec['N'], aa_spec['dN'])
        raa = aa_spec['N']/pp_spec['N']
        draa = raa*np.sqrt(aa_spec['dN']*aa_spec['dN']/(aa_spec['N']*aa_spec['N'])+\
                           pp_spec['dN']*pp_spec['dN']/(pp_spec['N']*pp_spec['N']) )
        ax.plot(aa_spec['pT'], raa, color=color)
        ax.fill_between(aa_spec['pT'], raa+draa, raa-draa, color=color, alpha=0.3)
        ax.text(0.4,0.05, cent + r'$\%$', transform=ax.transAxes)
        #raa = aay/ppy
        #draa = raa*np.sqrt(aady*aady/(aay*aay) + ppdy*ppdy/(ppy*ppy))
        #ax.plot(aax, raa, color=color)
        #ax.fill_between(aax, raa+draa, raa-draa, color=color, alpha=0.3)

theory_labels = [Line2D([],[],label='MATTER+'+eloss.upper(), color=module_colors['MATTER+'+eloss.upper()]) for eloss in ['cujet','martini']]
artist = axes[1].legend(loc='upper left', handles=theory_labels, ncols=2,fontsize=20)
axes[1].add_artist(artist)
for ax in axes:
    ax.set_ylabel(r'$R^{h^{\pm}}_{\mathrm{AA}}$')
axes[2].set_xlabel(r'$p_T$ (GeV)')

for iax, ax in enumerate(axes):
    ax.legend(loc='best', handles=exp_labels[iax], bbox_to_anchor=(0.98,0.8))
axes[0].text(0.01, 0.7, r'Pb-Pb $\sqrt{s}=5.02$ ATeV'+'\n'+'$|\eta|<1.0$',transform=axes[0].transAxes)

axes[0].set_xscale('log')
plt.show()
