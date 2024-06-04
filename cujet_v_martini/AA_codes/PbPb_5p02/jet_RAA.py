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
exp_loc="../../../expt/PbPb_5p02/Jets/"

atlas_data = {}
ATLAS_centralities = {'00-10':19,'10-20':20,'30-40':22,'40-50':23}
ATLAS_axes = {'00-10':0, '10-20':1, '30-40':2, '40-50':2}
ATLAS_markers = {'00-10':'P', '10-20':'o', '30-40':'<', '40-50':'>'}
for centrality in ATLAS_centralities:
    table_num = ATLAS_centralities[centrality]
    tmp = read_csv(exp_loc + f"ATLAS/Table{table_num}.csv", comment='#').rename(columns=ddicts.colnames_ATLAS_RAA)
    atlas_data[centrality] = tmp

cms_data = {}
CMS_axes = {'00-10':0, '10-30':1, '30-50':2}
CMS_markers = {'00-10':'s', '10-30':'*', '30-50':'D'}
for centrality in CMS_markers:
    cms_data[centrality] = {}
    for R in ['0p2','0p3','0p4','0p6', '0p8']:
        fname = exp_loc + f"CMS/AA_R_{R}_{centrality}.csv"
        tmp = read_csv(fname, comment='#').rename(columns=ddicts.colnames_CMS_RAA)
        tmp = tmp[tmp['y'] > 0]
        cms_data[centrality][R]=tmp


## pp baseline:
maxcuts = {'00-10':850,'10-20':850, '30-50':600}
pp_specs = {}
for R in ['0p2','0p3','0p4','0p6', '0p8']:
    fname = "../../jetscape_data/max_time/maxT_200_highstat/pp_5020_jet_spec_jet_rad_{r}_2.00.csv"
    r = R.replace('p','.')
    tmp = read_csv(fname.format(r=r), comment='#')
    tmp['pT'] = 0.5*(tmp['ptmin'] + tmp['ptmax'])
    tmp = tmp[tmp['pT'] > 50 ]
    pp_specs[R] = tmp


AA_specs = {}
fname_template = '../../jetscape_data/sqrt_s_5020/maxt_200/{eloss}/PbPb5020_{cent}_jet_spec_jet_rad_{r}_2.00.csv'
for eloss in ['cujet','martini']:
    AA_specs[eloss] = {}
    for cent in ['00-10','10-20','30-50']:
        AA_specs[eloss][cent] = {}
        for R in ['0p2','0p3','0p4','0p6', '0p8']:
            r = R.replace('p','.')
            tmp = read_csv(fname_template.format(eloss=eloss, cent=cent, r=r),comment='#')
            tmp['pT'] = 0.5*(tmp['ptmin'] + tmp['ptmax'])
            AA_specs[eloss][cent][R] = tmp

## figure for R = 0.4
fig, axes = plt.subplots(3,1,sharex=True, sharey=True, figsize=(16,9))
labels_exp = {}
for cent, axIndx in ATLAS_axes.items():
    ax = axes[axIndx]
    data = atlas_data[cent]
    marker = ATLAS_markers[cent]
    util.plot_expr_data_on_axis(ax, data, marker,'black','gray', s=40)
    if axIndx not in labels_exp:
        labels_exp[axIndx] = []
    labels_exp[axIndx].append(Line2D([],[],color='black',marker=marker,markersize=10,label=r'ATLAS ' + cent+r"$\%$"))

r = '0p4'
for cent, axIndx in CMS_axes.items():
    ax = axes[axIndx]
    marker = CMS_markers[cent]
    data = cms_data[cent][r]
    util.plot_expr_data_on_axis(ax, data, marker,'black','gray', s=40)
    labels_exp[axIndx].append(Line2D([],[],color='black',marker=marker,markersize=10,label=r'CMS ' + cent+r"$\%$"))

for eloss in AA_specs:
    color = module_colors['MATTER+'+eloss.upper()] 
    for indx, cent in enumerate(AA_specs[eloss]):
        ax = axes[indx]
        pp = pp_specs[r]
        pp = pp[pp['pT'].between(50, maxcuts[cent])]
        aa = AA_specs[eloss][cent][r]
        aa = aa[aa['pT'].between(50, maxcuts[cent])]
        raa = aa['Ncut']/pp['Ncut']
        draa = raa*np.sqrt(aa['dNcut']*aa['dNcut']/(aa['Ncut']*aa['Ncut']) + \
                           pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']))

        ax.plot(pp['pT'], raa, color=color, zorder=0)
        ax.fill_between(pp['pT'], raa+draa, raa-draa, color=color, alpha=0.3, zorder=0)
        ax.text(0.1, 0.8, f'{cent}'+ r'$\%$', transform=ax.transAxes)

theory_labels = [Line2D([],[],label='MATTER+'+eloss.upper(), color=module_colors['MATTER+'+eloss.upper()]) for eloss in ['cujet','martini']]
artist = axes[2].legend(handles=theory_labels, loc='lower right', ncols=1)

axes[2].add_artist(artist)
axes[0].set_ylim(bottom=-0.01, top=1.01)

for iax, ax in enumerate(axes):
    ax.legend(handles=labels_exp[iax], bbox_to_anchor=(1.,0.7), loc='upper left')

axes[2].set_xlabel(r'$p_T$ (GeV)', fontsize=35)
for ax in axes:
    ax.set_ylabel(r'$R^{\mathrm{jet}}_{\mathrm{AA}}$')

axes[0].text(0.05, 0.1, r'Pb-Pb $\sqrt{s}=5.02$ ATeV, $|\eta^{\mathrm{jet}}|<2.0$', transform=axes[0].transAxes)
axes[2].text(0.05, 0.1, r'CMS: $|\eta^{\mathrm{jet}}|<2.0$' + '\n'+r'ATLAS: $|y^{\mathrm{jet}}|<2.8$', transform=axes[2].transAxes)
axes[2].text(0.35, 0.1, r'Anti-$k_T$, R=0.4', transform=axes[2].transAxes)
axes[0].set_xscale('log')
## Now compute all the RAAs:
raa_calcs = {}
for eloss in ['martini', 'cujet']:
    raa_calcs[eloss] = {}
    for cent in ['00-10','10-20','30-50']:
        raa_calcs[eloss][cent] = {}
        for R in ['0p2','0p3','0p4','0p6', '0p8']:
            aa = AA_specs[eloss][cent][R]
            aa = aa[aa['pT'].between(200, 800)]
            pp = pp_specs[R]
            pp = pp[pp['pT'].between(200, 800)]
            (xa, ay, day) = aa['pT'],aa['Ncut'],aa['dNcut']
            (xp, py, dpy) = pp['pT'],pp['Ncut'],pp['dNcut']
            raa = ay/py
            draa = raa*np.sqrt(day*day/(ay*ay) + dpy*dpy/(py*py))
            raa_calcs[eloss][cent][R] = (xa, raa, draa)

fig1, axes = plt.subplots(3,5,sharex=True, sharey=True, figsize=(16,9))
#axes[0][0].set_xscale('log')
markers = {0:'*', 1:'P', 2:'D'}
cms_handles = {0:[]}
for eloss in ['cujet','martini']:
    cms_handles[0].append(Line2D([],[],color=module_colors['MATTER+'+eloss.upper()], label='MATTER+'+eloss.upper()))
for cent, icent in CMS_axes.items():
    marker = markers[icent]
    if icent not in cms_handles:
        cms_handles[icent] = []
    cms_handles[icent].append(Line2D([],[],marker=marker,label=cent+r'$\%$',color='black',markersize=10))
    for iR, R in enumerate(['0p2','0p3','0p4','0p6', '0p8']):
        ax = axes[icent][iR]
        data = cms_data[cent][R]
        util.plot_expr_data_on_axis(ax, data, marker, 'black','gray',s=50)

for icent, cent in enumerate(['00-10','10-20','30-50']):
    for iR, R in enumerate(['0p2','0p3','0p4','0p6', '0p8']):
        ax = axes[icent][iR]
        for eloss in ['martini','cujet']:
            color = module_colors['MATTER+'+eloss.upper()]
            x, raa, draa = raa_calcs[eloss][cent][R]
            ax.plot(x, raa, color=color)
            ax.fill_between(x, raa-draa, raa+draa, color=color, alpha=0.2)

first_row = axes[0]
for iR, R in enumerate(['0p2','0p3','0p4','0p6', '0p8']):
    ax = first_row[iR]
    r = R.replace('p','.')
    txt = r'$R=$'+f'{r}'
    ax.text(0.3,0.8,txt,transform=ax.transAxes)

for icent, cent in enumerate(['00-10','10-20','30-50']):
    axes[icent][0].text(0.05, 0.1, cent + r"$\%$", transform=axes[icent][0].transAxes)
    axes[icent][4].legend(loc='upper left', bbox_to_anchor=(0.9,0.8), handles=cms_handles[icent], fontsize=20)

axes[1][4].text(0.4,0.08,r'$|\eta|<2.0$', transform=axes[1][4].transAxes)
axes[2][4].text(0.4,0.08,r'Anti-$k_{T}$', transform=axes[2][4].transAxes)
fig1.suptitle(r'Pb-Pb $\sqrt{s}=5.02$ ATeV')
for ax in [ax[0] for ax in axes]:
    ax.set_ylabel(r'$R^{\mathrm{jet}}_{\mathrm{AA}}$')
for ax in axes[-1]:
    ax.set_xlabel(r'$p_T$ (GeV)')
plt.show()


## Now plot the ratio of RAA's and compare to the available CMS data

