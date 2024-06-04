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
import jetDicts as ddicts

my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 30,
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
    "axes.spines.top"  : False,
    "legend.frameon"   : False
}
plt.rcParams.update(my_rcParams)

def construct_new_centrality(df1, df2, xname):
    xmin = df1[f'{xname}min']
    xmax = df1[f'{xname}max']
    ncut = 0.5*(df1['Ncut'] + df2['Ncut'])
    dncut = np.sqrt(df1['dNcut']*df1['dNcut'] + df2['dNcut']*df2['dNcut'])
    wcut = 0.5*(df1['Wcut'] + df2['Wcut'])
    dwcut = np.sqrt(df1['dWcut']*df1['dWcut'] + df2['dWcut']*df2['dWcut'])
    return pd.DataFrame({f'{xname}min':xmin, f'{xname}max':xmax,'Ncut':ncut,'dNcut':dncut,'Wcut':wcut,'dWcut':dwcut})

centralities = {0:'00-10',1:'20-30',2:'30-40'}

fnames_DPT = {'00-10':'Table9','20-30':'Table10','30-40':'Table11'}
fnames_DZ = {'00-10':'Table25','20-30':'Table26','30-40':'Table27'}

expt_loc = "../../../expt/PbPb_2p76/jet_fragmentation/ATLAS/HEPData-ins1511869-v1-csv/"
experiment_pT = {cent:pd.read_csv(expt_loc+f"{fname}.csv",comment='#').rename(columns=ddicts.col_names_FF_ATLAS_DPT) for cent, fname in fnames_DPT.items()}
experiment_z = {cent:pd.read_csv(expt_loc+f"{fname}.csv",comment='#').rename(columns=ddicts.col_names_FF_ATLAS_DZ) for cent, fname in fnames_DZ.items()}

experiments = {0:experiment_pT}
from COLORS import module_colors
maxpT = 120
minpT = 0.1
zmin, zmax = 0.005,2
## Read in the pp data:
pp_DPT = pd.read_csv('../../jetscape_data/max_time/maxT_200_highstat/pp_2760_LHC_FF_pT_jet_rad_0.4_2.10.csv',comment='#')
pp_DPT = pp_DPT[pp_DPT['pTmax'] < maxpT]
pp_DPT = pp_DPT[pp_DPT['pTmin'] > minpT]
pp_jetspec = pd.read_csv('../../jetscape_data/max_time/maxT_200_highstat/pp_2760_LHC_jet_spec_FF_jet_rad_0.4_2.10.csv',comment='#')
njets_pp =  sum(pp_jetspec['Ncut'].tolist())

## read in the AA data:
elosses = ['martini','cujet']
AA_DPT = {}
AA_Spec = {}
AA_centralities = ['00-05','05-10','20-30','30-40']
AA_LOCATION_FF = "../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_LHC_FF_{obs}_jet_rad_0.4_2.10.csv"
AA_LOCATION_SP = "../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_LHC_jet_spec_FF_jet_rad_0.4_2.10.csv" 
for eloss in elosses:
    AA_DPT[eloss] = {}
    AA_Spec[eloss] = {}
    for cent in AA_centralities:
        tmp = pd.read_csv(AA_LOCATION_FF.format(eloss=eloss,cent=cent,obs='pT'),comment='#')
        tmp = tmp[tmp['pTmax'] < maxpT] 
        tmp = tmp[tmp['pTmin'] > minpT]
        AA_DPT[eloss][cent] = tmp
        AA_Spec[eloss][cent] = pd.read_csv(AA_LOCATION_SP.format(eloss=eloss,cent=cent),comment='#')

    new_cent_DPT = construct_new_centrality(AA_DPT[eloss]['00-05'],AA_DPT[eloss]['05-10'],'pT')
    new_cent_spec = construct_new_centrality(AA_Spec[eloss]['00-05'],AA_Spec[eloss]['05-10'],'pT')
    AA_DPT[eloss]['00-10']=new_cent_DPT
    AA_Spec[eloss]['00-10']=new_cent_spec

observables_aa = {0:AA_DPT}
observables_pp = {0:pp_DPT}
x_axis = {0:'pT'}
fig, axes = plt.subplots(1, 3, figsize=(16,9), sharey=True,sharex='row')

for experiment in experiments:
    data = experiments[experiment]
    pp = observables_pp[experiment]
    for eloss in elosses:
        for icent, cent in centralities.items():
            ax = axes[icent]
            centrality_data = data[cent]
            util.plot_expr_data_on_axis(ax, centrality_data, marker='*')
            color = module_colors['MATTER+'+eloss.upper()]
            aa = observables_aa[experiment][eloss][cent]
            x = 0.5*(aa[f'{x_axis[experiment]}max']+aa[f'{x_axis[experiment]}min'])
            dx = aa[f'{x_axis[experiment]}max']-aa[f'{x_axis[experiment]}min']
            njets_AA = sum(AA_Spec[eloss][cent]['Ncut'].tolist())
            ratio = (aa['Ncut']/pp['Ncut']) * (njets_pp/njets_AA)
            dratio = ratio*np.sqrt(aa['dNcut']*aa['dNcut']/(aa['Ncut']*aa['Ncut']) + pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']))
            errorboxes = [Rectangle((x-delx, y - dy), width=2*delx, height=2*dy, zorder=2)
                    for x, delx, y, dy in zip(x, 0.5*dx, ratio, dratio)]
            pc = PatchCollection(errorboxes, facecolor=color, edgecolor=color, alpha=0.6)
            ax.add_collection(pc)
            ax.scatter(x, ratio, color=color,marker=ddicts.eloss_marker[eloss], zorder=2)

labels = [Line2D([],[],label=r'ATLAS $|y|<2.1$ (2017)', marker='*',color='black')]
for eloss in elosses:
    labels.append(Line2D([],[],label='MATTER+'+eloss.upper(), color=module_colors['MATTER+'+eloss.upper()],marker=ddicts.eloss_marker[eloss]))

for icent, cent in centralities.items():
    axes[icent].text(0.05,0.1, f"{cent}$\%$", transform=axes[icent].transAxes)
    axes[icent].text(0.05,0.1, f"{cent}$\%$", transform=axes[icent].transAxes)
axes[0].set_xscale('log')
axes[0].set_xticks([1,5,10,20,50,100])
axes[0].set_xticklabels(['1','5', '10','20','50','100'])
axes[1].set_xlabel(r'$p^{h^{\pm}}_{T}$ (GeV)',fontsize=25)
axes[0].set_ylabel(r'$R_{D\left(p_T\right)}$', fontsize=25)
axes[2].legend(loc='upper left', handles=labels, fontsize=30)

#axes[1][0].set_xscale('log')
#axes[1][0].set_xticks([0.01,0.04,0.1,0.4,1])
#axes[1][0].set_xticklabels(['0.01','0.04','0.1','0.4','1'])
for ax in axes:
    ax.set_xlabel(r'$p^{h^{\pm}}_T$ (GeV)',fontsize=30)
axes[0].text(0.2,0.9,s=r'$p^{j}_{T} \in (100,398)$ GeV',transform=axes[0].transAxes)
axes[0].text(0.2,0.8,s=r'$|\eta|<2$',transform=axes[0].transAxes)
axes[0].set_ylim(top=2.)
plt.show()
