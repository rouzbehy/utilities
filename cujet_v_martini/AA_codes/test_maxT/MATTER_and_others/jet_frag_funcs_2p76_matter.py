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

centralities = {0:'00-05'}

maxpT = 120
minpT = 0.1
zmin, zmax = 0.005,2
## Read in the pp data:
pp_DPT = pd.read_csv('../../../jetscape_data/max_time/maxT_200_highstat/pp_2760_LHC_FF_pT_jet_rad_0.4_2.10.csv',comment='#')
pp_DPT = pp_DPT[pp_DPT['pTmax'] < maxpT]
pp_DPT = pp_DPT[pp_DPT['pTmin'] > minpT]
pp_DZ = pd.read_csv('../../../jetscape_data/max_time/maxT_200_highstat/pp_2760_LHC_FF_z_jet_rad_0.4_2.10.csv',comment='#')
pp_DZ = pp_DZ[pp_DZ['zmax'] < zmax]
pp_DZ = pp_DZ[pp_DZ['zmin'] > zmin]
pp_jetspec = pd.read_csv('../../../jetscape_data/max_time/maxT_200_highstat/pp_2760_LHC_jet_spec_FF_jet_rad_0.4_2.10.csv',comment='#')
njets_pp =  sum(pp_jetspec['Ncut'].tolist())

## read in the AA data:
elosses = ['martini','cujet']
AA_DPT = {}
AA_DZ = {}
AA_Spec = {}
AA_centralities = ['00-05']
AA_LOCATION_FF = "../../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_LHC_FF_{obs}_jet_rad_0.4_2.10.csv"
AA_LOCATION_SP = "../../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_LHC_jet_spec_FF_jet_rad_0.4_2.10.csv" 
for eloss in elosses:
    AA_DPT[eloss] = {}
    AA_DZ[eloss] = {}
    AA_Spec[eloss] = {}
    for cent in AA_centralities:
        tmp = pd.read_csv(AA_LOCATION_FF.format(eloss=eloss,cent=cent,obs='pT'),comment='#')
        tmp = tmp[tmp['pTmax'] < maxpT] 
        tmp = tmp[tmp['pTmin'] > minpT]
        AA_DPT[eloss][cent] = tmp
        tmp1 = pd.read_csv(AA_LOCATION_FF.format( eloss=eloss,cent=cent,obs='z'),comment='#')
        tmp1 = tmp1[tmp1['zmax'] < zmax] 
        tmp1 = tmp1[tmp1['zmin'] > zmin]
        AA_DZ[eloss][cent] = tmp1
        AA_Spec[eloss][cent] = pd.read_csv(AA_LOCATION_SP.format(eloss=eloss,cent=cent),comment='#')

matter_fname_FF = "../../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760_maxT_{maxT}/PbPb2760_{cent}_LHC_FF_{obs}_jet_rad_0.4_2.10.csv"
matter_fname_SP = "../../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760_maxT_{maxT}/PbPb2760_{cent}_LHC_jet_spec_FF_jet_rad_0.4_2.10.csv" 
eloss = 'matter'
options = [20, 200]
AA_DPT[eloss] = {}
AA_DZ[eloss] = {}
AA_Spec[eloss] = {}
for opt in options:
    AA_DPT[eloss][opt] = {}
    AA_DZ[eloss][opt] = {}
    AA_Spec[eloss][opt] = {}
    for cent in AA_centralities:
        tmp = pd.read_csv(matter_fname_FF.format(eloss=eloss,cent=cent,obs='pT',maxT=opt),comment='#')
        tmp = tmp[tmp['pTmax'] < maxpT] 
        tmp = tmp[tmp['pTmin'] > minpT]
        AA_DPT[eloss][opt][cent] = tmp
        tmp1 = pd.read_csv(matter_fname_FF.format(eloss=eloss,cent=cent,obs='z',maxT=opt),comment='#')
        tmp1 = tmp1[tmp1['zmax'] < zmax] 
        tmp1 = tmp1[tmp1['zmin'] > zmin]
        AA_DZ[eloss][opt][cent] = tmp1
        AA_Spec[eloss][opt][cent] = pd.read_csv(matter_fname_SP.format(eloss=eloss,cent=cent,maxT=opt),comment='#')

observables_aa = {0:AA_DPT,1:AA_DZ}
observables_pp = {0:pp_DPT,1:pp_DZ}
x_axis = {0:'pT',1:'z'}

fig, axes = plt.subplots(1, 2, figsize=(16,9), 
                    gridspec_kw={'top':0.958,'bottom':0.126,'left':0.097,
                                 'right':0.977,'hspace':0.35,'wspace':0.15},
                    sharey=True,sharex=False)

axes = axes.flatten()
for i in [0,1]:
    pp = observables_pp[i]
    for eloss in ['martini','cujet']:
        for icent, cent in centralities.items():
            ax = axes[i]
            aa = observables_aa[i][eloss][cent]
            x = 0.5*(aa[f'{x_axis[i]}max']+aa[f'{x_axis[i]}min'])
            dx = aa[f'{x_axis[i]}max']-aa[f'{x_axis[i]}min']
            njets_AA = sum(AA_Spec[eloss][cent]['Ncut'].tolist())
            ratio = (aa['Ncut']/pp['Ncut']) * (njets_pp/njets_AA)
            dratio = ratio*np.sqrt(aa['dNcut']*aa['dNcut']/(aa['Ncut']*aa['Ncut']) + pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']))
            #errorboxes = [Rectangle((x-delx, y - dy), width=2*delx, height=2*dy)
            #        for x, delx, y, dy in zip(x, 0.5*dx, ratio, dratio)]
            #pc = PatchCollection(errorboxes, facecolor=ddicts.eloss_colours[eloss], edgecolor=ddicts.eloss_colours[eloss], alpha=0.5)
            #ax.add_collection(pc)
            #ax.scatter(x, ratio, color=ddicts.eloss_colours[eloss],marker=ddicts.eloss_marker[eloss])
            ax.plot(x, ratio, color=ddicts.eloss_colours[eloss])
            ax.fill_between(x, ratio+dratio, ratio-dratio, color=ddicts.eloss_colours[eloss], alpha=0.3)
    eloss = 'matter'
    for opt in options:
        color = ddicts.eloss_colours[f'{eloss}_{opt}']
        for icent, cent in centralities.items():
            ax = axes[i]
            aa = observables_aa[i][eloss][opt][cent]
            x = 0.5*(aa[f'{x_axis[i]}max']+aa[f'{x_axis[i]}min'])
            dx = aa[f'{x_axis[i]}max']-aa[f'{x_axis[i]}min']
            njets_AA = sum(AA_Spec[eloss][opt][cent]['Ncut'].tolist())
            ratio = (aa['Ncut']/pp['Ncut']) * (njets_pp/njets_AA)
            dratio = ratio*np.sqrt(aa['dNcut']*aa['dNcut']/(aa['Ncut']*aa['Ncut']) + pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']))
            #errorboxes = [Rectangle((x-delx, y - dy), width=2*delx, height=2*dy)
            #        for x, delx, y, dy in zip(x, 0.5*dx, ratio, dratio)]
            #pc = PatchCollection(errorboxes, facecolor=color, edgecolor=color, alpha=0.5)
            #ax.add_collection(pc)
            #ax.scatter(x, ratio, color=color,marker=ddicts.eloss_marker[eloss])

            ax.plot(x, ratio, color=color)
            ax.fill_between(x, ratio+dratio, ratio-dratio, color=color, alpha=0.3)
labels = []
for eloss in elosses:
    labels.append(Line2D([],[],label=eloss.upper(), color=ddicts.eloss_colours[eloss],marker=ddicts.eloss_marker[eloss]))
for opt in options:
    labels.append(Line2D([],[],label=f'MATTER (maxT={opt})',color=ddicts.eloss_colours[f'matter_{opt}'],marker=ddicts.eloss_marker[f'matter']))

axes[0].set_xscale('log')
axes[0].set_xticks([1,5,10,20,50,100])
axes[0].set_xticklabels(['1','5', '10','20','50','100'])
axes[0].set_xlabel(r'$p^{h^{\pm}}_{T}$ (GeV)',fontsize=30)
axes[0].set_ylabel(r'$R_{D\left(p_T\right)}$')
axes[0].set_ylim(top=1.5,bottom=0.2)

axes[1].set_xscale('log')
axes[1].set_xticks([0.01,0.04,0.1,0.4,1])
axes[1].set_xticklabels(['0.01','0.04','0.1','0.4','1'])
axes[1].set_xlabel(r'$z$',fontsize=30)
axes[1].set_ylabel(r'$R_{D\left(z\right)}$')
axes[1].legend(loc='upper right', handles=labels, fontsize=18)

axes[0].text(0.1,0.2,s=r'$100$ GeV $< p^{j}_{T}< 398$ GeV',transform=axes[0].transAxes)
axes[0].text(0.1,0.15,s=r'$|\eta|<2$',transform=axes[0].transAxes)
axes[0].text(0.1,0.1,s=r'Anti-$k_{T}$', transform=axes[0].transAxes)
axes[0].text(0.1,0.05,s=r'$0-5\%$', transform=axes[0].transAxes)

plt.show()
