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



pT_lower_lim = 60
pT_upper_lim = 310
inel_Xsec = 62.03948 

mother_path ='/Users/rmyazdi/Documents/research/jetscape_project/expt/PbPb_2p76/jets/'
jet_radii = {0:'0.2',1:'0.3',2:'0.4', 3:'0.6',4:'0.8'}
linestyles = {'0.3':'solid','0.4':'dashed','0.6':'dotted','0.8':'dashdot'}

cents_jetscape = {0:'00-05'}#,1:'05-10',2:'10-20',3:'20-30',4:'30-40',5:'40-50'}

## CMS data for R=0.2, 0.3, 0.4 (ALICE also has 10-30 for R=0.2)
cents_R_dependence = {0:'00-05'}#, 1:'05-10',2:'10-30',3:'30-50'}

## ATLAS
cents_medium_size  = {0:'00-10'}#,1:'10-20',2:'20-30',3:'30-40',4:'40-50'}

CMS_DATA = {}
for icent, cent in cents_R_dependence.items():
    CMS_DATA[cent] = {}
    for iR, radius in enumerate(['0.2','0.3','0.4']):
        R = radius.replace('.','p')
        tmp = pd.read_csv(f"{mother_path}CMS/jetscape/cent_{cent}_R_{R}.csv",comment='#').rename(columns=ddicts.colnames_CMS_RAA)
        tmp = tmp[tmp['xhigh'].between(pT_lower_lim,pT_upper_lim)]
        CMS_DATA[cent][radius] = tmp

jetscape = {}
elosses = ['matter']
matter_options = [20,200]
for eloss in elosses:
    jetscape[eloss] = {}
    for opt in matter_options:
        jetscape[eloss][opt] = {}
        for ir, radius in jet_radii.items():
            jetscape[eloss][opt][radius] = {}
            for icent, cent in cents_jetscape.items():
                tmp = pd.read_csv(f"../../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760_maxT_{opt}/PbPb2760_{cent}_jet_spec_jet_rad_{radius}_2.00.csv", comment='#')
                tmp['pT'] = 0.5*(tmp['ptmax']+tmp['ptmin'])
                tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
                jetscape[eloss][opt][radius][cent] = tmp

for eloss in ['martini', 'cujet']:
    jetscape[eloss] = {}
    for ir, radius in jet_radii.items():
        jetscape[eloss][radius] = {}
        for icent, cent in cents_jetscape.items():
            tmp = pd.read_csv(f"../../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_jet_spec_jet_rad_{radius}_2.00.csv", comment='#')
            tmp['pT'] = 0.5*(tmp['ptmax']+tmp['ptmin'])
            tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
            jetscape[eloss][radius][cent] = tmp
fname = '../../../jetscape_data/sqrt_s_2760/test_MaxT_200/PbPb_2760/PbPb2760_{cent}_jet_spec_jet_rad_{radius}_2.00.csv'
tmpeloss = 'martini+matter'
jetscape[tmpeloss] = {}
for ir, radius in jet_radii.items():
    jetscape[tmpeloss][radius] = {}
    for icent, cent in cents_jetscape.items():
        tmp = pd.read_csv(fname.format(cent=cent,radius=radius, comment='#'))
        tmp['pT'] = 0.5*(tmp['ptmax']+tmp['ptmin'])
        tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
        jetscape[tmpeloss][radius][cent] = tmp
        
pp_baseline =  {}
for ir, radius in jet_radii.items():
    tmp = pd.read_csv(f'../../../jetscape_data/max_time/maxT_200_highstat/pp_2760_jet_spec_jet_rad_{radius}_2.00.csv',comment='#')
    tmp['pT'] = 0.5 * (tmp['ptmin']+tmp['ptmax'])
    tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
    pp_baseline[radius] = tmp

fig, axes = plt.subplots(1,4,gridspec_kw={'top':0.985,'bottom':0.12,
                            'left':0.11,'right':0.95,
                            'hspace':0.15,'wspace':0.1},
                           sharex='col', sharey='row', figsize=(20,10))

axes = axes.flatten()
cent = '00-05'
eloss = 'matter'
for ir, radius in enumerate(['0.2','0.3','0.4']):
    pp = pp_baseline[radius]
    util.plot_expr_data_on_axis(axes[ir],CMS_DATA[cent][radius],marker=ddicts.markers_experiment['cms'])
    for opt in matter_options:
        AA = jetscape[eloss][opt][radius][cent]
        raa = AA['Ncut']/pp['Ncut']
        draa = raa*np.sqrt(AA['dNcut']*AA['dNcut']/(AA['Ncut']*AA['Ncut']) + pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']) ) 
        color = ddicts.eloss_colours[f'{eloss}_{opt}']
        axes[ir].plot(AA['pT'], raa, color=color)
        axes[ir].fill_between(AA['pT'], raa+draa,raa-draa,color=color,alpha=0.2)
for ir, radius in enumerate(['0.2','0.3','0.4']):
    pp = pp_baseline[radius]
    for eloss in ['martini', 'cujet']:
        AA = jetscape[eloss][radius][cent]
        raa = AA['Ncut']/pp['Ncut']
        draa = raa*np.sqrt(AA['dNcut']*AA['dNcut']/(AA['Ncut']*AA['Ncut']) + pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']) ) 
        color = ddicts.eloss_colours[f'{eloss}']
        axes[ir].plot(AA['pT'], raa, color=color)
        axes[ir].fill_between(AA['pT'], raa+draa,raa-draa,color=color,alpha=0.2)
    eloss = 'maartini+matter'

for i in range(3):
    axes[i].text(0.2, 0.1, r'$R=$'+f"{jet_radii[i]}",transform=axes[i].transAxes)
    axes[i].set_xlabel(r'$p_T$ (GeV)')

label  = r"$R^{\mathrm{jet}}_{\mathrm{AA}}$"
fig.supylabel(label, fontsize=30)
#fig.supxlabel(r'$p_T$ (GeV)', fontsize=30)

handles = [Line2D([],[],label='MATTER (maxT = 20)', color=ddicts.eloss_colours['matter_20'])]
handles.append(Line2D([],[],label='MATTER (maxT = 200)', color=ddicts.eloss_colours['matter_200']))
handles.append(Line2D([],[],label='MARTINI+MATTER', color=ddicts.eloss_colours['martini']))
handles.append(Line2D([],[],label='CUJET+MATTER', color=ddicts.eloss_colours['cujet']))
handle_exp = [Line2D([],[],marker=ddicts.markers_experiment['cms'],color='black',label='CMS $|\eta|<2$ (2017)')]


axes[3].legend(handles=handle_exp, loc='lower right',fontsize=18)
axes[3].text(0.1,0.15,r'Pb-Pb, $\sqrt{s}=2.76$ ATeV' + "\n" + "$|\eta|<2$" + "\n" + r"$0-5\%$", fontsize=20, transform=axes[3].transAxes)
axes[3].legend(handles=handles, bbox_to_anchor=(-0.09,0.75), loc='upper left',prop={'size': 18})
axes[3].axis('off')
plt.show()
