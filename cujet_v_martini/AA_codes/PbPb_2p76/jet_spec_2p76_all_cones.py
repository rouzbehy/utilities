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

from COLORS import module_colors
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 20,
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

elosses = ['cujet', 'martini']

pT_lower_lim = 60
pT_upper_lim = 310
inel_Xsec = 62.03948 

mother_path ='/Users/rmyazdi/Documents/research/jetscape_project/expt/PbPb_2p76/jets/'
jet_radii = {0:'0.2',1:'0.3',2:'0.4', 3:'0.6',4:'0.8'}
linestyles = {'0.3':'solid','0.4':'dashed','0.6':'dotted','0.8':'dashdot'}

cents_jetscape = {0:'00-05',1:'05-10',2:'10-20',3:'20-30',4:'30-40',5:'40-50'}

## CMS data for R=0.2, 0.3, 0.4 (ALICE also has 10-30 for R=0.2)
cents_R_dependence = {0:'00-05', 1:'05-10',2:'10-30',3:'30-50'}

## ATLAS
cents_medium_size  = {0:'00-10',1:'10-20',2:'20-30',3:'30-40',4:'40-50'}

CMS_DATA = {}
for icent, cent in cents_R_dependence.items():
    CMS_DATA[cent] = {}
    for iR, radius in enumerate(['0.2','0.3','0.4']):
        R = radius.replace('.','p')
        tmp = pd.read_csv(f"{mother_path}CMS/jetscape/cent_{cent}_R_{R}.csv",comment='#').rename(columns=ddicts.colnames_CMS_RAA)
        tmp = tmp[tmp['xhigh'].between(pT_lower_lim,pT_upper_lim)]
        CMS_DATA[cent][radius] = tmp

jetscape = {}
elosses = ['cujet', 'martini']
for eloss in elosses:
    jetscape[eloss] = {}
    for ir, radius in jet_radii.items():
        jetscape[eloss][radius] = {}
        for icent, cent in cents_jetscape.items():
            tmp = pd.read_csv(f"../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_jet_spec_jet_rad_{radius}_2.00.csv", comment='#')
            tmp['pT'] = 0.5*(tmp['ptmax']+tmp['ptmin'])
            #print(tmp['pT'])
            #tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
            jetscape[eloss][radius][cent] = tmp
        ## compute 10-30 and 30-50 centralities here:
        tmp1, tmp2 = jetscape[eloss][radius]['10-20'], jetscape[eloss][radius]['20-30']
        tmp3 = {}
        tmp3['ptmin'] = tmp1['ptmin']
        tmp3['ptmax'] = tmp1['ptmax']
        tmp3['pT'] = tmp1['pT']
        tmp3['Ncut'] = 0.5*(tmp1['Ncut'] + tmp2['Ncut'])
        tmp3['dNcut'] = np.sqrt(tmp1['dNcut']*tmp1['dNcut'] + tmp2['dNcut']*tmp2['dNcut'])
        jetscape[eloss][radius]['10-30'] = pd.DataFrame(tmp3)
        tmp4, tmp5 = jetscape[eloss][radius]['30-40'], jetscape[eloss][radius]['40-50'] 
        tmp6 = {}
        tmp6['ptmin'] = tmp4['ptmin']
        tmp6['ptmax'] = tmp4['ptmax']
        tmp6['pT']    = tmp4['pT']
        tmp6['Ncut']  = 0.5*(tmp4['Ncut'] + tmp5['Ncut'])
        tmp6['dNcut'] = np.sqrt(tmp4['dNcut']*tmp4['dNcut'] + tmp5['dNcut']*tmp5['dNcut'])
        jetscape[eloss][radius]['30-50'] = pd.DataFrame(tmp6)

pp_baseline =  {}
for ir, radius in jet_radii.items():
    tmp = pd.read_csv(f'../../jetscape_data/max_time/maxT_200_highstat/pp_2760_jet_spec_jet_rad_{radius}_2.00.csv',comment='#')
    tmp['pT'] = 0.5 * (tmp['ptmin']+tmp['ptmax'])
    #tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
    pp_baseline[radius] = tmp

fig, axes = plt.subplots(5,4,sharex='col', sharey='row', figsize=(16,9))

for ir, radius in jet_radii.items():
    pp = pp_baseline[radius]
    pp = pp[pp['pT'].between(pT_lower_lim,pT_upper_lim)]
    for icent, cent in cents_R_dependence.items():
        if ir < 3:
            util.plot_expr_data_on_axis(axes[ir][icent],CMS_DATA[cent][radius],marker=ddicts.markers_experiment['cms'])
        for eloss in elosses:
            color = module_colors['MATTER+'+eloss.upper()]
            AA = jetscape[eloss][radius][cent]
            AA = AA[AA['pT'].between(pT_lower_lim,pT_upper_lim)]
            raa = AA['Ncut']/pp['Ncut']
            draa = raa*np.sqrt(AA['dNcut']*AA['dNcut']/(AA['Ncut']*AA['Ncut']) + pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']) ) 
            axes[ir][icent].plot(AA['pT'], raa, color=color)
            axes[ir][icent].fill_between(AA['pT'], raa+draa,raa-draa,color=color,alpha=0.2)

#axes[0][0].set_ylim(top=1.2,bottom=-0.01)
for i in range(5):
        axes[i][0].text(0.2, 0.1, r'$R=$'+f"{jet_radii[i]}",transform=axes[i][0].transAxes)
        axes[i][0].set_ylim(bottom=-0.01, top=1.01)
for iax, ax in enumerate(axes[0]):
    centrality = cents_R_dependence[iax]
    centrality = r'$' + centrality + '$'
    ax.text(0.1, 0.9, f"{centrality}" + r"$\%$", transform=ax.transAxes)

#for ax in axes[-1]:
#    ax.set_xlabel(r'$p_T$ (GeV)', fontsize=22)

label = r"$R^{\mathrm{jet}}_{\mathrm{AA}}$"
fig.supylabel(label, fontsize=30)
fig.supxlabel(r'$p_T$ (GeV)', fontsize=30)

handles = [Line2D([],[],label='MATTER+'+eloss.upper(), color=module_colors['MATTER+'+eloss.upper()]) for eloss in elosses]
handle_exp = [Line2D([],[],marker=ddicts.markers_experiment['cms'],color='black',label='CMS (2017)')]

axes[0][3].legend(handles=handle_exp, loc='lower right')
axes[1][3].text(0.1,0.25,r'Pb-Pb, $\sqrt{s}=2.76$ ATeV' + "\n" + "$|\eta|<2$", fontsize=20, transform=axes[1][3].transAxes)
axes[2][3].legend(handles=handles, loc='lower right')

## ax0.set_ylabel(r'$E\frac{\mathrm{d}N^{h^{\pm}}}{\mathrm{d}^3p}$ (GeV${}^{-2}$)', fontsize=22)
## ax1.set_ylabel('Data/Calc.')
#axes[0].set_ylim(bottom=-0.01, top=1.41)
#axes[0].text(0.1,0.7,r'$|\eta|<1.0$', transform=axes[0].transAxes)
#fig.savefig("../../Plots/jet_RAA_PbPb_2p76_CUJET_vs_MARTINI.pdf",dpi=100)
#plt.close(fig)

fig, ax = plt.subplots(1,1,sharey=True, sharex=True,figsize=(20,10))
icent, cent = 0, '00-05'

for eloss in elosses:
    AA_0p2 = jetscape[eloss]['0.2'][cent]
    AA_0p2 = AA_0p2[AA_0p2['pT'].between(pT_lower_lim,pT_upper_lim)]
    pT, y2 = AA_0p2['pT'], AA_0p2['Ncut']
    pp_0p2 = pp_baseline['0.2']
    pp_0p2 = pp_0p2[pp_0p2['pT'].between(pT_lower_lim,pT_upper_lim)]
    raa1   = AA_0p2['Ncut']/pp_0p2['Ncut']
    draa1 = raa1*np.sqrt(AA_0p2['dNcut']**2/AA_0p2['Ncut']**2 + pp_0p2['dNcut']**2/pp_0p2['Ncut']**2)
    for ir, radius in jet_radii.items():
        if radius == '0.2':
            continue
        color = module_colors['MATTER+'+eloss.upper()]
        AA = jetscape[eloss][radius][cent]
        AA = AA[AA['pT'].between(pT_lower_lim,pT_upper_lim)]
        x1, y1 = AA['pT'], AA['Ncut']
        pp = pp_baseline[radius]
        pp = pp[pp['pT'].between(pT_lower_lim,pT_upper_lim)]
        raa2 = AA['Ncut']/pp['Ncut']
        draa2 = raa2*np.sqrt(AA['dNcut']**2/AA['Ncut']**2 + pp['dNcut']**2/pp['Ncut']**2)

        raa = raa2/raa1
        draa = raa*np.sqrt(draa1**2/raa1**2 + draa2**2/raa**2) 
        ax.plot(pT, raa, color=color, linestyle=linestyles[radius])
        ax.fill_between(pT, raa+draa,raa-draa,color=color,alpha=0.2)

for item in linestyles:
    handles.append(Line2D([],[],color='black', linestyle=linestyles[item], label=r'$R=$' + item))
#handles.append(Line2D([],[],color='black',linestyle=linestyles[radius],label=f'R={radius}'))        
label  = r"$R^{\mathrm{jet}}_{\mathrm{AA}}(R)/R^{\mathrm{jet}}_{\mathrm{AA}}(R=0.2)$"
ax.set_ylabel(label, fontsize=30)
ax.legend(handles=handles, loc='upper center')
ax.text(0.85,0.75,r'Pb-Pb' + "\n" + r'$\sqrt{s}=2.76$ ATeV' + "\n" + "$|\eta|<2$" + "\n" + "$00-05\%$", fontsize=30, transform=ax.transAxes)
ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30)
#fig.savefig("../../Plots/Jet_RAA_Cone_Size_Dependence_CUJET_vs_MARTINI.pdf",dpi=200)
plt.show()