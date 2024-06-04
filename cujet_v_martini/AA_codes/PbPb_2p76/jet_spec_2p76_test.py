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

mpl.use('MacOSX')
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

elosses = ['cujet', 'martini']

pT_lower_lim = 60
pT_upper_lim = 310
inel_Xsec = 62.03948 

mother_path ='/Users/rmyazdi/Documents/research/jetscape_project/expt/PbPb_2p76/jets/'
jet_radii = {0:'0.2',1:'0.3',2:'0.4'}#, 3:'0.6',4:'0.8'}
linestyles = {'0.3':'solid','0.4':'dashed','0.6':'dotted','0.8':'dashdot'}

cents_jetscape = {0:'00-05'}#,1:'05-10',2:'10-20',3:'20-30',4:'30-40',5:'40-50'}

## CMS data for R=0.2, 0.3, 0.4 (ALICE also has 10-30 for R=0.2)
cents_R_dependence = {0:'00-05'}#, 1:'05-10',2:'10-30',3:'30-50'}

## ATLAS
cents_medium_size  = {0:'00-10'}

CMS_DATA = {}
for icent, cent in cents_R_dependence.items():
    CMS_DATA[cent] = {}
    for iR, radius in enumerate(['0.2','0.3','0.4']):
        R = radius.replace('.','p')
        tmp = pd.read_csv(f"{mother_path}CMS/jetscape/cent_{cent}_R_{R}.csv",comment='#').rename(columns=ddicts.colnames_CMS_RAA)
        tmp = tmp[tmp['xhigh'].between(pT_lower_lim,pT_upper_lim)]
        CMS_DATA[cent][radius] = tmp

jetscape_old = {}
elosses = ['martini']
for eloss in elosses:
    jetscape_old[eloss] = {}
    for ir, radius in jet_radii.items():
        jetscape_old[eloss][radius] = {}
        for icent, cent in cents_jetscape.items():
            tmp = pd.read_csv(f"../jetscape_data/{eloss}/PbPb_2760/PbPb2760_{cent}_jet_spec_jet_rad_{radius}_2.00.csv", comment='#')
            tmp['pT'] = 0.5*(tmp['ptmax']+tmp['ptmin'])
            tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
            jetscape_old[eloss][radius][cent] = tmp

jetscape_new = {}
for eloss in elosses:
    jetscape_new[eloss] = {}
    for ir, radius in jet_radii.items():
        jetscape_new[eloss][radius] = {}
        for icent, cent in cents_jetscape.items():
            tmp = pd.read_csv(f"../jetscape_data/martini/test_max_alphas_0p37/PbPb2760_{cent}_jet_spec_jet_rad_{radius}_2.00.csv", comment='#')
            tmp['pT'] = 0.5*(tmp['ptmax']+tmp['ptmin'])
            tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
            jetscape_new[eloss][radius][cent] = tmp

pp_baseline =  {}
for ir, radius in jet_radii.items():
    tmp = pd.read_csv(f'../jetscape_data/pp/2760/pp_2760_jet_spec_jet_rad_{radius}_2.00.csv',comment='#')
    tmp['pT'] = 0.5 * (tmp['ptmin']+tmp['ptmax'])
    tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
    pp_baseline[radius] = tmp

fig, axes = plt.subplots(3,1,gridspec_kw={'top':0.985,'bottom':0.1,
                            'left':0.12,'right':0.95,
                            'hspace':0.08,'wspace':0.1},
                            sharex='col', sharey='row', figsize=(16,9))
#axes = axes.flatten()

for ir, radius in enumerate(['0.2','0.3','0.4']):
    pp = pp_baseline[radius]
    for icent, cent in cents_R_dependence.items():
        util.plot_expr_data_on_axis(axes[ir],CMS_DATA[cent][radius],marker=ddicts.markers_experiment['cms'])
        AA = jetscape_old[eloss][radius][cent]
        raa = AA['Ncut']/pp['Ncut']
        draa = raa*np.sqrt(AA['dNcut']*AA['dNcut']/(AA['Ncut']*AA['Ncut']) + pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']) ) 
        axes[ir].plot(AA['pT'], raa, color=ddicts.eloss_colours['cujet'])
        axes[ir].fill_between(AA['pT'], raa+draa,raa-draa,color=ddicts.eloss_colours['cujet'],alpha=0.2)

        AA = jetscape_new[eloss][radius][cent]
        raa = AA['Ncut']/pp['Ncut']
        draa = raa*np.sqrt(AA['dNcut']*AA['dNcut']/(AA['Ncut']*AA['Ncut']) + pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']) ) 
        axes[ir].plot(AA['pT'], raa, color=ddicts.eloss_colours[eloss])
        axes[ir].fill_between(AA['pT'], raa+draa,raa-draa,color=ddicts.eloss_colours[eloss],alpha=0.2)


#axes[0][0].set_ylim(top=1.2,bottom=-0.01)
for i in range(3):
        label  = r"$R^{\mathrm{jet}}_{\mathrm{AA}}$"
        axes[i].set_ylabel(label, fontsize=22)
        axes[i].text(0.2, 0.1, r'$R=$'+f"{jet_radii[i]}",transform=axes[i].transAxes)
# for iax, ax in enumerate(axes):
#     centrality = cents_R_dependence[iax]
#     centrality = r'$' + centrality + '$'
#     ax.text(0.1, 0.9, f"{centrality}" + r"$\%$", transform=ax.transAxes)
axes[2].set_xlabel(r'$p_T$ (GeV)', fontsize=22)
handles_theory = [Line2D([],[],color=ddicts.eloss_colours['martini'],label=r'Martini-$\alpha_s^{\mathrm{max}}=0.37$'),
                  Line2D([],[],color=ddicts.eloss_colours['cujet'],label=r'Martini-$\alpha_s^{\mathrm{max}}=0.42$')]
handle_exp = [Line2D([],[],marker=ddicts.markers_experiment['cms'],color='black',label='CMS $|\eta|<2$ (2017)')]


axes[0].legend(handles=handle_exp, loc='lower right',fontsize=18)
axes[1].text(0.1,0.25,r'Pb-Pb, $\sqrt{s}=2.76$ ATeV' + "\n" + "$|\eta|<2$", fontsize=20, transform=axes[1].transAxes)
axes[2].legend(handles=handles_theory, loc='lower right')


plt.show()