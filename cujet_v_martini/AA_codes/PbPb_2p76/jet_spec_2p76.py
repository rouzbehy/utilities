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

plt.rcParams.update(util.my_rcParams)

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
            tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
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
    tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
    pp_baseline[radius] = tmp

fig, axes = plt.subplots(3,4,sharex='col', sharey='row', figsize=(16,9))

for ir, radius in enumerate(['0.2','0.3','0.4']):
    pp = pp_baseline[radius]
    for icent, cent in cents_R_dependence.items():
        util.plot_expr_data_on_axis(axes[ir][icent],CMS_DATA[cent][radius],marker=ddicts.markers_experiment['cms'])
        for eloss in elosses:
            color = module_colors['MATTER+'+eloss.upper()]
            AA = jetscape[eloss][radius][cent]
            raa = AA['Ncut']/pp['Ncut']
            draa = raa*np.sqrt(AA['dNcut']*AA['dNcut']/(AA['Ncut']*AA['Ncut']) + pp['dNcut']*pp['dNcut']/(pp['Ncut']*pp['Ncut']) ) 
            axes[ir][icent].plot(AA['pT'], raa, color=color)
            axes[ir][icent].fill_between(AA['pT'], raa+draa,raa-draa,color=color,alpha=0.2)

#axes[0][0].set_ylim(top=1.2,bottom=-0.01)
for i in range(3):
        #label  = r"$R^{\mathrm{jet}}_{\mathrm{AA}}$"
        #axes[i][0].set_ylabel(label, fontsize=22)
        axes[i][0].text(0.2, 0.1, r'$R=$'+f"{jet_radii[i]}",transform=axes[i][0].transAxes)
        axes[i][0].set_ylabel(r"$R^{\mathrm{jet}}_{\mathrm{AA}}$", fontsize=30)
for iax, ax in enumerate(axes[0]):
    centrality = cents_R_dependence[iax]
    centrality = r'$' + centrality + '$'
    ax.text(0.1, 0.9, f"{centrality}" + r"$\%$", transform=ax.transAxes)
for ax in axes[2]:
    ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30)

handles = [Line2D([],[],label='MATTER+'+eloss.upper(), color=module_colors['MATTER+'+eloss.upper()]) for eloss in elosses]
handle_exp = [Line2D([],[],marker=ddicts.markers_experiment['cms'],color='black',label='CMS $|\eta|<2$ (2017)')]


axes[0][3].legend(handles=handle_exp, loc='lower right',fontsize=25)
axes[1][3].text(0.1,0.25,r'Pb-Pb, $\sqrt{s}=2.76$ ATeV' + "\n" + "$|\eta|<2$", fontsize=30, transform=axes[1][3].transAxes)
axes[2][3].legend(handles=handles, loc='lower right',fontsize=25)

## ax0.set_ylabel(r'$E\frac{\mathrm{d}N^{h^{\pm}}}{\mathrm{d}^3p}$ (GeV${}^{-2}$)', fontsize=22)
## ax1.set_ylabel('Data/Calc.')
#axes[0].set_ylim(bottom=-0.01, top=1.41)
#axes[0].text(0.1,0.7,r'$|\eta|<1.0$', transform=axes[0].transAxes)
# fig.savefig("../../Plots/jet_RAA_PbPb_2p76_CUJET_vs_MARTINI.pdf",dpi=100)
# plt.close(fig)

#fig, ax = plt.subplots(1,1,
#                        gridspec_kw={'top':0.985,'bottom':0.1,
#                        'left':0.12,'right':0.95,
#                        'hspace':0.08,'wspace':0.1},
#                        sharey=True, sharex=True,figsize=(20,10))
##axes = ax.flatten()
#icent, cent = 0, '00-05'
#for ir, radius in jet_radii.items():
#    #ax = axes[ir]
#    #pp = pp_baseline[radius]
#    #(pp_xmin, pp_xmax, pp_y, pp_dy) = util.combine_bins(pp['ptmin'].tolist(), pp['ptmax'].tolist(),\
#    #                                                    pp['Ncut'].tolist(), pp['dNcut'].tolist())
#    #pT = 0.5*(pp_xmin+pp_xmax)
#    
#    for eloss in elosses:
#        if radius == '0.2':
#            continue
#        AA = jetscape[eloss][radius][cent]
#        
#        (aa_xmin, aa_xmax, aa_y, aa_dy) = util.combine_bins(AA['ptmin'].tolist(), AA['ptmax'].tolist(),\
#                                                        AA['Ncut'].tolist(), AA['dNcut'].tolist())
#        AA_0p2 = jetscape[eloss]['0.2'][cent]
#        (_, _, aa_y_0p2, aa_dy_0p2) = util.combine_bins(AA_0p2['ptmin'].tolist(), AA_0p2['ptmax'].tolist(),\
#                                                        AA_0p2['Ncut'].tolist(), AA_0p2['dNcut'].tolist())
#        pT = 0.5*(aa_xmin+aa_xmax)
#        raa = aa_y/aa_y_0p2
#        draa = raa*np.sqrt(aa_dy*aa_dy/(aa_y*aa_y) + (aa_dy_0p2*aa_dy_0p2)/(aa_y_0p2*aa_y_0p2)) 
#        ax.plot(pT, raa, color=ddicts.eloss_colours[eloss], linestyle=linestyles[radius])
#        ax.fill_between(pT, raa+draa,raa-draa,color=ddicts.eloss_colours[eloss],alpha=0.2)
#    #if ir > 0:
#    #    ax.spines['left'].set_visible(False)
#    #    ax.yaxis.set_visible(False)
#    #ax.text(0.2, 0.9, f"R={radius}", transform=ax.transAxes)
#
#ax.set_ylabel(r'$R/(R=0.2)$')
#
#label  = r"$R^{\mathrm{jet}}_{\mathrm{AA}}$"
##axes[0].set_ylabel(label, fontsize=40)
##axes[1].legend(handles=handles, loc='lower center')
##axes[2].text(0.05,0.1,r'Pb-Pb' + "\n" + r'$\sqrt{s}=2.76$ ATeV' + "\n" + "$|\eta|<2$" + "\n" + "$00-05\%$", fontsize=30, transform=axes[2].transAxes)
##fig.savefig("../../Plots/Jet_RAA_Cone_Size_Dependence_CUJET_vs_MARTINI.pdf",dpi=200)
plt.show()
