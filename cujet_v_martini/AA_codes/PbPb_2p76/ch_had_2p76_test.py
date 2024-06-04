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
    "font.size": 22,
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

pT_lower_lim = 5
pT_upper_lim = 210
inel_Xsec= 62.03948 

absol_path_expt='/Users/rmyazdi/Documents/research/jetscape_project/expt/PbPb_2p76/charged'
centralities = {0:'00-05'}#,1:'05-10',2:'10-20',
#                3:'20-30',4:'30-40',5:'40-50'}

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
for icent, cent in enumerate(['00-05']):
    tmp = pd.read_csv(f"{absol_path_expt}/CMS/Table5_{cent}.csv",comment='#').rename(columns=ddicts.colnames_raa_cms)
    tmp = tmp[tmp['xhigh'].between(pT_lower_lim,pT_upper_lim)]
    cms_experimental_data[cent] = tmp

pp_spec = pd.read_csv("../../jetscape_data/pp/2760/pp_2760_charged_hadrons_eta_cut_1.csv") ## eta < 1, ptmin, ptmax, N, dN^2
pp_spec['pT']  = 0.5*(pp_spec['pTmax']+pp_spec['pTmin'])
pp = pp_spec[pp_spec['pT'].between(pT_lower_lim,pT_upper_lim)]


tmp = pd.read_csv(f"../../jetscape_data/martini/test_max_alphas_0p37/PbPb2760_00-05_charged_hadrons_eta_cut_1.csv",comment='#')
tmp['pT']  = 0.5*(tmp['pTmax']+tmp['pTmin'])
tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
AA_spec_new = tmp


tmp = pd.read_csv(f"../../jetscape_data/martini/PbPb_2760/PbPb2760_00-05_charged_hadrons_eta_cut_1.csv",comment='#')
tmp['pT']  = 0.5*(tmp['pTmax']+tmp['pTmin'])
tmp = tmp[tmp['pT'].between(pT_lower_lim,pT_upper_lim)]
AA_spec_old = tmp

fig, axes = plt.subplots(1,1,gridspec_kw={'top':0.985,'bottom':0.1,
                                          'left':0.10,'right':0.95,
                                          'hspace':0.1,'wspace':0.0},
                                          sharex=True, sharey=True, figsize=(16,9))
#axes = axes.flatten()

for icent, cent in centralities.items():
    alice = alice_experimental_data[cent]
    util.plot_expr_data_on_axis(axes, alice, marker=ddicts.markers['alice'], s=50)

    atlas2 = atlas2_experimental_data[cent]
    util.plot_expr_data_on_axis(axes, atlas2, marker=ddicts.markers['atlas2'], s=50)

    if cent in cms_experimental_data:
        cms = cms_experimental_data[cent]
        util.plot_expr_data_on_axis(axes, cms, marker=ddicts.markers['cms'], s=50) 
    
    if cent in atlas1_experimental_data:
        atlas = atlas1_experimental_data[cent]
        util.plot_expr_data_on_axis(axes, atlas, marker=ddicts.markers['atlas1'], s=50) 
(ptmin, ptmax, proton, dproton) = util.combine_bins(pp['pTmin'].to_list(), pp['pTmax'].to_list(),
                                                    pp['N'].to_list(), pp['dN'].to_list())
pT = 0.5*(ptmax + ptmin)
#pT = pp['pT']#0.5*(pp['pTmin'] + pp['pTmax'])
aa_new = AA_spec_new
(_, _, PbPb, dPbPb) = util.combine_bins(aa_new['pTmin'].to_list(), \
                                        aa_new['pTmax'].to_list(),\
                                        aa_new['N'].to_list(),\
                                        aa_new['dN'].to_list())
raa_new = PbPb/proton
err_new = raa_new*np.sqrt(dPbPb*dPbPb/(PbPb*PbPb) + dproton*dproton/(proton*proton))
axes.plot(pT, raa_new, color=ddicts.eloss_colours['martini'])
axes.fill_between(pT, raa_new+err_new, raa_new-err_new, color=ddicts.eloss_colours['martini'], alpha=0.3)
axes.text(0.65, 0.1, '00-05%', transform=axes.transAxes)

aa_old = AA_spec_old
(_, _, PbPb, dPbPb) = util.combine_bins(aa_old['pTmin'].to_list(), \
                                        aa_old['pTmax'].to_list(),\
                                        aa_old['N'].to_list(),\
                                        aa_old['dN'].to_list())
raa_old = PbPb/proton
err_old = raa_old*np.sqrt(dPbPb*dPbPb/(PbPb*PbPb) + dproton*dproton/(proton*proton))
axes.plot(pT, raa_old, color=ddicts.eloss_colours['cujet'])
axes.fill_between(pT, raa_old+err_old, raa_old-err_old, color=ddicts.eloss_colours['cujet'], alpha=0.3)
axes.text(0.65, 0.1, '00-05%', transform=axes.transAxes)





handles_experiment = [Line2D([],[],marker=ddicts.markers['alice'] ,markersize=10,label=r'ALICE $|\eta|<1.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['atlas2'],markersize=10,label=r'ATLAS $|\eta|<2.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['atlas1'],markersize=10,label=r'ATLAS $|\eta|<1.0$',color='black'),
                      Line2D([],[],marker=ddicts.markers['cms'],markersize=10, label=r'CMS $|\eta|<1.0$',color='black')]
handles_theory = [Line2D([],[],color=ddicts.eloss_colours['martini'],label=r'Martini-$\alpha_s^{\mathrm{max}}=0.37$'),
                  Line2D([],[],color=ddicts.eloss_colours['cujet'],label=r'Martini-$\alpha_s^{\mathrm{max}}=0.42$')]

axes.legend(handles=handles_experiment,loc='best')
axes.legend(handles=handles_theory,loc='upper left')
axes.set_ylim(bottom=-0.01, top=1.51)
axes.text(0.1,0.6,r'$|\eta|<1.0$', transform=axes.transAxes)
axes.text(0.07,0.8,r'Pb-Pb, $\sqrt{s}=2.76$ ATeV', fontsize=30, transform=axes.transAxes)
plt.show()
