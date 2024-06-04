## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
#from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.interpolate import interp1d as interpolate
from matplotlib.colors import CSS4_COLORS as css
## my custom modules
import util
import dictionaries as my_dicts
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
    "axes.spines.top" : False,
    "legend.frameon":False
}
plt.rcParams.update(my_rcParams)
from COLORS import module_colors as modules
## Read the data:
jetscape_data_dir = '../../jetscape_data/sqrt_s_200/'
master_pp_ch = pd.read_csv(jetscape_data_dir+"pp/charged_hadrons_eta_cut_0.5.csv",comment='#')
master_pp_ch['pT'] = 0.5*(master_pp_ch['pTmax'] + master_pp_ch['pTmin'])
master_pp_k = pd.read_csv(jetscape_data_dir+"pp/kaon_eta_cut_0.5.csv", comment='#')
master_pp_k['pT'] = 0.5*(master_pp_k['pTmax']+master_pp_k['pTmin'])
master_pp_pich = pd.read_csv(jetscape_data_dir+"pp/pion_eta_cut_0.5.csv", comment='#')
master_pp_pich['pT'] = 0.5*(master_pp_k['pTmax']+master_pp_k['pTmin'])

centralities = ['00-10','10-20']
models = ['martini','cujet']
tmpl_loc = jetscape_data_dir+"AuAu/{model}/{cent}/charged_hadrons_eta_cut_0.5.csv"
AA_ch = {m: {c: pd.read_csv(tmpl_loc.format(model=m,cent=c), comment='#') 
            for c in centralities} for m in models}
tmpl_loc = jetscape_data_dir+"AuAu/{model}/{cent}/kaon_eta_cut_0.5.csv"
AA_k = {m: {c: pd.read_csv(tmpl_loc.format(model=m,cent=c), comment='#') 
            for c in centralities} for m in models}
tmpl_loc = jetscape_data_dir+"AuAu/{model}/{cent}/pion_eta_cut_0.5.csv"
AA_pich = {m: {c: pd.read_csv(tmpl_loc.format(model=m,cent=c), comment='#') 
            for c in centralities} for m in models}
## EXPERIMENT
exp_loc = '../../../../KERNELS_NLO_NP_PART2/exp_data/sqrts_200GeV'
star_charged = exp_loc + '/Charged/RAA_{c}.csv'
data = {c: pd.read_csv(star_charged.format(c=c),comment='#').rename(columns=my_dicts.colnames_STAR_RAA_charged) for c in ['0_5','10_20']}
data_pi_pm = pd.read_csv(exp_loc+'/IdenHad/STAR_2012/Figure3.1_RAAvsp_Tin0-12%Au+Aucollisions..csv',comment='#').rename(columns=my_dicts.colnames_STAR_RAA_identified)
data_kp = pd.read_csv(exp_loc+'/IdenHad/STAR_2012/Figure3.2_RAAvsp_Tin0-12%Au+Aucollisions..csv',comment='#').rename(columns=my_dicts.colnames_STAR_RAA_identified) 

## neutral pion:
data_pi0 = {c:pd.read_csv(exp_loc + f'/IdenHad/PHENIX_2013/pi0_{c}.csv',comment='#').rename(columns=my_dicts.colnames_PHENIX_RAA_pi0) for c in ['0_10','10_20']}
for item in data_pi0:
    spec = data_pi0[item]
    spec['dy_syst+'] = np.sqrt(spec['dy_syst1+']**2 + spec['dy_syst2+']**2)
    spec['dy_syst-'] = np.sqrt(spec['dy_syst1-']**2 + spec['dy_syst2-']**2)

fig, axes = plt.subplots(2,1,sharex=True, figsize=(16,9), sharey=True)

pTmin, pTmax = 5, 15
for icent, cent in zip([0,1],['0_5','10_20']):
    util.plot_expr_data_on_axis(axes[icent], data[cent], 's')
    tag = cent.replace('_','-')
pp = master_pp_ch[master_pp_ch['pT'].between(pTmin,pTmax)]
pp_err =  pp['dN']*pp['dN']/(pp['N']*pp['N'])
ppy = pp['N']
axes = axes.flatten()
for model in models:
    color = modules['MATTER+'+model.upper()]
    AA_spec = AA_ch[model]
    for icent, cent in enumerate(centralities):
        ax = axes[icent]
        tag = cent.replace('_','-')
        axes[icent].text(0.7,0.5, f'{tag}'+r'$\%$', transform=axes[icent].transAxes)
        aa = AA_spec[cent]
        aa['pT'] = 0.5*(aa['pTmax'] + aa['pTmin']) 
        aa = aa[aa['pT'].between(pTmin,pTmax)]
        y1, dy1 = aa['N'], aa['dN']
        raa1 = y1/ppy
        draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + pp_err)
        ax.plot(aa['pT'], raa1, color=color) 
        ax.fill_between(aa['pT'], raa1+draa1, raa1-draa1, color=color, alpha=0.2)

axes[1].set_xlabel(r'$p_T$ (GeV)', fontsize=30)
axes[0].set_ylabel(r'$R^{\mathrm{h}^{\pm}}_{\mathrm{AA}}$', fontsize=30)
handels = [Line2D([],[],label=l,color=c) for l, c in modules.items() if l !='MATTER']
handels.append(Line2D([],[],label=r'STAR(2003), $|\eta|<0.5$', marker='s', color='black'))
axes[0].legend(handles=handels, loc='best')
## Reference for charged hadron RAA: Phys.Rev.Lett. 91 (2003) 172302, 2003.

fig, axes = plt.subplots(nrows=2,ncols=1,sharex=True, figsize=(16,9), sharey=True)

pTmin, pTmax = 5, 16
data_pi_pm = data_pi_pm[data_pi_pm['x'].between(pTmin, pTmax)]
data_kp = data_kp[data_kp['x'].between(pTmin,pTmax)]
util.plot_expr_data_on_axis(axes[0], data_pi_pm, '^')
util.plot_expr_data_on_axis(axes[1], data_kp, 'v')

pp_ch = master_pp_ch[master_pp_ch['pT'].between(pTmin,pTmax)]
pp_pich = master_pp_pich[master_pp_pich['pT'].between(pTmin,pTmax)]
pp_kp  = np.array([v1-v2 for v1,v2 in zip(pp_ch['N'],pp_pich['N'])]) 
pp_dkp = np.sqrt([v1**2 + v2**2 for v1,v2 in zip(pp_ch['dN'],pp_pich['dN'])])

for model in models:
    color = modules['MATTER+'+model.upper()]
    AA_spec_ch = AA_ch[model]
    AA_spec_pich = AA_pich[model]
    cent = '00-10'
    for icent in range(2):
        aa_ch = AA_spec_ch[cent]
        aa_pich = AA_spec_pich[cent]
        aa_ch['pT'] = 0.5*(aa_ch['pTmax'] + aa_ch['pTmin']) 
        aa_pich['pT'] = 0.5*(aa_pich['pTmax'] + aa_pich['pTmin']) 
        aa_ch_cut = aa_ch[aa_ch['pT'].between(pTmin,pTmax)]
        aa_pich_cut = aa_pich[aa_pich['pT'].between(pTmin,pTmax)]
        x = aa_ch_cut['pT']
        aa_kp  = np.array([v1-v2 for v1,v2 in zip(aa_ch_cut['N'],aa_pich_cut['N'])]) 
        aa_dkp = np.sqrt([v1**2 + v2**2 for v1,v2 in zip(aa_ch_cut['dN'],aa_pich_cut['dN'])])
        ax = axes[icent]
        y1, dy1 = (aa_pich_cut['N'], aa_pich_cut['dN']) if icent == 0 else (aa_kp, aa_dkp)
        yp, dyp = (pp_pich['N'], pp_pich['dN']) if icent == 0 else (pp_kp, pp_dkp)
        #if icent == 0:
        #    ax.text(0.4, 0.5, r'$\pi^{\pm}$', transform=ax.transAxes)
        #else:
        #    ax.text(0.4, 0.5, r'$h^{\pm}-\pi^{\pm}$', transform=ax.transAxes)
        raa1 = y1/yp
        draa1 = raa1*np.sqrt(dy1**2/y1**2 + dyp**2/yp**2)
        ax.plot(aa['pT'], raa1, color=color) 
        ax.fill_between(aa['pT'], raa1+draa1, raa1-draa1, color=color, alpha=0.2)

handles = [Line2D([],[],label=l,color=c) for l, c in modules.items() if l !='MATTER']
handles.append(Line2D([],[],color='black', label=r'$\pi^{\pm}$', marker='^' ))
handles.append(Line2D([],[],color='black', label=r'$h^{\pm}-\pi^{\pm}$', marker='v' ))
axes[1].text(0.1,0.9,'Data: STAR(2012), $|\eta|<0.5$, 0-12$\%$', transform=axes[1].transAxes)
axes[1].text(0.1,0.75,'Simulation: 0-10$\%$', transform=axes[1].transAxes)
axes[0].legend(loc='upper left', handles=handles, ncols=2)
axes[0].set_ylabel(r'$R_{\mathrm{AA}}^{\pi^{\pm}}$')
axes[1].set_ylabel(r'$R_{\mathrm{AA}}\left(h^{\pm}-\pi^{\pm}\right)$')
axes[1].set_xlabel(r'$p_T$ (GeV)')
plt.show()
# pTmin, pTmax = 4.5, 20.1
# pp = master_pp[master_pp['pT'].between(pTmin,pTmax)]
# fig, axes = plt.subplots(nrows=1,ncols=2,sharex=True, figsize=(16,9), sharey=True)

# util.plot_expr_data_on_axis(axes[0], data_pi0['0_10'], 's' )
# util.plot_expr_data_on_axis(axes[1], data_pi0['10_20'], 's')
# #util.plot_expr_data_on_axis(axes[0], data_pi_pm, '^')
# for icent, cent in enumerate([r'0-10$\%$',r'10-20$\%$']):
#     axes[icent].text(0.8, 0.9, cent, transform=axes[icent].transAxes)
#     axes[icent].set_xlabel(r'$p_T$ (GeV)')

# axes[0].set_ylabel(r'$R_{\mathrm{AA}}\left(\pi^{0}\right)$')
# for rate_set in rate_names:
#     color = rate_colours[rate_names[rate_set]]
#     AA_spec = AA[rate_set]
#     for icent, cent in enumerate(['0_10','10_20']):
#         aa = AA_spec[cent]
#         ax = axes[icent]
#         aa['pT'] = 0.5*(aa['pTmax'] + aa['pTmin']) 
#         aa = aa[aa['pT'].between(pTmin,pTmax)]
#         tag = 'pi0'
#         y1, dy1 = aa[tag], aa[f'd{tag}']
#         yp, dyp = pp[tag], pp[f'd{tag}']
#         raa1 = y1/yp
#         draa1 = raa1*np.sqrt(dy1**2/y1**2 + dyp**2/yp**2)
#         ax.plot(aa['pT'], raa1, color=color) 
#         ax.fill_between(aa['pT'], raa1+draa1, raa1-draa1, color=color, alpha=0.2)
#         #tag = 'pich'
#         #y1, dy1 = aa[tag], aa[f'd{tag}']
#         #yp, dyp = pp[tag], pp[f'd{tag}']
#         #raa1 = y1/yp
#         #draa1 = raa1*np.sqrt(dy1**2/y1**2 + dyp**2/yp**2)
#         #ax.plot(aa['pT'], raa1, color=color,linestyle='dotted') 
#         #ax.fill_between(aa['pT'], raa1+draa1, raa1-draa1, color=color, alpha=0.2)
# handles = [Line2D([],[],label=l, color=c) for l, c in rate_colours.items()]
# handles.append(Line2D([],[],label='PHENIX (2013) $|\eta|<0.35$', marker='s', color='black'))
# axes[0].legend(handles=handles, loc='upper left')
# axes[0].text(0.05, 0.7, r'Simulation: $|\eta|<0.5$', transform=axes[0].transAxes)
# plt.show()