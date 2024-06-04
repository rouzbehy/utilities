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
from COLORS import rate_set_colors as rate_colours
rate_names = {1:'LO',2:'NLO',3:'NP'}
## Read the data:
master_pp = pd.read_csv("../martini_results/final_AuAu_200/pp/hadron_spectra.csv",comment='#')
master_pp['pT'] = 0.5*(master_pp['pTmax'] + master_pp['pTmin'])
tmpl_loc = "../martini_results/final_AuAu_200/AuAu/rset_{r}/cent_{cent}/hadron_spectra.csv"
centralities = ['0_5','5_10','10_20']
AA = {r: {c: pd.read_csv(tmpl_loc.format(r=r,cent=c), comment='#') 
            for c in centralities} for r in rate_names}
## Combine the 0-5 and 5-10
for rate_set in rate_names:
    AA_spec = AA[rate_set]
    y1 = 0.5*(AA_spec['0_5']['pich'] + AA_spec['5_10']['pich'])
    dy1 = 0.5*np.sqrt(AA_spec['0_5']['dpich']**2 + AA_spec['5_10']['dpich']**2)
    y2 = 0.5*(AA_spec['0_5']['Nch'] + AA_spec['5_10']['Nch'])
    dy2 = 0.5*np.sqrt(AA_spec['0_5']['dNch']**2 + AA_spec['5_10']['dNch']**2)
    y4  = 0.5*(AA_spec['0_5']['pi0'] + AA_spec['5_10']['pi0']) 
    dy4 = 0.5*np.sqrt(AA_spec['0_5']['dpi0']**2 + AA_spec['5_10']['dpi0']**2) 
    y3 = y2 - y1
    dy3 = np.sqrt(dy2**2 + dy1**2)
    AA_spec['0_10'] = pd.DataFrame({'pTmin':AA_spec['0_5']['pTmin'].to_list(), 
                                    'pTmax':AA_spec['0_5']['pTmax'].to_list(),
                                    'pich':y1.to_list(), 'dpich':dy1.to_list(),
                                    'kp':y3.to_list(), 'dkp':dy3.to_list(),'pi0':y4.to_list(),'dpi0':dy4.to_list()})
    AA_spec['10_20']['kp']  = AA_spec['10_20']['Nch'] -AA_spec['10_20']['pich'] 
    AA_spec['10_20']['dkp'] = np.sqrt(AA_spec['10_20']['dNch']**2 -AA_spec['10_20']['dpich']**2)

## EXPERIMENT
exp_loc = '../../../exp_data/sqrts_200GeV'
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

fig, axes = plt.subplots(1,3,gridspec_kw={'top':0.985,'bottom':0.11,
                                      'left':0.08,'right':0.98,
                                      'hspace':0.11,'wspace':0.09},
                              sharex=True, figsize=(16,9), sharey=True)

pTmin, pTmax = 3, 15
for icent, cent in zip([0,2],['0_5','10_20']):
    util.plot_expr_data_on_axis(axes[icent], data[cent], 's')
    tag = cent.replace('_','-')


pp = master_pp[master_pp['pT'].between(pTmin,pTmax)]
pp_err =  pp['dNch']*pp['dNch']/(pp['Nch']*pp['Nch'])
ppy = pp['Nch']
axes = axes.flatten()
for rate_set in rate_names:
    color = rate_colours[rate_names[rate_set]]
    AA_spec = AA[rate_set]
    for icent, cent in enumerate(centralities):
        ax = axes[icent]
        tag = cent.replace('_','-')
        axes[icent].text(0.7,0.5, f'{tag}'+r'$\%$', transform=axes[icent].transAxes)
        aa = AA_spec[cent]
        aa['pT'] = 0.5*(aa['pTmax'] + aa['pTmin']) 
        aa = aa[aa['pT'].between(pTmin,pTmax)]
        y1, dy1 = aa['Nch'], aa['dNch']
        raa1 = y1/ppy
        draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + pp_err)
        ax.plot(aa['pT'], raa1, color=color) 
        ax.fill_between(aa['pT'], raa1+draa1, raa1-draa1, color=color, alpha=0.2)

for ax in axes:
    ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30)
axes[0].set_ylabel(r'$R^{\mathrm{h}^{\pm}}_{\mathrm{AA}}$', fontsize=30)
handels = [Line2D([],[],label=l,color=c) for l, c in rate_colours.items()]
handels.append(Line2D([],[],label=r'STAR(2003), $|\eta|<0.5$', marker='s', color='black'))
axes[1].legend(handles=handels, loc='best')
## Reference for charged hadron RAA: Phys.Rev.Lett. 91 (2003) 172302, 2003.

fig, axes = plt.subplots(nrows=1,ncols=2,sharex=True, figsize=(16,9), sharey=True)

util.plot_expr_data_on_axis(axes[0], data_pi_pm, '^')
util.plot_expr_data_on_axis(axes[1], data_kp, 'v')

pTmin, pTmax = 3, 16
pp = master_pp[master_pp['pT'].between(pTmin,pTmax)]
pp['kp']  = np.array([v1-v2 for v1,v2 in zip(pp['Nch'],pp['pich'])]) 
pp['dkp'] = np.sqrt([v1**2 + v2**2 for v1,v2 in zip(pp['dNch'],pp['dpich'])])

for rate_set in rate_names:
    color = rate_colours[rate_names[rate_set]]
    AA_spec = AA[rate_set]
    cent='0_10'
    for icent in range(2):
        aa = AA_spec[cent]
        ax = axes[icent]
        aa['pT'] = 0.5*(aa['pTmax'] + aa['pTmin']) 
        aa = aa[aa['pT'].between(pTmin,pTmax)]
        tag = 'pich' if icent == 0 else 'kp'
        if tag == 'pich':
            ax.text(0.4, 0.5, r'$\pi^{\pm}$', transform=ax.transAxes)
        else:
            ax.text(0.4, 0.5, r'$h^{\pm}-\pi^{\pm}$', transform=ax.transAxes)
        y1, dy1 = aa[tag], aa[f'd{tag}']
        yp, dyp = pp[tag], pp[f'd{tag}']
        raa1 = y1/yp
        draa1 = raa1*np.sqrt(dy1**2/y1**2 + dyp**2/yp**2)
        ax.plot(aa['pT'], raa1, color=color) 
        ax.fill_between(aa['pT'], raa1+draa1, raa1-draa1, color=color, alpha=0.2)
handles = [Line2D([],[],label=l,color=c) for l, c in rate_colours.items()]
handles.append(Line2D([],[],color='black', label=r'$\pi^{\pm}$', marker='^' ))
handles.append(Line2D([],[],color='black', label=r'$h^{\pm}-\pi^{\pm}$', marker='v' ))
axes[1].text(0.32,0.9,'Data: STAR(2012), $|\eta|<0.5$, 0-12$\%$', transform=axes[1].transAxes)
axes[1].text(0.32,0.8,'Simulation: 0-10$\%$', transform=axes[1].transAxes)
axes[0].legend(loc='upper left', handles=handles)
axes[0].set_ylabel(r'$R_{\mathrm{AA}}^{\pi^{\pm}}$'+' or ' + r'$R_{\mathrm{AA}}\left(h^{\pm}-\pi^{\pm}\right)$')
for ax in axes:
    ax.set_xlabel(r'$p_T$ (GeV)')

pTmin, pTmax = 4.5, 20.1
pp = master_pp[master_pp['pT'].between(pTmin,pTmax)]
fig, axes = plt.subplots(nrows=1,ncols=2,sharex=True, figsize=(16,9), sharey=True)

util.plot_expr_data_on_axis(axes[0], data_pi0['0_10'], 's' )
util.plot_expr_data_on_axis(axes[1], data_pi0['10_20'], 's')
#util.plot_expr_data_on_axis(axes[0], data_pi_pm, '^')
for icent, cent in enumerate([r'0-10$\%$',r'10-20$\%$']):
    axes[icent].text(0.8, 0.9, cent, transform=axes[icent].transAxes)
    axes[icent].set_xlabel(r'$p_T$ (GeV)')

axes[0].set_ylabel(r'$R_{\mathrm{AA}}^{\pi^{0}}$')
for rate_set in rate_names:
    color = rate_colours[rate_names[rate_set]]
    AA_spec = AA[rate_set]
    cent = ['0_10']
    for icent, cent in enumerate(['0_10','10_20']):
        aa = AA_spec[cent]
        ax = axes[icent]
        aa['pT'] = 0.5*(aa['pTmax'] + aa['pTmin']) 
        aa = aa[aa['pT'].between(pTmin,pTmax)]
        tag = 'pi0'
        y1, dy1 = aa[tag], aa[f'd{tag}']
        yp, dyp = pp[tag], pp[f'd{tag}']
        raa1 = y1/yp
        draa1 = raa1*np.sqrt(dy1**2/y1**2 + dyp**2/yp**2)
        ax.plot(aa['pT'], raa1, color=color) 
        ax.fill_between(aa['pT'], raa1+draa1, raa1-draa1, color=color, alpha=0.2)
        #tag = 'pich'
        #y1, dy1 = aa[tag], aa[f'd{tag}']
        #yp, dyp = pp[tag], pp[f'd{tag}']
        #raa1 = y1/yp
        #draa1 = raa1*np.sqrt(dy1**2/y1**2 + dyp**2/yp**2)
        #ax.plot(aa['pT'], raa1, color=color,linestyle='dotted') 
        #ax.fill_between(aa['pT'], raa1+draa1, raa1-draa1, color=color, alpha=0.2)
handles = [Line2D([],[],label=l, color=c) for l, c in rate_colours.items()]
handles.append(Line2D([],[],label='PHENIX (2013) $|\eta|<0.35$', marker='s', color='black'))
axes[0].legend(handles=handles, loc='upper left')
axes[0].text(0.05, 0.7, r'Simulation: $|\eta|<0.5$', transform=axes[0].transAxes)
plt.show()