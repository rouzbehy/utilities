#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d
import util 
import jetDicts as ddicts
from dictionaries import multiplicity_PbPb_2p76 as multiplicity
from dictionaries import channel_colors, channel_linestyles, oversampling_factors
from COLORS import module_colors, colors
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
inel_Xsec = 62.03948 
elosses = ['cujet', 'martini']
## read in experimental data
exp_loc = "../../../expt/PbPb_2p76/photons/"
centralities = {'00-20':exp_loc+'HEPData-ins1394677-v1-Table_1.csv',
                '20-40':exp_loc+'HEPData-ins1394677-v1-Table_2.csv'}

expdata = {}

pT_lower_lim = {'00-20':0,'20-40':0}
pT_upper_lim = {'00-20':0,'20-40':0}
for cent in centralities:
    tmp = pd.read_csv(centralities[cent],comment='#')
    pTmin, pTmax = tmp['x'].min(), tmp['x'].max()
    pT_lower_lim[cent] = pTmin
    pT_upper_lim[cent] = pTmax
    tmp = tmp[tmp['x'].between(pT_lower_lim[cent],pT_upper_lim[cent])]
    expdata[cent] = tmp


elosses = ['martini','cujet']
jetscape_cents = ['00-05','05-10','10-20','20-30','30-40']
jetscape = {}
fname_tmpl = '../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_photon_spec_0.80.csv'
for eloss in elosses:
    jetscape[eloss] = {}
    for cent in jetscape_cents:
        fname = fname_tmpl.format(eloss=eloss,cent=cent)
        factor = oversampling_factors[f'{eloss}-{cent}']
        Nbin = multiplicity[cent]
        tmp = pd.read_csv(fname,comment='#')
        xmin = tmp['ptmin']
        xmax = tmp['ptmax']
        x = 0.5*(xmin+xmax)
        dx = xmax - xmin
        tmp['conv']  *= Nbin/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        tmp['dconv'] *= Nbin/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        tmp['brem']  *= Nbin/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        tmp['dbrem'] *= Nbin/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        jetscape[eloss][cent] = tmp
    ## construct the 0-20 and 20-40 centralities from the above
    tmp1 = util.combine_centralities(jetscape[eloss]['00-05'],jetscape[eloss]['05-10'])
    tmp2 = util.combine_centralities(tmp1, jetscape[eloss]['10-20'])
    jetscape[eloss]['00-20'] = tmp2[tmp2['pT'].between(pT_lower_lim['00-20'],pT_upper_lim['00-20'])]
    tmp3 = util.combine_centralities(jetscape[eloss]['20-30'],jetscape[eloss]['30-40'])
    jetscape[eloss]['20-40'] = tmp3[tmp3['pT'].between(pT_lower_lim['20-40'],pT_upper_lim['20-40'])]
centralities = ['00-20','20-40']
channels = ['thermal','preEq', 'prompt']
nice_channels = {'prompt':'Prompt','thermal':'Thermal','preEq':'Pre-Equilibrium'}
other_photons = {}
tmpl_jf = "../../other_data/JF_MultiMessenger/PbPb2760_{cent}_{chan}.csv"

for cent in centralities:
    other_photons[cent] = {}
    for ch in channels:
        if ch == 'prompt':
            continue
        tmp = np.loadtxt(tmpl_jf.format(cent=cent,chan=ch),unpack=True,delimiter=',')
        x = tmp[0]
        y = tmp[1]
        other_photons[cent][ch] = pd.DataFrame({'pT':x,'N':y})

kfactor = -1
tmp = pd.read_csv('../../jetscape_data/prompt_photons/PbPb_2760/gamma_spectra.csv', comment='#')
tmp_x   = 0.5*(tmp['pTmin'] + tmp['pTmax'])
tmp_dx  = tmp['pTmax'] - tmp['pTmin']
tmp_y   = tmp['prompt']/(2*np.pi*tmp_x*tmp_dx*inel_Xsec*2*0.8)
tmp_dy  = tmp['dprompt']/(2*np.pi*tmp_x*tmp_dx*inel_Xsec*2*0.8)
for cent in centralities:
    yy  = tmp_y*multiplicity[cent]
    dyy = tmp_dy*np.sqrt(multiplicity[cent])
    if kfactor < 0 and cent =='00-20':
        f = interp1d(tmp_x, np.log(yy), kind='linear', fill_value='extrapolate')
        dat_x, dat_y = expdata['00-20']['x'].to_list(), expdata['00-20']['y'].to_list()
        print(dat_x[-1])
        kfactor = dat_y[-1]/np.exp(f(dat_x[-1]))
    print(kfactor)
    yy = kfactor*yy
    dyy = np.sqrt(kfactor)*dyy
    other_photons[cent]['prompt'] = pd.DataFrame({'pT':tmp_x.to_list(), 'N':yy.to_list(), 'dN':dyy.to_list()})
    
x = np.linspace(1, 22, 15)
for cent in other_photons:
    specs = other_photons[cent]
    f_prompt  = interp1d(specs['prompt']['pT'] , np.log(specs['prompt']['N'])  , kind='linear', fill_value='extrapolate')
    f_thermal = interp1d(specs['thermal']['pT'], np.log(specs['thermal']['N']) , kind='linear', fill_value='extrapolate')
    f_preEq   = interp1d(specs['preEq']['pT']  , np.log(specs['preEq']['N'])   , kind='linear', fill_value='extrapolate')
    df_prompt = interp1d(specs['prompt']['pT'] , np.log(specs['prompt']['dN']) , kind='linear', fill_value='extrapolate')
    ytot      = np.exp(f_preEq(x)) + np.exp(f_prompt(x)) + np.exp(f_thermal(x)) 
    dytot     = np.exp(df_prompt(x))
    other_photons[cent]['total'] = (x,ytot,dytot)

total_with_jet_medium = {}
for eloss in elosses:
    total_with_jet_medium[eloss] = {}
    for cent in centralities:
        tmpx  = jetscape[eloss][cent]['pT']
        dx    = jetscape[eloss][cent]['dpT']
        spec  = jetscape[eloss][cent]['conv']+ jetscape[eloss][cent]['brem']
        dspec = np.sqrt(jetscape[eloss][cent]['dconv']**2 +jetscape[eloss][cent]['dbrem']**2)
        fjetmed     = interp1d(tmpx, np.log(spec), kind='linear', fill_value='extrapolate')
        dfjetmed    = interp1d(tmpx, np.log(dspec), kind='linear', fill_value='extrapolate')
        oth         = other_photons[cent]['total']
        total_spec  = oth[1] + np.exp(fjetmed(x))
        dtotal_spec = np.sqrt(oth[2]**2 + np.exp(dfjetmed(x))**2)
        total_with_jet_medium[eloss][cent] = (x, total_spec, dtotal_spec)


fig, axes = plt.subplots(2,2,figsize=(16,12),
                        gridspec_kw={'height_ratios':(1,1)},
                        sharey=True, sharex=True)

for icent, cent in enumerate(centralities):
    data = expdata[cent]
    for ieloss, eloss in enumerate(elosses):
        color = module_colors['MATTER+'+eloss.upper()]
        spec_ax = axes[ieloss][icent]
        util.plot_expr_data_on_axis(spec_ax, data, marker='o')
        if ieloss == 0:
            spec_ax.text(0.1,0.2,f'{cent}\%', transform=spec_ax.transAxes)
        pT,spec_with_jet_medium,dspec_with_jet_medium = total_with_jet_medium[eloss][cent]
        pT,spec_without_jet_medium,dspec_without_jet_medium = other_photons[cent]['total']
        f_jm = interp1d(pT, np.log(spec_with_jet_medium), kind='linear', fill_value='extrapolate')
        df_jm = interp1d(pT, np.log(dspec_with_jet_medium), kind='linear', fill_value='extrapolate')
        f_njm = interp1d(pT,  np.log(spec_without_jet_medium), kind='linear', fill_value='extrapolate')
        df_njm = interp1d(pT,  np.log(dspec_without_jet_medium), kind='linear', fill_value='extrapolate')
        spec_jm = np.exp(f_jm(data['x']))
        dspec_jm = np.exp(df_jm(data['x']))
        spec_njm = np.exp(f_njm(data['x']))
        dspec_njm = np.exp(df_njm(data['x']))
        pT = data['x']
        spec_ax.plot(pT, spec_jm, color=color)
        spec_ax.fill_between(pT, spec_jm-dspec_jm, spec_jm+dspec_jm, color=color, alpha=0.3)
        spec_ax.plot(pT, spec_njm, color=colors[2], linestyle='dashed')
        spec_ax.fill_between(pT, spec_njm-dspec_njm, spec_njm+dspec_njm, color=colors[2], alpha=0.3)

for ax in axes[1]:
    ax.set_yscale('log')
    ax.set_xlabel(r'$p_T$ (GeV)')
for ax in [axes[0][0],axes[1][0]]:
    ax.set_ylabel(r'$\frac{1}{N_{\mathrm{evt}}\,2\pi\,p_T}\frac{\mathrm{d}N^{\gamma}}{\mathrm{d}p_T\mathrm{d}\eta}$, (GeV$^{-2}$)')

handles = [Line2D([],[],label=r'ALICE (2016) $|\eta|<0.8$', marker='o',color='black')]
for eloss in elosses:
    handles.append(Line2D([],[],color=module_colors['MATTER+'+eloss.upper()],label='MATTER+'+eloss.upper()))
handles.append(Line2D([],[],color=colors[2],label='Without Jet-Medium', linestyle='dashed'))
axes[0][1].legend(loc='upper right', handles=handles, fontsize=25)
axes[0][0].text(0.3,0.8, r'Pb-Pb @ $\sqrt{s}=2.76$ ATeV', transform=axes[0][0].transAxes)

plt.show()