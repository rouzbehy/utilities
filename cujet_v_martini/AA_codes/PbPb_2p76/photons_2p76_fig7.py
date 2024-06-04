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
from dictionaries import channel_linestyles, oversampling_factors
from COLORS import module_colors, channel_colors
plt.rcParams.update(util.my_rcParams)
inel_Xsec = 62.03948 

## read in experimental data
exp_loc = "../../../expt/PbPb_2p76/photons/"
centralities = {'00-20':exp_loc+'HEPData-ins1394677-v1-Table_1.csv',
                '20-40':exp_loc+'HEPData-ins1394677-v1-Table_2.csv'}

expdata = {}
pT_lower_lim = {'00-20':0,'20-40':0}
pT_upper_lim = {'00-20':0,'20-40':0}
for cent in centralities:
    tmp = pd.read_csv(centralities[cent],comment='#')
    min_pT = min(tmp['x'])
    max_pT = max(tmp['x'])
    pT_lower_lim[cent] = 2#min_pT
    pT_upper_lim[cent] = 14#max_pT
    tmp = tmp[tmp['x'].between(pT_lower_lim[cent],pT_upper_lim[cent])]
    expdata[cent] = tmp

elosses = ['martini', 'cujet']
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
    jetscape[eloss]['00-20'] = tmp2
    tmp3 = util.combine_centralities(jetscape[eloss]['20-30'],jetscape[eloss]['30-40'])
    jetscape[eloss]['20-40'] = tmp3
centralities = ['00-20','20-40']
channels = ['prompt','thermal','preEq']
nice_channels = {'prompt':'Prompt','thermal':'Thermal','preEq':'Pre-Equilibrium'}
jfp_spec = {}
tmpl_jf = "../../other_data/JF_MultiMessenger/PbPb2760_{cent}_{chan}.csv"
for cent in centralities:
    jfp_spec[cent] = {}
    for ch in channels:
        if ch == 'prompt':
            continue
        tmp = np.loadtxt(tmpl_jf.format(cent=cent,chan=ch),unpack=True,delimiter=',')
        x = tmp[0]
        y = tmp[1]
        jfp_spec[cent][ch] = pd.DataFrame({'pT':x,'N':y})

kfactor = -1
tmp = pd.read_csv('../../jetscape_data/prompt_photons/PbPb_2760/gamma_spectra.csv', comment='#')

tmp['x'] = 0.5*(tmp['pTmin'] + tmp['pTmax'])
tmp = tmp[tmp['x'].between(3, 50)]
tmp_dx  = tmp['pTmax'] - tmp['pTmin']
tmp_y   = tmp['prompt']/(2*np.pi*tmp['x']*tmp_dx*inel_Xsec*2*0.8)
tmp_dy  = tmp['dprompt']/(2*np.pi*tmp['x']*tmp_dx*inel_Xsec*2*0.8)
for cent in centralities:
    yy  = tmp_y*multiplicity[cent]
    dyy = tmp_dy*np.sqrt(multiplicity[cent])
    if kfactor < 0 and cent =='00-20':
        f = interp1d(tmp['x'], np.log(yy), kind='linear', fill_value='extrapolate')
        dat_x, dat_y = expdata['00-20']['x'].to_list(), expdata['00-20']['y'].to_list()
        kfactor = dat_y[-1]/np.exp(f(dat_x[-1]))
    print(kfactor)
    yy = kfactor*yy
    dyy = np.sqrt(kfactor)*dyy
    jfp_spec[cent]['prompt'] = pd.DataFrame({'pT':tmp['x'].to_list(), 'N':yy.to_list(), 'dN':dyy.to_list()})



rebinned_jf = {}
for cent in centralities:
    rebinned_jf[cent] = {}
    pTmin, pTmax = jetscape['martini'][cent]['ptmin'], jetscape['martini'][cent]['ptmax']
    for ch in channels:
        master_data = jfp_spec[cent][ch]
        x, dx, xmins, xmaxes, y = [],[],[],[],[]
        f = interp1d(master_data['pT'],np.log(master_data['N']), kind='linear', fill_value='extrapolate') 
        for (xlow, xhigh) in zip(pTmin, pTmax):
            xval = 0.5*(xhigh + xlow)
            dxval = xhigh - xlow
            yval = np.exp(f(xval))
            x.append(xval)
            xmins.append(xlow)
            xmaxes.append(xhigh)
            dx.append(dxval)
            y.append(yval)
        new_master = pd.DataFrame({'ptmin':xmins,'ptmax':xmaxes,'pT':x,'dpT':dx,'N':y})
        tmp = new_master
        rebinned_jf[cent][ch] = tmp

total_v1 = {}
for eloss in elosses:
    total_v1[eloss] = {}
    for cent in centralities:
        num_binary_coll = multiplicity[cent]
        x = jetscape[eloss][cent]['pT'].tolist()
        dx = jetscape[eloss][cent]['dpT'].tolist()
        spec = jetscape[eloss][cent]['conv']+ jetscape[eloss][cent]['brem']
        dspec = np.sqrt(jetscape[eloss][cent]['dconv']*jetscape[eloss][cent]['dconv']+\
                        jetscape[eloss][cent]['dbrem']*jetscape[eloss][cent]['dbrem'])
        spec = spec.tolist() 
        for ch in channels:
            new = rebinned_jf[cent][ch]['N'].tolist()
            spec = [v1 + v2 for (v1,v2) in zip(spec, new)]
        total_v1[eloss][cent] = (np.array(x), np.array(dx), np.array(spec), np.array(dspec.tolist()))

ratios_v1 = {}
jet_medium_photons = {}
locc = '/Users/rmyazdi/Documents/research/jetscape_project/v2/jetscape_data/photon_sp/'
for eloss in elosses:
    ratios_v1[eloss] = {}
    jet_medium_photons[eloss]= {}
    for cent in centralities:
        ratios_v1[eloss][cent] = {}
        pT, dpT, spec, _ = total_v1[eloss][cent]
        num_binary_coll = multiplicity[cent]
        jet_medium_sources = jetscape[eloss][cent]['conv']+jetscape[eloss][cent]['brem']
        jet_medium_sources = np.array(jet_medium_sources.tolist())
        jet_medium_ratio = jet_medium_sources/spec
        if cent == '00-20':
            with open(locc+f'photons_00-20_2p76_{eloss}_jetscape.csv', 'w') as f:
                f.write('pT,jmed,othr,total\n')
                for item in zip(pT, jet_medium_sources, spec):
                    x, yjet, tot = item
                    ynonjet = tot - yjet
                    f.write(f'{x:0.3e},{yjet:0.6e},{ynonjet:0.6e},{tot:0.6e}\n')
        ratios_v1[eloss][cent]['jet_medium'] = jet_medium_ratio
        for ch in channels:
            ratios_v1[eloss][cent][ch]= np.array([v1/v2 for (v1,v2) in zip(rebinned_jf[cent][ch]['N'],spec)])

        jet_medium_photons[eloss][cent] = {'pT':pT.tolist(),'dpT':dpT,
                                           'N':jet_medium_sources}
fig, axes = plt.subplots(2,2,figsize=(16,12),
                        gridspec_kw={'height_ratios':(1,1)},
                        sharey=True, sharex=True)
subfigure_labels = [['(a)', '(b)'], ['(c)', '(d)']]
for icent, cent in enumerate(centralities):
    for ieloss, eloss in enumerate(elosses): 
        spec_ax = axes[ieloss][icent]
        spec_ax.text(0.1, 0.9, subfigure_labels[icent][ieloss], transform=spec_ax.transAxes)
        if ieloss == 0:
            spec_ax.text(0.1,0.7,f'{cent}\%', transform=spec_ax.transAxes)
        pT,_,_,_ = total_v1[eloss][cent]
        ratios = ratios_v1[eloss][cent]
        for channel in ratios:
            col = ''
            if channel != 'jet_medium':
                col = channel_colors[channel]
            else:
                col = module_colors[f'MATTER+{eloss.upper()}']
            lstyle = channel_linestyles[channel]
            ratio= ratios_v1[eloss][cent]['jet_medium'].tolist()
            mark='s'
            spec_ax.plot(pT, ratios[channel], color=col,linestyle=lstyle)

for ax in axes[1]:
    ax.set_xlabel(r'$p_T$ (GeV)')
for ax in [axes[0][0],axes[1][0]]:
    ax.set_ylabel(r'Channel/Total')
#fig.supxlabel(r'$p_T$ (GeV)')
#fig.supylabel(r'Channel/Total')

handles = []
for eloss in elosses:
    handles.append(Line2D([],[],color=module_colors[f'MATTER+{eloss.upper()}'],label=f'MATTER+{eloss.upper()}'))
jf_handles = []
for ch in channels:
    jf_handles.append(Line2D([],[],color=channel_colors[ch], linestyle=channel_linestyles[ch],label=nice_channels[ch]))

axes[1][1].legend(loc='center right', handles=handles)#, fontsize=20)
axes[0][1].legend(loc='center right', handles=jf_handles)#, fontsize=20)
axes[0][0].text(0.4,0.6, r'Pb-Pb @ $\sqrt{s}=2.76$ ATeV', transform=axes[0][0].transAxes)#, fontsize=20)
axes[0][0].text(0.6,0.5, r'$|\eta|<0.8$', transform=axes[0][0].transAxes)#, fontsize=20)
#fig.savefig("../../Plots/Photon_Plots/fig7_photon_channel_breakdown_ratio_total.pdf", dpi=200)
plt.show()
