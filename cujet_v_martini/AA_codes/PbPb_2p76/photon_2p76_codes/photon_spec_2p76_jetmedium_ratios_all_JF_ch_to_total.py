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
from dictionaries import multiplicity, channel_colors, channel_linestyles, oversampling_factors

mpl.use('Qt5Agg')
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
elosses = ['martini', 'cujet']
jetscape_cents = ['00-05','05-10','10-20','20-30','30-40']
jetscape = {}
fname_tmpl = '../jetscape_data/{eloss}/PbPb_2760/PbPb2760_{cent}_photon_spec_0.80.csv'
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
tmpl_jf = "../other_data/JF_MultiMessenger/PbPb2760_{cent}_{chan}.csv"
for cent in centralities:
    jfp_spec[cent] = {}
    for ch in channels:
        tmp = np.loadtxt(tmpl_jf.format(cent=cent,chan=ch),unpack=True,delimiter=',')
        x = tmp[0]
        y = tmp[1]
        jfp_spec[cent][ch] = pd.DataFrame({'pT':x,'N':y})

rebinned_jf = {}
for cent in centralities:
    rebinned_jf[cent] = {}
    pTmin, pTmax = jetscape['martini'][cent]['ptmin'], jetscape['martini'][cent]['ptmax']
    for ch in channels:
        master_data = jfp_spec[cent][ch]
        x, dx, xmins, xmaxes, y = [],[],[],[],[]
        f = InterpolatedUnivariateSpline(master_data['pT'],np.log(master_data['N'])) 
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
        ratios_v1[eloss][cent]['jet_medium'] = jet_medium_ratio
        for ch in channels:
            ratios_v1[eloss][cent][ch]= np.array([v1/v2 for (v1,v2) in zip(rebinned_jf[cent][ch]['N'],spec)])

        jet_medium_photons[eloss][cent] = {'pT':pT.tolist(),'dpT':dpT,
                                           'N':jet_medium_sources}
fig, axes = plt.subplots(2,2,figsize=(16,12),
                        gridspec_kw={'height_ratios':(1,1),
                                     'top':0.945,'bottom':0.105,
                                     'left':0.12,'right':0.965,
                                     'hspace':0.05,'wspace':0.01},
                        sharey=True, sharex=True)

for icent, cent in enumerate(centralities):
    for ieloss, eloss in enumerate(elosses): 
        spec_ax = axes[ieloss][icent]
        if ieloss == 0:
            spec_ax.text(0.1,0.7,f'{cent}\%', transform=spec_ax.transAxes)
        pT,_,_,_ = total_v1[eloss][cent]
        ratios = ratios_v1[eloss][cent]
        for channel in ratios:
            col = ''
            if channel != 'jet_medium':
                col = channel_colors[channel]
            else:
                col = channel_colors[f'jet_medium-{eloss}']
            lstyle = channel_linestyles[channel]
            ratio= ratios_v1[eloss][cent]['jet_medium'].tolist()
            mark='s'
            spec_ax.plot(pT, ratios[channel], color=col,linestyle=lstyle)

for ax in axes[1]:
    ax.set_xlabel(r'$p_T$ (GeV)')
for ax in [axes[0][0],axes[1][0]]:
    ax.set_ylabel(r'Channel/Total')

handles = []
for eloss in elosses:
    handles.append(Line2D([],[],color=ddicts.eloss_colours[eloss],label="Jet-Medium: "+eloss.upper()))
for ch in channels:
    handles.append(Line2D([],[],color=channel_colors[ch], linestyle=channel_linestyles[ch],label=nice_channels[ch]))

axes[1][1].legend(loc='center right', handles=handles, fontsize=25)
plt.show()