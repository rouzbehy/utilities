#!/usr/bin/env python3
from calendar import c
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
elosses = ['cujet', 'martini']
## read in experimental data
exp_loc = "../../expt/PbPb_2p76/photons/"
centralities = {'00-20':exp_loc+'HEPData-ins1394677-v1-Table_1.csv',
                '20-40':exp_loc+'HEPData-ins1394677-v1-Table_2.csv'}

expdata = {}

pT_lower_lim = {'00-20':0,'20-40':0}
pT_upper_lim = {'00-20':13,'20-40':13}
for cent in centralities:
    tmp = pd.read_csv(centralities[cent],comment='#')
    tmp = tmp[tmp['x'].between(pT_lower_lim[cent],pT_upper_lim[cent])]
    expdata[cent] = tmp


elosses = ['martini','cujet']
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
    jetscape[eloss]['00-20'] = tmp2[tmp2['pT'].between(pT_lower_lim['00-20'],pT_upper_lim['00-20'])]
    tmp3 = util.combine_centralities(jetscape[eloss]['20-30'],jetscape[eloss]['30-40'])
    jetscape[eloss]['20-40'] = tmp3[tmp3['pT'].between(pT_lower_lim['20-40'],pT_upper_lim['20-40'])]
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
        tmp = pd.DataFrame({'ptmin':xmins,'ptmax':xmaxes,'pT':x,'dpT':dx,'N':y})
        rebinned_jf[cent][ch] = tmp

total_with_jet_medium = {}
for eloss in elosses:
    total_with_jet_medium[eloss] = {}
    for cent in centralities:
        x = jetscape[eloss][cent]['pT'].tolist()
        dx = jetscape[eloss][cent]['dpT'].tolist()
        spec = jetscape[eloss][cent]['conv']+ jetscape[eloss][cent]['brem']
        dspec =np.sqrt(jetscape[eloss][cent]['dconv']*jetscape[eloss][cent]['dconv']+\
                        jetscape[eloss][cent]['dbrem']*jetscape[eloss][cent]['dbrem'])
        spec = spec.tolist() 
        for ch in channels:
            new = rebinned_jf[cent][ch]['N'].tolist()
            spec = [v1 + v2 for (v1,v2) in zip(spec, new)]
        total_with_jet_medium[eloss][cent] = (np.array(x), np.array(dx), np.array(spec), np.array(dspec.tolist()))

total_without_jet_medium = {}
for eloss in elosses:
    total_without_jet_medium[eloss] = {}
    for cent in centralities:
        num_binary_coll = multiplicity[cent]
        x = jetscape[eloss][cent]['pT'].tolist()
        dx = jetscape[eloss][cent]['dpT'].tolist()
        spec = [0 for i in x]
        for ch in channels:
            new = rebinned_jf[cent][ch]['N'].tolist()
            spec = [v1 + v2 for (v1,v2) in zip(spec, new)]
        total_without_jet_medium[eloss][cent] = (np.array(x), np.array(dx), np.array(spec), np.array(dspec.tolist()))

fig, axes = plt.subplots(2,2,figsize=(16,12),
                        gridspec_kw={'height_ratios':(1,1),
                                     'top':0.945,'bottom':0.105,
                                     'left':0.12,'right':0.965,
                                     'hspace':0.05,'wspace':0.01},
                        sharey=True, sharex=True)



for icent, cent in enumerate(centralities):
    data = expdata[cent]
    for ieloss, eloss in enumerate(elosses): 
        spec_ax = axes[ieloss][icent]
        util.plot_expr_data_on_axis(spec_ax, data, marker='o')
        if ieloss == 0:
            spec_ax.text(0.1,0.2,f'{cent}\%', transform=spec_ax.transAxes)
        pT,_,spec_with_jet_medium,_ = total_with_jet_medium[eloss][cent]
        pT,_,spec_without_jet_medium,_ = total_without_jet_medium[eloss][cent]
        spec_ax.plot(pT, spec_with_jet_medium, color=ddicts.eloss_colours[eloss])
        spec_ax.plot(pT, spec_without_jet_medium, color='gray')


for ax in axes[1]:
    ax.set_yscale('log')

fig.supylabel(r'$\frac{1}{N_{\mathrm{evt}}\,2\pi\,p_T}\frac{\mathrm{d}N^{\gamma}}{\mathrm{d}p_T\mathrm{d}\eta}$, (GeV$^{-2}$)')
fig.supxlabel(r'$p_T$ (GeV)')

handles = [Line2D([],[],label='ALICE (2016)', marker='o',color='black')]
for eloss in elosses:
    handles.append(Line2D([],[],color=ddicts.eloss_colours[eloss],label=eloss.upper()))
handles.append(Line2D([],[],color='gray',label='Without Jet-Medium'))


axes[0][1].legend(loc='upper right', handles=handles, fontsize=25)
plt.show()