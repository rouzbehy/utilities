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
def combine_centralities(df1, df2):
    xmin = df1['ptmin']
    xmax = df1['ptmax']
    x = 0.5*(xmin+xmax)
    dx = xmax - xmin
    conv =0.5*(df1['conv']+ df2['conv'])/(x*dx*inel_Xsec*2*np.pi*2*0.8) 
    brem =0.5*(df1['brem']+ df2['brem'])/(x*dx*inel_Xsec*2*np.pi*2*0.8) 
    othr =0.5*(df1['othr']+ df2['othr'])/(x*dx*inel_Xsec*2*np.pi*2*0.8) 

    dconv =np.sqrt(df1['dconv']* df1['dconv']+ df2['dconv']* df2['dconv'])/(x*dx*inel_Xsec*2*np.pi*2*0.8) 
    dbrem =np.sqrt(df1['dbrem']* df1['dbrem']+ df2['dbrem']* df2['dbrem'])/(x*dx*inel_Xsec*2*np.pi*2*0.8)  
    dothr =np.sqrt(df1['dothr']* df1['dothr']+ df2['dothr']* df2['dothr'])/(x*dx*inel_Xsec*2*np.pi*2*0.8)  
    return pd.DataFrame({'pT':x,'dpT':dx,'ptmin':xmin,'ptmax':xmax,
                          'conv':conv,'dconv':dconv,'brem':brem,'dbrem':dbrem,
                          'othr':othr,'dothr':dothr})

elosses = ['cujet', 'martini']


## read in experimental data
exp_loc = "../../expt/PbPb_2p76/photons/"
centralities = {'00-20':exp_loc+'HEPData-ins1394677-v1-Table_1.csv',
                '20-40':exp_loc+'HEPData-ins1394677-v1-Table_2.csv'}

expdata = {}
# pT_lower_lim = {'00-20':0,'20-40':1.2}
# pT_upper_lim = {'00-20':16,'20-40':16}
pT_lower_lim = {'00-20':0,'20-40':0}
pT_upper_lim = {'00-20':0,'20-40':0}
for cent in centralities:
    tmp = pd.read_csv(centralities[cent],comment='#')
    min_pT = 2.5#min(tmp['x'])
    max_pT = max(tmp['x'])
    pT_lower_lim[cent] = min_pT
    pT_upper_lim[cent]=max_pT
    tmp = tmp[tmp['x'].between(pT_lower_lim[cent],pT_upper_lim[cent])]
    expdata[cent] = tmp
## read jetscape:
elosses = ['martini','cujet']
jetscape_cents = ['00-05','05-10','10-20','20-30','30-40']
jetscape = {}
fname_tmpl = '../jetscape_data/{eloss}/PbPb_2760/PbPb2760_{cent}_photon_spec_0.80.csv'
for eloss in elosses:
    jetscape[eloss] = {}
    for cent in jetscape_cents:
        fname = fname_tmpl.format(eloss=eloss,cent=cent)
        factor = oversampling_factors[f'{eloss}-{cent}']
        tmp = pd.read_csv(fname,comment='#')
        tmp['conv'] /= factor
        tmp['dconv'] /= factor
        tmp['brem'] /= factor
        tmp['dbrem'] /= factor
        jetscape[eloss][cent] = tmp
    ## construct the 0-20 and 20-40 centralities from the above
    tmp1 = combine_centralities(jetscape[eloss]['00-05'],jetscape[eloss]['05-10'])
    tmp2 = combine_centralities(tmp1, jetscape[eloss]['10-20'])
    jetscape[eloss]['00-20'] = tmp2[tmp2['pT'].between(pT_lower_lim['00-20'],pT_upper_lim['00-20'])]
    print(jetscape[eloss]['00-20']['pT'])
    tmp3 = combine_centralities(jetscape[eloss]['20-30'],jetscape[eloss]['30-40'])
    jetscape[eloss]['20-40'] = tmp3[tmp3['pT'].between(pT_lower_lim['20-40'],pT_upper_lim['20-40'])]
    print(jetscape[eloss]['20-40']['pT'])

## Read in the prompt and thermal calculations of JF Paquet
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

## Rebin the JF results to match my bins (and the experimental data)
## as a sample. The x axis is shared among all of them:
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
        rebinned_jf[cent][ch] = new_master

## now form the total spectra
total_v1 = {}

for eloss in elosses:
    total_v1[eloss] = {}
    for cent in centralities:
        num_binary_coll = multiplicity[cent]
        x = jetscape[eloss][cent]['pT']
        dx = jetscape[eloss][cent]['dpT']
        spec = num_binary_coll*(jetscape[eloss][cent]['conv']+ jetscape[eloss][cent]['brem'])
        dspec =num_binary_coll*( np.sqrt(jetscape[eloss][cent]['dconv']*jetscape[eloss][cent]['dconv']+\
                        jetscape[eloss][cent]['dbrem']*jetscape[eloss][cent]['dbrem']))
        spec = spec.tolist() 
        for ch in channels:
            new = rebinned_jf[cent][ch]['N'].tolist()
            spec = [v1 + v2 for (v1,v2) in zip(spec, new)]
        total_v1[eloss][cent] = (x, dx, np.array(spec), np.array(dspec.tolist()))


## Now look at a channel by channel breakdown
## use total_v1

ratios_v1 = {}
jet_medium_photons = {}
for eloss in elosses:
    ratios_v1[eloss] = {}
    jet_medium_photons[eloss]= {}
    for cent in centralities:
        data = expdata[cent]
        pT, dpT, spec, dspec = total_v1[eloss][cent]
        num_binary_coll = multiplicity[cent]
        x = jetscape[eloss][cent]['pT']
        dx = jetscape[eloss][cent]['dpT']
        ## now make individual ratios:
        jet_medium_sources = num_binary_coll*(jetscape[eloss][cent]['conv']+jetscape[eloss][cent]['brem'])
        print(jet_medium_sources)
        djet_medium_sources = num_binary_coll*np.sqrt(jetscape[eloss][cent]['dconv']*jetscape[eloss][cent]['dconv']+\
                                      jetscape[eloss][cent]['dbrem']*jetscape[eloss][cent]['dbrem'])
        spec = data['y']

        jet_medium_ratio = np.array([v1/v2 for v1,v2 in zip(jet_medium_sources,spec)])
        djet_medium_ratio = djet_medium_sources/spec
        prmt_ratio = rebinned_jf[cent]['prompt']['N']/spec
        thrm_ratio = rebinned_jf[cent]['thermal']['N']/spec
        preq_ratio = rebinned_jf[cent]['preEq']['N']/spec
        jet_medium_photons[eloss][cent] = {'pT':pT,'dpT':dpT,
                                           'N':jet_medium_sources,
                                           'dN':djet_medium_sources}
        ratios_v1[eloss][cent] = {'jet_medium':jet_medium_ratio , 
                                  'prompt':prmt_ratio ,
                                  'thermal':thrm_ratio,
                                  'preEq':preq_ratio   }

fig, axes = plt.subplots(2,2,figsize=(16,12),
                        gridspec_kw={'height_ratios':(3,1),
                                     'top':0.945,'bottom':0.105,
                                     'left':0.12,'right':0.965,
                                     'hspace':0.05,'wspace':0.01},
                        sharey='row', sharex=True)

for icent, cent in enumerate(expdata.keys()):
    spec_ax = axes[0][icent]
    ratio_ax = axes[1][icent]
    data = expdata[cent]
    util.plot_expr_data_on_axis(spec_ax, data, marker='o')
    spec_ax.text(0.1,0.1,f'{cent}\%', transform=spec_ax.transAxes)
    #pT, dpT, _,_ = total_v1['martini'][cent]
    pT = data['x']
    for channel in channels:
        col = channel_colors[channel]
        lstyle = channel_linestyles[channel]
        spec = rebinned_jf[cent][channel]
        spec_ax.plot(pT, spec['N'], color=col,linestyle=lstyle)
        ratios = ratios_v1['cujet'][cent]
        for e in ratios[channel].tolist():
            print(e)
        ratio_ax.plot(pT, ratios[channel], color=col,linestyle=lstyle) 
        
    for eloss in elosses:
        key = f"jet_medium-{eloss}"
        col = channel_colors[key]
        lstyle = channel_linestyles['jet_medium']
        spec = jet_medium_photons[eloss][cent]
        ratio= ratios_v1[eloss][cent]['jet_medium']
        mark='s'
        spec_ax.plot(pT, spec['N'], color=col,linestyle=lstyle)
        ratio_ax.plot(pT, ratio, color=col,linestyle=lstyle) 

axes[0][0].set_yscale('log')
axes[0][0].set_ylabel(r'$\frac{1}{N_{\mathrm{evt}}\,2\pi\,p_T}\frac{\mathrm{d}N^{\gamma}}{\mathrm{d}p_T\mathrm{d}\eta}$, (GeV$^{-2}$)')
axes[1][0].set_ylabel(r'Data/Theory')
handles = [Line2D([],[],label='ALICE (2016)', marker='o',color='black')]

for eloss in elosses:
    handles.append(Line2D([],[],color=ddicts.eloss_colours[eloss],label="Jet-Medium: "+eloss.upper()))
jf_handles = []
for ch in channels:
    jf_handles.append(Line2D([],[],color=channel_colors[ch], linestyle=channel_linestyles[ch],label=nice_channels[ch]))

for ax in axes[1]:
    ax.axhline(1, linestyle='dotted',color='black')
    ax.set_xlabel(r'$p_T$ (GeV)')
axes[0][1].legend(loc='upper right', handles=handles, fontsize=20)
axes[0][0].legend(loc='upper right', handles=jf_handles, fontsize=20)
axes[0][1].text(0.7, 0.55, r'$|\eta|<0.8$',transform=axes[0][1].transAxes)
plt.show()