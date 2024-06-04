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


elosses = ['cujet', 'martini']


## read in experimental data
exp_loc = "../../expt/PbPb_2p76/photons/"
centralities = {'00-20':exp_loc+'HEPData-ins1394677-v1-Table_1.csv',
                '20-40':exp_loc+'HEPData-ins1394677-v1-Table_2.csv'}

expdata = {}
pT_lower_lim = {'00-20':0,'20-40':0}
pT_upper_lim = {'00-20':0,'20-40':0}
for cent in centralities:
    tmp = pd.read_csv(centralities[cent],comment='#')
    min_pT = min(tmp['x'])
    max_pT = max(tmp['x'])
    pT_lower_lim[cent]= 2.5#min_pT
    pT_upper_lim[cent]= max_pT
    tmp = tmp[tmp['x'].between(pT_lower_lim[cent],pT_upper_lim[cent])]
    expdata[cent] = tmp
## read jetscape:
elosses = ['martini','cujet']
jetscape_cents = ['00-05','05-10','10-20','20-30','30-40']
jetscape = {}
fname_tmpl = '../jetscape_data/{eloss}/PbPb_2760/PbPb2760_{cent}_photon_spec_0.80.csv'
inel_Xsec = 62.03948 
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

## my old Pythia prompts: 
prompt = pd.read_csv("/Users/rmyazdi/Documents/research/jetscape_project/v2/other_data/my_prompt_nPDF.csv", comment='#')
## rebin this prompt calculation to match the rest:
rebinned_prompt = {}
f = InterpolatedUnivariateSpline (np.log(prompt['pT']),np.log(prompt['N'])) 
df = InterpolatedUnivariateSpline(np.log(prompt['pT']),np.log(prompt['dN']))
for cent in centralities:
    pTmin, pTmax = jetscape['martini'][cent]['ptmin'].tolist(), jetscape['martini'][cent]['ptmax'].tolist()
    x, dx, xmins, xmaxes, y, dy = [],[],[],[],[],[]
    for (xlow, xhigh) in zip(pTmin, pTmax):
        xval = 0.5*(xhigh + xlow)
        dxval = xhigh - xlow
        yval  = np.exp(f(np.log(xval)))
        dyval = np.exp(df(np.log(xval)))
        x.append(xval)
        xmins.append(xlow)
        xmaxes.append(xhigh)
        dx.append(dxval)
        y.append(yval)
        dy.append(dyval)
    new_master = pd.DataFrame({'ptmin':xmins,'ptmax':xmaxes,'pT':x,'dpT':dx,'N':y,'dN':dy})
    rebinned_prompt[cent] = new_master

total_v3 = {}
for eloss in elosses:
    total_v3[eloss] = {}
    for cent in centralities:
        nbin = multiplicity[cent]
        x = jetscape[eloss][cent]['pT']
        dx = jetscape[eloss][cent]['dpT']
        conv = jetscape[eloss][cent]['conv'].tolist()
        brem = jetscape[eloss][cent]['brem'].tolist()
        prmpt = nbin*np.array(rebinned_prompt[cent]['N'])
        spec = np.array([sum(e) for e in zip(conv,brem, prmpt)])

        dconv = jetscape[eloss][cent]['dconv'].tolist() 
        dbrem = jetscape[eloss][cent]['dbrem'].tolist() 
        dprmp = nbin*np.array(rebinned_prompt[cent]['dN'].tolist())
        dspec = np.sqrt([sum(map(lambda x:x*x, e)) for e in zip(dconv, dbrem, dprmp)])
        for ch in ['thermal','preEq']:
            new = rebinned_jf[cent][ch]['N'].tolist()
            spec = [v1 + v2 for (v1,v2) in zip(spec, new)]
        total_v3[eloss][cent] = (x, dx, np.array(spec), np.array(dspec.tolist()))

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
    for eloss in elosses:
        pT, dpT, spec, dspec = total_v3[eloss][cent]
        col = ddicts.eloss_colours[eloss]
        mark = ddicts.eloss_marker[eloss]
        errorboxes_1 = [Rectangle((x-delx, y - dy), width=2*delx, height=2*dy)
                    for x, delx, y, dy in zip(pT, 0.5*dpT, spec, dspec)]
        pc1 = PatchCollection(errorboxes_1, facecolor=col, edgecolor=col, alpha=0.5)
        spec_ax.add_collection(pc1)
        spec_ax.scatter(pT, spec, color=col,marker=mark)
        ## take ratio
        ratio = data['y']/spec
        ## error on ratio: ignore the statistical uncertaintiy 
        # of the jetscape results: it's too small to matter.
        dratio_stat_plus = ratio*(data['dy_stat+']/data['y'])
        dratio_stat_minus = ratio*(data['dy_stat+']/data['y'])
        dratio_syst_plus = ratio*(data['dy_syst+']/data['y'])
        dratio_syst_minus = ratio*(data['dy_syst-']/data['y'])

        errorboxes_2 = [Rectangle((x-delx, y - abs(dy_minus)), width=2*delx, height=(abs(dy_minus)+abs(dy_plus)))
                    for x, delx, y, dy_minus, dy_plus in zip(pT, 0.5*dpT, ratio, dratio_syst_minus, dratio_syst_plus)]
        pc2 = PatchCollection(errorboxes_2, facecolor=col, edgecolor=col, alpha=0.65)
        ratio_ax.add_collection(pc2)
        ratio_ax.errorbar(pT, ratio, xerr=0.5*dpT,yerr=[dratio_stat_minus, dratio_syst_plus], color=col,marker=mark,linewidth=1,fmt='none')
        ratio_ax.scatter(pT, ratio, color=col,marker=mark)
handles = [Line2D([],[],label='ALICE (2016)', marker='o',color='black')]
for eloss in elosses:
    handles.append(Line2D([],[],color=ddicts.eloss_colours[eloss],label=eloss.upper(),marker=ddicts.eloss_marker[eloss]))

for ax in axes[1]:
    ax.axhline(1,color='black', linestyle='dotted')

axes[0][0].set_yscale('log')
axes[0][0].set_ylabel(r'$\frac{1}{N_{\mathrm{evt}}\,2\pi\,p_T}\frac{\mathrm{d}N^{\gamma}}{\mathrm{d}p_T\mathrm{d}\eta}$, (GeV$^{-2}$)')
axes[1][0].set_ylabel(r'Data/Theory')
fig.supxlabel(r'$p_T$ (GeV)')
axes[0][1].legend(loc='upper right', handles=handles)
axes[0][1].text(0.7, 0.55, r'$|\eta|<0.8$',transform=axes[0][1].transAxes)
plt.show()