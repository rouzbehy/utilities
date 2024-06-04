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
    conv =0.5*(df1['conv']+ df2['conv'])
    brem =0.5*(df1['brem']+ df2['brem']) 
    othr =0.5*(df1['othr']+ df2['othr']) 

    dconv =np.sqrt(df1['dconv']* df1['dconv']+ df2['dconv']* df2['dconv'])
    dbrem =np.sqrt(df1['dbrem']* df1['dbrem']+ df2['dbrem']* df2['dbrem'])  
    dothr =np.sqrt(df1['dothr']* df1['dothr']+ df2['dothr']* df2['dothr'])  
    return pd.DataFrame({'pT':x,'dpT':dx,'ptmin':xmin,'ptmax':xmax,
                          'conv':conv,'dconv':dconv,'brem':brem,'dbrem':dbrem,
                          'othr':othr,'dothr':dothr})

elosses = ['cujet', 'martini']

expdata = {}
# pT_lower_lim = {'00-20':0,'20-40':1.2}
# pT_upper_lim = {'00-20':16,'20-40':16}
pT_lower_lim = {'00-20':0,'20-40':0}
pT_upper_lim = {'00-20':0,'20-40':0}

## read jetscape:
elosses = ['martini','cujet']
centralities = ['00-05','05-10','10-20','20-30','30-40', '40-50']
jetscape = {}
fname_tmpl = '../jetscape_data/{eloss}/PbPb_2760/PbPb2760_{cent}_photon_spec_0.80.csv'
for eloss in elosses:
    jetscape[eloss] = {}
    for cent in centralities:
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

fig, axes = plt.subplots(3,2,figsize=(16,12),
                        gridspec_kw={'top':0.945,'bottom':0.105,
                                     'left':0.12,'right':0.965,
                                     'hspace':0.05,'wspace':0.1},
                        sharey='row', sharex=True)
axes = axes.flatten()

linestyles= {'brem':'solid', 'conv':'dashed'}
channels = ['brem', 'conv']
marker_style = {'brem':'v', 'conv':'P'}

for icent, cent in enumerate(centralities):
    ax = axes[icent]
    ax.text(0.7,0.7,f'{cent}\%', transform=ax.transAxes)

    for eloss in elosses:
        col = ddicts.eloss_colours[eloss]
        conv, brem = jetscape[eloss][cent]['conv'], jetscape[eloss][cent]['brem']
        dconv, dbrem = jetscape[eloss][cent]['dconv'], jetscape[eloss][cent]['dbrem']
        total = conv + brem
        dtotal = np.sqrt(dconv * dconv + dbrem * dbrem)
        conv_ratio = conv/total
        dconv_ratio = conv_ratio* np.sqrt( dtotal*dtotal / (total*total) + dconv*dconv/(conv*conv))
        brem_ratio = brem/total
        dbrem_ratio = brem_ratio*np.sqrt( dtotal*dtotal / (total*total) + dbrem*dbrem/(brem*brem))
        ptmin, ptmax = jetscape[eloss][cent]['ptmin'], jetscape[eloss][cent]['ptmax']

        pT, dpT = 0.5*(ptmin+ptmax), (ptmax-ptmin)
        mark = ddicts.eloss_marker[eloss]
        errorboxes_1 = [Rectangle((x-delx, y - dy), width=2*delx, height=2*dy)
                    for x, delx, y, dy in zip(pT, 0.5*dpT, conv_ratio, dconv_ratio)]
        pc1 = PatchCollection(errorboxes_1, facecolor=col, edgecolor=col, alpha=0.5)
        ax.add_collection(pc1)

        errorboxes_2 = [Rectangle((x-delx, y - dy), width=2*delx, height=2*dy)
                    for x, delx, y, dy in zip(pT, 0.5*dpT, brem_ratio, dbrem_ratio)]
        pc2 = PatchCollection(errorboxes_2, facecolor=col, edgecolor=col, alpha=0.5)
        ax.add_collection(pc2)



        ax.scatter(pT, conv_ratio, color=col, marker=marker_style['conv'])
        ax.scatter(pT, brem_ratio, color=col, marker=marker_style['brem'])

fig.supxlabel(r'$p_T$ (GeV)')
fig.supylabel(r'Channel/[Jet-Medium Total]')
fig.title('Pb-Pb @')
plt.show()