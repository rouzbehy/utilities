#!/usr/bin/env python3
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from COLORS import colors
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
    "axes.spines.top"  : False,
    "legend.frameon"   : False,
    "axes.labelsize" : 25}
plt.rcParams.update(my_rcParams)


tmpl = '~/Documents/research/hydro_plots/VISHNU/HydroProfiles/{name}-Avg-{cent}/event-0/photon_total_Spvn.dat'
collSystems = {'AuAu200':['00-10','10-20'],'PbPb2760':['00-05','05-10','10-20','20-30','30-40','40-50'], 'PbPb5020':['00-10','10-20','30-50']}

tmpl_JF = '/Users/rmyazdi/Documents/research/jetscape_project/v2/other_data/JF_Photon_Calc/Spectra/thermal_from_hydro/{name}/C{cent}/average_sp.dat'
#tmpl_JF = '/Users/rmyazdi/Documents/research/jetscape_project/v2/other_data/JF_MultiMessenger/{collsys}_{cent}_thermal.csv'

pbpb_5020_C0_5 = np.loadtxt(tmpl_JF.format(name='PbPb5020', cent='0-5'), unpack=True, delimiter=' ')
pbpb_5020_C5_10 = np.loadtxt(tmpl_JF.format(name='PbPb5020', cent='5-10'), unpack=True, delimiter=' ')

pbpb_5020_C0_10 = 0.5 * (pbpb_5020_C0_5[1] + pbpb_5020_C5_10[1])

fig, axes = plt.subplots(3, 1, sharex=True)

for iax, name in enumerate(collSystems):
    ax = axes[iax]
    legend = []
    ax.set_yscale('log')
    for icent, cent in enumerate(collSystems[name]):
        col = colors[icent]
        data = pd.read_csv(tmpl.format(name=name,cent=cent))
        data.columns = ['pT', 'sp', 'N/A']
        ax.plot(data['pT'], data['sp'], color=col)
        ## now try JF:
        try:
            if cent == '00-10' and name == 'PbPb5020':
                ax.plot(pbpb_5020_C0_5[0], pbpb_5020_C0_10, color=col)
            else:
                cent = cent if cent!='00-10' else '0-10'
                cent = cent if cent!='00-05' else '0-5'
                fname = tmpl_JF.format(name=name, cent=cent)
                jf_data = np.loadtxt(fname, unpack=True, delimiter=' ')
                ax.plot(jf_data[0], jf_data[1], color=col, linestyle='dotted')
            legend.append(Line2D([],[],color=col,label=f'{cent}%'))
        except:
            # fname = tmpl_JF.format(name=name, cent=cent)
            # if not os.path.exists(fname):
            #     print(f"{fname} doesn't exist")
            continue
    ax.legend(loc='best', handles=legend)

plt.show()        