#!/usr/bin/env python3
import numpy as np
from pandas import read_csv, DataFrame
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from matplotlib import ticker
from COLORS import colors
import helper
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 30,
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

## Read in the charged hadron multiplicity file
ch_mult = read_csv('ch_mult.csv')
ch_mult.sort_values(by=['nch'], inplace=True)

def func(x, *params):
    return params[0] + params[1] * x

pcuts = [1, 2, 3, 4, 5]
for i in range(1,12):
    pcuts.append(5 + i*0.5)

alphas = {'martini':[],'cujet':[],'none':[]}
betas  = {'martini':[],'cujet':[],'none':[]}
results = {}
for pcut in pcuts:
    data = [[] for i in range(7)]
    for index, row in ch_mult.iterrows():
        (nuclei, energy, cent, nch) = row['nuclei'], row['energy'], row['centrality'], row['nch']
        martini = helper.construct_photon_spec(nuclei, energy, cent, 'martini', pcut)
        cujet = helper.construct_photon_spec(nuclei, energy, cent, 'cujet', pcut)
        no_jetmedium = helper.construct_photon_spec(nuclei, energy, cent, 'none', pcut)
        data[0].append(nch)
        for i in range(2):
            data[i+1].append(martini[i])
            data[i+3].append(cujet[i])
            data[i+5].append(no_jetmedium[i])


    ps_martini,_ = curve_fit(func,np.log(data[0]), np.log(data[1]), p0=(1,1))
    ps_martini,_ = curve_fit(func,np.log(data[0]), np.log(data[1]), p0=ps_martini)
    ps_cujet,_   = curve_fit(func,np.log(data[0]), np.log(data[3]), p0=ps_martini)
    ps_none,_    = curve_fit(func,np.log(data[0]), np.log(data[5]), p0=ps_martini)
    
    alphas['martini'].append(np.exp(ps_martini[0]))
    alphas['cujet']  .append(np.exp(ps_cujet[0]))
    alphas['none']   .append(np.exp(ps_none[0]))

    betas['martini'].append(ps_martini[1])
    betas['cujet']  .append(ps_cujet[1])
    betas['none']   .append(ps_none[1])
    
    results[pcut] = data

for item in zip(pcuts, betas['martini'], betas['cujet'], betas['none']):
    p, a1, a2, a3 = item
    print(f"{p:0.2f},{a1:0.2f},{a2:0.2f},{a3:0.2f}")
fig, axes = plt.subplots(2,1, figsize=(16,9), sharex=True)
handles = [Line2D([],[],label=eloss, color=c) for (eloss, c) in zip(['MATTER+MARTINI','MATTER+CUJET','No Jet-Medium'],
                                                                    [colors[0],colors[1],colors[2]])]
axes[0].plot(pcuts, betas['martini'] , color=colors[0])
axes[0].plot(pcuts, betas['cujet']   , color=colors[1])
axes[0].plot(pcuts, betas['none']    , color=colors[2])
axes[1].plot(pcuts, alphas['martini'], color=colors[0])
axes[1].plot(pcuts, alphas['cujet']  , color=colors[1])
axes[1].plot(pcuts, alphas['none']   , color=colors[2])
axes[0].legend(loc='upper left', handles=handles, bbox_to_anchor=(1.,0.3))
axes[1].set_yscale('log')
axes[0].set_ylim(bottom=0.1, top=3.2)
axes[1].set_xlabel(r'$p^{\mathrm{cut}}_T$ (GeV)')
axes[0].set_ylabel(r'$\alpha$')
axes[1].set_ylabel(r'$A$')

fig1, ax1 = plt.subplots(1,1,figsize=(16,9))

cols = {'martini':colors[0], 'cujet':colors[1], 'none':colors[2]}
markers = ['*','>','<', '+']
indices = {'martini':1, 'cujet':3, 'none':5}
xvals = np.linspace(240, 2800, 100)
linestyles = {1: 'solid', 3:'dashed', 5.5:'dotted'}#, 5.5:'dashdot'}
pcut_index = {1:0, 3:2, 5:4, 5.5:5}
for indx, item in enumerate(linestyles):
   handles.append(Line2D([],[],color='black', 
                   linestyle=linestyles[item], 
                   label=r'$p^{\mathrm{cut}}_T=$'+f'{item} GeV', 
                   marker=markers[indx],
                   markersize=15))

for ip, pcut in enumerate([1, 3, 5.5]):
   data = results[pcut]
   lstyle = linestyles[pcut]
   for eloss in cols:
       color = cols[eloss]
       marker = markers[ip]
       parms = [np.log(alphas[eloss][pcut_index[pcut]]), betas[eloss][pcut_index[pcut]]]
       ax1.scatter(data[0], data[indices[eloss]] , color=color,  marker=marker)
       ax1.plot(xvals, [np.exp(func(np.log(x), *parms)) for x in xvals], color=color, linestyle=lstyle)
       
ax1.legend(loc='upper right', handles=handles)#, bbox_to_anchor=(1.,0.8))
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel(r'$\frac{dN^{\gamma}}{d\eta}$')
ax1.set_xlabel(r'$\frac{dN^{h^{\pm}}}{d\eta}$')

plt.show()
