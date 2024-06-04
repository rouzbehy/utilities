#!/usr/bin/env python3
"""
    Plot the computed spectra for jets 
"""
## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
## my custom modules
import util 
import jetDicts as ddicts

#mpl.use('Qt5Agg')
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

## read in jetscape calculation
maxpT = 120
minpT = 0.1
zmin, zmax = 0.005,2
## Read in the pp data:
pp_DPT     = pd.read_csv('../martini_results/pp/jet_FF_pT.csv',comment='#')
pp_DZ      = pd.read_csv('../martini_results/pp/jet_FF_z.csv',comment='#')

pp_DPT = pp_DPT[pp_DPT['pTmax'] < maxpT]
pp_DPT = pp_DPT[pp_DPT['pTmin'] > minpT]
pp_DZ = pp_DZ[pp_DZ['zmax'] < zmax]
pp_DZ = pp_DZ[pp_DZ['zmin'] > zmin]
njets_pp, nevts_pp = 0 , 0
with open('../martini_results/pp/jet_FF_pT.csv','r') as f:
    line = f.readline()
    line = line.split(' ')
    njets_pp = float(line[1])
    nevts_pp = float(line[-1])

## read in ATLAS results
data_pT = pd.read_csv('../../../exp_data/jet_fragmentation/ATLAS/Table_pT_100-398_pp.csv',comment='#').rename(columns=ddicts.colnames_ATLAS_FF_pT)
data_z = pd.read_csv("../../../exp_data/jet_fragmentation/ATLAS/Table_z_100-398_pp.csv",comment='#').rename(columns=ddicts.colnames_ATLAS_FF_z)

## make the plots:
fig, axes = plt.subplots(2, 2, figsize=(16,9), sharex='col',
                        gridspec_kw={'height_ratios':(3,1)})

ax = axes[0][0]## dNch/dpT
util.plot_expr_data_on_axis(ax, data_pT, 's')
dpT = pp_DPT['pTmax'] - pp_DPT['pTmin']
pT = 0.5*(pp_DPT['pTmax'] + pp_DPT['pTmin'])
scale = njets_pp
y  = pp_DPT['N'] /(dpT*scale)
dy = pp_DPT['dN']/(dpT*scale)
errorboxes = [Rectangle((x-delx, yval - yerrval), width=2*delx, height=yerrval) 
                for x, delx, yval, yerrval in zip(pT, 0.5*dpT, y,dy) ]
pc = PatchCollection(errorboxes, facecolor="#d64045", alpha=1)
ax.add_collection(pc)
ax.scatter(pT, y, color="#d64045", marker='p',s=50)

## take ratio:
datay = data_pT['y'].tolist()
ratio = np.array([v1/v2 for (v1,v2) in zip(datay,y)])

dratio_stat_pos = ratio*data_pT['dy_stat+']/data_pT['y']
dratio_stat_neg = ratio*data_pT['dy_stat-']/data_pT['y']
dratio_syst_pos = ratio*data_pT['dy_syst+']/data_pT['y']
dratio_syst_neg = ratio*data_pT['dy_syst-']/data_pT['y']
ratio_frame = pd.DataFrame({'x':pT.tolist(),
                            'xlow':pp_DPT['pTmin'].to_list(),
                            'xhigh':pp_DPT['pTmax'].tolist(),
                            'y':ratio,
                            'dy_stat+':dratio_stat_pos.tolist(),
                            'dy_stat-':dratio_stat_neg.tolist(),
                            'dy_syst+':dratio_syst_pos.tolist(),
                            'dy_syst-':dratio_syst_neg.tolist()})
ax = axes[1][0]
util.plot_expr_data_on_axis(ax, ratio_frame, marker='p',color="#d64045",face="#d64045")

axes[0][0].set_xscale('log')
axes[0][0].set_yscale('log')

ax = axes[0][1]## dNch/dz
util.plot_expr_data_on_axis(ax, data_z, 's')

dz = pp_DZ['zmax'] - pp_DZ['zmin']
z = 0.5*(pp_DZ['zmax'] + pp_DZ['zmin'])
y  = pp_DZ['N'] /(dz*scale)
dy = pp_DZ['dN']/(dz*scale)
errorboxes = [Rectangle((x-delx, yval - yerrval), width=2*delx, height=yerrval) 
                for x, delx, yval, yerrval in zip(z, 0.5*dz, y,dy) ]
pc = PatchCollection(errorboxes, facecolor="#d64045", alpha=1)
ax.add_collection(pc)
ax.scatter(z, y, color="#d64045", marker='p',s=50)

## take ratio:
datay = data_z['y'].tolist()
ratio = np.array([v1/v2 for (v1,v2) in zip(datay,y)])

dratio_stat_pos = ratio*data_z['dy_stat+']/data_z['y']
dratio_stat_neg = ratio*data_z['dy_stat-']/data_z['y']
dratio_syst_pos = ratio*data_z['dy_syst+']/data_z['y']
dratio_syst_neg = ratio*data_z['dy_syst-']/data_z['y']
ratio_frame = pd.DataFrame({'x':z.tolist(),
                            'xlow':pp_DZ['zmin'].to_list(),
                            'xhigh':pp_DZ['zmax'].tolist(),
                            'y':ratio,
                            'dy_stat+':dratio_stat_pos.tolist(),
                            'dy_stat-':dratio_stat_neg.tolist(),
                            'dy_syst+':dratio_syst_pos.tolist(),
                            'dy_syst-':dratio_syst_neg.tolist()})
ax = axes[1][1]
util.plot_expr_data_on_axis(ax, ratio_frame, marker='p',color="#d64045",face="#d64045")
axes[0][1].set_xscale('log')
axes[0][1].set_yscale('log')

for ax in axes[1]:
    ax.set_ylim(bottom=0.7, top=1.1)
    ax.axhline(1, color='black', linestyle='dotted')

axes[0][0].set_ylabel(r'$\frac{1}{N_{\mathrm{jet}}} \frac{\mathrm{d}N_{\mathrm{ch}}}{\mathrm{d}p_T}$ (GeV$^{-1}$)')
axes[0][1].set_ylabel(r'$\frac{1}{N_{\mathrm{jet}}} \frac{\mathrm{d}N_{\mathrm{ch}}}{\mathrm{d}z}$')

axes[1][0].set_xlabel(r'$p^{h^{\pm}}_T$ (GeV)')
axes[1][1].set_xlabel(r'$z$')
axes[1][0].set_ylabel("Data/Theory")

labels = [Line2D([],[],label=r'ATLAS $|y|<2.1$(2017)', color='black',marker='s'),
          Line2D([],[],label=r'Theory $|\eta|<2$', color='#d64045', marker='p')]
axes[0][1].legend(loc='lower left', handles=labels)

axes[0][0].text(0.05, 0.2, r'$100< p^{\mathrm{jet}}_T < 398$ GeV', transform=axes[0][0].transAxes)
axes[0][0].text(0.05,0.3,s=r'p-p, $\sqrt{s}=2.76$ TeV', transform=axes[0][0].transAxes)

plt.show()
