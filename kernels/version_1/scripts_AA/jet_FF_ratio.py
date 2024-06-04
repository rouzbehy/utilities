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
import dictionaries as my_dicts
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

def get_martini_results(rateset, cent):
    fname_pt = f'../martini_results/final_PbPb_2p76/rset_{rateset}/cent_{cent}/jet_FF_pT.csv'
    fname_z  = f'../martini_results/final_PbPb_2p76/rset_{rateset}/cent_{cent}/jet_FF_z.csv'
    tmp_DPT  = pd.read_csv(fname_pt,comment='#')
    tmp_DZ   = pd.read_csv(fname_z,comment='#')
    tmp_DPT  = tmp_DPT[tmp_DPT['pTmax'] < maxpT]
    tmp_DPT  = tmp_DPT[tmp_DPT['pTmin'] > minpT]
    tmp_DZ   = tmp_DZ[tmp_DZ['zmax'] < zmax]
    tmp_DZ   = tmp_DZ[tmp_DZ['zmin'] > zmin]
    njets_pp, nevts_pp = 0 , 0
    with open(fname_pt.format(r=rset),'r') as f:
        line = f.readline()
        line = line.split(' ')
        njets_pp = float(line[1])
        nevts_pp = float(line[-1])
    
    dpT = tmp_DPT['pTmax'] - tmp_DPT['pTmin']
    pT = 0.5*(tmp_DPT['pTmax'] + tmp_DPT['pTmin'])
    scale = njets_pp
    tmpy  = tmp_DPT['N'] /(dpT*scale)
    dtmpy = tmp_DPT['dN']/(dpT*scale)
    

    dz = tmp_DZ['zmax'] - tmp_DZ['zmin']
    z = 0.5*(tmp_DZ['zmax'] + tmp_DZ['zmin'])
    tmpyz  = tmp_DZ['N'] /(dz*scale)
    dtmpyz = tmp_DZ['dN']/(dz*scale)

    return pT, dpT, tmpy, dtmpy, z, dz, tmpyz, dtmpyz

## read in jetscape calculation
minpT, maxpT = 0.1, 120
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

dpT = pp_DPT['pTmax'] - pp_DPT['pTmin']
pT = 0.5*(pp_DPT['pTmax'] + pp_DPT['pTmin'])
scale = njets_pp
y  = pp_DPT['N'] /(dpT*scale)
dy = pp_DPT['dN']/(dpT*scale)
pp_dpt = pd.DataFrame({'pT':pT,'dpT':dpT, 'y':y, 'dy':dy})

dz = pp_DZ['zmax'] - pp_DZ['zmin']
z = 0.5*(pp_DZ['zmax'] + pp_DZ['zmin'])
y  = pp_DZ['N'] /(dz*scale)
dy = pp_DZ['dN']/(dz*scale)
pp_dz = pd.DataFrame({'z':z,'dz':dz,'y':y,'dy':dy})

## Read in the AA data:
AA_dpt = {}
AA_dz = {}
from COLORS import rate_set_colors as rate_colours
# = my_dicts.rate_colours
rate_names = {1:'LO', 2:'NLO',3:'NP'}

for rset in rate_names:
    #color = rate_colours[rate_names[rset]]
    pT1, dpT1, tmpy1, dtmpy1, z1, dz1, tmpyz1, dtmpyz1 = get_martini_results(rset, '0_5')
    pT2, dpT2, tmpy2, dtmpy2, z2, dz2, tmpyz2, dtmpyz2 = get_martini_results(rset, '5_10')
    tmpy = 0.5*(tmpy1+tmpy2)
    dtmpy = 0.5*np.sqrt(dtmpy1**2 + dtmpy2**2)
    tmpyz = 0.5*(tmpyz1+tmpyz2)
    dtmpyz = 0.5*np.sqrt(dtmpyz1**2 + dtmpyz1**2)
    AA_dpt[rset] = pd.DataFrame({'pT':pT,'dpT':dpT, 'y':tmpy, 'dy':dtmpy})
    AA_dz[rset] = pd.DataFrame({'z':z,'dz':dz,'y':tmpyz,'dy':dtmpyz})

exp_loc = '../../../exp_data/jet_fragmentation/ATLAS/'
## read in ATLAS results
data_pT = pd.read_csv(exp_loc+'HEPData-ins1511869-v1-csv/Table9.csv',comment='#').rename(columns=ddicts.colnames_ATLAS_FF_DPT)
data_z = pd.read_csv(exp_loc+"HEPData-ins1511869-v1-csv/Table25.csv",comment='#').rename(columns=ddicts.colnames_ATLAS_FF_DZ)
## make the plots:
fig, axes = plt.subplots(1, 2, figsize=(16,9), sharex='col')

ax = axes[0]
util.plot_expr_data_on_axis(ax, data_pT, 's')
## cycle through the rate sets:
yp, dyp = pp_dpt['y'], pp_dpt['dy']
rate_names = {1:'LO', 2:'NLO',3:'NP'}

for rate_set in rate_names:
    color = rate_colours[rate_names[rate_set]]
    aa = AA_dpt[rate_set]
    ya, dya = aa['y'], aa['dy']
    pT = aa['pT']
    dpT = aa['dpT']
    ratio = ya/yp
    dratio = np.array([r * np.sqrt(erra*erra/(va*va) + errp*errp/(vp*vp)) \
                        for (r, erra,va, errp, vp) in zip(ratio, dya,ya,dyp,yp)])

    # errorboxes = [Rectangle((x-delx, yval - yerrval), width=2*delx, height=2*yerrval) 
    #             for x, delx, yval, yerrval in zip(pT, 0.5*dpT, ratio,dratio) ]

    # pc = PatchCollection(errorboxes, facecolor=color, alpha=0.2)
    # ax.add_collection(pc)
    # ax.scatter(pT, ratio, color=color, marker='p',s=50)
    ax.plot(pT, ratio, color=color)
    ax.fill_between(pT, ratio-dratio, ratio+dratio, color=color, alpha=0.2)
ax = axes[1]
util.plot_expr_data_on_axis(ax, data_z, 's')
## cycle through the rate sets:
yp, dyp = pp_dz['y'], pp_dz['dy']
for rate_set in rate_names:
    color = rate_colours[rate_names[rate_set]]
    aa = AA_dz[rate_set]
    ya, dya = aa['y'], aa['dy']
    z = aa['z']
    dz = aa['dz']
    ratio = ya/yp
    dratio = np.array([r * np.sqrt(erra*erra/(va*va) + errp*errp/(vp*vp)) for (r, erra,va, errp, vp) in zip(ratio, dya,ya,dyp,yp)])

    # errorboxes = [Rectangle((x-delx, yval - yerrval), width=2*delx, height=2*yerrval) 
    #            for x, delx, yval, yerrval in zip(z, 0.5*dz, ratio,dratio) ]

    # pc = PatchCollection(errorboxes, facecolor=color, alpha=0.2)
    # ax.add_collection(pc)
    #ax.scatter(z, ratio, color=color, marker='p',s=50)
    ax.plot(z, ratio, color=color)
    ax.fill_between(z, ratio-dratio, ratio+dratio, color=color, alpha=0.2)

labels = [Line2D([],[],color=c,label=l) for l, c in rate_colours.items()]
labels.append(Line2D([],[],color='black', marker='s', label='ATLAS $0$-$10\%$'))
axes[1].legend(loc='upper left', handles=labels)
axes[0].set_xscale('log')
axes[0].set_ylabel(r'$R_{D(p_T)}$')
axes[1].set_ylabel(r'$R_{D(z)}$')
axes[0].set_xlabel(r'$p^{h^{\pm}}_T$ (GeV)')
axes[1].set_xlabel(r'$z$')
axes[0].text(0.28,0.85,s=r'Pb-Pb, $\sqrt{s}=2.76$ ATeV', transform=axes[0].transAxes)
axes[0].text(0.28,0.80,s=r'$100< p^{\mathrm{jet}}_T < 398$ GeV', transform=axes[0].transAxes)
axes[0].text(0.28,0.75,s=r'$0$-$10\%$', transform=axes[0].transAxes)
plt.show()
