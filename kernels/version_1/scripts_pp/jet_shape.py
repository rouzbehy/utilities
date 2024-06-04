#!/usr/bin/env python3
from re import A
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import trapezoid, simpson
import util 
import dictionaries as ddicts

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


pT_lower_lim = 20
pT_upper_lim = 96
inel_Xsec= 62.03948 

## read in the experimental results
data_fname = "/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/exp_data/jet_shape/CMS_ppsmeared_2760_Jet_Shape_70-100_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
CMS_DATA = np.loadtxt(data_fname,comments='#',unpack=True, delimiter='\t')

CMS_DATA = pd.DataFrame({'x':CMS_DATA[0],'y':CMS_DATA[1],'dx':CMS_DATA[2],'dy':CMS_DATA[3]})
## Plot the data:
fig, (ax,ax1) = plt.subplots(2, 1,
                            gridspec_kw={'height_ratios':(3,1),'top':0.985,
                                'bottom':0.1, 'left':0.1,
                                'right':0.95, 'hspace':0.02,
                                'wspace':0.18}, figsize=(16,9), sharex=True)

color="#d64045"
path = '../martini_results/pp/jet_shape.csv'
calcs = pd.read_csv(path,comment='#')
errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                for x, delx, y, yerrlow, yerrhigh in
                zip(CMS_DATA["x"], CMS_DATA['dx'], CMS_DATA["y"], CMS_DATA["dy"], CMS_DATA["dy"])] 
# Create patch collection with specified colour/alpha
pc = PatchCollection(errorboxes, facecolor='grey', alpha=1)
ax.add_collection(pc)
ax.scatter(CMS_DATA['x'],CMS_DATA['y'],color='black',marker='s',s=30)

nevt_N, njet_N = 1, 1
with open(path, 'r') as f:
    line = f.readline()
    line = line.split(' ')[-1]
    njet = float(line)

dat = calcs[calcs['rmax'] < 0.31]
delta_r = dat['rmax'] - dat['rmin']
r = 0.5*(dat['rmax'] + dat['rmin'])
rho  = dat['N']  /(njet_N )
drho = dat['dN'] /(njet_N )
norm = sum(rho.to_list())#
rho_normed = rho  /(delta_r * norm)
drho_normed = drho/(delta_r * norm)

errorboxes = [Rectangle((x-delx, y - dy), width=2*delx, height=2*dy)
            for x, delx, y, dy in zip(r, 0.5*delta_r, rho_normed, drho_normed)]
pc = PatchCollection(errorboxes, facecolor=color, alpha=1.)
ax.add_collection(pc)

ratio = CMS_DATA['y']/rho_normed
dratio = ratio*np.sqrt(CMS_DATA["dy"]*CMS_DATA["dy"]/(CMS_DATA["y"]*CMS_DATA["y"]) + drho_normed*drho_normed/(rho_normed*rho_normed))
errorboxes_rat = [Rectangle((x-delx, y - yerrlow), 
            width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
            for x, delx, y, yerrlow, yerrhigh in
            zip(r, 0.5*delta_r, ratio, dratio, -1*dratio)]
pc_rat = PatchCollection(errorboxes_rat, facecolor=color, alpha=1)
ax1.add_collection(pc_rat)

labels = [Line2D([],[],color=color, label='Theory'),\
          Line2D([],[],color='black', marker='s', markersize=10, label='CMS (2014)')]

ax.text(0.02,0.1,s=r'$0.3<|\eta|<2.0$'+'\n'+r'$100$ GeV $< p^{\mathrm{jet}}_T$', transform=ax.transAxes)
ax.text(0.02,0.25,s=r'$R=0.3$, Anti-$k_{T}$', transform=ax.transAxes)
ax.text(0.02,0.35,s=r'p-p, $\sqrt{s}=2.76$ TeV', fontsize=30, transform=ax.transAxes)

ax.legend(handles=labels, loc='upper right')
ax.set_yscale('log')
ax.set_ylabel(r'$\rho(r)$')
ax1.set_xlabel(r'$r$')
ax1.set_ylim(top=1.36,bottom=0.65)
ax1.axhline(1,linestyle='dotted',color='black')
ax1.set_ylabel('Data/Theory')
plt.show()
