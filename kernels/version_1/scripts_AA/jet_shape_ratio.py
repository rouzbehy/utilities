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
# from scipy.interpolate import InterpolatedUnivariateSpline
# from scipy.integrate import trapezoid, simpson
import util 
import dictionaries as my_dicts
# my_rcParams={
#     "text.usetex": True,
#     "font.family": "Georgia",
#     "font.size": 25,
#     "lines.linewidth": 2,
#     "xtick.direction": "in",
#     "ytick.direction": "in",
#     "xtick.minor.visible": True,
#     "ytick.minor.visible": True,
#     "xtick.major.size" : 8,
#     "ytick.major.size" : 8,
#     "xtick.minor.size" : 4,
#     "ytick.minor.size" : 4,
#     "axes.spines.right": False,
#     "axes.spines.top" : False,
#     "legend.frameon":False
# }
plt.rcParams.update(util.my_rcParams)
def read_martini_jet_shape(fname):
    tmp = pd.read_csv(fname,comment='#')
    njet = 1
    with open(fname, 'r') as f:
        line = f.readline()
        line = line.split(' ')[-1]
        njet = float(line)
    dat     = tmp[tmp['rmax'] < 0.31]
    delta_r = tmp['rmax'] - tmp['rmin']
    r       = 0.5*(tmp['rmax'] + tmp['rmin'])
    rho  = tmp['N']  /(njet)
    drho = tmp['dN'] /(njet)
    norm = sum(rho.to_list())#
    rho_normed  = rho #/(delta_r * norm)
    drho_normed = drho#/(delta_r * norm)
    return rho_normed, drho_normed, r, delta_r
pT_lower_lim = 20
pT_upper_lim = 96
inel_Xsec= 62.03948 

## read in the experimental results
data_fname = "../../../exp_data/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_00-10_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
tmp = np.loadtxt(data_fname,comments='#',unpack=True, delimiter='\t')
cms_data = pd.DataFrame({'x':tmp[0],'y':tmp[1],'dx':tmp[2],'dy':tmp[3]})

## Get the Proton-Proton Baseline
path = '../martini_results/pp/jet_shape.csv'
calcs = pd.read_csv(path,comment='#')
njet = 1
with open(path, 'r') as f:
    line = f.readline()
    line = line.split(' ')[-1]
    njet = float(line)

dat = calcs[calcs['rmax'] < 0.31]
delta_r = dat['rmax'] - dat['rmin']
r = 0.5*(dat['rmax'] + dat['rmin'])
rho  = dat['N']  /(njet)
drho = dat['dN'] /(njet)
norm = sum(rho.to_list())#
rho_normed = rho  /(delta_r * norm)
drho_normed = drho/(delta_r * norm)

pp_data = pd.DataFrame({'r':r, 'dr':delta_r, 
                        'rho_normed':rho_normed, 
                        'drho_normed':drho_normed})
## Plot the data:
fig, ax = plt.subplots(1, 1, gridspec_kw={'top':0.985,
                                'bottom':0.1, 'left':0.08,
                                'right':0.95, 'hspace':0.02,
                                'wspace':0.18}, figsize=(16,9), sharex=True, sharey=True)

errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                for x, delx, y, yerrlow, yerrhigh in
                zip(cms_data["x"], cms_data['dx'], cms_data["y"], cms_data["dy"], cms_data["dy"])]
# Create patch collection with specified colour/alpha
pc = PatchCollection(errorboxes, facecolor='black', edgecolor="black", alpha=0.4)##eca400
ax.add_collection(pc)
scatter = ax.scatter(cms_data['x'],cms_data['y'],color='black',marker='s',s=30, label='CMS (2014)')
from COLORS import rate_set_colors as rate_colours
#rate_colours = #my_dicts.rate_colours
rate_names = {1:'LO', 2:'NLO',3:'NP'}

for rate_set in rate_names:
    color = rate_colours[rate_names[rate_set]]
    fname = f'../martini_results/final_PbPb_2p76/rset_{rate_set}/cent_0_5/jet_shape.csv'
    rho_normed1, drho_normed1, r, delta_r = read_martini_jet_shape(fname)
    fname = f'../martini_results/final_PbPb_2p76/rset_{rate_set}/cent_5_10/jet_shape.csv'
    rho_normed2, drho_normed2, r, delta_r = read_martini_jet_shape(fname)
    rho_normed = 0.5*(rho_normed1+rho_normed2)
    scale = sum(rho_normed)
    rho_normed = rho_normed/(delta_r*scale)
    drho_normed = 0.5*np.sqrt(drho_normed1**2 + drho_normed2**2)
    drho_normed = drho_normed/(delta_r*scale)
    ratio = rho_normed/pp_data['rho_normed']
    dr = ratio*np.sqrt((drho_normed/rho_normed)**2 + (pp_data['drho_normed']/pp_data['rho_normed'])**2)
    ax.plot(r, ratio, color=color)
    ax.fill_between(r, ratio-dr, ratio+dr, alpha=0.2, color=color)

ax.text(0.02,0.8,s=r'$0.3<|\eta|<2.0$'+'\n'+r'$100$ GeV $< p^{\mathrm{jet}}_T$', transform=ax.transAxes)
ax.text(0.02,0.7,s=r'$R=0.3$, Anti-$k_{T}$', transform=ax.transAxes)
ax.text(0.02,0.6, s=r'Centrality: $0$-$10\%$', transform=ax.transAxes)
# labels = [Line2D([],[],color=c,label=l) for (l,c) in momentum_cut_and_colors.items()]
# ax.legend(handles=labels, loc='upper center')
labels = [Line2D([],[],c=c,label=l) for l, c in rate_colours.items()]
labels.append(Line2D([],[],label='CMS',marker='s',color='black', markersize=10, linewidth=0))
ax.legend(loc='lower right', handles=labels)
ax.set_ylabel(r'$\rho(r)_{PbPb}/\rho(r)_{pp}$')
ax.set_xlabel(r'$r$')
#ax.axhline(1,linestyle='dotted',color='black')
plt.show()
