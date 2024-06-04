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
import jetDicts as ddicts
from COLORS import module_colors as colors
mpl.use('MacOSX')
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

jet_shape_cents = {0:'00-05'}


## Read in the PP jet shape and massage it:
path_pp = '../../../jetscape_data/max_time/maxT_200_highstat/pp_2760_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv'
tmp = pd.read_csv(path_pp,comment='#')
nevt_wcut, njet_wcut = 1, 1
with open(path_pp, 'r') as f:
    line = f.readline()
    line = line.split()
    (nevt_wcut, njet_wcut) = float(line[-3]), float(line[-1])
tmp = tmp[tmp['rmax'] < 0.31]
delta_r = tmp['rmax'] - tmp['rmin']
r = 0.5*(tmp['rmax'] + tmp['rmin'])
rho  = tmp['wcut']  /(nevt_wcut * njet_wcut )
drho = tmp['dwcut'] /(nevt_wcut * njet_wcut )
norm = sum(rho.to_list())#
rho_normed = rho  /(delta_r * norm)
drho_normed = drho/(delta_r * norm)
pp = pd.DataFrame({'r':r, 'dr':delta_r, 'rho':rho_normed, 'drho':drho_normed})

## Read in the AA jet shape results:
path_AA = '../../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv'
path_AA_matter = '../../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760_maxT_{maxT}/PbPb2760_{cent}_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv'
jetscape = {}
cents_jetscape = {0:'00-05'}
eloss = 'matter'
options = [200]
jetscape[eloss] = {}
for opt in options:
    jetscape[eloss][opt] = {}
    for icent, cent in cents_jetscape.items():
        path_jetscape = path_AA_matter.format(cent=cent, eloss=eloss,maxT=opt)
        tmp = pd.read_csv(path_jetscape, comment='#')
        nevt_wcut, njet_wcut = 1, 1
        with open(path_jetscape, 'r') as f:
            line = f.readline()
            line = line.split()
            (nevt_wcut, njet_wcut) = float(line[-3]), float(line[-1])
        tmp = tmp[tmp['rmax'] < 0.31]
        delta_r = tmp['rmax'] - tmp['rmin']
        r = 0.5*(tmp['rmax'] + tmp['rmin'])
        rho  = tmp['wcut']  /(nevt_wcut * njet_wcut )
        drho = tmp['dwcut'] /(nevt_wcut * njet_wcut )
        norm = sum(rho.to_list())
        #norm = simpson(rho.to_list(), r.to_list())
        rho_normed = rho  /(delta_r * norm)
        drho_normed = drho/(delta_r * norm)
        jetscape[eloss][opt][cent] = pd.DataFrame({'r':r, 'dr':delta_r, 'rho':rho_normed, 'drho':drho_normed})

elosses = ['cujet', 'martini']
for eloss in elosses:
    jetscape[eloss] = {}
    for icent, cent in cents_jetscape.items():
        path_jetscape = path_AA.format(cent=cent, eloss=eloss)
        tmp = pd.read_csv(path_jetscape, comment='#')
        nevt_wcut, njet_wcut = 1, 1
        with open(path_jetscape, 'r') as f:
            line = f.readline()
            line = line.split()
            (nevt_wcut, njet_wcut) = float(line[-3]), float(line[-1])
        tmp = tmp[tmp['rmax'] < 0.31]
        delta_r = tmp['rmax'] - tmp['rmin']
        r = 0.5*(tmp['rmax'] + tmp['rmin'])
        rho  = tmp['wcut']  /(nevt_wcut * njet_wcut )
        drho = tmp['dwcut'] /(nevt_wcut * njet_wcut )
        norm = sum(rho.to_list())
        #norm = simpson(rho.to_list(), r.to_list())
        rho_normed = rho  /(delta_r * norm)
        drho_normed = drho/(delta_r * norm)
        jetscape[eloss][cent] = pd.DataFrame({'r':r, 'dr':delta_r, 'rho':rho_normed, 'drho':drho_normed})

fig, ax= plt.subplots(1, 1,gridspec_kw={'top':0.98,'bottom':0.12, 'left':0.08,
                                'right':0.95, 'hspace':0.02,'wspace':0.18}, figsize=(14,10),sharex=True, sharey=True)
cent = '00-05'
for eloss in elosses:
    AA = jetscape[eloss][cent]
    ratio = AA['rho']/pp['rho']
    dratio = ratio*np.sqrt(AA['drho']*AA['drho']/(AA['rho']*AA['rho']) + pp['drho']*pp['drho']/(pp['rho']*pp['rho']))
    errorboxes = [Rectangle((x-delx, y - dy), width=2*delx, height=2*dy)
                for x, delx, y, dy in zip(AA['r'], 0.5*AA['dr'], ratio, dratio)]
    #pc = PatchCollection(errorboxes, facecolor=ddicts.eloss_colours[eloss], edgecolor=ddicts.eloss_colours[eloss], alpha=0.5)
    pc = PatchCollection(errorboxes, facecolor=colors[f'MATTER+{eloss}'.upper()], edgecolor=colors[f'MATTER+{eloss}'.upper()], alpha=0.5)
    ax.add_collection(pc)
    #ax.scatter(AA[['r']], ratio, color=ddicts.eloss_colours[eloss],marker=ddicts.eloss_marker[eloss])
    ax.scatter(AA[['r']], ratio, color=colors[f'MATTER+{eloss}'.upper()],marker=ddicts.eloss_marker[eloss])

eloss = 'matter'
for opt in options:
    AA = jetscape[eloss][opt][cent]
    ratio = AA['rho']/pp['rho']
    dratio = ratio*np.sqrt(AA['drho']*AA['drho']/(AA['rho']*AA['rho']) + pp['drho']*pp['drho']/(pp['rho']*pp['rho']))
    errorboxes = [Rectangle((x-delx, y - dy), width=2*delx, height=2*dy)
                for x, delx, y, dy in zip(AA['r'], 0.5*AA['dr'], ratio, dratio)]
    #color = ddicts.eloss_colours[f'{eloss}_{opt}']
    color = colors[f'{eloss}'.upper()]
    pc = PatchCollection(errorboxes, facecolor=color, edgecolor=color, alpha=0.5)
    ax.add_collection(pc)
    ax.scatter(AA[['r']], ratio, color=color,marker=ddicts.eloss_marker[eloss])

# handles_theory = [Line2D([],[],color=ddicts.eloss_colours[eloss],label=label,marker=ddicts.eloss_marker['matter'])\
#                      for eloss,label in zip(['matter_20','matter_200'],
#                                             ['MATTER ($t_{\mathrm{max}}$ = 20 $\mathrm{fm}/c$)',\
#                                              'MATTER ($t_{\mathrm{max}}$ = 200 $\mathrm{fm}/c$)'])]
# handles_theory.append(Line2D([],[],color=ddicts.eloss_colours['martini'], label='MATTER+MARTINI',marker=ddicts.eloss_marker['martini']))
# handles_theory.append(Line2D([],[],color=ddicts.eloss_colours['cujet'], label='MATTER+CUJET',marker=ddicts.eloss_marker['cujet']))
handles_theory = [Line2D([],[],color=c, label=l) for l,c in colors.items()]
ax.legend(loc='upper right', handles=handles_theory)

ax.set_ylabel(r'$\rho_{\mathrm{PbPb}}\left(r\right)/\rho_{\mathrm{pp}}\left(r\right)$')
ax.set_ylim(top=1.55, bottom=0.5)
ax.set_xticks([0,0.1,0.2,0.3])
ax.set_xticklabels(['0.0','0.1','0.2','0.3'])

ax.text(0.05, 0.8, 'Pb-Pb, $\sqrt{s}=2.76$ ATeV' + "\n" + "$0.3<|\eta|<2.0$" + "\n" + r"$R=0.3$, $0$-$5\%$", fontsize=30, transform=ax.transAxes)
ax.text(0.05, 0.65, r'Anti-$k_{T}$', transform=ax.transAxes)
ax.set_xlabel(r'$r$')

ax.axhline(1,color='black',linestyle='dotted')

plt.show()
