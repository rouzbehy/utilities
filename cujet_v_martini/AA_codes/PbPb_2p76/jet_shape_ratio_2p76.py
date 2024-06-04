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
from COLORS import module_colors

plt.rcParams.update(util.my_rcParams)

def construct_new_centrality(cent1, cent2):
   tmp = {}
   tmp['r'] = cent1['r']
   tmp['dr'] = cent1['dr']
   tmp['rho'] = 0.5*(cent1['rho'] + cent2['rho'])
   tmp['drho'] = np.sqrt(cent1['drho']*cent1['drho'] + cent2['drho']*cent2['drho'])
   return pd.DataFrame(tmp) 

jet_shape_cents = {0:'00-10',1:'10-30',2:'30-50'}

## read in the experimental results
data_fname = "../../../expt/PbPb_2p76/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_{cent}_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
CMS_DATA = {}
for _,cent in jet_shape_cents.items():
    tmp = np.loadtxt(data_fname.format(cent=cent),comments='#',unpack=True, delimiter='\t')
    CMS_DATA[cent] = pd.DataFrame({'x':tmp[0],'y':tmp[1],'dx':tmp[2],'dy':tmp[3]})

## Plot the data:
fig, axes = plt.subplots(1, 3, figsize=(18,6), sharex=True, sharey=True)

## Read in the PP jet shape and massage it:
path_pp = '../../jetscape_data/max_time/maxT_200_highstat/pp_2760_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv'
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
path_AA = '../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv'
jetscape = {}
elosses = ['cujet', 'martini']
cents_jetscape = {0:'00-05',1:'05-10',2:'10-20',3:'20-30',4:'30-40',5:'40-50'}
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

    ## compute 00-10, 10-30 and 30-50 centralities here:
    tmp1, tmp2 = jetscape[eloss]['00-05'], jetscape[eloss]['05-10']
    jetscape[eloss]['00-10'] = construct_new_centrality(tmp1, tmp2)
    tmp3, tmp4 = jetscape[eloss]['10-20'], jetscape[eloss]['20-30']
    jetscape[eloss]['10-30'] = construct_new_centrality(tmp3, tmp4)
    tmp5, tmp6 = jetscape[eloss]['30-40'], jetscape[eloss]['40-50'] 
    jetscape[eloss]['30-50'] = construct_new_centrality(tmp5, tmp6)

for icent, cent in jet_shape_cents.items():
    data = CMS_DATA[cent]
    errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                    for x, delx, y, yerrlow, yerrhigh in
                    zip(data["x"], data['dx'], data["y"], data["dy"], data["dy"])]
    
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor='black', edgecolor="black", alpha=0.4)##eca400
    axes[icent].add_collection(pc)
    scatter = axes[icent].scatter(data['x'],data['y'],color='black',marker='s',s=30, label='CMS (2014)')
    centrality_label = "$"+f"{cent}" + "\%$"
    axes[icent].text(0.2,0.1,centrality_label, transform=axes[icent].transAxes)
    for eloss in elosses:
        color = module_colors['MATTER+'+eloss.upper()]
        AA = jetscape[eloss][cent]
        ratio = AA['rho']/pp['rho']
        dratio = ratio*np.sqrt(AA['drho']*AA['drho']/(AA['rho']*AA['rho']) + pp['drho']*pp['drho']/(pp['rho']*pp['rho']))
        errorboxes = [Rectangle((x-delx, y - dy), width=2*delx, height=2*dy)
                    for x, delx, y, dy in zip(AA['r'], 0.5*AA['dr'], ratio, dratio)]
        pc = PatchCollection(errorboxes, facecolor=color, edgecolor=color, alpha=0.5)
        axes[icent].add_collection(pc)
        axes[icent].scatter(AA[['r']], ratio, color=color,marker=ddicts.eloss_marker[eloss])

labels = [Line2D([],[],label='MATTER+'+eloss.upper(),color=module_colors['MATTER+'+eloss.upper()], marker=ddicts.eloss_marker[eloss],markersize=8) for eloss in elosses]
labels.append(Line2D([],[], marker='s', color='black', label='CMS (2014)',markerfacecolor='black', markersize=8))
axes[2].legend(loc='upper right', handles=labels)

for ax in axes:
    ax.set_xlabel(r'$r$')

#axes[0].set_ylabel(r'$\rho_{\mathrm{PbPb}}\left(r\right)/\rho_{\mathrm{pp}}\left(r\right)$')
axes[0].set_ylabel(r'$R_{\rho}\left(r\right)$')
axes[0].set_ylim(top=1.55, bottom=0.5)
axes[1].set_xticks([0,0.1,0.2,0.3])
axes[1].set_xticklabels(['0.0','0.1','0.2','0.3'])

axes[1].text(0.05, 0.8, 'Pb-Pb, $\sqrt{s}=2.76$ ATeV' + "\n" + "$0.3<|\eta|<2.0$", fontsize=30, transform=axes[1].transAxes)
plt.show()
