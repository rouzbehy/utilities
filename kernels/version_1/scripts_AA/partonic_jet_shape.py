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

def get_njet(fname):
    with open(fname, 'r') as f:
        line = f.readline()
        line = line.split(' ')[-1]
        njet = float(line)
        return njet

pT_lower_lim = 20
pT_upper_lim = 96
inel_Xsec= 62.03948 
rmax = 0.8
successful_shapes = {}
failed_shapes = {}
momentum_cut_and_colors = {#'0p2':'orange',#'0p5':'green',
                           #'1p0':'red','2p0':'blue',
                           '4T':'cyan','8T':'purple','10T':'black',
                           '4T_old':'green', '8T_old':'red', '10T_old':'blue'}
for pcut in momentum_cut_and_colors:
    iset = 9 if 'old' not in pcut else 8
    pcut_tmp = pcut.replace('_old','')
    #iset = 4 if 'T' not in pcut else 6
    loc = "./processed/data_set{iset}/pcut_{pcut}/"
    loc = loc.format(iset=iset, pcut=pcut_tmp)
    fname = loc + "parton_jet_shape_failed.csv"
    tmp = pd.read_csv(fname, comment='#')
    njets = get_njet(fname)
    dat = tmp[tmp['rmax'] < rmax]
    delta_r = dat['rmax'] - dat['rmin']
    r = 0.5*(dat['rmax'] + dat['rmin'])
    failed_dict = {}
    failed_dict['r'] = r
    failed_dict['dr'] = delta_r
    norm_factor_1 = 0.
    for item in ['fermion','gluon']:
        rho  = dat[item]  /(njets)
        drho = dat[f'd{item}'] /(njets)
        norm = sum(rho.to_list())
        rho_normed  = rho #/(delta_r)# * norm)
        drho_normed = drho#/(delta_r)# * norm)
        failed_dict[item] = rho_normed
        failed_dict[f'd{item}'] = drho_normed
        norm_factor_1 += norm 
    failed_dict['norm'] = norm_factor_1
    failed_shapes[pcut] = pd.DataFrame(failed_dict)

    fname = loc + "parton_jet_shape_success.csv"
    tmp = pd.read_csv(fname, comment='#')
    njets = get_njet(fname)
    dat = tmp[tmp['rmax'] < rmax]
    delta_r = dat['rmax'] - dat['rmin']
    r = 0.5*(dat['rmax'] + dat['rmin'])
    success_dict = {}
    success_dict['r'] = r
    success_dict['dr'] = delta_r
    norm_factor_2 = 0.
    for item in ['fermion','gluon']:
        rho  = dat[item]  /(njets)
        drho = dat[f'd{item}'] /(njets)
        norm = sum(rho.to_list())
        rho_normed  = rho  #/(delta_r)# * norm)
        drho_normed = drho #/(delta_r)# * norm)
        norm_factor_2 += norm
        success_dict[item] = rho_normed
        success_dict[f'd{item}'] = drho_normed
    success_dict['norm'] = norm_factor_2
    successful_shapes[pcut] = pd.DataFrame(success_dict)

fig, axes = plt.subplots(1, 2, gridspec_kw={'top':0.985,
                      'bottom':0.1, 'left':0.1,
                      'right':0.95, 'hspace':0.02,
                      'wspace':0.18}, figsize=(16,9), sharex=True, sharey=True)
axes = axes.flatten()

#colors = {'0p2':'orange','0p5':'green','1p0':'red','2p0':'blue'}
parton_lines = {'fermion':'solid', 'gluon':'dashed'}

for pcut in momentum_cut_and_colors:
    color = momentum_cut_and_colors[pcut]
    tmp = successful_shapes[pcut]
    ax = axes[0]
    #ax.text(0.02,0.1,s=r'$0.3<|\eta|<2.0$'+'\n'+r'$100$ GeV $< p^{\mathrm{jet}}_T$', transform=ax.transAxes)
    ax.text(0.6,0.8,'Successful Events', transform=ax.transAxes)
    for parton in ['fermion','gluon']:
        lstyle = parton_lines[parton]
        r, y, dy = tmp['r'], tmp[parton], tmp[f'd{parton}']
        normalization = tmp['norm']
        y = y/normalization
        dy = dy/normalization
        ax.plot(r, y, color=color, linestyle=lstyle)
        #ax.fill_between(r, y-dy, y+dy, color=color, alpha=0.2)

    ax = axes[1]
    tmp = failed_shapes[pcut]
    ax.text(0.6,0.8,'Failed Events', transform=ax.transAxes)
    for parton in ['fermion','gluon']:
        lstyle = parton_lines[parton]
        r, y, dy = tmp['r'], tmp[parton], tmp[f'd{parton}']
        normalization = tmp['norm']
        y = y/normalization
        dy = dy/normalization
        ax.plot(r, y, color=color, linestyle=lstyle)
        #ax.plot(r, y, color='red', linestyle=lstyle)
        #ax.fill_between(r, y-dy, y+dy, color=color, alpha=0.2)

for ax in axes:
    ax.set_xlabel(r'$r$')
    ax.set_xlim(left=-0.01,right=0.41)

line_labels = [Line2D([],[],label='pCut='+l.replace('p','.')+" GeV",color=c) for l,c in momentum_cut_and_colors.items()]
axes[1].legend(loc='center right', handles=line_labels)
axes[0].set_ylabel(r'$\rho(r)$')
axes[0].text(0.5,0.7,r'$p_{T}>20$ GeV, $|\eta|<1.0$'+"\n"+"R=0.4", transform=axes[0].transAxes)

fig1, ax1 = plt.subplots(1, 1, gridspec_kw={'top':0.985,
                      'bottom':0.1, 'left':0.1,
                      'right':0.95, 'hspace':0.02,
                      'wspace':0.18}, figsize=(16,9), sharex=True, sharey=True)
for pcut in momentum_cut_and_colors:
    #color = momentum_cut_and_colors[pcut]
    color = 'blue'
    tmp = successful_shapes[pcut]
    #ax.text(0.02,0.1,s=r'$0.3<|\eta|<2.0$'+'\n'+r'$100$ GeV $< p^{\mathrm{jet}}_T$', transform=ax.transAxes)
    #ax.text(0.6,0.8,'Successful Events', transform=ax.transAxes)
    for parton in ['fermion','gluon']:
        lstyle = parton_lines[parton]
        r, y, dy = tmp['r'], tmp[parton], tmp[f'd{parton}']
        normalization = tmp['norm']
        y = y/normalization
        dy = dy/normalization
        ax1.plot(r, y, color=color, linestyle=lstyle)
        #ax.fill_between(r, y-dy, y+dy, color=color, alpha=0.2)

    tmp = failed_shapes[pcut]
    #ax.text(0.6,0.8,'Failed Events', transform=ax.transAxes)
    color = 'red'
    for parton in ['fermion','gluon']:
        lstyle = parton_lines[parton]
        r, y, dy = tmp['r'], tmp[parton], tmp[f'd{parton}']
        normalization = tmp['norm']
        y = y/normalization
        dy = dy/normalization
        ax1.plot(r, y, color=color, linestyle=lstyle)
        #ax.plot(r, y, color='red', linestyle=lstyle)
        #ax.fill_between(r, y-dy, y+dy, color=color, alpha=0.2)

ax1.set_xlabel(r'$r$')
ax1.set_xlim(left=-0.01,right=0.41)
label_dict = {'successful' : 'blue', 'failed':'red'}
parton_dict = {'gluon':'dashed','fermion':'solid'}
#line_labels = [Line2D([],[],label='pCut='+l.replace('p','.')+" GeV",color=c) for l,c in momentum_cut_and_colors.items()]
line_labels = [Line2D([],[],label=l, color=c) for l, c in label_dict.items()]
for p, lstyle in parton_dict.items():
    line_labels.append(Line2D([],[],label=p,color='black',linestyle=lstyle))

ax1.legend(loc='center right', handles=line_labels)
ax1.set_ylabel(r'$\rho(r)$')
ax1.text(0.5,0.7,r'$p_{T}>20$ GeV, $|\eta|<1.0$'+"\n"+"R=0.4", transform=ax1.transAxes)
plt.show()
