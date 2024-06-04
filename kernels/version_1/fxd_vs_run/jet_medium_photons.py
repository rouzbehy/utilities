## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
import util
import dictionaries as my_hdicts
import jetDicts as my_jdicts
import getSpecs as spc
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
pTlow, pThigh = 1, 200
from COLORS import rate_set_colors as rate_colours
## read in the calculations
rate_names = {1:'LO',2:'NLO',3:'NP'}

## JET MEDIUM:
## now read in fixed alphas where all rate sets have alphas=0.3
fname_tmp = '../martini_results/fixed_alphas_const/rset_{r}/gamma_spectra.csv'
AA_fxd_c = {r: spc.get_charged_spec(fname_tmp.format(r=r), pTlow, pThigh ) for r in rate_names}
fname_tmp = '../martini_results/fixed_alphas_fitted/rset_{r}/gamma_spectra.csv'
AA_fxd_f = {r: spc.get_charged_spec(fname_tmp.format(r=r), pTlow, pThigh ) for r in rate_names}
fname_tmp = '../martini_results/run_alphas_fit/rset_{r}/gamma_spectra.csv'
AA_run_f = {r: spc.get_charged_spec(fname_tmp.format(r=r), pTlow, pThigh ) for r in rate_names}

## One figure only: since fixed alpha_s runs are only 0-5% 
## Compare the jet-medium channels only
fig, axes = plt.subplots(3, 3, sharex=True, sharey='row', height_ratios=[2,1,1])
column = 0
for run_collection, label in zip([AA_fxd_c, AA_fxd_f, AA_run_f],
                                 [r'$\alpha_s=0.3$',
                                  r'$\alpha_s=\alpha_{s,i}$',
                                  r'$\alpha_s=\alpha_{s}(p,T)$']):
    axes[0][column].text(0.05, 0.05, label, fontsize=25, transform=axes[0][column].transAxes)
    for rateset in rate_names:
        color = rate_colours[rate_names[rateset]]
        run = run_collection[rateset]
        run['pT']  = 0.5*(run['pTmin']+run['pTmax'])
        run['dpT'] = run['pTmax'] - run['pTmin']
        #tmp = run[run['pT'].between(pTlow,pThigh)]
        x, dx = run['pT'], run['dpT']
        yconv, dyconv = run['conv'], run['dconv']
        ybrem, dybrem = run['brem'], run['dbrem']
        yconv  = yconv/(2*np.pi*x*dx*2*0.8*1000)
        dyconv = dyconv/(2*np.pi*x*dx*2*0.8*1000)
        ybrem  = ybrem/(2*np.pi*x*dx*2*0.8*1000)
        dybrem = dybrem/(2*np.pi*x*dx*2*0.8*1000)
        total = yconv + ybrem
        dtotal = np.sqrt(dybrem*dybrem+dyconv*dyconv)
        run['total'] = total.to_list()
        run['dtotal'] = dtotal.to_list()
        
        row = 0
        axes[row][column].plot(x, total, color=color, linestyle='solid')
        axes[row][column].fill_between(x, total-dtotal, total+dtotal, color=color, alpha=0.2)
        row = 1
        ## ratios of channel to total
        total_err = dtotal*dtotal/(total*total)
        rconv = yconv/total
        drconv = rconv*np.sqrt(dyconv*dyconv/(yconv*yconv)+total_err)
        rbrem = ybrem/total
        drbrem = rbrem*np.sqrt(dybrem*dybrem/(ybrem*ybrem)+total_err)
        axes[row][column].plot(x, rbrem, color=color, linestyle='solid')
        axes[row][column].fill_between(x, rbrem-drbrem, rbrem+drbrem, color=color, alpha=0.2)
        axes[row][column].plot(x, rconv, color=color, linestyle='dashed')
        axes[row][column].fill_between(x, rconv-drconv, rconv+drconv, color=color, alpha=0.2)
    column+=1
## ratios of rate sets to LO
column = 0
for run_collection in [AA_fxd_c, AA_fxd_f, AA_run_f]:
    y1, dy1 = run_collection[1]['total'], run_collection[1]['dtotal']
    errfac = dy1*dy1/(y1*y1)
    for rateset in [2,3]:
        color = rate_colours[rate_names[rateset]]
        run = run_collection[rateset]
        y2, dy2 = run['total'], run['dtotal']
        ratio = y2/y1
        dratio = ratio*np.sqrt(dy2*dy2/(y2*y2) + errfac)
        axes[2][column].plot(run['pT'], ratio, color=color)
        axes[2][column].fill_between(run['pT'], ratio+dratio, ratio-dratio, alpha=0.2, color=color)
    column += 1

for ax in axes[2]:
    ax.set_xlabel(r'$p^{\gamma}_T$ (GeV)')

axes[0][0].set_yscale('log')
axes[0][0].set_ylabel(r'$\frac{1}{N_{\mathrm{bin}}2\pi p_T}\frac{d\sigma^{\gamma}}{d\eta dp_T}$ (mb.GeV${}^{-2}$)')
axes[1][0].set_ylabel('Channel Ratios', fontsize=18)
axes[2][0].set_ylabel('Ratio to LO'   , fontsize=18)
labels_1 = [Line2D([],[],color=c,label=l) for l, c in rate_colours.items()]
labels_2 = []
labels_1.append(Line2D([],[],color='black',linestyle='solid',label='Brem/Total'))
labels_1.append(Line2D([],[],color='black',linestyle='dashed',label='Conv/Total'))
#labels_3 = []
#labels_3.append(Line2D([],[],color=my_hdicts.rate_colours['NLO'],linestyle='solid',label='NLO/LO'))
#labels_3.append(Line2D([],[],color=my_hdicts.rate_colours['NP'],linestyle='dashed',label='NP/LO'))
axes[0][2].legend(handles=labels_1, bbox_to_anchor=(0.6,1.0), loc='upper left')
#axes[1][2].legend(handles=labels_2, bbox_to_anchor=(0.9,0.5))
#axes[2][2].legend(handles=labels_3, bbox_to_anchor=(1.01,0.8))

## Simple figure to take ratio of fixed and fitted alphas run
## to the fitted and running alphas run

## ratios of rate sets to LO
fig, ax = plt.subplots(1, 1)
for rateset in [1,2,3]:
    color = rate_colours[rate_names[rateset]]
    x, dx = AA_fxd_f[rateset]['pT'], AA_fxd_f[rateset]['dpT'] 
    y1, dy1 = AA_fxd_f[rateset]['total'], AA_fxd_f[rateset]['dtotal']
    y2, dy2 = AA_run_f[rateset]['total'], AA_run_f[rateset]['dtotal']
    ratio  = y2/y1
    dratio = ratio*np.sqrt(dy2*dy2/(y2*y2) + dy1*dy1/(y1*y1))
    ax.plot(x, ratio, color=color)
    ax.fill_between(x, ratio-dratio, ratio+dratio, color=color,alpha=0.2)

labels = [Line2D([],[],color=c,label=l) for l, c in rate_colours.items()]
ax.legend(loc='best', handles=labels)
ax.set_xlabel(r'$p^{\gamma}_T$ (GeV)')
ax.set_ylabel(r'Ratio $\alpha_{s}(p,T)/\alpha_{s,i}$')
plt.show()