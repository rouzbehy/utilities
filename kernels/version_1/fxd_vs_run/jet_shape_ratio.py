#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
import getSpecs as spc
import util 
import dictionaries as my_dicts
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

#pT_lower_lim = 20
#pT_upper_lim = 96
rmax = 0.31
inel_Xsec= 62.03948 

## read in the experimental results
data_fname = "../../../exp_data/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_00-10_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
tmp = np.loadtxt(data_fname,comments='#',unpack=True, delimiter='\t')
cms_data = pd.DataFrame({'x':tmp[0],'y':tmp[1],'dx':tmp[2],'dy':tmp[3]})

rate_colours = my_dicts.rate_colours
rate_names = {1:'LO', 2:'NLO',3:'NP'}

fname_tmp = '../martini_results/pp/jet_shape.csv'
pp_data = spc.get_jet_shape(fname_tmp, rmax)
fname_tmp = '../martini_results/fixed_alphas_const/rset_{r}/jet_shape.csv'
#fname_tmp = '../martini_results/final_PbPb_2p76/rset_{r}/cent_0_5/jet_shape.csv'
AA_fxd_c = {r: spc.get_jet_shape(fname_tmp.format(r=r), rmax) for r in rate_names}
fname_tmp = '../martini_results/fixed_alphas_fitted/rset_{r}/jet_shape.csv'
AA_fxd_f = {r: spc.get_jet_shape(fname_tmp.format(r=r), rmax) for r in rate_names}
fname_tmp = '../martini_results/run_alphas_fit/rset_{r}/jet_shape.csv'
AA_run_f_v2 = {r: spc.get_jet_shape(fname_tmp.format(r=r), rmax) for r in rate_names}
fname_tmp = '../martini_results/rset_{r}/jet_shape.csv'
AA_run_f_v1 = {r: spc.get_jet_shape(fname_tmp.format(r=r), rmax) for r in rate_names}


fig, ax = plt.subplots(1, 1, gridspec_kw={'top':0.985,
                                'bottom':0.1, 'left':0.1,
                                'right':0.95, 'hspace':0.02,
                                'wspace':0.18}, figsize=(16,9), sharex=True, sharey=True)

errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                for x, delx, y, yerrlow, yerrhigh in
                zip(cms_data["x"], cms_data['dx'], cms_data["y"], cms_data["dy"], cms_data["dy"])]
# Create patch collection with specified colour/alpha
pc = PatchCollection(errorboxes, facecolor='black', edgecolor="black", alpha=0.4)##eca400
ax.add_collection(pc)
scatter = ax.scatter(cms_data['x'],cms_data['y'],color='black',marker='s',s=30, label='CMS (2014)')

rate_colours = my_dicts.rate_colours
rate_names = {1:'LO', 2:'NLO',3:'NP'}
yp,  dyp = pp_data['rho'], pp_data['drho']
err_fac = dyp*dyp/(yp*yp)
r = pp_data['r']
for rate_set in rate_names:
    color = rate_colours[rate_names[rate_set]]
    y1, dy1 = AA_run_f_v1[rate_set]['rho'], AA_run_f_v1[rate_set]['drho']
    y2, dy2 = AA_run_f_v2[rate_set]['rho'], AA_run_f_v2[rate_set]['drho']
    y3, dy3 = AA_fxd_f[rate_set]['rho'], AA_fxd_f[rate_set]['drho']
    y4, dy4 = AA_fxd_c[rate_set]['rho'], AA_fxd_c[rate_set]['drho']
    ratio_1 = y1/yp
    dration_1 = ratio_1*np.sqrt(dy1*dy1/(y1*y1) + err_fac)
    ratio_2   = y2/yp
    dration_2 = ratio_2*np.sqrt(dy2*dy2/(y2*y2) + err_fac)
    ratio_3   = y3/yp
    dration_3 = ratio_3*np.sqrt(dy3*dy3/(y3*y3) + err_fac)
    ratio_4   = y4/yp
    dration_4 = ratio_4*np.sqrt(dy4*dy4/(y4*y4) + err_fac)

    #ax.plot(r, ratio_1, color=color, linestyle='solid')
    ax.plot(r, ratio_2, color=color, linestyle='solid')
    ax.plot(r, ratio_3, color=color, linestyle='dashed')
    ax.plot(r, ratio_4, color=color, linestyle='dotted')

ax.text(0.02,0.60,s=r'$0.3<|\eta|<2.0$'+'\n'+r'$100$ GeV $< p^{\mathrm{jet}}_T$', transform=ax.transAxes)
ax.text(0.02,0.55,s=r'$R=0.3$, Anti-$k_{T}$', transform=ax.transAxes)
ax.text(0.02,0.50,s=r'$0$-$5$\%', transform=ax.transAxes)

extra_labels = {r'$\alpha_s=\alpha_s(\mu)$':'solid',
                r'$\alpha_s=\alpha_{s,i}$':'dashed',
                r'$\alpha_s=0.3$':'dotted'}
theor_hands = [Line2D([],[],label=l,color=c) for (l,c) in rate_colours.items()]
for l, st in extra_labels.items():
    theor_hands.append(Line2D([],[],label=l,linestyle=st,color='black'))
theor_hands.append(Line2D([],[], marker='s', color='black', label='CMS $0$-$10\%$ (2014)',markerfacecolor='black', markersize=8))
ax.legend(handles=theor_hands, loc='upper center')
ax.set_ylabel(r'$\rho(r)_{PbPb}/\rho(r)_{pp}$')
ax.set_xlabel(r'$r$')
#ax.axhline(1,linestyle='dotted',color='black')
plt.show()
