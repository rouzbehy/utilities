#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import util
from COLORS import rate_set_colors as rate_colours

plt.rcParams.update(util.my_rcParams)
rate_names = {1:'LO', 2:'NLO',3:'NP'}

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
    rho_normed  = rho /(delta_r * norm)
    drho_normed = drho/(delta_r * norm)
    return rho_normed, drho_normed, r, delta_r

pT_lower_lim, pT_upper_lim, inel_Xsec = 20, 96, 62.03948

## read in the experimental results
data_fname = "../../../exp_data/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_{c}_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
cms_data = {}
for c in ['00-10', '10-30']:
    tmp = np.loadtxt(data_fname.format(c=c),comments='#',unpack=True, delimiter='\t')
    cms_data[c] = pd.DataFrame({'x':tmp[0],'y':tmp[1],'dx':tmp[2],'dy':tmp[3]})

## Get the Proton-Proton Baseline
path = "../../calcs/pp_2p76/jet_shape.csv"

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

## Compute the jet shape ratio for each rate set and the two
## centralities of interest:
jet_shape_ratios = {}
for rate_set in rate_names:
    fname = "../../calcs/final_runs/rset_{r}/2p76/{c}/jet_shape.csv"
    rho_normed1, drho_normed1, r, delta_r = read_martini_jet_shape(fname.format(r=rate_set, c='0_5'))
    rho_normed2, drho_normed2, r, delta_r = read_martini_jet_shape(fname.format(r=rate_set, c='5_10'))
    rho_normed = 0.5*(rho_normed1+rho_normed2)
    scale = sum(rho_normed)
    rho_normed = rho_normed/(delta_r*scale)
    drho_normed = 0.5*np.sqrt(drho_normed1**2 + drho_normed2**2)
    drho_normed = drho_normed/(delta_r*scale)
    ratio = rho_normed/pp_data['rho_normed']
    dr = ratio*np.sqrt((drho_normed/rho_normed)**2 + (pp_data['drho_normed']/pp_data['rho_normed'])**2)
    jet_shape_ratios[rate_set] = {'0_10':(r, delta_r, ratio, dr)}


    rho_normed, drho_normed, r, delta_r = read_martini_jet_shape(fname.format(r=rate_set, c='10_20'))
    scale = sum(rho_normed)
    rho_normed = rho_normed/(delta_r*scale)
    drho_normed = 0.5*np.sqrt(drho_normed1**2 + drho_normed2**2)
    drho_normed = drho_normed/(delta_r*scale)
    ratio = rho_normed/pp_data['rho_normed']
    dr = ratio*np.sqrt((drho_normed/rho_normed)**2 + (pp_data['drho_normed']/pp_data['rho_normed'])**2)
    jet_shape_ratios[rate_set]['10_20'] = (r, delta_r, ratio, dr)

## Plot the data:
fig, axes = plt.subplots(1, 2, gridspec_kw={'top':0.985,
                                'bottom':0.1, 'left':0.08,
                                'right':0.95, 'hspace':0.02,
                                'wspace':0.1}, figsize=(16,9), sharex=True, sharey=True)

for c, ax in zip(cms_data, axes):
    data = cms_data[c]
    errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                    for x, delx, y, yerrlow, yerrhigh in
                    zip(data["x"], data['dx'], data["y"], data["dy"], data["dy"])]
    pc = PatchCollection(errorboxes, facecolor='black', edgecolor="black", alpha=0.4)
    ax.add_collection(pc)
    scatter = ax.scatter(data['x'], data['y'],color='black',marker='s',s=30, label='CMS (2014)')


for rate_set in rate_names:
    for iax, c in enumerate(['0_10','10_20']):
        ax = axes[iax]
        clo, chi = c.split('_')
        txt: str = f"${clo}$-${chi}$\%"
        ax.text(0.4,0.4, s=txt, transform=ax.transAxes)
        color = rate_colours[rate_names[rate_set]]
        r, dr, ratio, dratio = jet_shape_ratios[rate_set][c]

        ax.plot(r, ratio, color=color)
        ax.fill_between(r, ratio-dratio, ratio+dratio, alpha=0.2, color=color)

ax = axes[0]
ax.text(0.02,0.9,s=r'Pb-Pb, $\sqrt{s}=2.76$ ATeV', transform=ax.transAxes)
ax.text(0.02,0.75,s=r'$0.3<|\eta|<2.0$'+'\n'+r'$100$ GeV $< p^{\mathrm{jet}}_T$', transform=ax.transAxes)
ax.text(0.02,0.65,s=r'$R=0.3$, Anti-$k_{T}$', transform=ax.transAxes)
for ax in axes:
    ax.set_xlabel(r'$r$')

labels = [Line2D([],[],c=c,label=l) for l, c in rate_colours.items()]
labels.append(Line2D([],[],label='CMS',marker='s',color='black', markersize=10, linewidth=0))
axes[1].legend(loc='upper left', handles=labels, handletextpad=0.05)
axes[1].text(0.02, 0.05, 'CMS: $10$-$30$\%', transform=axes[1].transAxes)
axes[0].set_ylabel(r'$R_{\rho}$')
plt.savefig("../../plots/compare_coll_kernels/jet_shape_ratio.png", dpi=200)
plt.show()
