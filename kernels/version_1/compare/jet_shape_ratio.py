#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
import dictionaries as my_dicts
import util
plt.rcParams.update(util.my_rcParams)
inel_Xsec= 62.03948 

def get_shape(fname):
    calcs = pd.read_csv(fname,comment='#')
    njet = 1
    with open(fname, 'r') as f:
        line = f.readline()
        line = line.split(' ')[-1]
        njet = float(line)

    dat = calcs[calcs['rmax'] < 0.31]
    delta_r = dat['rmax'] - dat['rmin']
    dr = delta_r.to_list()
    r = 0.5*(dat['rmax'] + dat['rmin'])
    if 'jetscape' not in fname:
        rho  = dat['N']  /(njet)
        drho = dat['dN'] /(njet)
    else:
        rho = dat['wcut']/(njet)
        drho = dat['dwcut']/(njet)
    norm = sum(rho.to_list())#
    rho_normed = rho  /(delta_r * norm)
    drho_normed = drho/(delta_r * norm)
    r = r.to_list()
    return pd.DataFrame({'y':rho_normed,'dy':drho_normed,'r':r,'dr':dr})

from COLORS import colors, shower_effects

## read in the experimental results
data_fname = "../../../exp_data/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_00-10_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
tmp = np.loadtxt(data_fname,comments='#',unpack=True, delimiter='\t')
cms_data = pd.DataFrame({'x':tmp[0],'y':tmp[1],'dx':tmp[2],'dy':tmp[3]})

## Get the shapes from the stand-alone runs
pp_fname = '../martini_results/pp/jet_shape.csv'
AA_fname = '../martini_results/final_PbPb_2p76/rset_1/cent_0_5/jet_shape.csv'
pp_shape_stdalone = get_shape(pp_fname)
AA_shape_stdalone = get_shape(AA_fname)

## Get shapes from JETSCAPE runs
prefix = '../../../../jetscape_project/v2/jetscape_data/'
fname_AA = 'PbPb2760_00-05_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv'
fname_pp = 'pp_2760_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv'
spec_dir_1 = 'sqrt_s_2760/martini/PbPb_2760/'
spec_dir_2 = 'sqrt_s_2760/martini_new_kappas/'
spec_dir_3 = 'sqrt_s_2760/martini_new_kap_matter_vac/'
spec_dir_pp = '/max_time/maxT_200_highstat/'

jetscape_old_kappas_AA = get_shape(prefix+spec_dir_1+fname_AA)
jetscape_new_kappas_AA = get_shape(prefix+spec_dir_2+fname_AA)
jetscape_vac_matter_AA = get_shape(prefix+spec_dir_3+fname_AA)
jetscape_pp = get_shape(prefix+spec_dir_pp+fname_pp)

## Plot the data:
fig, ax = plt.subplots(1, 1, gridspec_kw={'top':0.985,
                                'bottom':0.1, 'left':0.07,
                                'right':0.99, 'hspace':0.02,
                                'wspace':0.18}, figsize=(16,9), sharex=True, sharey=True)

errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                for x, delx, y, yerrlow, yerrhigh in
                zip(cms_data["x"], cms_data['dx'], cms_data["y"], cms_data["dy"], cms_data["dy"])]
# Create patch collection with specified colour/alpha
pc = PatchCollection(errorboxes, facecolor='black', edgecolor="black", alpha=0.4)##eca400
ax.add_collection(pc)
scatter = ax.scatter(cms_data['x'],cms_data['y'],color='black',marker='s',s=30, label='CMS (2014)')


## Now do the RAA calculation
y_aa, dy_aa = AA_shape_stdalone['y'], AA_shape_stdalone['dy']
y_pp, dy_pp = pp_shape_stdalone['y'], pp_shape_stdalone['dy']
raa = y_aa/y_pp
draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
pt_1 = pp_shape_stdalone['r']
ax.plot(pt_1, raa, color=colors[0])
ax.fill_between(pt_1, raa+draa, raa-draa, color=colors[0], alpha=0.2)


y_aa, dy_aa = jetscape_old_kappas_AA['y'], jetscape_old_kappas_AA['dy']
y_pp, dy_pp = jetscape_pp['y'], jetscape_pp['dy']
raa = y_aa/y_pp
draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
pt_1 = jetscape_pp['r']
ax.plot(pt_1, raa, color=colors[1])
ax.fill_between(pt_1, raa+draa, raa-draa, color=colors[1], alpha=0.2)

## Now do the RAA calculation
y_aa, dy_aa = jetscape_new_kappas_AA['y'], jetscape_new_kappas_AA['dy']
y_pp, dy_pp = jetscape_pp['y'], jetscape_pp['dy']
raa = y_aa/y_pp
draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
pt_1 = jetscape_pp['r']
ax.plot(pt_1, raa, color=colors[2])
ax.fill_between(pt_1, raa+draa, raa-draa, color=colors[2], alpha=0.2)

## Now do the RAA calculation
y_aa, dy_aa = jetscape_vac_matter_AA['y'], jetscape_vac_matter_AA['dy']
y_pp, dy_pp = jetscape_pp['y'], jetscape_pp['dy']
raa = y_aa/y_pp
draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
pt_1 = jetscape_pp['r']
ax.plot(pt_1, raa, color=colors[3])
ax.fill_between(pt_1, raa+draa, raa-draa, color=colors[3], alpha=0.2)


theor_hands = [Line2D([],[],label=l,color=c) for (l,c) in shower_effects.items()]
artist = ax.legend(loc='lower right', handles=theor_hands, ncol=1, fontsize=15)
ax.add_artist(artist)
ax.text(0.02,0.55,s=r'$0.3<|\eta|<2.0$'+'\n'+r'$100$ GeV $< p^{\mathrm{jet}}_T$' + '\n' + r'$0$-$5\%$', transform=ax.transAxes)
ax.text(0.02,0.75,s=r'$R=0.3$, Anti-$k_{T}$', transform=ax.transAxes)

#labels = [Line2D([],[],color=c,label=l) for (l,c) in momentum_cut_and_colors.items()]
expt_labels = [Line2D([],[],color='black',label=r'CMS (2014) $0$-$10\%$', marker='s')]
ax.legend(handles=expt_labels, loc='upper left')
ax.set_ylabel(r'$\rho(r)_{PbPb}/\rho(r)_{pp}$')
ax.set_xlabel(r'$r$')
#ax.axhline(1,linestyle='dotted',color='black')
plt.show()