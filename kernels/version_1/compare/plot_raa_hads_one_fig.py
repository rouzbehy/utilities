## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
#from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.interpolate import interp1d as interpolate
## my custom modules
import util
import dictionaries as my_dicts

plt.rcParams.update(util.my_rcParams)
from COLORS import shower_effects, colors
pT_low_cut, pT_high_cut = 10, 100
def get_spectrum(fname):
    tmp = pd.read_csv(fname, comment='#')
    print(tmp.head(2))
    if 'pTmin' in tmp:
        tmp['pT'] = 0.5 * (tmp['pTmin'] + tmp['pTmax'])
        tmp['dpT'] = tmp['pTmax'] - tmp['pTmin']
    else:
        tmp['pT'] = 0.5 * (tmp['ptmin'] + tmp['ptmax'])
        tmp['dpT'] = tmp['ptmax'] - tmp['ptmin']
    result = tmp[tmp['pT'].between(pT_low_cut, pT_high_cut)]
    return result


## read in the experimental RAA:
cms     = pd.read_csv('../../../exp_data/charged/CMS/Table5_1.csv',comment='#').rename(columns=my_dicts.colnames_raa_cms)
alice   = pd.read_csv('../../../exp_data/charged/ALICE/PbPb_RAA_cent_00_05_eta_0p8.csv',comment='#').rename(columns=my_dicts.colnames_raa_alice)
atlas_1 = pd.read_csv("../../../exp_data/charged/ATLAS/Table49.csv",comment='#').rename(columns=my_dicts.colnames_raa_atlas)
atlas_2 = pd.read_csv("../../../exp_data/charged/ATLAS/Table33.csv",comment='#').rename(columns=my_dicts.colnames_raa_atlas) 
atlas_1 = atlas_1[atlas_1['x'].between(8, 104)]
atlas_2 = atlas_2[atlas_2['x'].between(8, 104)]


## READ IN THE CALCULATION:

pp_loc  = '../martini_results/pp/hadron_spectra.csv'
AA_loc = '../martini_results/final_PbPb_2p76/rset_1/cent_0_5/hadron_spectra.csv'
### read in the calculated pp charged hadron spec
# Stand alone martini
standalone_pp = get_spectrum(pp_loc) 
standalone_AA = get_spectrum(AA_loc)
# martini in jetscape
prefix = '../../../../jetscape_project/v2/jetscape_data/'
fname_AA = 'PbPb2760_00-05_charged_hadrons_eta_cut_1.csv'
fname_pp = 'pp_2760_charged_hadrons_eta_cut_1.csv'
spec_dir_1 = 'sqrt_s_2760/martini/PbPb_2760/'
spec_dir_2 = 'sqrt_s_2760/martini_new_kappas/'
spec_dir_3 = '/sqrt_s_2760/martini_new_kap_matter_vac/'
spec_dir_pp = '/max_time/maxT_200_highstat/'

jetscape_old_kappas_AA = get_spectrum(prefix+spec_dir_1+fname_AA)
jetscape_new_kappas_AA = get_spectrum(prefix+spec_dir_2+fname_AA)
jetscape_vac_matter_AA = get_spectrum(prefix+spec_dir_3+fname_AA) 
jetscape_pp = get_spectrum(prefix+spec_dir_pp+fname_pp)

fig, ax = plt.subplots(1,1,gridspec_kw={'top':0.985,'bottom':0.11,
                                      'left':0.08,'right':0.98,
                                      'hspace':0.11,'wspace':0.09},
                              sharex=True, figsize=(16,9), sharey=True)


cms = cms[cms['x'].between(pT_low_cut, pT_high_cut)]
alice = alice[alice['x'].between(pT_low_cut, pT_high_cut)]

#for ax in [ax1, ax2]:
util.plot_expr_data_on_axis(ax,cms,'*')
util.plot_expr_data_on_axis(ax,alice,'o')
util.plot_expr_data_on_axis(ax,atlas_1,'P')
util.plot_expr_data_on_axis(ax,atlas_2,'s')

## Now do the RAA calculation
y_aa, dy_aa = standalone_AA['N'], standalone_AA['dN']
y_pp, dy_pp = standalone_pp['N'], standalone_pp['dN']
raa = y_aa/y_pp
draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
pt_1 = standalone_AA['pT']
ax.plot(pt_1, raa, color=colors[0])
ax.fill_between(pt_1, raa+draa, raa-draa, color=colors[0], alpha=0.2)

y_aa, dy_aa = jetscape_old_kappas_AA['N'], jetscape_old_kappas_AA['dN']
y_pp, dy_pp = jetscape_pp['N'], jetscape_pp['dN']
raa = y_aa/y_pp
draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
pt_1 = jetscape_pp['pT']
ax.plot(pt_1, raa, color=colors[1])
ax.fill_between(pt_1, raa+draa, raa-draa, color=colors[1], alpha=0.2)

## Now do the RAA calculation
y_aa, dy_aa = jetscape_new_kappas_AA['N'], jetscape_new_kappas_AA['dN']
y_pp, dy_pp = jetscape_pp['N'], jetscape_pp['dN']
raa = y_aa/y_pp
draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
pt_1 = jetscape_pp['pT']
ax.plot(pt_1, raa, color=colors[2])
ax.fill_between(pt_1, raa+draa, raa-draa, color=colors[2], alpha=0.2)

## Now do the RAA calculation
y_aa, dy_aa = jetscape_vac_matter_AA['N'], jetscape_vac_matter_AA['dN']
y_pp, dy_pp = jetscape_pp['N'], jetscape_pp['dN']
raa = y_aa/y_pp
draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
pt_1 = jetscape_pp['pT']
ax.plot(pt_1, raa, color=colors[3])
ax.fill_between(pt_1, raa+draa, raa-draa, color=colors[3], alpha=0.2)

theor_hands = [Line2D([],[],label=l,color=c) for (l,c) in shower_effects.items()]
expt_hands = [Line2D([],[],color=css['black'], label='CMS $|\eta|<1.0$', marker='*', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ALICE $|\eta|<0.8$', marker='o', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ATLAS $|\eta|<1.0$', marker='P', markersize=10, linestyle="None"),
Line2D([],[],color=css['black'], label='ATLAS $|\eta|<2.0$', marker='s', markersize=10, linestyle="None")]

artist = ax.legend(loc='upper right', handles=expt_hands, ncol=2)#, bbox_to_anchor=(0.95, 0.2))#,prop={'size':20})
ax.add_artist(artist)

ax.legend(loc='lower right', handles=theor_hands, ncol=1)#,fontsize=15)

ax.set_ylim(top=1.21,bottom=-0.01)
ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30)
ax.set_ylabel(r'$R^{\mathrm{h}^{\pm}}_{\mathrm{AA}}$', fontsize=30)
ax.text(0.05, 0.9, r'Pb-Pb \@ $\sqrt{s}=2.76$ ATeV' +'\n'+r'$|\eta|<1$, $0$-$5\%$', transform=ax.transAxes) 
plt.show()
