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
import dictionaries as my_hdicts
import jetDicts as my_jdicts

plt.rcParams.update(util.my_rcParams)
from COLORS import colors, shower_effects
pTlow  = 60
pThigh = 300

def get_spectrum(fname):
    tmp = pd.read_csv(fname, comment='#')
    if 'pTmin' in tmp:
        tmp['pT'] = 0.5 * (tmp['pTmin'] + tmp['pTmax'])
        tmp['dpT'] = tmp['pTmax'] - tmp['pTmin']
    else:
        tmp['pT'] = 0.5 * (tmp['ptmin'] + tmp['ptmax'])
        tmp['dpT'] = tmp['ptmax'] - tmp['ptmin']
    result = tmp[tmp['pT'].between(pTlow, pThigh)]
    return result

## read in the experimental RAA:
cms = {}
for R in ['0p2','0p3','0p4']:
    cms[R] = pd.read_csv(f'../../../exp_data/jets/CMS/Table_JETRAA_{R}.csv',comment='#').rename(columns=my_jdicts.colnames_CMS_RAA)

rate_colours = my_hdicts.rate_colours
rate_names = {1:'LO', 2:'NLO',3:'NP'}

#### read in standalone stuff
pp_loc  = '../martini_results/pp/jet_spectra.csv'
standalone_pp = get_spectrum(pp_loc)
AA_loc= '../martini_results/final_PbPb_2p76/rset_1/cent_0_5//jet_spectra.csv'
standalone_AA = get_spectrum(AA_loc)
#### read in the jetscape runs
r_vs_R = {'0p2':'0.2','0p3':'0.3','0p4':'0.4'}
prefix = '../../../../jetscape_project/v2/jetscape_data'
fname_AA = 'PbPb2760_00-05_jet_spec_jet_rad_{r}_2.00.csv'
fname_pp = 'pp_2760_jet_spec_jet_rad_{r}_2.00.csv'
spec_dir_1 = '/sqrt_s_2760/martini/PbPb_2760/'
spec_dir_2 = '/sqrt_s_2760/martini_new_kappas/'
spec_dir_3 = '/sqrt_s_2760/martini_new_kap_matter_vac/'
spec_dir_pp = '/max_time/maxT_200_highstat/'
jetscape_pp =          {r: get_spectrum(prefix+spec_dir_pp+fname_pp.format(r=R)) for (r,R) in r_vs_R.items()}
jetscape_AA_old_kset = {r: get_spectrum(prefix+spec_dir_1+fname_AA.format(r=R)) for (r,R) in r_vs_R.items()}
jetscape_AA_new_kset = {r: get_spectrum(prefix+spec_dir_2+fname_AA.format(r=R)) for (r,R) in r_vs_R.items()}
jetscape_AA_vac_matter = {r: get_spectrum(prefix+spec_dir_3+fname_AA.format(r=R)) for (r,R) in r_vs_R.items()}


fig, axes = plt.subplots(3,1,gridspec_kw={'top':0.985,'bottom':0.11,
                                      'left':0.08,'right':0.7,
                                      'hspace':0.1,'wspace':0.09},
                              sharex=True, figsize=(16,9), sharey=True)

axes = axes.flatten()
column='N{R}'
dcolumn='dN{R}'
for ir, r in enumerate(cms):
    ax = axes[ir]
    util.plot_expr_data_on_axis(ax,cms[r],'*')
    rtxt = r'$R=$' + r.replace('p','.')
    ax.text(0.45, 0.1, rtxt, transform=ax.transAxes)


for ir, R in enumerate(['0p2','0p3','0p4']):
    ax = axes[ir]
    y_aa, dy_aa = standalone_AA[column.format(R=R)], standalone_AA[dcolumn.format(R=R)]
    y_pp, dy_pp = standalone_pp[column.format(R=R)], standalone_pp[dcolumn.format(R=R)]
    raa = y_aa/y_pp
    draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
    pt_1 = standalone_AA['pT']
    ax.plot(pt_1, raa, color=colors[0])
    ax.fill_between(pt_1, raa+draa, raa-draa, color=colors[0], alpha=0.2)

    y_aa, dy_aa = jetscape_AA_old_kset[R]['Ncut'], jetscape_AA_old_kset[R]['dNcut']
    y_pp, dy_pp = jetscape_pp[R]['Ncut'], jetscape_pp[R]['dNcut']
    raa = y_aa/y_pp
    draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
    pt_2 = jetscape_AA_old_kset[R]['pT']
    ax.plot(pt_2, raa, color=colors[1])
    ax.fill_between(pt_2, raa+draa, raa-draa, color=colors[1], alpha=0.2)

    y_aa, dy_aa = jetscape_AA_new_kset[R]['Ncut'], jetscape_AA_new_kset[R]['dNcut']
    y_pp, dy_pp = jetscape_pp[R]['Ncut'], jetscape_pp[R]['dNcut']
    raa = y_aa/y_pp
    draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
    pt_2 = jetscape_AA_old_kset[R]['pT']
    ax.plot(pt_2, raa, color=colors[2])
    ax.fill_between(pt_2, raa+draa, raa-draa, color=colors[2], alpha=0.2)

    y_aa, dy_aa = jetscape_AA_vac_matter[R]['Ncut'], jetscape_AA_vac_matter[R]['dNcut']
    y_pp, dy_pp = jetscape_pp[R]['Ncut'], jetscape_pp[R]['dNcut']
    raa = y_aa/y_pp
    draa = raa*np.sqrt(dy_aa*dy_aa/(y_aa*y_aa) + dy_pp*dy_pp/(y_pp*y_pp))
    pt_2 = jetscape_AA_old_kset[R]['pT']
    ax.plot(pt_2, raa, color=colors[3])
    ax.fill_between(pt_2, raa+draa, raa-draa, color=colors[3], alpha=0.2)

theor_hands = [Line2D([],[],label=l,color=c) for (l,c) in shower_effects.items()]

theor_hands.append(Line2D([],[],color=css['black'], label='CMS, $|\eta^{\mathrm{jet}}|<2.0$ (2017)', marker='*', markersize=10, linestyle="None"))

axes[0].legend(loc='lower right', handles=theor_hands, ncol=1, bbox_to_anchor=(1.5,-0.35), fontsize=20)
axes[1].text(1.21, 0.65, r'$0$-$5\%$', transform=axes[1].transAxes)
ax = axes[2]
ax.set_ylim(bottom=-0.01, top=1.01)
ax.set_xlabel(r'$p^{\mathrm{jet}}_T$ (GeV)')#, fontsize=30)
for ax in axes:
    ax.set_ylabel(r'$R^{\mathrm{jet}}_{\mathrm{AA}}$')#, fontsize=30)

plt.show()
