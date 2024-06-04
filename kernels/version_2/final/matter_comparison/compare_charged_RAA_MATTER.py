## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.lines import Line2D
## my custom modules
import util
import dictionaries as my_dicts
from COLORS import colors_comparison as colors
from COLORS import experiment_color
plt.rcParams.update(util.my_rcParams)
pT_low_cut, pT_high_cut = 33, 100
"""
    PLOT CHARGED HADRON RAA

    Plot 0-5\% centrality stuff.
    Compare:
        1. MARTINI new fit kappas
        2. MARTINI old fit kappas
        3. MARTINI (old kappas) + MATTER (vac)
        4. MARTINI (old kappas) + MATTER (medium)
        5. MARTINI (Chanwook kappas) + MATTER (medium)
"""
## read in the experimental RAA:
exp_loc = '../../../exp_data/charged/'
cms_cents = {1:'0_5'}
cms     = {c:pd.read_csv(exp_loc+f'CMS/Table5_{i}.csv',comment='#').rename(columns=my_dicts.colnames_raa_cms) 
           for i,c in cms_cents.items()}

alice_cents = {'00_05':'0_5'}
alice   = {c:pd.read_csv(exp_loc+f'ALICE/PbPb_RAA_cent_{i}_eta_0p8.csv',comment='#').rename(columns=my_dicts.colnames_raa_alice)
           for i,c in alice_cents.items()}

atlas_1 = pd.read_csv(exp_loc+"ATLAS/Table49.csv",comment='#').rename(columns=my_dicts.colnames_raa_atlas)
atlas_2_cents = {33:'0_5'}
atlas_2 = {c:pd.read_csv(exp_loc+f"ATLAS/Table{i}.csv",comment='#').rename(columns=my_dicts.colnames_raa_atlas)
           for i, c in atlas_2_cents.items()}

cms = {c: itm[itm['x'].between(pT_low_cut,pT_high_cut)] for c, itm in cms.items()}
alice = {c:itm[itm['x'].between(pT_low_cut,pT_high_cut)] for c, itm in alice.items()}

atlas_1 = atlas_1[atlas_1['x'].between(pT_low_cut,pT_high_cut)]
atlas_2 = {c:itm[itm['x'].between(pT_low_cut,pT_high_cut)] for c, itm in atlas_2.items()}

#### read in the calculated pp hadron spec for the new formation calculations MARTINI
pp_loc  = "../../calcs/pp_2p76/chgd_spec.csv"
pp_data_ftime = pd.read_csv(pp_loc, comment='#')
pp_data_ftime['pT'] = 0.5 * (pp_data_ftime['pTmin'] + pp_data_ftime['pTmax'])
pp_data_ftime['dpT'] = pp_data_ftime['pTmax'] - pp_data_ftime['pTmin']
pp_data_ftime = pp_data_ftime[pp_data_ftime['pT'].between(pT_low_cut, pT_high_cut)]

## read in the calculated pp hadron spectrum for the non-formation time MARTINI
pp_loc  = "../../../v2/proton_proton/data/hadron_spectra.csv"
pp_data_no_ftime = pd.read_csv(pp_loc, comment='#')
pp_data_no_ftime['pT'] = 0.5 * (pp_data_no_ftime['pTmin'] + pp_data_no_ftime['pTmax'])
pp_data_no_ftime['dpT'] = pp_data_no_ftime['pTmax'] - pp_data_no_ftime['pTmin']
pp_data_no_ftime = pp_data_no_ftime[pp_data_no_ftime['pT'].between(pT_low_cut, pT_high_cut)]

## read in the pp jet spectrum for the jetscape work
pp_loc = "../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/"
pp_loc += "pp/pp_2760_charged_hadrons_eta_cut_1.csv"
pp_data_jetscape = {}
dat = pd.read_csv(pp_loc, comment='#')
dat['pT'] = 0.5*(dat['pTmin'] + dat['pTmax'])
dat['dpT'] = dat['pTmax'] - dat['pTmin']
dat = dat[dat['pT'].between(pT_low_cut, pT_high_cut)]
pp_data_jetscape = dat


## Read in the AA charged hadron results. First the new calculations with formation time
loc = "../../calcs/final_runs/rset_1/2p76/{cent}/hadron_spectra.csv"
martini_form_time = {}
for c in ['0_5']:
    dat = pd.read_csv(loc.format(cent=c), comment='#')
    dat['pT'] = 0.5 * (dat['pTmin'] + dat['pTmax'])
    dat['dpT'] = dat['pTmax'] - dat['pTmin']
    martini_form_time[c] = dat[dat['pT'].between(pT_low_cut, pT_high_cut)]

# Read in the AA jet results for MARTINI without formation time
loc = "../../../v2/production_v2/martini_results/final_PbPb_2p76"
loc += "/rset_1/cent_{c}/hadron_spectra.csv"
martini_nform_time = {}
for c in ['0_5']:
    dat = pd.read_csv(loc.format(c=c), comment='#')
    dat['pT']  = 0.5 * (dat['pTmin'] + dat['pTmax'])
    dat['dpT'] =        dat['pTmax'] - dat['pTmin']
    martini_nform_time[c] = dat[dat['pT'].between(pT_low_cut, pT_high_cut)]

# Read in the AA hadron results for MARTINI in JETSCAPE:
aa_loc = "../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/martini/"
aa_loc += "PbPb_2760/PbPb2760_00-05_charged_hadrons_eta_cut_1.csv"
dat = pd.read_csv(aa_loc, comment='#')
dat['pT'] = 0.5*(dat['pTmin'] + dat['pTmax'])
dat['dpT'] = dat['pTmax'] - dat['pTmin']
aa_data_jetscape_chanwook = dat[dat['pT'].between(pT_low_cut, pT_high_cut)]

# Read in the AA hadron results for MARTINI in JETSCAPE, old fit kappa:
aa_loc = "../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/"
aa_loc += "martini_new_kappas/PbPb2760_00-05_charged_hadrons_eta_cut_1.csv"
dat = pd.read_csv(aa_loc, comment='#')
dat['pT'] = 0.5*(dat['pTmin'] + dat['pTmax'])
dat['dpT'] = dat['pTmax'] - dat['pTmin']
aa_data_jetscape_k1 = dat[dat['pT'].between(pT_low_cut, pT_high_cut)]

# Read in the AA hadron results for MARTINI in JETSCAPE, old fit kappa:
aa_loc = "../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/"
aa_loc += "martini_new_kap_matter_vac/PbPb2760_00-05_charged_hadrons_eta_cut_1.csv"
dat = pd.read_csv(aa_loc, comment='#')
dat['pT'] = 0.5*(dat['pTmin'] + dat['pTmax'])
dat['dpT'] = dat['pTmax'] - dat['pTmin']
aa_data_jetscape_k1_vac = dat[dat['pT'].between(pT_low_cut, pT_high_cut)]
    
fig, ax = plt.subplots(1,1,gridspec_kw={'top':0.985,'bottom':0.11,
                                      'left':0.08,'right':0.98,
                                      'hspace':0.11,'wspace':0.02},
                            figsize=(16,9))


for cent in alice:
    alice_spectrum = alice[cent]
    atlas_spectrum = atlas_2[cent]
    util.plot_expr_data_on_axis(ax, alice_spectrum, 'P', s=150, color=experiment_color, face=experiment_color)
    util.plot_expr_data_on_axis(ax, atlas_spectrum, 's', s=150, color=experiment_color, face=experiment_color)

for cent in cms:
    cms_spectrum = cms[cent]
    util.plot_expr_data_on_axis(ax,cms_spectrum,'*', s=150, color=experiment_color, face=experiment_color)

util.plot_expr_data_on_axis(ax,atlas_1,'o', color=experiment_color, face=experiment_color)

for cent in alice:
    ## old kappa set
    pp, AA = pp_data_no_ftime, martini_nform_time[cent]
    x = pp['pT']
    y1, dy1 = AA['N'], AA['dN']
    raa1 = y1/pp['N'] 
    draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + pp['dN']*pp['dN']/(pp['N']*pp['N']))
    ax.plot(x, raa1, color=colors[0]) 
    ax.fill_between(x, raa1+draa1, raa1-draa1, color=colors[0], alpha=0.2)
    ## new kappa set
    pp, AA = pp_data_ftime, martini_form_time[cent]
    x = pp['pT']
    y1, dy1 = AA['N'], AA['dN']
    raa1 = y1/pp['N'] 
    draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + pp['dN']*pp['dN']/(pp['N']*pp['N']))
    ax.plot(x, raa1, color=colors[1]) 
    ax.fill_between(x, raa1+draa1, raa1-draa1, color=colors[1], alpha=0.2)
    ## jetscape, Chanwook kappas
    pp, AA = pp_data_jetscape, aa_data_jetscape_chanwook
    x = pp['pT']
    y1, dy1 = AA['N'], AA['dN']
    raa1 = y1/pp['N'] 
    draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + pp['dN']*pp['dN']/(pp['N']*pp['N']))
    ax.plot(x, raa1, color=colors[2], linestyle='dashed') 
    ax.fill_between(x, raa1+draa1, raa1-draa1, color=colors[2], alpha=0.2)
    ## jetscape, MARTINI old fit, MATTER in vacuum
    pp, AA = pp_data_jetscape, aa_data_jetscape_k1_vac
    x = pp['pT']
    y1, dy1 = AA['N'], AA['dN']
    raa1 = y1/pp['N'] 
    draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + pp['dN']*pp['dN']/(pp['N']*pp['N']))
    ax.plot(x, raa1, color=colors[3], linestyle='dotted') 
    ax.fill_between(x, raa1+draa1, raa1-draa1, color=colors[3], alpha=0.2)
    ## jetscape, MARTINI old fit, MATTER in medium
    pp, AA = pp_data_jetscape, aa_data_jetscape_k1
    x = pp['pT']
    y1, dy1 = AA['N'], AA['dN']
    raa1 = y1/pp['N'] 
    draa1 = raa1*np.sqrt(dy1*dy1/(y1*y1) + pp['dN']*pp['dN']/(pp['N']*pp['N']))
    ax.plot(x, raa1, color=colors[4]) 
    ax.fill_between(x, raa1+draa1, raa1-draa1, color=colors[4], alpha=0.2)


line_colors = {r'MARTINI (no $\tau_f$)'              : (colors[0], 'solid'), 
               r'MARTINI (with $\tau_f$)'            : (colors[1], 'solid'),
               r'MATTER + MARTINI ($\kappa_0$)'      : (colors[2], 'dashed'),
               r'MATTER (vac) + MARTINI ($\kappa_1$)': (colors[3], 'dotted'),
               r'MATTER + MARTINI ($\kappa_1$)'      : (colors[4], 'solid')}

theo_hands = [Line2D([],[],label=l,color=c[0], linestyle=c[1]) for (l,c) in line_colors.items()]
expt_hands = [Line2D([],[],color=css['black'], label='CMS $|\eta|<1.0$', marker='*', markersize=10, linestyle="None"),
              Line2D([],[],color=css['black'], label='ALICE $|\eta|<0.8$', marker='P', markersize=10, linestyle="None"),
              Line2D([],[],color=css['black'], label='ATLAS $|\eta|<1.0$', marker='o', markersize=10, linestyle="None"),
              Line2D([],[],color=css['black'], label='ATLAS $|\eta|<2.0$', marker='s', markersize=10, linestyle="None")]

artist = ax.legend(loc='lower right', handles=expt_hands, 
                        ncol=1, bbox_to_anchor=(0.98, 0.01),handletextpad=0.1)
ax.add_artist(artist)

ax.legend(loc='lower left', handles=theo_hands, handletextpad=0.1, fontsize=28)

ax.set_ylim(top=1.21,bottom=-0.19)
for cent in alice:
    ax.set_xlabel(r'$p_T$ (GeV)')
    low, high = cent.split('_')
    centrality = f'${low}$-${high}$\%'
    ax.text(0.8,0.8,centrality, transform=ax.transAxes)

ax.text(0.05, 0.85, r'Pb-Pb,$\sqrt{s}=2.76$ ATeV' +'\n'+r'$|\eta|<1$', transform=ax.transAxes)
ax.set_ylabel(r'$R^{\mathrm{h}^{\pm}}_{\mathrm{AA}}$')
#plt.savefig("../../plots/formation_time_effect/charged_hadron_raa_formation_time_effect.png", dpi=200)
plt.show()
