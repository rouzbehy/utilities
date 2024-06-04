## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.lines import Line2D
## my custom modules
import util
import jetDicts as my_jdicts
from COLORS import colors_comparison as colors
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 30,
    "lines.linewidth": 4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size" : 8,
    "ytick.major.size" : 8,
    "xtick.minor.size" : 4,
    "ytick.minor.size" : 4,
    "axes.spines.right": True,
    "axes.spines.top" : True,
    "legend.frameon":False
}
plt.rcParams.update(my_rcParams)
pTlow, pThigh  = 60, 300

## read in the experimental RAA:
cms = {}
cents = {'00-05':'0_5'}
radii = ['0p4']
for R in radii:
    cms[R] = pd.read_csv(f'../../../exp_data/jets/CMS/jetscape/cent_00-05_R_{R}.csv',comment='#').rename(columns=my_jdicts.colnames_CMS_RAA)

"""
    Proton-Proton baselines
"""
#### read in the calculated pp jet spectrum for the new formation calculations MARTINI
pp_loc  = "../../calcs/pp_2p76/jet_spec.csv"
pp_data_ftime = pd.read_csv(pp_loc, comment='#')
pp_data_ftime['pT'] = 0.5 * (pp_data_ftime['pTmin'] + pp_data_ftime['pTmax'])
pp_data_ftime['dpT'] = pp_data_ftime['pTmax'] - pp_data_ftime['pTmin']
pp_data_ftime = pp_data_ftime[pp_data_ftime['pT'].between(pTlow, pThigh)]

## read in the calculated pp jet spectrum for the non-formation time MARTINI
pp_loc  = "../../../v2/proton_proton/data/jet_spectra.csv"
pp_data_no_ftime = pd.read_csv(pp_loc, comment='#')
pp_data_no_ftime['pT'] = 0.5 * (pp_data_no_ftime['pTmin'] + pp_data_no_ftime['pTmax'])
pp_data_no_ftime['dpT'] = pp_data_no_ftime['pTmax'] - pp_data_no_ftime['pTmin']
pp_data_no_ftime = pp_data_no_ftime[pp_data_no_ftime['pT'].between(pTlow, pThigh)]

## read in the pp jet spectrum for the jetscape work
pp_loc = "../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/"
pp_loc += "pp/pp_2760_jet_spec_jet_rad_{R}_2.00.csv"
pp_data_jetscape = {}
for rval in radii:
    r = rval.replace('p','.')
    loc = pp_loc.format(R=r)
    dat = pd.read_csv(loc, comment='#')
    dat['pT'] = 0.5*(dat['ptmin'] + dat['ptmax'])
    dat['dpT'] = dat['ptmax'] - dat['ptmin']
    dat = dat[dat['pT'].between(pTlow, pThigh)]
    pp_data_jetscape[r] = dat

"""
    Heavy ion collision results
"""
## Read in the AA jet results. First the new calculations with formation time
loc = "../../calcs/final_runs/rset_1/2p76/0_5/jet_spectra.csv"
martini_form_time = pd.read_csv(loc, comment='#')
martini_form_time['pT'] = 0.5 * (martini_form_time['pTmin'] + martini_form_time['pTmax'])
martini_form_time['dpT'] = martini_form_time['pTmax'] - martini_form_time['pTmin']
martini_form_time = martini_form_time[martini_form_time['pT'].between(pTlow, pThigh)]

# Read in the AA jet results for MARTINI without formation time
loc = "../../../v2/production_v2/martini_results/final_PbPb_2p76"
loc += "/rset_1/cent_0_5/jet_spectra.csv"
martini_nform_time = pd.read_csv(loc, comment='#')
martini_nform_time['pT'] = 0.5 * (martini_nform_time['pTmin'] + martini_nform_time['pTmax'])
martini_nform_time['dpT'] = martini_nform_time['pTmax'] - martini_nform_time['pTmin']
martini_nform_time = martini_nform_time[martini_nform_time['pT'].between(pTlow, pThigh)]

# Read in the AA jet results for MARTINI in JETSCAPE:
aa_loc = "../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/martini/"
aa_loc += "PbPb_2760/PbPb2760_00-05_jet_spec_jet_rad_{R}_2.00.csv"
aa_data_jetscape = {}
for rval in radii:
    r = rval.replace('p','.')
    loc = aa_loc.format(R=r)
    dat = pd.read_csv(loc, comment='#')
    dat['pT'] = 0.5*(dat['ptmin'] + dat['ptmax'])
    dat['dpT'] = dat['ptmax'] - dat['ptmin']
    dat = dat[dat['pT'].between(pTlow, pThigh)]
    aa_data_jetscape[r] = dat

# Read in the AA jet results for MARTINI (with new kappa0) in jetscape
aa_loc = "../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/"
aa_loc += "martini_new_kappas/PbPb2760_00-05_jet_spec_jet_rad_{R}_2.00.csv"
aa_data_jetscape_k1 = {}
for rval in radii:
    r = rval.replace('p','.')
    loc = aa_loc.format(R=r)
    dat = pd.read_csv(loc, comment='#')
    dat['pT'] = 0.5*(dat['ptmin'] + dat['ptmax'])
    dat['dpT'] = dat['ptmax'] - dat['ptmin']
    dat = dat[dat['pT'].between(pTlow, pThigh)]
    aa_data_jetscape_k1[r] = dat 
# Read in the AA jet results for MARTINI (with new kappa0) in jetscape
aa_loc = "../../../../jetscape_project/v2/jetscape_data/sqrt_s_2760/"
aa_loc += "martini_new_kap_matter_vac/PbPb2760_00-05_jet_spec_jet_rad_{R}_2.00.csv"
aa_data_jetscape_k1_vac = {}
for rval in radii:
    r = rval.replace('p','.')
    loc = aa_loc.format(R=r)
    dat = pd.read_csv(loc, comment='#')
    dat['pT'] = 0.5*(dat['ptmin'] + dat['ptmax'])
    dat['dpT'] = dat['ptmax'] - dat['ptmin']
    dat = dat[dat['pT'].between(pTlow, pThigh)]
    aa_data_jetscape_k1_vac[r] = dat 

"""
    Now plot them all
"""
fig, axes = plt.subplots(1,len(radii), gridspec_kw={'top':0.985,'bottom':0.1,
                                      'left':0.08,'right':0.98,
                                      'hspace':0.05,'wspace':0.07},
                         sharex=True, figsize=(16,9), sharey=True)


for i, r in enumerate(radii):
    ax = axes[i] if len(radii) > 1 else axes
    rtxt = r'$R=$' + r.replace('p','.')
    ax.text(0.05, 0.8, rtxt, transform=ax.transAxes)
    util.plot_expr_data_on_axis(ax,cms[r],'*')

    ## Plot the martini without formation time
    pp, AA = pp_data_no_ftime, martini_nform_time
    x = pp['pT']
    Rtag = f"N{r}"
    DRtag = f"dN{r}"
    y, dy = AA[Rtag], AA[DRtag]
    raa = y/pp[Rtag]
    draa = raa*np.sqrt(dy*dy/(y*y) + pp[DRtag]*pp[DRtag]/(pp[Rtag]*pp[Rtag]))
    ax.plot(x, raa, color=colors[0]) 
    ax.fill_between(x, raa+draa, raa-draa, color=colors[0], alpha=0.2)
    ## Plot the martini _with_ formation time
    pp, AA = pp_data_ftime, martini_form_time
    x = pp['pT']
    Rtag = f"N{r}"
    DRtag = f"dN{r}"
    y, dy = AA[Rtag], AA[DRtag]
    raa = y/pp[Rtag]
    draa = raa*np.sqrt(dy*dy/(y*y) + pp[DRtag]*pp[DRtag]/(pp[Rtag]*pp[Rtag]))
    ax.plot(x, raa, color=colors[1]) 
    ax.fill_between(x, raa+draa, raa-draa, color=colors[1], alpha=0.2)
    ## Plot MATTER + MARTINI (with Chanwook Kappa set)
    jetscape_r = r.replace('p','.')
    pp, AA = pp_data_jetscape[jetscape_r], aa_data_jetscape[jetscape_r]
    x = pp['pT']
    Rtag = f"Ncut"
    DRtag = f"dNcut"
    y, dy = AA[Rtag], AA[DRtag]
    raa = y/pp[Rtag]
    draa = raa*np.sqrt(dy*dy/(y*y) + pp[DRtag]*pp[DRtag]/(pp[Rtag]*pp[Rtag]))
    ax.plot(x, raa, color=colors[2], linestyle='dashed') 
    ax.fill_between(x, raa+draa, raa-draa, color=colors[2], alpha=0.2)
    ## Plot MATTER + MARTINI (with old fitted kappa set, matter in vacuum)
    AA = aa_data_jetscape_k1_vac[jetscape_r]
    y, dy = AA[Rtag], AA[DRtag]
    raa = y/pp[Rtag]
    draa = raa*np.sqrt(dy*dy/(y*y) + pp[DRtag]*pp[DRtag]/(pp[Rtag]*pp[Rtag]))
    ax.plot(x, raa, color=colors[3], linestyle='dotted') 
    ax.fill_between(x, raa+draa, raa-draa, color=colors[3], alpha=0.2)
    ## Plot MATTER + MARTINI (with old fitted kappa set)
    AA = aa_data_jetscape_k1[jetscape_r]
    y, dy = AA[Rtag], AA[DRtag]
    raa = y/pp[Rtag]
    draa = raa*np.sqrt(dy*dy/(y*y) + pp[DRtag]*pp[DRtag]/(pp[Rtag]*pp[Rtag]))
    ax.plot(x, raa, color=colors[4]) 
    ax.fill_between(x, raa+draa, raa-draa, color=colors[4], alpha=0.2)



line_colors = {r'MARTINI (no $\tau_f$)'              : (colors[0], 'solid'), 
               r'MARTINI (with $\tau_f$)'            : (colors[1], 'solid'),
               r'MATTER + MARTINI ($\kappa_0$)'      : (colors[2], 'dashed'),
               r'MATTER (vac) + MARTINI ($\kappa_1$)': (colors[3], 'dotted'),
               r'MATTER + MARTINI ($\kappa_1$)'      : (colors[4], 'solid')}

theor_hands = [Line2D([],[],color=css['black'], label='CMS $|\eta|<1.0$', marker='*', markersize=10, linestyle="None")]
for (l, c) in line_colors.items():
    theor_hands.append(Line2D([],[],label=l,color=c[0],linestyle=c[1]))

if len(radii) > 1:
    axes[0].set_ylim(bottom=-0.09, top=1.2)
    for ax in axes:
        ax.set_xlabel(r'$p_T$ (GeV)')
    axes[0].set_ylabel(r'$R^{\mathrm{jet}}_{\mathrm{AA}}$', fontsize=30)
    axes[0].text(0.08, 0.1, 'PbPb, $\sqrt{s}=2.76$ATeV'+'\n'+'$|\eta|<2.0$', transform=axes[0].transAxes)
    axes[2].legend(handles=theor_hands, loc='lower center', ncols=1, handletextpad=1, fontsize=20)
else:
    axes.set_ylim(bottom=-0.09, top=1.2)
    axes.set_xlabel(r'$p_T$ (GeV)')
    axes.set_ylabel(r'$R^{\mathrm{jet}}_{\mathrm{AA}}$', fontsize=30)
    axes.text(0.05, 0.85, 'PbPb, $\sqrt{s}=2.76$ ATeV $|\eta|<2.0$', transform=axes.transAxes)
    axes.legend(handles=theor_hands, loc='upper right', ncols=1, handletextpad=1, fontsize=30)

#plt.savefig("../../plots/formation_time_effect/jet_raa_formation_time_effect.png", dpi=200)
plt.show()