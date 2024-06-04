import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pandas import read_csv
from numpy import pi, log, arange,exp, loadtxt,sqrt, array
from scipy.interpolate import interp1d
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import ticker as mticker
import dictionaries as ddicts
import jetDicts as jdicts
from util import plot_expr_data_on_axis
from COLORS import colors
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
    "xtick.minor.size" : 6,
    "ytick.minor.size" : 6,
    "axes.spines.right": False,
    "axes.spines.top" : False,
    "legend.frameon":False
}
plt.rcParams.update(my_rcParams)
colour = '#d64045'
inelXsec = 64.029
ptmin, ptmax = 2, 120
cms_hads = read_csv("../../../exp_data/charged/CMS/Table1.csv",
        comment='#').rename(columns=ddicts.colnames_CMS_spectra_had)
cms_hads = cms_hads[cms_hads['x'].between(ptmin,ptmax)]
cms_hads['y'] *= inelXsec

atlas_hads_eta_1 = read_csv('../../../exp_data/charged/ATLAS/ATLAS_ch_had_spec_eta_1.csv', comment='#')
atlas_hads_eta_1 = atlas_hads_eta_1[atlas_hads_eta_1['x'].between(ptmin,ptmax)]

atlas_hads_eta_2 = read_csv("../../../exp_data/charged/ATLAS/Table1.csv", comment='#').rename(columns=ddicts.colnames_ATLAS_pp_chHad_spec)
atlas_hads_eta_2 = atlas_hads_eta_2[atlas_hads_eta_2['x'].between(ptmin,ptmax)]

alice_hads = read_csv("../../../exp_data/charged/ALICE/ALICE_pp_charged_hadron_yield_eta_0p8.csv", comment='#').rename(columns=ddicts.colnames_ALICE_pp_chHad_spec)
alice_hads = alice_hads[alice_hads['x'].between(ptmin,ptmax)]


martini_hads = read_csv("../martini_results/pp/hadron_spectra.csv",comment='#')


exp_labels = [r'CMS $|\eta|<1.0$', r'ALICE $|\eta|<0.8$', r'ATLAS $|\eta|<1.0$', r'ATLAS $|\eta|<2.0$']

exp_data = {r'CMS $|\eta|<1.0$':cms_hads, 
            r'ALICE $|\eta|<0.8$':alice_hads,
            r'ATLAS $|\eta|<1.0$':atlas_hads_eta_1,
            r'ATLAS $|\eta|<2.0$':atlas_hads_eta_2}

markers = {r'CMS $|\eta|<1.0$'  :'o', 
           r'ALICE $|\eta|<0.8$':'s',
           r'ATLAS $|\eta|<1.0$':'^',
           r'ATLAS $|\eta|<2.0$':'v'}

exp_cols = {label : colors[i] for i, label in enumerate(exp_labels)}

fig, (top_ax,bottom_ax) = plt.subplots(2,1,figsize=(16,9), 
                            #gridspec_kw={"right":0.96,
                            #"top":0.98, "bottom":0.12,
                            #"left":0.12,"hspace":0.05,
                            #"wspace":0.1,
                            height_ratios=(3,1),sharex=True)
pT  = 0.5*(martini_hads['pTmax']+martini_hads['pTmin'])
dpT = martini_hads['pTmax']-martini_hads['pTmin']

spec = martini_hads['N']/(pT*dpT*2*pi*2)
err  = martini_hads['dN']/(pT*dpT*2*pi*2)
f_theo = interp1d(log(pT), log(spec), kind='linear')
df_theo = interp1d(log(pT), log(err), kind='linear')

theory = exp(f_theo(log(cms_hads['x'])))
dtheory = exp(df_theo(log(cms_hads['x'])))


top_ax.plot(cms_hads['x'], theory, color='black')
top_ax.fill_between(cms_hads['x'], theory-dtheory, theory+dtheory, color='black', alpha=0.2)

for experiment in exp_labels:
    edata = exp_data[experiment]
    col = exp_cols[experiment]
    marker = markers[experiment]
    plot_expr_data_on_axis(top_ax, edata, marker, color=col,face=col)

    theory = exp(f_theo(log(edata['x'])))
    ratio = edata['y']/theory

    err_ratio_stat_pos = ratio*edata['dy_stat+']/edata['y']
    err_ratio_stat_neg = ratio*edata['dy_stat-']/edata['y']
    
    err_ratio_syst_pos = ratio*edata['dy_syst+']/edata['y']
    err_ratio_syst_neg = ratio*edata['dy_syst-']/edata['y']

    dx_exp = (edata['xhigh']-edata['xlow'])
    # Create patch collection with specified colour/alpha
    errorboxes = [Rectangle((x-0.5*delx, y - yerrlow), 
                            width=delx, 
                            height=abs(yerrlow)+abs(yerrhigh))
                for x, delx, y, yerrlow, yerrhigh in
                zip(edata['x'], dx_exp, ratio, err_ratio_syst_pos, -1*err_ratio_syst_neg)]
    pc = PatchCollection(errorboxes, facecolor=col, alpha=0.1)
    bottom_ax.add_collection(pc)
    bottom_ax.errorbar(edata['x'], ratio, xerr=0.5*(edata['xhigh']-edata['xlow']),
                                          yerr=[abs(err_ratio_stat_pos),abs(err_ratio_stat_neg)], 
                                          color=col, ls='none', marker=marker,markersize=8)
    


handles2 = [Line2D([],[],label=l, color=c, marker=markers[l], markersize=10) for l, c in exp_cols.items()]
handles2.append(Line2D([],[],label='MARTINI', color='black'))
top_ax.legend(loc='upper right', handles=handles2)
top_ax.text(0.05,0.4,r'p-p \@ $\sqrt{s}=2.76$ TeV' +'\n'+ r'$|\eta^{h^{\pm}}|<1.0$', transform=top_ax.transAxes)
bottom_ax.axhline(1.,linestyle='dotted',color='black')

top_ax.set_xlim(left=1,right=114)
top_ax.set_ylim(top=3e-1, bottom=1e-14)
top_ax.yaxis.set_major_locator(mticker.LogLocator(numticks=999))
top_ax.yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
#bottom_ax.set_ylim(top=1.8,bottom=0.2)
top_ax.set_ylabel(r'$E \frac{\mathrm{d}\sigma^{h^{\pm}}}{d^3p}$ (GeV$^{-2}$)')
bottom_ax.axhline(y=1,color='black',linestyle='dashed')
bottom_ax.set_xlabel(r'$p^{h^{\pm}}_T$ (GeV)')
bottom_ax.set_ylabel('Data over' + '\nTheory')
top_ax.set_yscale('log')
bottom_ax.set_xscale('log')

plt.show()
