import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pandas import read_csv
from numpy import pi, log, arange,exp, loadtxt,sqrt, array
from scipy.interpolate import interp1d
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import dictionaries as ddicts
import jetDicts as jdicts
from util import plot_expr_data_on_axis
#mpl.use('pdf')
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

fig, (top_ax,bottom_ax) = plt.subplots(2,1,figsize=(16,9), 
                            gridspec_kw={"right":0.96,
                            "top":0.98, "bottom":0.12,
                            "left":0.12,"hspace":0.05,
                            "wspace":0.1,"height_ratios":(3,1)},sharex=True)
cols = { 2 : "#422040",
         3 : "#d64045",
         4 : "#0092cc",
         5 : "#eca400",
         6 : "#7cae7a"}
handles = []
factors = {2: 1, 3:5, 4:16}
markers = ['o','P','*']

hadron_jets = read_csv("../martini_results/pp/jet_spectra.csv", comment='#')
parton_jets = read_csv("../martini_results/pp/partonic_jet_spectra.csv", comment='#')
for (i, marker) in zip([2,3,4],['o','o','o']):
    cms_jets = read_csv(f"../../../exp_data/jets/CMS/Table4_0p{i}.csv",comment='#').rename(columns=jdicts.colnames_CMS_spec_pp)
    pT = log(cms_jets['x'])

    fac = factors[i]
    plot_expr_data_on_axis(top_ax,cms_jets,marker,color=cols[i],face=cols[i], factor=fac) 
    handles.append(Line2D([],[],label=rf'R=0.{i} $\times$ {fac:d}',marker=marker,color=cols[i],markersize=10))

    x  = 0.5*(hadron_jets['pTmin']+hadron_jets['pTmax'])
    dx = hadron_jets['pTmax']-hadron_jets['pTmin']

    y  = 1e6*hadron_jets[f"N0p{i}"] /(dx*4)
    dy = 1e6*hadron_jets[f"dN0p{i}"]/(dx*4)

    f_had = interp1d(log(x), log(y))
    df_had = interp1d(log(x), log(dy))

    hjet = exp(f_had(log(cms_jets['x'])))
    dhjet = exp(df_had(log(cms_jets['x'])))

    top_ax.plot(cms_jets['x'],fac*hjet,color=cols[i])
    top_ax.fill_between(cms_jets['x'],fac*(hjet+dhjet),fac*(hjet-dhjet),color=cols[i],alpha=0.2)

    x1  = 0.5*(parton_jets['pTmin']+parton_jets['pTmax'])
    dx1 = parton_jets['pTmax']-parton_jets['pTmin']

    y1  = 1e6*parton_jets[f"N0p{i}"] /(dx*4)
    dy1 = 1e6*parton_jets[f"dN0p{i}"]/(dx*4)

    f_part = interp1d(log(x1), log(y1))
    df_part = interp1d(log(x1), log(dy1))

    pjet = exp(f_part(log(cms_jets['x'])))
    dpjet = exp(df_part(log(cms_jets['x'])))

    top_ax.plot(cms_jets['x'],fac*pjet,color=cols[i], linestyle='dotted')
    top_ax.fill_between(cms_jets['x'],fac*(pjet+dpjet),fac*(pjet-dpjet),color=cols[i],alpha=0.2)

    ratio = pjet/hjet
    bottom_ax.plot(cms_jets['x'], ratio, color=cols[i])

#handles += [Line2D([],[],label=f"R = 0.{i}",color=cols[i]) for i in [2,3,4]]
artist = top_ax.legend(loc='lower left', handles=handles)
top_ax.add_artist(artist)
#handles2 = [Line2D([],[],label='CMS: Phys.Rev.C 96 (2017) 015202, 2017.',color='black',marker='o',markersize=10)]
handles2 = [Line2D([],[],label='CMS(2017)',color='black',marker='o',markersize=10)]

top_ax.legend(loc='upper right', handles=handles2)
top_ax.text(0.6,0.7,r'p-p \@ $\sqrt{s}=2.76$ TeV' +'\n'+ r'$|\eta_{\mathrm{jet}}|<2.0$', transform=top_ax.transAxes)
top_ax.set_yscale('log')
top_ax.set_xlim(left=62,right=303)
top_ax.set_ylim(bottom=7e-5, top=2e2)
bottom_ax.set_ylim(bottom=0.7,top=1.3)
top_ax.set_ylabel(r'$\frac{\mathrm{d}\sigma^{\mathrm{jet}}}{\mathrm{d}p_T\mathrm{d}\eta}$ (nb GeV$^{-1}$)')
bottom_ax.axhline(y=1,color='black',linestyle='dashed')
bottom_ax.set_xlabel(r'$p^{\mathrm{jet}}_T$ (GeV)')
bottom_ax.set_ylabel('Partonic/Hadronic', fontsize=20)
plt.show()
