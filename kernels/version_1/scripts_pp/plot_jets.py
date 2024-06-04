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

fig, (top_ax,bottom_ax) = plt.subplots(2,1,figsize=(16,9), 
                            #gridspec_kw={"right":0.96,
                            #"top":0.98, "bottom":0.12,
                            #"left":0.12,"hspace":0.05,
                            #"wspace":0.1,
                            height_ratios=(3,1),sharex=True)
#cols = { 2: "#CC3311",
#        3: "#33BBEE",
#        4: "#EE7733"}

#d05525, #ad8c34, #384d56, #955f90, #00abb8
#cols = {2:"#d05525", 3:"#ad8c34", 4:"#384d56", 5:"#955f90", 6:"#00abb8"}
cols = { 2 : "#422040",
         3 : "#d64045",
         4 : "#0092cc",
         5 : "#eca400",
         6 : "#7cae7a"}
cols = {ix : colors[ix] for ix in range(2,7)}
handles = []
factors = {2: 1, 3:5, 4:16}
markers = ['o','P','*']

#pp_jets = read_csv("./data/jet_spectra.csv", comment='#')
pp_jets = read_csv("../martini_results/pp/jet_spectra.csv", comment='#')
for (i, marker) in zip([2,3,4],['o','o','o']):
    cms_jets = read_csv(f"../../../exp_data/jets/CMS/Table4_0p{i}.csv",comment='#').rename(columns=jdicts.colnames_CMS_spec_pp)
    pT = log(cms_jets['x'])

    fac = factors[i]
    plot_expr_data_on_axis(top_ax,cms_jets,marker,color=cols[i],face=cols[i], factor=fac) 
    handles.append(Line2D([],[],label=rf'R=0.{i} $\times$ {fac:d}',marker=marker,color=cols[i],markersize=10))

    x  = 0.5*(pp_jets['pTmin']+pp_jets['pTmax'])
    dx = pp_jets['pTmax']-pp_jets['pTmin']

    y  = 1e6*pp_jets[f"N0p{i}"] /(dx*4)
    dy = 1e6*pp_jets[f"dN0p{i}"]/(dx*4)

    f_theo = interp1d(log(x), log(y))
    df_theo = interp1d(log(x), log(dy))

    clean_y  = exp(f_theo(pT))
    clean_dy = exp(df_theo(pT))
    top_ax.plot(cms_jets['x'],fac*clean_y,color=cols[i])
    top_ax.fill_between(cms_jets['x'],fac*(clean_y+clean_dy),fac*(clean_y-clean_dy),color=cols[i],alpha=0.2)


    ratio = cms_jets['y']/exp(f_theo(pT))
    err_ratio_stat_pos = ratio*cms_jets['dy_stat+']/cms_jets['y']
    err_ratio_stat_neg = ratio*cms_jets['dy_stat-']/cms_jets['y']

    err_ratio_syst_pos = ratio*cms_jets['dy_syst+']/cms_jets['y']
    err_ratio_syst_neg = ratio*cms_jets['dy_syst-']/cms_jets['y']

    dx_exp = (cms_jets['xhigh']-cms_jets['xlow'])
    # Create patch collection with specified colour/alpha
    errorboxes = [Rectangle((x-0.5*delx, y - yerrlow), 
                            width=delx, 
                            height=abs(yerrlow)+abs(yerrhigh))
                    for x, delx, y, yerrlow, yerrhigh in
                    zip(cms_jets['x'], dx_exp, ratio, err_ratio_syst_pos, -1*err_ratio_syst_neg)]
    pc = PatchCollection(errorboxes, facecolor=cols[i], alpha=0.1)
    bottom_ax.add_collection(pc)
    bottom_ax.errorbar(cms_jets['x'], ratio, xerr=0.5*(cms_jets['xhigh']-cms_jets['xlow']), \
        yerr=[err_ratio_stat_pos,-1*err_ratio_stat_neg], color=cols[i], ls='none',markersize=8,marker='o')
    #bottom_ax.plot(cms_jets['x'], ratio, color=cols[i])
    #bottom_ax.fill_between(cms_jets['x'], ratio+err_pos,ratio-err_neg, alpha=0.1, color=cols[i])
handles.append(Line2D([],[],color='black', label='MARTINI'))

top_ax.legend(loc='upper right', handles=handles, title=r'CMS $|\eta|<2.0$')
top_ax.text(0.01,0.1,r'p-p \@ $\sqrt{s}=2.76$ TeV' +'\n'+ r'$|\eta_{\mathrm{jet}}|<2.0$', transform=top_ax.transAxes)
top_ax.set_yscale('log')
top_ax.set_xlim(left=62,right=303)
top_ax.set_ylim(bottom=7e-5, top=2e2)
bottom_ax.set_ylim(bottom=0.7,top=1.3)
top_ax.set_ylabel(r'$\frac{\mathrm{d}\sigma^{\mathrm{jet}}}{\mathrm{d}p_T\mathrm{d}\eta}$ (nb GeV$^{-1}$)')
bottom_ax.axhline(y=1,color='black',linestyle='dashed')
bottom_ax.set_xlabel(r'$p^{\mathrm{jet}}_T$ (GeV)')
bottom_ax.set_ylabel('Data/Theory')
plt.show()
