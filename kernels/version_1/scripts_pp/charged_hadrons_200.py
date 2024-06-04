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
inelXsec = 42

loc = '../../../exp_data/sqrts_200GeV/Charged/pp_NSD.csv'
star_data    = read_csv(loc,comment='#').rename(columns=ddicts.colnames_STAR_pp_ch_spec)
martini_hads = read_csv("../martini_results/final_AuAu_200/pp/hadron_spectra.csv",comment='#')
fig, (top_ax,bottom_ax) = plt.subplots(2,1,figsize=(16,9), 
                            gridspec_kw={"right":0.96,
                            "top":0.98, "bottom":0.12,
                            "left":0.12,"hspace":0.05,
                            "wspace":0.1,"height_ratios":(3,1)},sharex=True)

plot_expr_data_on_axis(top_ax, star_data, 'o')
pT  = 0.5*(martini_hads['pTmax']+martini_hads['pTmin'])
dpT = martini_hads['pTmax']-martini_hads['pTmin']

spec = martini_hads['N']/(2*pT*dpT*2*pi*2*0.5*inelXsec)
err  = martini_hads['dN']/(2*pT*dpT*2*pi*2*0.5*inelXsec)
f_theo = interp1d(log(pT), log(spec), kind='linear')
df_theo = interp1d(log(pT), log(err), kind='linear')

theory = exp(f_theo(log(star_data['x'])))
dtheory = exp(df_theo(log(star_data['x'])))
top_ax.plot(star_data['x'], theory, color=colour)
top_ax.fill_between(star_data['x'], theory-dtheory, theory+dtheory, color=colour, alpha=0.2)

theory = exp(f_theo(log(star_data['x'])))
ratio = star_data['y']/theory

err_ratio_pos = ratio*star_data['dy_+']/star_data['y']
err_ratio_neg = ratio*star_data['dy_-']/star_data['y']

dx_exp = (star_data['xhigh']-star_data['xlow'])
# Create patch collection with specified colour/alpha
errorboxes = [Rectangle((x-0.5*delx, y - yerrlow), 
                        width=delx, 
                        height=abs(yerrlow)+abs(yerrhigh))
                for x, delx, y, yerrlow, yerrhigh in
                zip(star_data['x'], dx_exp, ratio, err_ratio_pos, -1*err_ratio_neg)]
pc = PatchCollection(errorboxes, facecolor=colour, alpha=0.3)
bottom_ax.add_collection(pc)
bottom_ax.scatter(star_data['x'], ratio, color=colour, marker='o',s=15)

handles2 = [Line2D([],[],label='STAR(2003)',color='black',marker='o',markersize=10)]
top_ax.legend(loc='upper right', handles=handles2)
top_ax.text(0.8,0.5,r'p-p \@ $\sqrt{s}=200$ GeV' +'\n'+ r'$|\eta^{h^{\pm}}|<0.5$', transform=top_ax.transAxes)
bottom_ax.axhline(1.,linestyle='dotted',color='black', zorder=0)

top_ax.set_xlim(left=0.2,right=11)
top_ax.set_ylim(top=2, bottom=5e-9)
top_ax.yaxis.set_major_locator(mticker.LogLocator(numticks=999))
top_ax.yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
bottom_ax.set_ylim(top=10.0,bottom=-0.1)
top_ax.set_ylabel(r'$E \frac{\mathrm{d}N^{h^{\pm}}}{d^3p}$ (GeV$^{-2}$)')
bottom_ax.axhline(y=1,color='black',linestyle='dashed')
bottom_ax.set_xlabel(r'$p^{h^{\pm}}_T$ (GeV)')
bottom_ax.set_ylabel('Data/Theory')
top_ax.set_yscale('log')
#bottom_ax.set_xscale('log')

plt.show()
