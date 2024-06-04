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
    "xtick.minor.size" : 4,
    "ytick.minor.size" : 4,
    "axes.spines.right": False,
    "axes.spines.top" : False,
    "legend.frameon":False
}
plt.rcParams.update(my_rcParams)

def get_spectrum(fname, oversampling_factor=1):
    pTmin, pTmax = (0,19) if 'final' in fname else (0, 20)
    tmp = pd.read_csv(fname, comment='#')
    print(fname,'\n')
    if 'pTmin' in tmp:
        tmp['pT'] = 0.5 * (tmp['pTmin'] + tmp['pTmax'])
        tmp['dpT'] = tmp['pTmax'] - tmp['pTmin']
    else:
        tmp['pT'] = 0.5 * (tmp['ptmin'] + tmp['ptmax'])
        tmp['dpT'] = tmp['ptmax'] - tmp['ptmin']
    tmp = tmp[tmp['pT'].between(pTmin, pTmax)]
    for col in ['conv','dconv','brem','dbrem']:
        tmp[col] /= oversampling_factor
    tmp['total']  = tmp['conv'] + tmp['brem']
    tmp['dtotal'] = np.sqrt(tmp['dconv']**2 + tmp['dbrem']**2)
    tmp['total'] /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT'])
    tmp['dtotal'] /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT']) 
    tmp['conv'] /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT'])
    tmp['dconv'] /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT']) 
    tmp['brem'] /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT'])
    tmp['dbrem'] /= (2*np.pi*2*0.8*tmp['pT']*tmp['dpT'])  
    return tmp

## READ IN THE CALCULATION:

AA_loc = '../martini_results/final_PbPb_2p76/rset_1/cent_0_5/gamma_spectra.csv'
### read in the calculated pp charged hadron spec
# Stand alone martini
standalone_AA = get_spectrum(AA_loc, 1000) 


## Read in the jetscape work
prefix = '../../../../jetscape_project/v2/jetscape_data'
fname_AA = 'PbPb2760_00-05_photon_spec_0.80.csv'
spec_dir_1 = '/sqrt_s_2760/martini/PbPb_2760/'
spec_dir_2 = '/sqrt_s_2760/martini_new_kappas/'
spec_dir_3 = '/sqrt_s_2760/martini_new_kap_matter_vac/'
jetscape_AA_old_kset   = get_spectrum(prefix+spec_dir_1+fname_AA, 1000000)
jetscape_AA_new_kset   = get_spectrum(prefix+spec_dir_2+fname_AA, 5000)
jetscape_AA_vac_matter = get_spectrum(prefix+spec_dir_3+fname_AA, 2000)

fig, ax = plt.subplots(1,1, figsize=(16,9),gridspec_kw={'left':0.079, 'bottom':0.09,'right':0.995,'top':0.995, 'hspace':0.038})

for run_type, color in zip([standalone_AA,jetscape_AA_old_kset, jetscape_AA_new_kset, jetscape_AA_vac_matter],
                           ['red','blue','green','orange']):
    y, dy = run_type['total'], run_type['dtotal']
    ax.plot(run_type['pT'], y, color=color)
    ax.fill_between(run_type['pT'], y+dy, y-dy, color=color, alpha=0.2)

ax.set_ylabel(r'$\frac{1}{2\pi p_T}\frac{d\sigma^{\gamma}}{dp_T d\eta}$ mb.GeV$^{-2}$')
ax.set_xlabel(r'$p^{\gamma}_T$ (GeV)')
ax.set_yscale('log')
from COLORS import rate_set_colors as colours
#colours = {'MARTINI(new $\kappa$-set)+PYTHIA':'red',
#           'MARTINI(old $\kappa$-set)+(MATTER:Med)':'blue',
#           'MARTINI(new $\kappa$-set)+(MATTER:Med)':'green',
#           'MARTINI(new $\kappa$-set)+(MATTER:Vac)':'orange'}
theor_hands = [Line2D([],[],label=l,color=c) for (l,c) in colours.items()]
ax.legend(loc='upper right', handles=theor_hands, ncol=1,fontsize=15)
ax.text(0.05,0.30, r'Pb-Pb @ $\sqrt{s}=2.76$ ATeV', transform=ax.transAxes)
ax.text(0.05,0.20, r'0$-$5$\%$, $|\eta|<0.8$',transform=ax.transAxes)

fig, axes = plt.subplots(2,1, height_ratios=[3,1], sharex=True, figsize=(16,9), gridspec_kw={'left':0.079, 'bottom':0.09,'right':0.995,'top':0.995, 'hspace':0.038})

for run_type, color in zip([standalone_AA,jetscape_AA_old_kset, jetscape_AA_new_kset, jetscape_AA_vac_matter],
                           ['red','blue','green','orange']):
    for channel, lstyle in zip(['brem','conv'],['solid','dashed']):
        ax = axes[0]
        y, dy = run_type[channel], run_type[f'd{channel}']
        ax.plot(run_type['pT'], y, color=color, linestyle=lstyle)
        ax.fill_between(run_type['pT'], y+dy, y-dy, color=color, alpha=0.2)
    ax = axes[1]
    y1, dy1 = run_type['brem'], run_type['dbrem'] 
    y2, dy2 = run_type['conv'], run_type['dconv'] 
    ratio = y2/y1
    dratio = ratio*np.sqrt(dy1*dy1/(y1*y1) + dy2*dy2/(y2*y2))
    ax.plot(run_type['pT'], ratio, color=color, linestyle='solid')
    ax.fill_between(run_type['pT'], ratio+dratio, ratio-dratio, color=color, alpha=0.2)

axes[0].set_ylabel(r'$\frac{1}{2\pi p_T}\frac{d\sigma^{\gamma}}{dp_T d\eta}$ mb.GeV$^{-2}$')
axes[1].set_xlabel(r'$p^{\gamma}_T$ (GeV)')
axes[1].set_ylabel('Conv. to Brem. Ratio', fontsize=18)
axes[0].set_yscale('log')
colours = {'MARTINI(new $\kappa$-set)+PYTHIA':'red',
           'MARTINI(old $\kappa$-set)+(MATTER:Med)':'blue',
           'MARTINI(new $\kappa$-set)+(MATTER:Med)':'green',
           'MARTINI(new $\kappa$-set)+(MATTER:Vac)':'orange'}
theor_hands = [Line2D([],[],label=l,color=c) for (l,c) in colours.items()]
axes[0].legend(loc='upper right', handles=theor_hands, ncol=1,fontsize=15)
axes[0].text(0.05,0.30, r'Pb-Pb @ $\sqrt{s}=2.76$ ATeV', transform=axes[0].transAxes)
axes[0].text(0.05,0.20, r'0$-$5$\%$, $|\eta|<0.8$',transform=axes[0].transAxes)
plt.show()
