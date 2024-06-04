## import scipy and matplotlib modules:
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle
from scipy.interpolate import interp1d
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
from COLORS import module_colors as modules
from COLORS import colors, channel_colors, channel_linestyles
xsec = 42 #mb
oversample = 2000
multiplicities = {'00-10':942.2, '10-20':591.55,'00-05':1053,'05-10':831.4}

## Read in MARTINI Calculations:
def get_spectrum(fname, osf=1, cent='00-05'):
    #pTmin, pTmax = (0,19) if 'final' in fname else (0, 20)
    tmp = pd.read_csv(fname, comment='#')
    xmin, xmax = '',''
    if 'pTmin' in tmp:
        xmin, xmax  = 'pTmin', 'pTmax'
        tmp['pT'] = 0.5 * (tmp['pTmin'] + tmp['pTmax'])
        tmp['dpT'] = tmp['pTmax'] - tmp['pTmin']
    else:
        xmin, xmax  = 'ptmin', 'ptmax'
        tmp['pT'] = 0.5 * (tmp['ptmin'] + tmp['ptmax'])
        tmp['dpT'] = tmp['ptmax'] - tmp['ptmin']
    #tmp = tmp[tmp['pT'].between(pTmin, pTmax)]
    tmp = tmp[tmp[xmax] < 21]
    for col in ['conv','dconv','brem','dbrem']:
        tmp[col] /= (xsec*osf)
        tmp[col] *= multiplicities[cent]
    tmp['total']  = tmp['conv'] + tmp['brem']
    tmp['dtotal'] = np.sqrt(tmp['dconv']**2 + tmp['dbrem']**2)
    tmp['total']  /= (2*np.pi*2*0.35*tmp['pT']*tmp['dpT'])
    tmp['dtotal'] /= (2*np.pi*2*0.35*tmp['pT']*tmp['dpT']) 
    tmp['conv']   /= (2*np.pi*2*0.35*tmp['pT']*tmp['dpT'])
    tmp['dconv']  /= (2*np.pi*2*0.35*tmp['pT']*tmp['dpT']) 
    tmp['brem']   /= (2*np.pi*2*0.35*tmp['pT']*tmp['dpT'])
    tmp['dbrem']  /= (2*np.pi*2*0.35*tmp['pT']*tmp['dpT'])  
    return tmp

centralities = ['00-10', '10-20']
channels = ['thermal','preEq']
models = ['cujet', 'martini']
jetscape_data_loc = '../../jetscape_data/sqrt_s_200/AuAu/'
fname = jetscape_data_loc + '{model}/{cent}/photon_spec_0.35.csv'
jet_medium = {m: {c: get_spectrum(fname.format(model=m,cent=c),\
                                                osf=oversample,cent=c) for c in centralities} for m in models}

# ## Read in non-jet medium sources, excluding prompts:
cents = ['00-05','05-10','10-20']
loc = '../../other_data/JF_MultiMessenger/'
tmpl  = 'AuAu200_{c}_{ch}.csv'
non_jet_med = {}
for cent in cents:
    non_jet_med[cent] = {}
    for ch in channels:
        tmp = np.loadtxt(loc+tmpl.format(c=cent,ch=ch), unpack=True, delimiter=',')
        non_jet_med[cent][ch] = pd.DataFrame({'pT':tmp[0], 'N':tmp[1], 'dN':np.zeros_like(tmp[0])})

## constrcut 0-20 % centrality from JF's work (non-jet medium photons)
non_jet_med['00-20'] = {}
for ch in ['thermal', 'preEq']:
    tmp_x = non_jet_med['00-05'][ch]['pT']
    tmp_y  = (non_jet_med['00-05'][ch]['N'] + non_jet_med['05-10'][ch]['N']+non_jet_med['10-20'][ch]['N'])/3.
    tmp_dy = np.sqrt(non_jet_med['00-05'][ch]['dN'] **2 + non_jet_med['05-10'][ch]['dN'] **2 +non_jet_med['10-20'][ch]['dN'] **2)/3.
    non_jet_med['00-20'][ch] = pd.DataFrame({'pT':tmp_x.to_list(), 'N':tmp_y.to_list(), 'dN':tmp_dy.to_list()}) 

## construct 0-20% jet medium spectrum:
for model in jet_medium:
    tmp_y  = 0.5*(jet_medium[model]['00-10']['total'] + jet_medium[model]['10-20']['total'] )
    tmp_dy = 0.5*np.sqrt(jet_medium[model]['00-10']['dtotal']**2 + jet_medium[model]['10-20']['dtotal']**2)
    jet_medium[model]['00-20'] = pd.DataFrame({'pT':jet_medium[model]['00-10']['pT'].to_list() ,'total':tmp_y.to_list(), 'dtotal':tmp_dy.to_list()})


fname = '/Users/rmyazdi/Documents/research/jetscape_project/expt/AuAu_200/NonPrompt/non_prompt_00-20.csv'
colnames = {'$p_{T}$':'x',
            '$p_{T}$ LOW':'xlow',
            '$p_{T}$ HIGH':'xhigh',
            'inv.yield':'y',
            'stat. +':'dy_stat+',
            'stat. -':'dy_stat-',
            'sys. +':'dy_syst+',
            'sys. -':'dy_syst-'}
data_phenix = pd.read_csv(fname, comment='#').rename(columns=colnames)
## finally, construct the total 0-20% non-prompt spectrum
x = data_phenix['x']
cent = '00-20'
totals = {}
f_therm = interp1d(non_jet_med[cent]['thermal']['pT'], np.log(non_jet_med[cent]['thermal']['N']), kind='linear', fill_value='extrapolate')
y_therm = np.exp(f_therm(x))
f_preEq = interp1d(non_jet_med[cent]['preEq']['pT'], np.log(non_jet_med[cent]['preEq']['N']), kind='linear', fill_value='extrapolate')
y_preEq = np.exp(f_preEq(x))
totals['No Jet-Medium'] = y_therm + y_preEq
for model in jet_medium:
    f = interp1d(jet_medium[model][cent]['pT'], np.log(jet_medium[model][cent]['total']), kind='linear', fill_value='extrapolate')
    y = np.exp(f(x))
    totals[model] = y_therm + y_preEq + y

fig, (ax1, ax2) = plt.subplots(2, 1, height_ratios=(3,1), figsize=(16,9), sharex=True)
util.plot_expr_data_on_axis(ax1, data_phenix, marker='*')
ax1.set_yscale('log')
delx = 0.5*(data_phenix['xhigh'] - data_phenix['xlow'])
syst_neg = np.sqrt(data_phenix['dy_syst-']**2/data_phenix['y']**2)
syst_pos = np.sqrt(data_phenix['dy_syst+']**2/data_phenix['y']**2)
stat_neg = np.sqrt(data_phenix['dy_stat-']**2/data_phenix['y']**2)
stat_pos = np.sqrt(data_phenix['dy_stat+']**2/data_phenix['y']**2)
zorders = {'martini':3, 'cujet':2, 'No Jet-Medium':0}
for model in ['martini', 'cujet', 'No Jet-Medium']:
    if model != 'No Jet-Medium':
        color = modules['MATTER+'+model.upper()]
    else:
        color = colors[2]
    ax1.plot(x, totals[model], color=color)

    ratio  = data_phenix['y']/totals[model]
    errorboxes = [Rectangle((xval-dx, y - yerrlow), width=2*dx, height=abs(yerrlow)+abs(yerrhigh), zorder=zorders[model])
                    for xval, dx, y, yerrlow, yerrhigh in
                    zip(x, delx, ratio, syst_pos, -1*syst_neg)]

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=color, alpha=0.3)
    ax2.add_collection(pc)
    ax2.scatter(x, ratio, color=color, marker='*', s=60, zorder=zorders[model])
    ax2.errorbar(x, ratio, xerr=delx, yerr=[stat_neg, stat_pos], fmt='none', color=color, zorder=zorders[model])

handles = []
nice_channels = {'prompt':'Prompt', 'preEq':'Pre-Eq.', 'thermal':'Thermal'}
for eloss in ['martini','cujet', 'No Jet-Medium']:
    if eloss != 'No Jet-Medium':
        handles.append(Line2D([],[],color=modules[f'MATTER+{eloss.upper()}'],label=f'MATTER+{eloss.upper()}'))
    else:
        handles.append(Line2D([],[],color=colors[2],label=eloss))

handles.append(Line2D([],[],color='black',marker='*', label='PHENIX (2022)', markersize=10))
ax1.legend(loc='upper right', handles=handles, fontsize=20, title='Non-prompt direct $\gamma$')
ax1.text(0.05,0.25, r'Au-Au @ $\sqrt{s}=200$ AGeV', transform=ax1.transAxes, fontsize=20)
ax1.text(0.05,0.1, r'$|\eta|<0.35$', transform=ax1.transAxes, fontsize=20)
ax1.set_ylabel(r'$\frac{1}{2\pi p_T}\frac{dN}{dp_T dy}$ (GeV${}^{-2}$)')
ax2.set_ylabel(r'Data over' + '\n'+ 'Theory')
ax2.set_xlabel(r'$p_T$ (GeV)')
plt.show()