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
    "font.size": 30,
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
    "axes.spines.top"  : False,
    "legend.frameon"   : False,
    "axes.labelsize" : 30}
plt.rcParams.update(my_rcParams)
from COLORS import module_colors as modules
from COLORS import colors, channel_colors, channel_linestyles
xsec = 42 #mb
oversample = 2000
multiplicities = {'00-10':942.2, '10-20':591.55,'00-05':1053,'05-10':831.4}

## Read in all the data:
cents_phenix = {'0_5':'00-05','5_10':'05-10','10_15':'10-15','15_20':'15-20'}
loc = '../../../../KERNELS_NLO_NP_PART2/exp_data/sqrts_200GeV/Photons/'
phenix_fname = loc + 'PHENIX/AuAu_direct_photon_{c}.csv'
PHENIX = {}
for c,v in cents_phenix.items():
    tmp = pd.read_csv(phenix_fname.format(c=c),comment='#').rename(columns=my_dicts.colnames_PHENIX_gamma)
    tmp = tmp[tmp['xhigh'] < 20]
    PHENIX[v] = tmp
          
y = 0.25*(PHENIX['00-05']['y'] + PHENIX['05-10']['y'] + PHENIX['10-15']['y'] + PHENIX['15-20']['y'])
dy_stat_p = 0.25*np.sqrt(PHENIX['00-05']['dy_stat+']**2 + \
                        PHENIX['05-10']['dy_stat+']**2 + PHENIX['10-15']['dy_stat+']**2 + PHENIX['15-20']['dy_stat+']**2)
dy_stat_n = -0.25*np.sqrt(PHENIX['00-05']['dy_stat-']**2 + \
                        PHENIX['05-10']['dy_stat-']**2 + PHENIX['10-15']['dy_stat-']**2 + PHENIX['15-20']['dy_stat-']**2)
dy_syst_p = 0.25*np.sqrt(PHENIX['00-05']['dy_syst+']**2 + \
                        PHENIX['05-10']['dy_syst+']**2 + PHENIX['10-15']['dy_syst+']**2 + PHENIX['15-20']['dy_syst+']**2)
dy_syst_n = -0.25*np.sqrt(PHENIX['00-05']['dy_syst-']**2 + \
                        PHENIX['05-10']['dy_syst-']**2 + PHENIX['10-15']['dy_syst-']**2 + PHENIX['15-20']['dy_syst-']**2)

tmpdata = pd.DataFrame({'xlow': PHENIX['00-05']['xlow'].to_list(),\
                    'xhigh':PHENIX['00-05']['xhigh'].to_list(),\
                    'x':PHENIX['00-05']['x'].to_list(),\
                    'y':y.to_list(),\
                    'dy_stat+':dy_stat_p.to_list(), 'dy_stat-':dy_stat_n.to_list(),\
                    'dy_syst+':dy_syst_p.to_list(), 'dy_syst-':dy_syst_n.to_list()})
PHENIX['00-20'] = tmpdata

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
        Nbin = multiplicities[cent] if 'd' not in col else np.sqrt(multiplicities[cent])
        tmp[col] *= Nbin
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
channels = ['prompt','thermal','preEq']
models = ['cujet', 'martini']
jetscape_data_loc = '../../jetscape_data/sqrt_s_200/AuAu/'
fname = jetscape_data_loc + '{model}/{cent}/photon_spec_0.35.csv'
jet_medium = {m: {c: get_spectrum(fname.format(model=m,cent=c),\
                                                osf=oversample,cent=c) for c in centralities} for m in models}
for model in jet_medium:
    specs = jet_medium[model]
    tmpy = 0.5*(specs['00-10']['total'] + specs['10-20']['total'])
    dtmpy = 0.5*np.sqrt(specs['00-10']['total']**2 + specs['10-20']['total']**2 )
    jet_medium[model]['00-20'] = pd.DataFrame({'pT':specs['00-10']['pT'].to_list(), 'total':tmpy.to_list(), 'dtotal':dtmpy.to_list()})
## Read in JF's work:
JFCents = ['00-05','05-10','10-20']
JFLOC = '../../other_data/JF_MultiMessenger/'
tmpl  = 'AuAu200_{c}_{ch}.csv'
other_photons = {}
for cent in JFCents:
    other_photons[cent] = {}
    for ch in channels:
        if ch == 'prompt':
            continue
        tmp = np.loadtxt(JFLOC+tmpl.format(c=cent,ch=ch), unpack=True, delimiter=',')
        other_photons[cent][ch] = pd.DataFrame({'pT':tmp[0], 'N':tmp[1], 'dN':np.zeros_like(tmp[0])})
other_photons['00-20'] = {}
## construct the 0-20% for the "other_photons"
tmp1 = 0.25*(other_photons['00-05']['thermal']['N'] + other_photons['05-10']['thermal']['N'] + other_photons['10-20']['thermal']['N'])
other_photons['00-20']['thermal'] = pd.DataFrame({'pT':other_photons['00-05']['thermal']['pT'], 'N':tmp1.to_list(), 'dN':np.zeros_like(tmp1)})

tmp1 = 0.25*(other_photons['00-05']['preEq']['N'] + other_photons['05-10']['preEq']['N'] + other_photons['10-20']['preEq']['N'])
other_photons['00-20']['preEq'] = pd.DataFrame({'pT':other_photons['00-05']['preEq']['pT'], 'N':tmp1.to_list(), 'dN':np.zeros_like(tmp1)})

tmp = pd.read_csv('../../jetscape_data/prompt_photons/AuAu_200/gamma_spectra.csv', comment='#')
tmp['pT'] = 0.5*(tmp['pTmin'] + tmp['pTmax'])
tmp = tmp[tmp['pT'].between(1, 20)]
tmp_x   = tmp['pT']
tmp_dx  = tmp['pTmax'] - tmp['pTmin']
tmp_y   = tmp['prompt']/(2*np.pi*tmp_x*tmp_dx*xsec*2*0.35)
tmp_dy  = tmp['dprompt']/(2*np.pi*tmp_x*tmp_dx*xsec*2*0.35)

mult = 0.25*(1053 + 831.4 + 660.1 + 523)
yy  = tmp_y *mult
dyy = tmp_dy*mult
## find k-factor:
f = interp1d(tmp_x, np.log(yy), kind='linear', fill_value='extrapolate')
dat_x, dat_y = PHENIX['00-20']['x'].to_list(), PHENIX['00-20']['y'].to_list()
kfactor = dat_y[-1]/np.exp(f(dat_x[-1]))
#kfactor -= 0.02
print(dat_x[-1])
print(kfactor)
yy = yy * kfactor
dyy = dyy * np.sqrt(kfactor)
other_photons['00-20']['prompt'] = pd.DataFrame({'pT':tmp_x, 'N':yy.to_list(), 'dN':dyy.to_list()})

# construct total 
x = np.linspace(1, 22, 15)
for cent in ['00-20']:
    specs = other_photons[cent]
    f_prompt  = interp1d(specs['prompt']['pT']   , np.log(specs['prompt']['N']) , kind='linear', fill_value='extrapolate')
    f_thermal = interp1d(specs['thermal']['pT']  , np.log(specs['thermal']['N']), kind='linear', fill_value='extrapolate')
    f_preEq   = interp1d(specs['preEq']['pT']    , np.log(specs['preEq']['N'])  , kind='linear', fill_value='extrapolate')
    df_prompt = interp1d(specs['prompt']['pT']   , np.log(specs['prompt']['dN']), kind='linear', fill_value='extrapolate')
    ytot      = np.exp(f_preEq(x)) + np.exp(f_prompt(x)) + np.exp(f_thermal(x)) 
    dytot     = np.exp(df_prompt(x))
    other_photons[cent]['total'] = (x,ytot,dytot)

ratios = {}
xx = np.linspace(2, 20, 25)
for cent in ['00-20']:
    ratios[cent] = {}
    ## construct ratios
    f_prmpt = interp1d(other_photons[cent]['prompt']['pT'], np.log(other_photons[cent]['prompt']['N']), kind='linear', fill_value='extrapolate')
    f_thmal = interp1d(other_photons[cent]['thermal']['pT'], np.log(other_photons[cent]['thermal']['N']), kind='linear', fill_value='extrapolate')
    f_preEq = interp1d(other_photons[cent]['preEq']['pT'], np.log(other_photons[cent]['preEq']['N']), kind='linear', fill_value='extrapolate')
    y_prompt = np.exp(f_prmpt(xx))
    y_thmal = np.exp(f_thmal(xx))
    y_preEq = np.exp(f_preEq(xx))
    for model in jet_medium:
        f_jet = interp1d(jet_medium[model][cent]['pT'], np.log(jet_medium[model][cent]['total']), kind='linear', fill_value='extrapolate')
        y_jet = np.exp(f_jet(xx))
        total = y_preEq + y_prompt + y_thmal + y_jet
        r_prompt  = y_prompt/total
        r_thermal = y_thmal/total
        r_preEq   = y_preEq/total
        r_jet     = y_jet/total
        ratios[cent][model] = {'prompt': r_prompt, 'preEq':r_preEq, 'thermal':r_thermal, 'jetmedium': r_jet}


"""
    * Figures:
        fig1 : compare to phenix, cent by cent
        fig2 : compare to STAR 0-20 and constructed 0-20 for phenix
"""
fig2, axes2 = plt.subplots(2, 1,sharex=True, figsize=(16,9), sharey=False, height_ratios=(3,1))

models.append('none')
## FIGURE 1 Work
# plot PHENIX data:
pTmin, pTmax = 2, 30
exp_handles = []
theory_handles = [Line2D([],[],color=modules['MATTER+'+m.upper()], label='MATTER+'+m.upper()) for m in models if m!='none']
theory_handles.append(Line2D([],[],color=colors[2], label='No Jet-Medium'))

# ## Figure 2: 0-20%
# # plot data first:
ax = axes2[0]
util.plot_expr_data_on_axis(ax, PHENIX['00-20'], marker='v')
f_oth  = interp1d(other_photons['00-20']['total'][0], 
                  np.log(other_photons['00-20']['total'][1]), 
                  kind='quadratic', fill_value='extrapolate')
df_oth = lambda x:-200

merged_x = PHENIX['00-20']['x'].to_list()
ax = axes2[0]
locc = '/Users/rmyazdi/Documents/research/jetscape_project/v2/jetscape_data/photon_sp/'
zorders = {'martini':0.1, 'cujet':0.4, 'none':0.6}
for model in models:
    if model != 'none':
        color  = modules['MATTER+'+model.upper()]
        spec   = jet_medium[model]
        f_jm   = interp1d(spec['00-20']['pT'], np.log(spec['00-20']['total']), kind='linear', fill_value='extrapolate')
        df_jm  = interp1d(spec['00-10']['pT'], np.log(spec['00-20']['dtotal']), kind='linear', fill_value='extrapolate')
        total_curve = np.exp(f_jm(merged_x)) + np.exp(f_oth(merged_x))
        dtotal_curve = np.sqrt(np.exp(df_jm(merged_x))**2 + np.exp(df_oth(merged_x))**2)
        total = np.exp(f_jm(merged_x)) + np.exp(f_oth(merged_x))
        ax.plot(merged_x, total, color=color, zorder=zorders[model])
        ax.fill_between(merged_x, total+dtotal_curve, total-dtotal_curve, color=color, alpha=0.2, zorder=zorders[model])
        ## Now do ratio:
        # constructed PHENIX :
        y_phenix = np.array([np.exp(f_jm(x)) + np.exp(f_oth(x)) for x in PHENIX['00-20']['x']])
        ratio_PHENIX = y_phenix/ PHENIX['00-20']['y']
        d_pos = PHENIX['00-20']['dy_stat+']**2 + PHENIX['00-20']['dy_syst+']**2
        d_neg = PHENIX['00-20']['dy_stat-']**2 + PHENIX['00-20']['dy_syst-']**2
        dy_sq = np.exp(df_jm(PHENIX['00-20']['x']))**2 + np.exp(df_oth(PHENIX['00-20']['x']))**2
        err_pos = ratio_PHENIX*np.sqrt(dy_sq/y_phenix**2 + d_pos/PHENIX['00-20']['y']**2)
        err_neg = ratio_PHENIX*np.sqrt(dy_sq/y_phenix**2 + d_neg/PHENIX['00-20']['y']**2)
        delx = 0.5*(PHENIX['00-20']['xhigh'] - PHENIX['00-20']['xlow'])
        axes2[1].scatter(PHENIX['00-20']['x'], ratio_PHENIX, marker='s', color=color, zorder=zorders[model])
        axes2[1].errorbar(PHENIX['00-20']['x'], ratio_PHENIX, xerr=delx, yerr=[err_neg, err_pos], fmt='none',color=color, zorder=zorders[model])

    else:
        color  = colors[2]
        total_curve = np.exp(f_oth(merged_x))
        dtotal_curve = np.exp(df_oth(merged_x))
        ax.plot(merged_x, total_curve, color=color,zorder=zorders[model])
        ax.fill_between(merged_x, total_curve+dtotal_curve, total_curve-dtotal_curve, color=color, zorder=zorders[model], alpha=0.4)
        ## Now do ratio:
        y_phenix = np.array([np.exp(f_oth(x)) for x in PHENIX['00-20']['x']])
        ratio_PHENIX = y_phenix/ PHENIX['00-20']['y']
        d_pos = PHENIX['00-20']['dy_stat+']**2 + PHENIX['00-20']['dy_syst+']**2
        d_neg = PHENIX['00-20']['dy_stat-']**2 + PHENIX['00-20']['dy_syst-']**2
        dy_sq = np.array([0 for nx in PHENIX['00-20']['x']])
        err_pos = ratio_PHENIX*np.sqrt(dy_sq/y_phenix**2 + d_pos/PHENIX['00-20']['y']**2)
        err_neg = ratio_PHENIX*np.sqrt(dy_sq/y_phenix**2 + d_neg/PHENIX['00-20']['y']**2)
        delx = 0.5*(PHENIX['00-20']['xhigh'] - PHENIX['00-20']['xlow'])
        axes2[1].scatter(PHENIX['00-20']['x'], ratio_PHENIX, marker='s', color=color,zorder=zorders[model])
        axes2[1].errorbar(PHENIX['00-20']['x'], ratio_PHENIX, xerr=delx, yerr=[err_neg, err_pos], fmt='none',color=color,zorder=zorders[model])

axes2[0].set_yscale('log')
#handles = [Line2D([],[],color=rate_colours[rate_names[r]], label=rate_names[r]) for r in rate_names]
#axes2[0].add_artist(artist)
#expt_handles = []
theory_handles.append(Line2D([],[],color='black', marker='s', label='PHENIX (2012) $0$-$20\%$'))
artist = axes2[0].legend(loc='best', handles=theory_handles)
#expt_handles.append(Line2D([],[],color='black', marker='^', label='STAR (2017) $0$-$20\%$'))
#axes2[0].legend(loc='center right', handles=expt_handles)
axes2[0].text(0.05, 0.15, 'Au-Au, $\sqrt{s}=200$ AGeV, $|\eta|<0.35$', transform=axes2[0].transAxes)
axes2[0].set_ylabel(r'$E_{\gamma} \frac{dN^{\gamma}}{d^3p^{\gamma}}$ (GeV${}^{-2}$)')
axes2[1].set_ylabel('Theory over' + '\n' + 'Data')
axes2[0].set_yscale('log')
axes2[1].set_xlabel('$p_T$ (GeV)')
axes2[1].axhline(1,color='black', linestyle='dotted')

fig3, axes3 = plt.subplots(1, 2, figsize=(16,9), sharex=True, sharey=True)
## plot channel ratios
cent = '00-20'
    
axes3[0].text(0.05,0.8, cent+r'$\%$', transform=axes3[0].transAxes)
spec = ratios[cent]
for imodel, model in enumerate(['martini','cujet']):
    ax = axes3[imodel]
    for ch in spec[model]:
        if ch == 'jetmedium':
            ax.plot(xx, spec[model][ch], color=modules['MATTER+'+model.upper()])
        else:
            ax.plot(xx, spec[model][ch], color=channel_colors[ch], linestyle=channel_linestyles[ch])

for ax in axes3:
    ax.set_xlabel(r'$p_T$ (GeV)')

axes3[0].set_ylabel(r'Channel/Total')

nice_channels = {'prompt':'Prompt', 'preEq':'Pre-Eq.', 'thermal':'Thermal'}

handles = [Line2D([],[],color=modules[f'MATTER+MARTINI'],label=f'MATTER+MARTINI')]
axes3[0].legend(loc='center right', handles=handles)
handles = [Line2D([],[],color=modules[f'MATTER+CUJET'],label=f'MATTER+CUJET')]
artist = axes3[1].legend(loc='center right', handles=handles)
axes3[1].add_artist(artist)
jf_handles = []
for ch in channels:
    jf_handles.append(Line2D([],[],color=channel_colors[ch], linestyle=channel_linestyles[ch],label=nice_channels[ch]))

axes3[1].legend(loc='upper right', handles=jf_handles, fontsize=30, bbox_to_anchor=(0.9,0.9))
fig3.suptitle('Au-Au @ $\sqrt{s}=200$ AGeV, $|\eta|<0.35$')
plt.show()
