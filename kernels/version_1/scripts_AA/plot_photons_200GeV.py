## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle
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
import COLORS 
rate_colours = COLORS.rate_set_colors
rate_names = {1:'LO',2:'NLO',3:'NP'}
xsec = 42 #mb
oversample_MARTINI = 1000
multiplicities = {'00-05':1053, '05-10':831.4, '10-20':591.55}
centralities = {'0_5':'00-05', '5_10':'05-10', '10_20':'10-20'}
channels = ['prompt','thermal','preEq']

## Read in MARTINI Calculations:
def get_spectrum(fname, osf=1, cent='00-05'):
    #pTmin, pTmax = (0,19) if 'final' in fname else (0, 20)
    tmp = pd.read_csv(fname, comment='#')
    if 'pTmin' in tmp:
        tmp['pT'] = 0.5 * (tmp['pTmin'] + tmp['pTmax'])
        tmp['dpT'] = tmp['pTmax'] - tmp['pTmin']
    else:
        tmp['pT'] = 0.5 * (tmp['ptmin'] + tmp['ptmax'])
        tmp['dpT'] = tmp['ptmax'] - tmp['ptmin']
    #tmp = tmp[tmp['pT'].between(pTmin, pTmax)]
    tmp = tmp[tmp['pTmax'] < 21]
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

fname = '../martini_results/final_AuAu_200/AuAu/rset_{r}/cent_{c}/gamma_spectra.csv'
jet_medium = {r: {centralities[c]: get_spectrum(fname.format(r=r,c=c),\
                                                osf=oversample_MARTINI,cent=centralities[c]) for c in centralities} for r in rate_names}

## Read in JF's work:
JFLOC='../../../../jetscape_project/v2/other_data/JF_MultiMessenger/'
tmpl = 'AuAu200_{c}_{ch}.csv'
other_photons = {v: {ch:np.loadtxt(JFLOC+tmpl.format(c=v,ch=ch), unpack=True, delimiter=',') for ch in channels if ch != 'prompt'} for v in centralities.values()}
tmp = pd.read_csv('~/Documents/research/jetscape_project/v2/jetscape_data/prompt_photons/AuAu_200/gamma_spectra.csv', comment='#')
tmp_x   = 0.5*(tmp['pTmin'] + tmp['pTmax'])
tmp_dx  = tmp['pTmax'] - tmp['pTmin']
tmp_y   = tmp['prompt']/(2*np.pi*tmp_x*tmp_dx*xsec*2*0.35)
tmp_dy  = tmp['dprompt']/(2*np.pi*tmp_x*tmp_dx*xsec*2*0.35)
for cent in other_photons:
    yy  = tmp_y*multiplicities[cent]
    dyy = tmp_dy*multiplicities[cent]
    other_photons[cent]['prompt'] = pd.DataFrame({'pT':tmp_x.to_list(), 'N':yy.to_list(), 'dN':dyy.to_list()})
## construct total 
x = np.linspace(1, 22, 15)
for cent in other_photons:
    specs = other_photons[cent]
    f_prompt  = interpolate(specs['prompt']['pT'], np.log(specs['prompt']['N']), kind='linear', fill_value='extrapolate')
    f_thermal = interpolate(specs['thermal'][0], np.log(specs['thermal'][1])   , kind='linear', fill_value='extrapolate')
    f_preEq   = interpolate(specs['preEq'][0], np.log(specs['preEq'][1])       , kind='linear', fill_value='extrapolate')
    ytot = np.exp(f_preEq(x)) + np.exp(f_prompt(x)) + np.exp(f_thermal(x)) 
    other_photons[cent]['total'] = (x,ytot)

## Read in all the data:
cents_phenix = {'0_5':'00-05','5_10':'05-10','10_15':'10-15','15_20':'15-20'}
loc = '../../../exp_data/sqrts_200GeV/Photons/'
phenix_fname = loc + 'PHENIX/AuAu_direct_photon_{c}.csv'
PHENIX = {v:pd.read_csv(phenix_fname.format(c=c),comment='#').rename(columns=my_dicts.colnames_PHENIX_gamma) \
          for c,v in cents_phenix.items()}
star_name = loc + 'STAR/HEPData-ins1474129-v1-Table_18.csv'
STAR = pd.read_csv(star_name, comment='#').rename(columns=my_dicts.colnames_STAR_gamma)

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

"""
    * Figures:
        fig1 : compare to phenix, cent by cent
        fig2 : compare to STAR 0-20 and constructed 0-20 for phenix
        fig3 : plot ratio of jet medium photons to the total photons for the 0-20% centrality
"""
fig1, axes1 = plt.subplots(2, 3,height_ratios=(2,1),
                          sharex='col', figsize=(16,9), sharey='row')
fig2, axes2 = plt.subplots(2, 1,sharex=True, figsize=(16,9), sharey=False, height_ratios=(2,1))

fig3, ax3 = plt.subplots(1, 1, figsize=(16,9))

## FIGURE 1 Work
# plot PHENIX data:
pTmin, pTmax = 2, 30
for icent, cent in enumerate(cents_phenix.values()):
    marker = 'x' if cent == '00-05' else 's'
    if cent =='15-20' or cent == '10-15':
        icent = 2
        marker = '<' if cent == '10-15' else '>'
    ax = axes1[0][icent]
    util.plot_expr_data_on_axis(ax, PHENIX[cent], marker=marker)


for rateset in rate_names:
    color = rate_colours[rate_names[rateset]]
    for icent, cent in enumerate(centralities.values()):
        ax = axes1[0][icent]
        ax.text(0.05, 0.1, f'{cent}'+'$\%$', transform=ax.transAxes)
        spec = jet_medium[rateset][cent]
        spec = spec[spec['pT'].between(pTmin,pTmax)]
        othr = other_photons[cent]['total']
        c = cent if cent != '10-20' else '15-20'
        dat = PHENIX[c]
        f_oth = interpolate(x=othr[0], y=np.log(othr[1]),kind='linear', fill_value='extrapolate')
        f_jetmed = interpolate(x=spec['pT'],y=np.log(spec['total']), kind='linear', fill_value='extrapolate')
        df = interpolate(x=spec['pT'], y=np.log(spec['dtotal']), kind='linear', fill_value='extrapolate')
        dy = np.exp(df(dat['x']))
        total = np.exp(f_oth(dat['x'])) + np.exp(f_jetmed(dat['x']))

        ax.plot(dat['x'], total, color=color)
        ax.fill_between(dat['x'], total+dy,total-dy, color=color, alpha=0.2)
        if cent == '10-20':
            continue
        ratio = total/dat['y']
        d_pos = dat['dy_stat+']**2 + dat['dy_syst+']**2
        d_neg = dat['dy_stat-']**2 + dat['dy_syst-']**2
        dy_sq = dy**2
        err_pos = ratio*np.sqrt(dy_sq/total**2 + d_pos/dat['y']**2)
        err_neg = ratio*np.sqrt(dy_sq/total**2 + d_neg/dat['y']**2)
        delx = 0.5*(dat['xhigh'] - dat['xlow'])
        ax = axes1[1][icent]
        marker = 'x'
        ax.scatter(dat["x"], ratio, color=color, marker=marker, s=60)
        ax.errorbar(dat["x"], ratio, xerr=delx, yerr=[err_neg, err_pos], fmt='none',color=color)

axes1[0][0].set_yscale('log')
axes1[0][0].set_ylabel(r'$E_{\gamma} \frac{dN^{\gamma}}{d^3p^{\gamma}}$ (GeV${}^{-2}$)')
axes1[1][0].set_ylabel('Theory/Data')
for ax in axes1[0]:
    ax.set_yscale('log')
for ax in axes1[1]:
    ax.set_xlabel(r'$p^{\gamma}_T$ (GeV)')
axes1[1][-1].axis('off')

labels = [Line2D([],[],color=rate_colours[rate_names[r]], label=rate_names[r]) for r in rate_names]
labels.append(Line2D([],[],marker='x',color='black',label='00-05$\%$'))
labels.append(Line2D([],[],marker='s',color='black',label='05-10$\%$'))
labels.append(Line2D([],[],marker='<',color='black',label='10-15$\%$'))
labels.append(Line2D([],[],marker='>',color='black',label='15-20$\%$'))
axes1[1][-1].legend(loc='best', handles=labels, ncols=2, fontsize=18)

## Figure 2: 0-20%
# plot data first:
ax = axes2[0]
util.plot_expr_data_on_axis(ax, PHENIX['00-20'], marker='s')
dx_STAR = np.array([x/20 for x in STAR['x']])
errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                for x, delx, y, yerrlow, yerrhigh in
                zip(STAR["x"], dx_STAR, STAR["y"], STAR["dy_syst+"], -1*STAR["dy_syst-"])]

pc = PatchCollection(errorboxes, facecolor='grey', alpha=0.5)
ax.add_collection(pc)
ax.scatter(STAR["x"], STAR["y"], color='black', marker='^', s=60)
ax.errorbar(STAR["x"], STAR["y"], yerr=[STAR["dy_stat+"], -1*STAR["dy_stat-"]], 
                                    fmt='none', color='black')
## add the jet medium photons

tmp_o = 0.25*(other_photons['00-05']['total'][1] + other_photons['05-10']['total'][1] +other_photons['10-20']['total'][1])
f_oth = interpolate(other_photons['00-05']['total'][0], np.log(tmp_o), kind='quadratic', fill_value='extrapolate')
merged_x = PHENIX['00-20']['x'].to_list()
merged_x.extend(STAR['x'].to_list())
merged_x = sorted(merged_x)
ax = axes2[0]
locc = '/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/v2/production_v2/compare/photon_sp/'
for rset in rate_names:
    color  = rate_colours[rate_names[rset]]
    spec   = jet_medium[rset]
    tmp_y  = 0.25*(spec['00-05']['total'] + spec['05-10']['total'] + spec['10-20']['total'])
    tmp_dy = 0.25*np.sqrt(spec['00-05']['dtotal']**2 + spec['05-10']['dtotal']**2 + spec['10-20']['dtotal']**2)
    f_jm   = interpolate(spec['00-05']['pT'], np.log(tmp_y), kind='quadratic', fill_value='extrapolate')
    df_jm  = interpolate(spec['00-05']['pT'], np.log(tmp_dy), kind='quadratic', fill_value='extrapolate')

    total = np.exp(f_jm(merged_x)) + np.exp(f_oth(merged_x))

    with open(locc + f'photons_00-20_rset_{rset}_200_martini_alone.csv', 'w') as f:
        f.write('pT,jmed,othr,total\n')
        for item in zip(merged_x, np.exp(f_jm(merged_x)) , np.exp(f_oth(merged_x)), total):
            line = [f'{v:0.5e}' for v in item]
            f.write(','.join(line)+'\n')
    ax.plot(merged_x, total, color=color)

    ## Now do ratio:
    y_phenix = np.array([np.exp(f_jm(x)) + np.exp(f_oth(x)) for x in PHENIX['00-20']['x']])
    ratio_PHENIX = y_phenix/PHENIX['00-20']['y']
    d_pos = PHENIX['00-20']['dy_stat+']**2 + PHENIX['00-20']['dy_syst+']**2
    d_neg = PHENIX['00-20']['dy_stat-']**2 + PHENIX['00-20']['dy_syst-']**2
    dy_sq = np.exp(df_jm(PHENIX['00-20']['x']))**2
    err_pos = ratio_PHENIX*np.sqrt(dy_sq/y_phenix**2 + d_pos/PHENIX['00-20']['y']**2)
    err_neg = ratio_PHENIX*np.sqrt(dy_sq/y_phenix**2 + d_neg/PHENIX['00-20']['y']**2)
    delx = 0.5*(PHENIX['00-20']['xhigh'] - PHENIX['00-20']['xlow'])
    axes2[1].scatter(PHENIX['00-20']['x'], ratio_PHENIX, marker='s', color=color)
    axes2[1].errorbar(PHENIX['00-20']['x'], ratio_PHENIX, xerr=delx, yerr=[err_neg, err_pos], fmt='none',color=color)

    y_STAR = np.exp(f_jm(STAR['x'])) + np.exp(f_oth(STAR['x']))
    ratio_STAR = y_STAR/STAR['y']
    d_pos = STAR["dy_stat+"]**2 + STAR["dy_syst+"]**2
    d_neg = STAR["dy_syst-"]**2 + STAR["dy_stat-"]**2
    dy_sq = np.exp(df_jm(STAR['x']))**2
    err_pos = ratio_STAR*np.sqrt(dy_sq/y_STAR**2 + d_pos/STAR['y']**2)
    err_neg = ratio_STAR*np.sqrt(dy_sq/y_STAR**2 + d_neg/STAR['y']**2)
    axes2[1].scatter (STAR['x'], ratio_STAR, marker='^', color=color)
    axes2[1].errorbar(STAR['x'], ratio_STAR, xerr=dx_STAR, yerr=[err_neg, err_pos], fmt='none',color=color)

handles = [Line2D([],[],color=rate_colours[rate_names[r]], label=rate_names[r]) for r in rate_names]
handles.append(Line2D([],[],color='black', marker='s', label='PHENIX (2012) $0$-$20\%$ (constructed)'))
handles.append(Line2D([],[],color='black', marker='^', label='STAR (2017) $0$-$20\%$'))
axes2[0].legend(loc='best', handles=handles)
axes2[0].set_ylabel(r'$E_{\gamma} \frac{dN^{\gamma}}{d^3p^{\gamma}}$ (GeV${}^{-2}$)')
axes2[1].set_ylabel('Theory/Data')
axes2[0].set_yscale('log')

plt.show()