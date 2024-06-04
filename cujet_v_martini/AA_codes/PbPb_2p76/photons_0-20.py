#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d
import util
import jetDicts as ddicts
from dictionaries import multiplicity_PbPb_2p76 as multiplicity
from dictionaries import oversampling_factors, channel_linestyles
from COLORS import module_colors, colors, channel_colors
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 35,
    "lines.linewidth": 4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size" : 12,
    "ytick.major.size" : 12,
    "xtick.minor.size" : 6,
    "ytick.minor.size" : 6,
    "axes.spines.right": False,
    "axes.spines.top"  : False,
    "legend.frameon"   : False,
    "axes.labelsize" : 40}
plt.rcParams.update(my_rcParams)

#do_incnlo = int(sys.argv[1])
do_incnlo = 1

inel_Xsec = 62.03948
elosses = ['cujet', 'martini']

## read in experimental data
exp_loc = "../../../expt/PbPb_2p76/photons/"
centralities = {'00-20':exp_loc+'HEPData-ins1394677-v1-Table_1.csv'}

expdata = {}
pT_lower_lim = {'00-20':2,'20-40':0}
pT_upper_lim = {'00-20':2,'20-40':0}
for cent in centralities:
    tmp = pd.read_csv(centralities[cent],comment='#')
    tmp['dx'] = tmp['xhigh']-tmp['xlow']
    #min_pT = min(tmp['x'])
    #max_pT = max(tmp['x'])
    #pT_lower_lim[cent] = min_pT
    #pT_upper_lim[cent] = max_pT
    #tmp = tmp#[tmp['x'].between(pT_lower_lim[cent],pT_upper_lim[cent])]
    expdata[cent] = tmp
## read jetscape:
elosses = ['martini','cujet']
jetscape_cents = ['00-05','05-10','10-20']
jetscape = {}
fname_tmpl = '../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_photon_spec_0.80.csv'
for eloss in elosses:
    jetscape[eloss] = {}
    for cent in jetscape_cents:
        fname = fname_tmpl.format(eloss=eloss,cent=cent)
        factor = oversampling_factors[f'{eloss}-{cent}']
        Nbin = multiplicity[cent]
        tmp = pd.read_csv(fname,comment='#')
        xmin = tmp['ptmin']
        xmax = tmp['ptmax']
        x = 0.5*(xmin+xmax)
        dx = xmax - xmin
        tmp['conv']  *= Nbin/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        tmp['dconv'] *= np.sqrt(Nbin)/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        tmp['brem']  *= Nbin/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        tmp['dbrem'] *= np.sqrt(Nbin)/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        jetscape[eloss][cent] = tmp
    ## construct the 0-20 and 20-40 centralities from the above
    tmp1 = util.combine_centralities(jetscape[eloss]['00-05'],jetscape[eloss]['05-10'])
    tmp2 = util.combine_centralities(tmp1, jetscape[eloss]['10-20'])
    tmp2['total'] = tmp2['conv'] + tmp2['brem']
    tmp2['dtotal'] = np.sqrt(tmp2['dconv']**2 + tmp2['dbrem']**2)
    jetscape[eloss]['00-20'] = tmp2

## Read in the prompt and thermal calculations of JF Paquet
channels = ['thermal','preEq','prompt']
nice_channels = {'prompt':'Prompt','thermal':'Thermal','preEq':'Pre-Equilibrium'}
other_photons = {}
tmpl_jf = "../../other_data/JF_MultiMessenger/PbPb2760_{cent}_{chan}.csv"
for cent in centralities:
    other_photons[cent] = {}
    for ch in channels:
        if ch == 'prompt':
            continue
        tmp = np.loadtxt(tmpl_jf.format(cent=cent,chan=ch),unpack=True,delimiter=',')
        x = tmp[0]
        y = tmp[1]
        other_photons[cent][ch] = pd.DataFrame({'pT':x,'N':y, 'dN':np.zeros_like(y)})

tmp = pd.read_csv('../../jetscape_data/prompt_photons/PbPb_2760/gamma_spectra.csv', comment='#')
tmp_x   = 0.5*(tmp['pTmin'] + tmp['pTmax'])
tmp_dx  = tmp['pTmax'] - tmp['pTmin']
tmp_y   = tmp['prompt']/(2*np.pi*tmp_x*tmp_dx*inel_Xsec*2*0.8)
tmp_dy  = tmp['dprompt']/(2*np.pi*tmp_x*tmp_dx*inel_Xsec*2*0.8)
yy  = tmp_y*multiplicity['00-20']
dyy = tmp_dy*np.sqrt(multiplicity['00-20'])
# find k-factor:
f = interp1d(tmp_x, np.log(yy), kind='linear', fill_value='extrapolate')
dat_x, dat_y = expdata['00-20']['x'].to_list(), expdata['00-20']['y'].to_list()
kfactor = dat_y[-1]/np.exp(f(dat_x[-1]))
print(kfactor)
yy = yy * kfactor
dyy = dyy * np.sqrt(kfactor)
other_photons['00-20']['prompt'] = pd.DataFrame({'pT':tmp_x, 'N':yy.to_list(), 'dN':dyy.to_list()})

ratios = {'martini':{},'cujet':{}}

x = np.linspace(1, 22, 15)
for cent in ['00-20']:
    specs = other_photons[cent]
    f_prompt  = interp1d(specs['prompt']['pT'] , np.log(specs['prompt']['N'])  , kind='linear', fill_value='extrapolate')
    f_thermal = interp1d(specs['thermal']['pT'], np.log(specs['thermal']['N']) , kind='linear', fill_value='extrapolate')
    f_preEq   = interp1d(specs['preEq']['pT']  , np.log(specs['preEq']['N'])   , kind='linear', fill_value='extrapolate')
    df_prompt = interp1d(specs['prompt']['pT'] , np.log(specs['prompt']['dN']) , kind='linear', fill_value='extrapolate')
    ypreq = np.exp(f_preEq(x))
    ytherm = np.exp(f_thermal(x))
    yprompt = np.exp(f_prompt(x))
    ytot    = ypreq + ytherm + yprompt
    dytot   = np.exp(df_prompt(x))
    other_photons[cent]['total'] = (x,ytot,dytot)
    for eloss in jetscape:
        ratios[eloss] = {}
        jspec = jetscape[eloss][cent]
        f_jetmed = interp1d(jspec['pT'], np.log(jspec['total']), kind='linear', fill_value='extrapolate')
        dfjetmed = interp1d(jspec['pT'], np.log(jspec['dtotal']), kind='linear', fill_value='extrapolate')
        yjetm = np.exp(f_jetmed(x))
        dyjetm = np.exp(dfjetmed(x))
        total = ytot + yjetm
        ratio_jmed = yjetm/total
        ratio_thermal = ytherm/total
        ratio_preq = ypreq/total
        ratio_prompt = yprompt/total
        ratios[eloss]['thermal'] = ratio_thermal
        ratios[eloss]['preEq'] = ratio_preq
        ratios[eloss]['prompt'] = ratio_prompt
        ratios[eloss]['jet_medium'] = ratio_jmed


total_v1 = {}
for eloss in elosses:
    total_v1[eloss] = {}
    for cent in centralities:
        tmpx        = jetscape[eloss][cent]['pT']
        dx          = jetscape[eloss][cent]['dpT']
        spec        = jetscape[eloss][cent]['conv']+ jetscape[eloss][cent]['brem']
        dspec       = np.sqrt(jetscape[eloss][cent]['dconv']**2 +jetscape[eloss][cent]['dbrem']**2)
        fjetmed     = interp1d(tmpx, np.log(spec), kind='linear', fill_value='extrapolate')
        dfjetmed    = interp1d(tmpx, np.log(dspec), kind='linear', fill_value='extrapolate')
        oth         = other_photons[cent]['total']
        total_spec  = oth[1] + np.exp(fjetmed(x))
        dtotal_spec = np.sqrt(oth[2]**2 + np.exp(dfjetmed(x))**2)
        total_v1[eloss][cent] = (x, dx, total_spec, dtotal_spec)

fig, axes = plt.subplots(2,1,figsize=(16,9),sharex=True, height_ratios=(2,1))

for icent, cent in enumerate(expdata.keys()):
    spec_ax = axes[0]
    ratio_ax = axes[1]
    data = expdata[cent]
    util.plot_expr_data_on_axis(spec_ax, data, marker='o')
    spec_ax.text(0.1,0.1,f'{cent}\%', transform=spec_ax.transAxes)
    for eloss in elosses:
        pT,_,spec, dspec = total_v1[eloss][cent]
        ftot = interp1d(pT, np.log(spec), kind='linear', fill_value='extrapolate')
        dftot = interp1d(pT, np.log(dspec), kind='linear', fill_value='extrapolate')
        spec = np.exp(ftot(data['x']))
        dspec = np.exp(dftot(data['x']))
        pT = data['x']
        dpT = data['dx']
        col  = module_colors['MATTER+'+eloss.upper()]
        mark = ddicts.eloss_marker[eloss]
        zorder = 0.5 if eloss=='martini' else 0.45
        spec_ax.plot(pT, spec, color=col, zorder=zorder)
        spec_ax.fill_between(pT, spec-dspec, spec+dspec, color=col, alpha=0.2, zorder=zorder)
        ## take ratio
        ratio = spec / data['y']
        dratio_stat_plus  = ratio*np.sqrt(data['dy_stat+']**2/data['y']**2 + dspec**2/spec**2)
        dratio_stat_minus = ratio*np.sqrt(data['dy_stat+']**2/data['y']**2 + dspec**2/spec**2)
        dratio_syst_plus  = ratio*data['dy_syst+']/data['y']
        dratio_syst_minus = ratio*data['dy_syst-']/data['y']

        errorboxes_2 = [Rectangle((x-delx, y - abs(dy_minus)), width=2*delx, height=(abs(dy_minus)+abs(dy_plus)), zorder=zorder)
                    for x, delx, y, dy_minus, dy_plus in zip(pT, 0.5*dpT, ratio, dratio_syst_minus, dratio_syst_plus)]
        pc2 = PatchCollection(errorboxes_2, facecolor=col, edgecolor=col, alpha=0.5)
        ratio_ax.add_collection(pc2)
        ratio_ax.errorbar(pT, ratio, xerr=0.5*dpT,yerr=[dratio_stat_minus, dratio_syst_plus], color=col,marker=mark,linewidth=1,fmt='none', zorder=zorder)
        ratio_ax.scatter(pT, ratio, color=col,marker=mark, zorder=zorder)
    ## Non jet medium ones:
    xspec, spec, dspec= other_photons[cent]['total']
    ftot  = interp1d(xspec, np.log(spec), kind='linear', fill_value='extrapolate')
    dftot = interp1d(xspec, np.log(dspec), kind='linear', fill_value='extrapolate')
    spec  = np.exp(ftot(data['x']))
    dspec = np.exp(dftot(data['x']))
    pT    = data['x']
    col   = colors[2]
    mark  = 's'#ddicts.eloss_marker[eloss]
    spec_ax.plot(pT, spec, color=col, zorder=0)
    spec_ax.fill_between(pT, spec-dspec, spec+dspec, color=col, alpha=0.2, zorder=0)
    ## take ratio
    ratio = spec / data['y']
    dratio_stat_plus  = ratio*np.sqrt(data['dy_stat+']**2/data['y']**2 + dspec**2/spec**2)
    dratio_stat_minus = ratio*np.sqrt(data['dy_stat+']**2/data['y']**2 + dspec**2/spec**2)
    dratio_syst_plus  = ratio*data['dy_syst+']/data['y']
    dratio_syst_minus = ratio*data['dy_syst-']/data['y']

    errorboxes_2 = [Rectangle((x-delx, y - abs(dy_minus)), width=2*delx, height=(abs(dy_minus)+abs(dy_plus)))
                for x, delx, y, dy_minus, dy_plus in zip(pT, 0.5*dpT, ratio, dratio_syst_minus, dratio_syst_plus)]
    pc2 = PatchCollection(errorboxes_2, facecolor=col, edgecolor=col, alpha=0.5)
    ratio_ax.add_collection(pc2)
    ratio_ax.errorbar(pT, ratio, xerr=0.5*dpT,yerr=[dratio_stat_minus, dratio_syst_plus], color=col,marker=mark,linewidth=1,fmt='none')
    ratio_ax.scatter(pT, ratio, color=col,marker=mark)


handles = [Line2D([],[],label=r'ALICE (2016) $|\eta|<0.8$', marker='o',color='black')]
for eloss in elosses:
    handles.append(Line2D([],[],color=module_colors['MATTER+'+eloss.upper()],label='MATTER+'+eloss.upper(),marker=ddicts.eloss_marker[eloss]))
handles.append(Line2D([],[],label='No Jet-Medium', color=colors[2], marker='s'))
axes[0].set_yscale('log')
axes[0].set_ylabel(r'$\frac{1}{N_{\mathrm{evt}}\,2\pi\,p_T}\frac{\mathrm{d}N^{\gamma}}{\mathrm{d}p_T\mathrm{d}\eta}$, (GeV$^{-2}$)', fontsize=30)
axes[1].set_ylabel(r'Theory over' + '\n' + 'Data', fontsize=30)
axes[1].axhline(1,color='black', linestyle='dotted')
axes[1].set_xlabel(r'$p_T$ (GeV)', fontsize=30)
axes[0].legend(loc='upper right', handles=handles)
axes[0].text(0.1,0.9, r'Pb-Pb @ $\sqrt{s}=2.76$ ATeV', transform=axes[0].transAxes)


nice_channels = {'prompt':'Prompt','thermal':'Thermal','preEq':'Pre-Equilibrium'}
fig, axes = plt.subplots(1,2,figsize=(16,9),
                        sharey=True, sharex=True)

for ieloss, eloss in enumerate(elosses):
    spec_ax = axes[ieloss]
    if ieloss == 0:
        #spec_ax.text(0.1,0.7,f'{cent}\%', transform=spec_ax.transAxes)
        spec_ax.set_ylabel('Channel to Total')
    spec_ax.set_xlabel('$p_T$ (GeV)')
    pT,_,_,_ = total_v1[eloss][cent]
    rs   = ratios[eloss]
    for channel in rs:
        col = ''
        if channel != 'jet_medium':
            col = channel_colors[channel]
        else:
            col = module_colors[f'MATTER+{eloss.upper()}']
        lstyle = channel_linestyles[channel]
        #ratio= rs[eloss]['jet_medium'].tolist()
        mark='s'
        spec_ax.plot(pT, rs[channel], color=col,linestyle=lstyle)

jf_handles = [Line2D([],[],color='white')]
for ch in channels:
    jf_handles.append(Line2D([],[],color=channel_colors[ch], linestyle=channel_linestyles[ch],label=nice_channels[ch]))

handles= [Line2D([],[],color=module_colors['MATTER+MARTINI'], label='MATTER+MARTINI')]
artist= axes[0].legend(loc='center right', handles=handles   , fontsize=30)
axes[0].add_artist(artist)
handles= [Line2D([],[],color=module_colors['MATTER+CUJET'], label='MATTER+CUJET')]
artist = axes[1].legend(loc='center right', handles=handles   , fontsize=30)
axes[1].add_artist(artist)
axes[0].legend(loc='upper left',bbox_to_anchor=(0.0,1.1) , handles=jf_handles, fontsize=30, )
fig.suptitle('Pb-Pb @ $\sqrt{s}=2.76$ ATeV, $|\eta|<0.8$, 0-20$\%$')
plt.show()
