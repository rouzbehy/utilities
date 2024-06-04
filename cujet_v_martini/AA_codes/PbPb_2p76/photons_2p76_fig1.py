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
from dictionaries import oversampling_factors
from COLORS import module_colors, colors

plt.rcParams.update(util.my_rcParams)

do_incnlo = int(sys.argv[1])
print("Are we doing INCNLO? ", bool(do_incnlo))

inel_Xsec = 62.03948 
elosses = ['cujet', 'martini']

## read in experimental data
exp_loc = "../../../expt/PbPb_2p76/photons/"
centralities = {'00-20':exp_loc+'HEPData-ins1394677-v1-Table_1.csv',
                '20-40':exp_loc+'HEPData-ins1394677-v1-Table_2.csv'}

expdata = {}
pT_lower_lim = {'00-20':0,'20-40':0}
pT_upper_lim = {'00-20':0,'20-40':0}
for cent in centralities:
    tmp = pd.read_csv(centralities[cent],comment='#')
    min_pT = min(tmp['x'])
    max_pT = max(tmp['x'])
    pT_lower_lim[cent] = 2#min_pT
    pT_upper_lim[cent] = 14#max_pT
    tmp = tmp[tmp['x'].between(pT_lower_lim[cent],pT_upper_lim[cent])]
    expdata[cent] = tmp
## read jetscape:
elosses = ['martini','cujet']
jetscape_cents = ['00-05','05-10','10-20','20-30','30-40']
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
        tmp['dconv'] *= np.sqrt(Nbin)/(np.sqrt(factor)*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        tmp['brem']  *= Nbin/(factor*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        tmp['dbrem'] *= np.sqrt(Nbin)/(np.sqrt(factor)*(x*dx*inel_Xsec*2*np.pi*2*0.8))
        jetscape[eloss][cent] = tmp
    ## construct the 0-20 and 20-40 centralities from the above
    tmp1 = util.combine_centralities(jetscape[eloss]['00-05'],jetscape[eloss]['05-10'])
    tmp2 = util.combine_centralities(tmp1, jetscape[eloss]['10-20'])
    jetscape[eloss]['00-20'] = tmp2[tmp2['pT'].between(pT_lower_lim['00-20'],pT_upper_lim['00-20'])]
    tmp3 = util.combine_centralities(jetscape[eloss]['20-30'],jetscape[eloss]['30-40'])
    jetscape[eloss]['20-40'] = tmp3[tmp3['pT'].between(pT_lower_lim['20-40'],pT_upper_lim['20-40'])]

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
if do_incnlo:
    for cent in other_photons:
        tmp = np.loadtxt(tmpl_jf.format(cent=cent,chan='prompt'), unpack=True, delimiter=',')
        tmp = pd.DataFrame({'pT':tmp[0], 'N':tmp[1], 'dN':np.zeros_like(tmp[1])})
        other_photons[cent]['prompt'] = tmp
kfactor = -1
if not do_incnlo:
    tmp = pd.read_csv('../../jetscape_data/prompt_photons/PbPb_2760/gamma_spectra.csv', comment='#')
    tmp_x   = 0.5*(tmp['pTmin'] + tmp['pTmax'])
    tmp_dx  = tmp['pTmax'] - tmp['pTmin']
    tmp_y   = tmp['prompt']/(2*np.pi*tmp_x*tmp_dx*inel_Xsec*2*0.8)
    tmp_dy  = tmp['dprompt']/(2*np.pi*tmp_x*tmp_dx*inel_Xsec*2*0.8)
    for cent in centralities:
        yy  = tmp_y*multiplicity[cent]
        dyy = tmp_dy*np.sqrt(multiplicity[cent])
        if kfactor < 0 and cent =='00-20':
            f = interp1d(tmp_x, np.log(yy), kind='linear', fill_value='extrapolate')
            dat_x, dat_y = expdata['00-20']['x'].to_list(), expdata['00-20']['y'].to_list()
            print(dat_x[-1])
            kfactor = dat_y[-1]/np.exp(f(dat_x[-1]))
        print(kfactor)
        yy = kfactor*yy
        dyy = np.sqrt(kfactor)*dyy
        other_photons[cent]['prompt'] = pd.DataFrame({'pT':tmp_x.to_list(), 'N':yy.to_list(), 'dN':dyy.to_list()})
 
x = np.linspace(1, 22, 15)
for cent in other_photons:
    specs = other_photons[cent]
    f_prompt  = interp1d(specs['prompt']['pT'] , np.log(specs['prompt']['N'])  , kind='linear', fill_value='extrapolate')
    f_thermal = interp1d(specs['thermal']['pT'], np.log(specs['thermal']['N']) , kind='linear', fill_value='extrapolate')
    f_preEq   = interp1d(specs['preEq']['pT']  , np.log(specs['preEq']['N'])   , kind='linear', fill_value='extrapolate')
    df_prompt = interp1d(specs['prompt']['pT'] , np.log(specs['prompt']['dN']) , kind='linear', fill_value='extrapolate')
    ytot      = np.exp(f_preEq(x)) + np.exp(f_prompt(x)) + np.exp(f_thermal(x)) 
    dytot     = np.exp(df_prompt(x))
    other_photons[cent]['total'] = (x,ytot,dytot)

total_v1 = {}
for eloss in elosses:
    total_v1[eloss] = {}
    for cent in centralities:
        tmpx  = jetscape[eloss][cent]['pT']
        dx    = jetscape[eloss][cent]['dpT']
        spec  = jetscape[eloss][cent]['conv']+ jetscape[eloss][cent]['brem']
        dspec = np.sqrt(jetscape[eloss][cent]['dconv']**2 +jetscape[eloss][cent]['dbrem']**2)
        fjetmed     = interp1d(tmpx, np.log(spec), kind='linear', fill_value='extrapolate')
        dfjetmed    = interp1d(tmpx, np.log(dspec), kind='linear', fill_value='extrapolate')
        oth         = other_photons[cent]['total']
        total_spec  = oth[1] + np.exp(fjetmed(x))
        dtotal_spec = np.sqrt(oth[2]**2 + np.exp(dfjetmed(x))**2)
        total_v1[eloss][cent] = (x, dx, total_spec, dtotal_spec)

fig, axes = plt.subplots(2,2,figsize=(16,12),sharey='row', sharex=True, height_ratios=(2,1))

for icent, cent in enumerate(expdata.keys()):
    spec_ax = axes[0][icent]
    ratio_ax = axes[1][icent]
    data = expdata[cent]
    util.plot_expr_data_on_axis(spec_ax, data, marker='o')
    spec_ax.text(0.1,0.1,f'{cent}\%', transform=spec_ax.transAxes)
    for eloss in elosses:
        pT, dpT, spec, dspec = total_v1[eloss][cent]
        ftot = interp1d(pT, np.log(spec), kind='linear', fill_value='extrapolate')
        dftot = interp1d(pT, np.log(dspec), kind='linear', fill_value='extrapolate')
        spec = np.exp(ftot(data['x']))
        dspec = np.exp(dftot(data['x']))
        pT = data['x']
        col  = module_colors['MATTER+'+eloss.upper()]
        mark = ddicts.eloss_marker[eloss]
        zorder = 5 if eloss=='martini' else 3
        spec_ax.plot(pT, spec, color=col, zorder=zorder)
        spec_ax.fill_between(pT, spec-dspec, spec+dspec, color=col, alpha=0.2, zorder=zorder)

        ## take ratio
        ratio = spec / data['y']
        dratio_stat_plus  = ratio*np.sqrt(data['dy_stat+']**2/data['y']**2 + dspec**2/spec**2)
        dratio_stat_minus = ratio*np.sqrt(data['dy_stat+']**2/data['y']**2 + dspec**2/spec**2)
        dratio_syst_plus  = ratio*data['dy_syst+']/data['y']#np.sqrt(data['dy_syst+']**2/data['y']**2 + dspec**2/spec**2)
        dratio_syst_minus = ratio*data['dy_syst-']/data['y']#np.sqrt(data['dy_syst-']**2/data['y']**2 + dspec**2/spec**2)

        errorboxes_2 = [Rectangle((x-delx, y - abs(dy_minus)), width=2*delx, height=(abs(dy_minus)+abs(dy_plus)), zorder=zorder)
                    for x, delx, y, dy_minus, dy_plus in zip(pT, 0.5*dpT, ratio, dratio_syst_minus, dratio_syst_plus)]
        pc2 = PatchCollection(errorboxes_2, facecolor=col, edgecolor=col, alpha=0.5)
        ratio_ax.add_collection(pc2)
        ratio_ax.errorbar(pT, ratio, xerr=0.5*dpT,yerr=[dratio_stat_minus, dratio_syst_plus], color=col,marker=mark,linewidth=1,fmt='none', zorder=zorder)
        ratio_ax.scatter(pT, ratio, color=col,marker=mark, zorder=zorder)
    ## Non jet medium ones:
    
    xspec, spec, dspec= other_photons[cent]['total']
    ftot = interp1d(xspec, np.log(spec), kind='linear', fill_value='extrapolate')
    dftot = interp1d(xspec, np.log(dspec), kind='linear', fill_value='extrapolate')
    spec = np.exp(ftot(data['x']))
    dspec = np.exp(dftot(data['x']))
    pT = data['x']
    col  = colors[2]
    mark = 's'#ddicts.eloss_marker[eloss]
    spec_ax.plot(pT, spec, color=col)
    spec_ax.fill_between(pT, spec-dspec, spec+dspec, color=col, alpha=0.2)
    ## take ratio
    ratio = spec / data['y']
    dratio_stat_plus  = ratio*np.sqrt(data['dy_stat+']**2/data['y']**2 + dspec**2/spec**2)
    dratio_stat_minus = ratio*np.sqrt(data['dy_stat+']**2/data['y']**2 + dspec**2/spec**2)
    dratio_syst_plus  = ratio*data['dy_syst+']/data['y']#np.sqrt(data['dy_syst+']**2/data['y']**2 + dspec**2/spec**2)
    dratio_syst_minus = ratio*data['dy_syst-']/data['y']#np.sqrt(data['dy_syst-']**2/data['y']**2 + dspec**2/spec**2)

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
axes[0][0].set_yscale('log')
axes[0][0].set_ylabel(r'$\frac{1}{N_{\mathrm{evt}}\,2\pi\,p_T}\frac{\mathrm{d}N^{\gamma}}{\mathrm{d}p_T\mathrm{d}\eta}$, (GeV$^{-2}$)', fontsize=30)
axes[1][0].set_ylabel(r'Theory over' + '\n' + 'Data', fontsize=30)
for ax in axes[1]:
    ax.axhline(1,color='black', linestyle='dotted')
    ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30)
axes[0][1].legend(loc='upper right', handles=handles)
#axes[0][1].text(0.7, 0.55, r'$|\eta|<0.8$',transform=axes[0][1].transAxes)
axes[0][0].text(0.2,0.85, r'Pb-Pb @ $\sqrt{s}=2.76$ ATeV', transform=axes[0][0].transAxes)
#fig.savefig("../../Plots/Photon_Plots/fig1_photon_spec_JF_prompts.pdf", dpi=200)
plt.show()
