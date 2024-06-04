#!/usr/bin/env python3
import sys
import numpy as np
from pandas import read_csv, DataFrame
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d
from matplotlib import ticker
import util 
import dictionaries as ddicts
from COLORS import module_colors, colors, channel_colors
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
sigma = 68
avail_cents = ['00-10','10-20','30-50']


## read in Experiment: 
colnames = {'PT [GEV/c]':'x',
            'PTLOW':'xlow',
            'PTHIGH':'xhigh',
            '1/(2*PI*Nevt)*(d**2N)/(PT*dPT*dY) [(GEV/c)**-2]':'y',
            'stat':'dy_stat',
            'sys':'dy_syst'}
fname_exp = "../../../expt/PbPb_5p02/Photons/fileDirPhotonSpectrum.txt"
alice = read_csv(fname_exp, comment='#').rename(columns=colnames)
alice['dy_stat+'] = alice['dy_stat']
alice['dy_stat-'] = alice['dy_stat']
alice['dy_syst+'] = alice['dy_syst']
alice['dy_syst-'] = alice['dy_syst']
alice = alice[alice['xlow'] > 2.6]
pT = alice['x']

## read in JETSCAPE events:
tmpl = "../../jetscape_data/sqrt_s_5020/maxt_200/{eloss}/PbPb5020_{cent}_photon_spec_0.80.csv"
data = {}
x = 0
for cent in avail_cents:
    data[cent] = {}
    for eloss in ['martini', 'cujet']:
        oversample = ddicts.exag_factor[eloss]
        ncoll = ddicts.multiplicity[cent]
        tmp = read_csv(tmpl.format(eloss=eloss, cent=cent), comment='#')
        tmp['pT'] = 0.5*(tmp['ptmin']+tmp['ptmax']) 
        x = tmp['pT']
        dx = tmp['ptmax'] - tmp['ptmin']
        tmp['dpT'] = dx
        for tag in ['conv','brem']:
            tmp['y'+tag]  = ncoll*tmp[tag]/(2*np.pi*2*0.8*tmp['pT']*oversample*dx*sigma)
            tmp['dy'+tag] = ncoll*tmp['d'+tag]/(2*np.pi*2*0.8*tmp['pT']*oversample*dx*sigma)
        tmp['N']  = (tmp['yconv'] + tmp['ybrem'])
        tmp['dN'] = np.sqrt(tmp['dyconv']*tmp['dyconv']+tmp['dybrem']*tmp['dybrem'])
        data[cent][eloss] = tmp

data['00-20'] = {}
 
## Read in JF's Photons:
for cent in ['00-20', '30-50']:
    if cent =='00-20':
        data[cent] = {}
        elossmods = {}
        for eloss in ['martini', 'cujet']:
            yconv1 = data['00-10'][eloss]['yconv']
            yconv2 = data['10-20'][eloss]['ybrem']
            ybrem1 = data['00-10'][eloss]['yconv']
            ybrem2 = data['10-20'][eloss]['ybrem']
        
            dyconv1 = data['00-10'][eloss]['dyconv']
            dyconv2 = data['10-20'][eloss]['dybrem']
            dybrem1 = data['00-10'][eloss]['dyconv']
            dybrem2 = data['10-20'][eloss]['dybrem']
        
            conv_0_20  = 0.5*(yconv1 + yconv2)
            dconv_0_20 = 0.5*np.sqrt(dyconv1*dyconv1 + dyconv2*dyconv2)
            brem_0_20  = 0.5*(ybrem1 + ybrem2)
            dbrem_0_20 = 0.5*np.sqrt(dybrem1*dybrem1 + dybrem2*dybrem2)
            total_jet_medium = conv_0_20 + brem_0_20
            dtotal_jet_medium = np.sqrt(dconv_0_20*dconv_0_20 + dbrem_0_20*dbrem_0_20)
            data[cent][eloss]= DataFrame({'pT':x, 'N':total_jet_medium, 'dN':dtotal_jet_medium})

    for ch in ['prompt','preEq','thermal']:
        if ch == "prompt":
            continue
        fname = f'../../other_data/JF_MultiMessenger/PbPb5020_{cent}_{ch}.csv'
        tmp = np.loadtxt(fname, delimiter=',',unpack=True)
        data[cent][ch] = DataFrame({'pT':tmp[0], 'N':tmp[1], 'dN':np.zeros_like(tmp[1])})

#kfactor = float(sys.argv[1])
tmp     = read_csv('../../jetscape_data/prompt_photons/PbPb_5020/gamma_spectra.csv', comment='#')
#tmp     = read_csv('/Users/rmyazdi/Documents/research/PYTHIA_CALCS/prompt_photons/results.csv', comment='#')
tmp_x   = 0.5*(tmp['pTmin'] + tmp['pTmax'])
tmp_dx  = tmp['pTmax'] - tmp['pTmin']
tmp_y   = tmp['prompt']/(2*np.pi*tmp_x*tmp_dx*sigma*2*0.8)
tmp_dy  = tmp['dprompt']/(2*np.pi*tmp_x*tmp_dx*sigma*2*0.8)

for cent in ['00-20','30-50']:
    yy  = tmp_y *ddicts.multiplicity[cent]
    dyy = tmp_dy*ddicts.multiplicity[cent]
    if cent == '00-20':
        f = interp1d(tmp_x, np.log(yy))
        y_alice = alice['y'].to_list()
        x_alice = pT.to_list()
        kfactor = y_alice[-1]/np.exp(f(x_alice[-1]))
        print(f"cent:{cent}, kfactor:{kfactor}")
        yy = kfactor*yy
        dyy = np.sqrt(kfactor)*dyy
        
    data[cent]['prompt'] = DataFrame({'pT':tmp_x.to_list(), 'N':yy.to_list(), 'dN':dyy.to_list()})



results = {}
martini = {}
cujet = {}
pT = np.linspace(2, 22, 23)
for cent in ['00-20','30-50']:
    results[cent] = {}
    martini[cent] = {}
    cujet[cent] = {}

    dat = data[cent]
    total = np.zeros_like(pT)

    f_jetmed_martini = interp1d(dat['martini']['pT'], np.log(dat['martini']['N']) , kind='linear', fill_value='extrapolate')
    f_jetmed_cujet   = interp1d(dat['martini']['pT'], np.log(dat['cujet']['N'])   , kind='linear', fill_value='extrapolate')
    f_prompt         = interp1d(dat['prompt']['pT'] , np.log(dat['prompt']['N'])  , kind='linear', fill_value='extrapolate')
    f_preEq          = interp1d(dat['preEq']['pT']  , np.log(dat['preEq']['N'])   , kind='linear', fill_value='extrapolate')
    f_thrml          = interp1d(dat['thermal']['pT'], np.log(dat['thermal']['N']) , kind='linear', fill_value='extrapolate')


    jmart = np.exp(f_jetmed_martini(pT))
    jcujt = np.exp(f_jetmed_cujet(pT))
    prmpt = np.exp(f_prompt(pT))
    preeq = np.exp(f_preEq(pT))
    thrml = np.exp(f_thrml(pT))

    results[cent]['martini']       = jmart + prmpt + preeq + thrml
    results[cent]['cujet']         = jcujt + prmpt + preeq + thrml
    results[cent]['no jet-medium'] = prmpt + preeq + thrml
    
    ## now compute ratios: 
    martini[cent]['jet-medium-martini'] = jmart/results[cent]['martini']
    martini[cent]['prompt']             = prmpt/results[cent]['martini']
    martini[cent]['thermal']            = thrml/results[cent]['martini']
    martini[cent]['preEq']              = preeq/results[cent]['martini']
    
    cujet[cent]['jet-medium-cujet'] = jcujt/results[cent]['cujet']
    cujet[cent]['prompt']           = prmpt/results[cent]['cujet']
    cujet[cent]['thermal']          = thrml/results[cent]['cujet']
    cujet[cent]['preEq']            = preeq/results[cent]['cujet']
    # if cent == '00-20':
    #     non_med = prmpt + preeq + thrml 
    #     locc = '/Users/rmyazdi/Documents/research/jetscape_project/v2/jetscape_data/photon_sp/'
    #     with open(locc+f'photons_00-20_5p02_martini_jetscape.csv', 'w') as f:
    #         f.write('pT,jmed,othr,total\n')
    #         for item in zip(pT, jmart, non_med, results[cent]['martini']):
    #             line = [f'{v:0.6e}' for v in item]
    #             f.write(','.join(line)+'\n')
    #     with open(locc+f'photons_00-20_5p02_cujet_jetscape.csv', 'w') as f:
    #         f.write('pT,jmed,othr,total\n')
    #         for item in zip(pT, jcujt, non_med, results[cent]['cujet']):
    #             line = [f'{v:0.6e}' for v in item]
    #             f.write(','.join(line)+'\n')

pT_alice = alice['x']
fig, axes = plt.subplots(1,2,sharex=True, sharey=True, figsize=(16,9))
axes = axes.flatten()
for icent, cent in enumerate(['00-20', '30-50']):
    ax = axes[icent]
    ax.text(0.05, 0.1, cent+r'$\%$', transform=ax.transAxes, fontsize=30)
    if icent == 0:
       util.plot_expr_data_on_axis(ax, alice, '*','black','gray', s=40, factor=1) 
    tmp = results[cent]
    for item in tmp:
        col = 'black'
        if item in ['cujet', 'martini']:
            col = module_colors['MATTER+'+item.upper()]
        else:
            col = colors[2]
        y = tmp[item]
        func = interp1d(pT, np.log(y), kind='linear', fill_value='extrapolate')
        yy = np.exp(func(pT_alice))
        ax.plot(pT_alice, yy, color=col)

labels = [Line2D([],[],color=module_colors['MATTER+'+f.upper()], label='MATTER+'+f.upper()) for f in ['cujet','martini']]
labels.append(Line2D([],[],color=colors[2], label='No Jet-Medium'))

axes[1].legend(loc='upper right', handles=labels)
for ax in axes:
    ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30)
axes[0].set_ylabel(r'$\frac{1}{2\pi p_T}\frac{\mathrm{d}^2N}{\mathrm{d}\eta \mathrm{d}p_T}$ (GeV${}^{-2}$)', fontsize=30)
axes[0].set_yscale('log')

labels_exp = [Line2D([],[],color='black', marker='*', label='ALICE (2022) 0-20$\%, |\eta|<0.8$', markersize=10)]
axes[0].legend(loc='best', handles=labels_exp)
axes[0].text(0.45, 0.75, 'Preliminary Data', transform=axes[0].transAxes, color='red')
axes[1].text(0.2,0.65,'Pb-Pb $\sqrt{s}=5.02$ ATeV, $|\eta|<0.8$', transform=axes[1].transAxes)


fig2, axes2 = plt.subplots(2,1,gridspec_kw={'height_ratios':(3,1)},
                                          sharex=True, sharey=False, figsize=(16,9))

ax = axes2[0]
util.plot_expr_data_on_axis(ax, alice, '*','black','gray', s=40, factor=1) 
tmp = results['00-20']
for item in tmp :
    col = 'black'
    if item in ['cujet', 'martini']:
        col = module_colors['MATTER+'+item.upper()]
    else:
        col = colors[2]
    y = tmp[item]
    func = interp1d(pT, np.log(y), kind='linear', fill_value='extrapolate')
    #if item == 'no jet-medium':
    #    print(f"For {item}: ")
    #    for yvalue, xvalue in zip(alice['y'], pT_alice):
    #        print(xvalue,"-->",yvalue/np.exp(func(xvalue)))
    #    exit(0)
    yy = np.exp(func(pT_alice))
    ax.plot(pT_alice, yy, color=col)
    r = yy/alice['y']
    dr_stat_plus  = r * alice['dy_stat+']/alice['y']
    dr_stat_minus = r* alice['dy_stat-']/alice['y']
    dr_syst_plus  = r*alice['dy_syst+']/alice['y']
    dr_syst_minus = r*alice['dy_syst-']/alice['y']
    ratio = DataFrame({'x':pT_alice, 'xhigh':alice['xhigh'], 'xlow':alice['xlow'],
                        'y':r, 'dy_stat+':dr_stat_plus, 'dy_stat-':dr_stat_minus,
                        'dy_syst+':dr_syst_plus, 'dy_syst-':dr_syst_minus})
    util.plot_expr_data_on_axis(axes2[1], ratio, '*',col,col, s=40, factor=1) 

labels2 = labels
labels2.append(labels_exp[0])
artist=axes2[0].legend(loc='best', handles=labels2)
axes2[1].set_xlabel(r'$p_T$ (GeV)', fontsize=30)
axes2[0].set_ylabel(r'$\frac{1}{2\pi p_T}\frac{\mathrm{d}^2N}{\mathrm{d}\eta \mathrm{d}p_T}$ (GeV${}^{-2}$)', fontsize=30)
axes2[1].set_ylabel('Theory over' + '\n' + 'Data', fontsize=30)
axes2[0].text(0.05, 0.1, 'Preliminary Data', transform=axes2[0].transAxes, color='red')
axes2[0].text(0.5,0.45,'Pb-Pb $\sqrt{s}=5.02$ ATeV, $|\eta|<0.8$', transform=axes2[0].transAxes)
ax.set_yscale('log')

fig3, axes3 = plt.subplots(2,2,sharex=True, sharey=True, figsize=(18,9))
for irow in range(2):
    cent = '00-20' if irow == 0 else '30-50'
    axes3[irow][0].text(0.05, 0.8, cent+r'$\%$', transform=axes3[irow][0].transAxes, fontsize=30)
    for icol in range(2):
        ratios = martini[cent] if icol == 0 else cujet[cent]
        ax = axes3[irow][icol]
        for ch in ratios:
            color = 'black'
            if ch in channel_colors:
                color = channel_colors[ch]
            else:
                color = module_colors['MATTER+MARTINI'] if icol == 0 else module_colors['MATTER+CUJET'] 
            ax.plot(pT, ratios[ch], color=color, linestyle=ddicts.channel_linestyles[ch])

pretty={'prompt':'Prompt', 'preEq':'Pre-Equil.', 'thermal':'Thermal',
        'jet-medium-cujet':'Jet-Med. (CUJET)', 'jet-medium-martini':'Jet-Med. (MARTINI)'}

labs = [Line2D([],[],color=module_colors['MATTER+'+eloss.upper()], label='MATTER+'+eloss.upper(), \
               linestyle=ddicts.channel_linestyles[f'jet-medium-{eloss}']) for eloss in ['cujet','martini'] ]

for channel in channel_colors:
    labs.append(Line2D([],[],color=channel_colors[channel], label=pretty[channel], linestyle=ddicts.channel_linestyles[channel]))

axes3[1][1].legend(handles=labs, loc='center right', fontsize=18)#, bbox_to_anchor=(1., 0.5))
for ax in [axes3[0][0], axes3[1][0]]:
    ax.set_ylabel('Channel/Total', fontsize=30)
for ax in axes3[1]:
    ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30)
axes3[0][1].text(0.05, 0.85, 'Pb-Pb $\sqrt{s}=5.02$ ATeV, $|\eta|<0.8$',transform=axes3[0][1].transAxes)


fig4, axes4 = plt.subplots(1,3, sharex=True, sharey=True, figsize=(18,9))

markers = {'martini':'P', 'cujet':'v'}
colours = {'conv-cujet'   : colors[0],
           'conv-martini' : colors[1],
           'brem-cujet'   : colors[2], 
           'brem-martini' : colors[3]}

for icent,  cent in enumerate(['00-10','10-20','30-50']):
    ax = axes4[icent]
    ax.set_xlabel(r'$p_T$ (GeV)', fontsize=30)
    ax.text(0.4, 0.2, cent+'$\%$', fontsize=30, transform=ax.transAxes)
    for eloss in ['martini', 'cujet']:
        tmp = data[cent][eloss]
        tmp = tmp[tmp['pT'].between(2., 20.)]
        ratio_conv = tmp['yconv']/tmp['N']
        ratio_brem = tmp['ybrem']/tmp['N']
        drat_conv = ratio_conv*np.sqrt(tmp['dyconv']**2/(tmp['yconv']**2)+ tmp['dN']**2/(tmp['N']**2))
        drat_brem = ratio_brem*np.sqrt(tmp['dyconv']**2/(tmp['yconv']**2)+ tmp['dN']**2/(tmp['N']**2))

        eboxes_conv = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                        for x, delx, y, yerrlow, yerrhigh in
                        zip(tmp['pT'], 0.5*tmp['dpT'], ratio_conv, drat_conv, drat_conv)]
        color_conv = f'conv-{eloss}'
        ax.scatter(tmp['pT'], ratio_conv, marker=markers[eloss], color=colours[color_conv], linewidths=0)
        pc = PatchCollection(eboxes_conv, facecolor=colours[color_conv], alpha=0.2)
        ax.add_collection(pc)

        color_brem = f'brem-{eloss}'
        eboxes_brem = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                        for x, delx, y, yerrlow, yerrhigh in
                        zip(tmp['pT'], 0.5*tmp['dpT'], ratio_brem, drat_brem, drat_brem)]
        ax.scatter(tmp['pT'], ratio_brem, marker=markers[eloss], color=colours[color_brem], linewidths=0)
        pc = PatchCollection(eboxes_brem, facecolor=colours[color_brem], alpha=0.2)
        ax.add_collection(pc)

labels_4 = [Line2D([],[],label=tag.split('-')[1].upper() + '-'+tag.split('-')[0].capitalize()+'.',color=colours[tag],marker=markers[tag.split('-')[1]]) for tag in colours]
axes4[-1].legend(loc='upper right', handles=labels_4)#,bbox_to_anchor=(1.01,0.8))
fig4.suptitle('Pb-Pb $\sqrt{s}=5.02$ ATeV, $|\eta|<0.8$', fontsize=30)
axes4[0].set_ylabel('Channel/[Jet-Medium Total]', fontsize=30)
plt.show()
