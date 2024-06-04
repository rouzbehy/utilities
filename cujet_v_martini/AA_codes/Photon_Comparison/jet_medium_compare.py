#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
from dictionaries import multiplicity
from COLORS import module_colors, colors
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
    "axes.spines.top"  : False,
    "legend.frameon"   : False,
    "axes.labelsize" : 25}
plt.rcParams.update(my_rcParams)

oversampling = {'AuAu200':2000,
                'PbPb2p76': {'martini-00-05':1000000,
                             'martini-05-10':10000,
                             'cujet-00-05':5000,
                             'cujet-05-10':5000},
                'PbPb5020': 2000}

mult_PbPb_5p02 = {'00-10': 0.5*(1762+1380)}
mult_PbPb_2p76 = {'00-05':1615,'05-10':1268}
mult_AuAu_200 = {'00-10': (1053+831.4)/2}

loc_5p02 = '../../jetscape_data/sqrt_s_5020/maxt_200/{eloss}/PbPb5020_{cent}_photon_spec_0.80.csv'
loc_2p76 = '../../jetscape_data/sqrt_s_2760/{eloss}/PbPb_2760/PbPb2760_{cent}_photon_spec_0.80.csv'
loc_200  = '../../jetscape_data/sqrt_s_200/AuAu/{eloss}/00-10/photon_spec_0.35.csv'
PbPb_5p02_spec = {}
PbPb_2p76_spec = {}
AuAu_spec = {}
for eloss in ['martini','cujet']:
    PbPb_5p02_spec[eloss] = {}
    PbPb_2p76_spec[eloss] = {}
    AuAu_spec[eloss] = {}
    for cent in ['00-10']:
        tmp = pd.read_csv(loc_5p02.format(eloss=eloss,cent=cent), comment='#')
        tmp = tmp[tmp['ptmin'] > 1.5]
        tmp = tmp[tmp['ptmax'] < 30]
        pT = 0.5*(tmp['ptmax']+tmp['ptmin'])
        dpT = tmp['ptmax']-tmp['ptmin']
        y_conv  = mult_PbPb_5p02[cent]*tmp['conv']/(pT*dpT*2*0.8*2*np.pi*oversampling['PbPb5020'] )
        dy_conv = mult_PbPb_5p02[cent]*tmp['dconv']/(pT*dpT*2*0.8*2*np.pi*oversampling['PbPb5020'])
        y_brem  = mult_PbPb_5p02[cent]*tmp['brem']/(pT*dpT*2*0.8*2*np.pi*oversampling['PbPb5020'] )
        dy_brem = mult_PbPb_5p02[cent]*tmp['dbrem']/(pT*dpT*2*0.8*2*np.pi*oversampling['PbPb5020'])
        PbPb_5p02_spec[eloss][cent] = pd.DataFrame({'pT':pT, 'conv':y_conv,
                                               'dconv':dy_conv,
                                               'brem':y_brem,'dbrem':dy_brem})

    for cent in ['00-05','05-10']:
        tmp = pd.read_csv(loc_2p76.format(eloss=eloss,cent=cent), comment='#')
        tmp = tmp[tmp['ptmin'] > 1.5]
        tmp = tmp[tmp['ptmax'] < 30]
        pT = 0.5*(tmp['ptmax']+tmp['ptmin'])
        dpT = tmp['ptmax']-tmp['ptmin']
        y_conv  = mult_PbPb_2p76[cent]*tmp['conv']/(pT*dpT*2*0.8*2*np.pi*oversampling['PbPb2p76'][f'{eloss}-{cent}'])
        dy_conv = mult_PbPb_2p76[cent]*tmp['dconv']/(pT*dpT*2*0.8*2*np.pi*oversampling['PbPb2p76'][f'{eloss}-{cent}'])
        y_brem  = mult_PbPb_2p76[cent]*tmp['brem']/(pT*dpT*2*0.8*2*np.pi*oversampling['PbPb2p76'][f'{eloss}-{cent}'])
        dy_brem = mult_PbPb_2p76[cent]*tmp['dbrem']/(pT*dpT*2*0.8*2*np.pi*oversampling['PbPb2p76'][f'{eloss}-{cent}'])
        PbPb_2p76_spec[eloss][cent] = pd.DataFrame({'pT':pT, 'conv':y_conv,
                                               'dconv':dy_conv,
                                               'brem':y_brem,'dbrem':dy_brem})
    ## Pb-Pb
    spec = PbPb_2p76_spec[eloss]
    tmp_conv = 0.5*(spec['00-05']['conv'] + spec['05-10']['conv'])
    tmp_dconv = 0.5*np.sqrt(spec['00-05']['dconv']**2 + spec['05-10']['dconv']**2)
    tmp_brem = 0.5*(spec['00-05']['brem'] + spec['05-10']['brem'])
    tmp_dbrem = 0.5*np.sqrt(spec['00-05']['dbrem']**2 + spec['05-10']['dbrem']**2)
    PbPb_2p76_spec[eloss]['00-10'] = pd.DataFrame({'pT':spec['00-05']['pT'],
                                              'conv':tmp_conv,
                                              'dconv':tmp_dconv,
                                              'brem':tmp_brem,
                                              'dbrem':tmp_dbrem})
    ## Au-Au
    cent = '00-10'
    tmp = pd.read_csv(loc_200.format(eloss=eloss), comment='#')
    tmp = tmp[tmp['ptmin'] > 1.5]
    tmp = tmp[tmp['ptmax'] < 25]
    pT = 0.5*(tmp['ptmax']+tmp['ptmin'])
    dpT = tmp['ptmax']-tmp['ptmin']
    y_conv  = mult_AuAu_200['00-10']*tmp['conv'] /(pT*dpT*2*0.35*2*np.pi*oversampling['AuAu200'])
    dy_conv = mult_AuAu_200['00-10']*tmp['dconv']/(pT*dpT*2*0.35*2*np.pi*oversampling['AuAu200'])
    y_brem  = mult_AuAu_200['00-10']*tmp['brem'] /(pT*dpT*2*0.35*2*np.pi*oversampling['AuAu200'])
    dy_brem = mult_AuAu_200['00-10']*tmp['dbrem']/(pT*dpT*2*0.35*2*np.pi*oversampling['AuAu200'])
    AuAu_spec[eloss][cent] = pd.DataFrame({'pT':pT, 'conv':y_conv,
                                           'dconv':dy_conv,
                                           'brem':y_brem,'dbrem':dy_brem})


fig, axes = plt.subplots(2, 2, figsize=(16,9), sharex=True, height_ratios=(3,1))

axes = axes.flatten()
pT_PBPB = PbPb_2p76_spec['martini']['00-10']['pT']
pT_AUAU = AuAu_spec['martini']['00-10']['pT']
print("Pb-Pb @ 5.02 ATeV:")
spec1 = PbPb_5p02_spec['martini']['00-10']['conv']
spec2 = PbPb_5p02_spec['cujet']['00-10']['conv']
spec3 = PbPb_5p02_spec['martini']['00-10']['brem']
spec4 = PbPb_5p02_spec['cujet']['00-10']['brem']
print('PbPb:',f'{sum((spec1+spec3)/(spec2+spec4))/len(spec1):0.3f}')
print("Conversion Ratios: ", f'{sum(spec1/spec2)/len(spec1):0.3f}')
print("Brem. Ratios: ", f'{sum(spec3/spec4)/len(spec3):0.3f}')
print("Pb-Pb @ 2.76 ATeV:")
spec1 = PbPb_2p76_spec['martini']['00-10']['conv']
spec2 = PbPb_2p76_spec['cujet']['00-10']['conv']
spec3 = PbPb_2p76_spec['martini']['00-10']['brem']
spec4 = PbPb_2p76_spec['cujet']['00-10']['brem']
print('PbPb:',f'{sum((spec1+spec3)/(spec2+spec4))/len(spec1):0.3f}')
print("Conversion Ratios: ", f'{sum(spec1/spec2)/len(spec1):0.3f}')
print("Brem. Ratios: ", f'{sum(spec3/spec4)/len(spec3):0.3f}')
print("Au-Au @ 200 AGeV:")
spec1 = AuAu_spec['martini']['00-10']['conv']
spec2 = AuAu_spec['cujet']['00-10']['conv']
spec3 = AuAu_spec['martini']['00-10']['brem']
spec4 = AuAu_spec['cujet']['00-10']['brem']
print('AuAu:',f'{sum((spec1+spec3)/(spec2+spec4))/len(spec1):0.3f}')
print("Conversion Ratios: ", f'{sum(spec1/spec2)/len(spec1):0.3f}')
print("Brem. Ratios: ", f'{sum(spec3/spec4)/len(spec3):0.3f}')


for eloss in ['cujet','martini']:
    lstyle = 'solid'
    color = module_colors['MATTER+'+eloss.upper()]
    ax = axes[0]
    spec = PbPb_5p02_spec[eloss]['00-10']['brem']
    ax.plot(pT_PBPB, spec, color=color, linestyle=lstyle)
    lstyle = 'dashed'
    spec = PbPb_2p76_spec[eloss]['00-10']['brem']
    ax.plot(pT_PBPB, spec, color=color, linestyle=lstyle)
    lstyle='dotted'
    spec = AuAu_spec[eloss]['00-10']['brem']
    ax.plot(pT_AUAU, spec, color=color, linestyle=lstyle)

    ax = axes[1]
    lstyle = 'solid'
    spec = PbPb_5p02_spec[eloss]['00-10']['conv']
    ax.plot(pT_PBPB, spec, color=color, linestyle=lstyle)
    lstyle = 'dashed'
    spec = PbPb_2p76_spec[eloss]['00-10']['conv']
    ax.plot(pT_PBPB, spec, color=color, linestyle=lstyle)
    lstyle='dotted'
    spec = AuAu_spec[eloss]['00-10']['conv']
    ax.plot(pT_AUAU, spec, color=color, linestyle=lstyle)

cent = '00-10'
color= 'black'
for iax, channel in zip([2,3],['brem','conv']):
    lstyle = 'solid'
    ax = axes[iax]
    spec_martini = PbPb_5p02_spec['martini'][cent][channel]
    spec_cujet   = PbPb_5p02_spec['cujet'][cent][channel]

    ratio = spec_martini/spec_cujet
    ax.plot(pT_PBPB, ratio, color=color, linestyle=lstyle)

    lstyle = 'dashed'
    spec_martini = PbPb_2p76_spec['martini'][cent][channel]
    spec_cujet   = PbPb_2p76_spec['cujet'][cent][channel]

    ratio = spec_martini/spec_cujet
    ax.plot(pT_PBPB, ratio, color=color, linestyle=lstyle)

    lstyle='dotted'
    spec_martini = AuAu_spec['martini'][cent][channel]
    spec_cujet   = AuAu_spec['cujet'][cent][channel]

    ratio = spec_martini/spec_cujet
    ax.plot(pT_AUAU, ratio, color=color, linestyle=lstyle)


#axes[0].set_yscale('log')
for ax in axes:
    ax.set_xlabel(r'$p_T$ (GeV)')
axes[0].set_ylabel(r'$\frac{1}{2\pi p_T}\frac{d\sigma^{\gamma}}{dp_Td\eta}$ (GeV)${}^{-2}$', fontsize=35)

for iax, ax in enumerate(axes):
    tag = 'Bremsstrahlung' if iax == 0 else 'Conversion'
    ax.text(0.4, 0.8, tag, transform=ax.transAxes)
labels = [Line2D([],[],color=c, label=l) for l, c in module_colors.items() if l != 'MATTER']
for name, lstyle in zip(['AuAu-200 AGeV, $|\eta|<0.35$','PbPb-2.76 ATeV, $|\eta|<0.8$', 'PbPb-5.02 ATeV, $|\eta|<0.8$ [Prelim.]'],['dotted','dashed','solid']):
    labels.append(Line2D([],[],color='black',linestyle=lstyle, label=name))
axes[1].legend(loc='lower left', handles=labels, fontsize=28)
axes[0].text(0.1, 0.1, r'Centrality: 0-10$\%$', transform=axes[0].transAxes)
plt.show()
