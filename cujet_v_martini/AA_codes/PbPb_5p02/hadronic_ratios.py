## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
#from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.interpolate import interp1d as interpolate
## my custom modules
import util
import dictionaries as my_hdicts
import jetDicts as my_jdicts
plt.rcParams.update(util.my_rcParams)
pTmax = 25
## Data:
data_centralities = ['00-05','05-10','10-20','30-40','40-50']
fAA_k_to_pi = '../../../expt/PbPb_5p02/Hadron_RAA/ALICE/PbPb_5p02_{c}_Kaon_to_Pi.csv'
fAA_p_to_pi = '../../../expt/PbPb_5p02/Hadron_RAA/ALICE/PbPb_5p02_{c}_Proton_to_Pi.csv'
fpp_k_to_pi = '../../../expt/PbPb_5p02/Hadron_RAA/ALICE/pp_5p02_Kaon_to_Pi.csv'
fpp_p_to_pi = '../../../expt/PbPb_5p02/Hadron_RAA/ALICE/pp_5p02_Proton_to_Pi.csv'

AA_data_kaon_to_pion = {cent: util.combine_systematics(pd.read_csv(
                                    fAA_k_to_pi.format(c=cent),comment='#')
                                    .rename(columns=my_hdicts.colnames_kaon_to_pion))
                                    for cent in data_centralities}
AA_data_proton_to_pion = {cent: util.combine_systematics(pd.read_csv(
                                fAA_p_to_pi.format(c=cent), comment='#')
                                .rename(columns=my_hdicts.colnames_proton_to_pion))
                                for cent in data_centralities}

pp_data_kaon_to_pion = pd.read_csv(fpp_k_to_pi,comment='#').rename(columns=my_hdicts.pp_colnames_kaon_to_pion) 
pp_data_proton_to_pion = pd.read_csv(fpp_p_to_pi, comment='#').rename(columns=my_hdicts.pp_colnames_proton_to_pion)

## JETSCAPE calculation:
direc = '../../jetscape_data/hadroChem/'
fname = direc + '{e}_PbPb5020_{c}_pi_k_proton_spec_0.5.csv'
jetscape_centralities = ['00-10','10-20','30-50']
pp = pd.read_csv(direc+'pp_5020_pi_k_proton_spec_0.5.csv',comment='#')
pp = pp[pp['pTmax'] < pTmax]
AA = {eloss: {cent: pd.read_csv(fname.format(e=eloss,c=cent,commen='#')) 
                for cent in jetscape_centralities} 
                for eloss in ['cujet','martini']}

## Plot the proton proton version:
fig1, ax1 = plt.subplots(1,1,figsize=(16,9),
                        gridspec_kw={'left':0.05,'right':0.95,'top':0.95,'bottom':0.05})

util.plot_expr_data_on_axis(ax1, pp_data_kaon_to_pion, marker='^')
util.plot_expr_data_on_axis(ax1, pp_data_proton_to_pion, marker='v')

## form the ratios using the calculation
pT = 0.5*(pp['pTmax'] + pp['pTmin'])
pi, dpi = pp['pions'], pp['dpions']
k, dk = pp['kaons'], pp['dkaons']
p, dp = pp['protons'], pp['dprotons']

r1 = k/pi
dr1 = r1*np.sqrt(dk*dk/(k*k) + dpi*dpi/(pi*pi))
r2 = p/pi
dr2 = r2*np.sqrt(dp*dpi/(p*p)+ dpi*dpi/(pi*pi))
ax1.plot(pT, r1, linestyle='solid', color='red')
ax1.fill_between(pT, r1-dr1, r1+dr1, color='red',alpha=0.2)
ax1.plot(pT, r2, linestyle='dashed', color='red')
ax1.fill_between(pT, r2-dr2, r2+dr2, color='red',alpha=0.2)



fig2, axes2 = plt.subplots(1,3,figsize=(16,9),
                        gridspec_kw={'left':0.05,'right':0.95,'top':0.95,'bottom':0.05}, sharex=True, sharey=True)

eloss_colors = my_hdicts.eloss_colours
for eloss, color in eloss_colors.items():
        if eloss == 'no jet-medium':
            continue
        curr_AA = AA[eloss]
        for icent, cent in enumerate(jetscape_centralities):
            ax = axes2[icent]
            if icent == 0:
                util.plot_expr_data_on_axis(ax, AA_data_kaon_to_pion['00-05'], marker='^')
                util.plot_expr_data_on_axis(ax, AA_data_kaon_to_pion['05-10'], marker='v')

                util.plot_expr_data_on_axis(ax, AA_data_proton_to_pion['00-05'], marker='s')
                util.plot_expr_data_on_axis(ax, AA_data_proton_to_pion['05-10'], marker='p')
            elif icent == 1:
                util.plot_expr_data_on_axis(ax, AA_data_kaon_to_pion['10-20'], marker='*')
                util.plot_expr_data_on_axis(ax, AA_data_proton_to_pion['10-20'], marker='o')
            else:
                util.plot_expr_data_on_axis(ax, AA_data_kaon_to_pion['30-40'], marker='<')
                util.plot_expr_data_on_axis(ax, AA_data_kaon_to_pion['40-50'], marker='>')
                util.plot_expr_data_on_axis(ax, AA_data_proton_to_pion['30-40'], marker='+')
                util.plot_expr_data_on_axis(ax, AA_data_proton_to_pion['40-50'], marker='x')

            spec = curr_AA[cent]
            spec = spec[spec['pTmax'] < pTmax]
            pT = 0.5*(spec['pTmax'] + spec['pTmin'])
            pi, dpi = spec['pions'], spec['dpions']
            k, dk = spec['kaons'], spec['dkaons']
            p, dp = spec['protons'], spec['dprotons']
            #pt1, pt2, pi, dpi = util.combine_bins(spec['pTmin'], spec['pTmax'], pi, dpi)
            #pt1, pt2, k, dk = util.combine_bins(spec['pTmin'], spec['pTmax'], k, dk)
            #pt1, pt2, p, dp = util.combine_bins(spec['pTmin'], spec['pTmax'], p, dp)
            #pt = 0.5*(pt1+pt2)
            r1 = k/pi
            dr1 = r1*np.sqrt(dk*dk/(k*k) + dpi*dpi/(pi*pi))
            r2 = p/pi
            dr2 = r2*np.sqrt(dp*dpi/(p*p)+ dpi*dpi/(pi*pi))
            ax.plot(pT, r1, linestyle='solid', color=color)
            ax.fill_between(pT, r1-dr1, r1+dr1, color=color,alpha=0.2)
            ax.plot(pT, r2, linestyle='dashed', color=color)
            ax.fill_between(pT, r2-dr2, r2+dr2, color=color,alpha=0.2)


plt.show()
