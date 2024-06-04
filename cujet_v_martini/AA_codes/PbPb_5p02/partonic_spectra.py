## import scipy and matplotlib modules:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
## my custom modules
import util
import dictionaries as my_hdicts
import jetDicts as my_jdicts
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

dirctory = '../../jetscape_data/hadroChem/'
fname = '{e}_PbPb5020_{c}_parton_spec_2.0.csv'
jetscape_centralities = ['00-10','10-20','30-50']
pTmin, pTmax = 1, 2000


parton_spectra_AA = {eloss : {cent: pd.read_csv(dirctory+fname.format(e=eloss,c=cent),comment='#')
                        for cent in jetscape_centralities} for eloss in ['martini','cujet']}
parton_spectra_pp = pd.read_csv('../../jetscape_data/hadroChem/pp_5020_parton_spec_2.0.csv', comment='#')
pT = 0.5*(parton_spectra_pp['pTmax'] + parton_spectra_pp['pTmin'])
parton_spectra_pp['pT'] = pT.to_list()
parton_spectra_pp = parton_spectra_pp[parton_spectra_pp['pT'].between(pTmin, pTmax)]


fig, axes = plt.subplots(1,3,figsize=(16,9),
                        gridspec_kw={'left':0.05,'right':0.95,
                                     'top':0.95,'bottom':0.1}, sharex=True,sharey=True)
q  = parton_spectra_pp['qqbar']
dq = parton_spectra_pp['dqqbar']
g  = parton_spectra_pp['glue']
dg = parton_spectra_pp['dglue']

pT = parton_spectra_pp['pT']
r1  = q/g
dr1 = r1*np.sqrt(dq*dq/(q*q) + dg*dg/(g*g))
colours = {'matter+cujet':"#7cae7a", 'matter+martini':"#d64045", 'matter(vac.)':'#eca400'}
for ax in axes:
    ax.plot(pT, r1, linestyle='solid', color='#eca400')
    ax.fill_between(pT, r1-dr1, r1+dr1, color='#eca400',alpha=0.2)

eloss_colors = my_hdicts.eloss_colours
for eloss, color in eloss_colors.items():
    if 'no' in eloss:
        continue
    for icent, cent in enumerate(jetscape_centralities):
        ax = axes[icent]
        ax.text(0.1, 0.8,s=f'{cent}'+r'$\%$', transform=ax.transAxes)
        spec = parton_spectra_AA[eloss][cent]
        #spec = spec[spec['qqbar']>0]
        pT = 0.5*(spec['pTmax'] + spec['pTmin'])
        pT = pT.to_list()
        spec['pT'] = pT
        spec = spec[spec['pT'].between(pTmin, pTmax)]
        q = spec['qqbar']
        dq = spec['dqqbar']
        g = spec['glue']
        dg = spec['dglue']

        r1  = q/g
        dr1 = r1*np.sqrt(dq*dq/(q*q) + dg*dg/(g*g))
        ax.plot(spec['pT'], r1, linestyle='solid', color=color)
        ax.fill_between(spec['pT'], r1-dr1, r1+dr1, color=color,alpha=0.2)

label_set_1 = [Line2D([],[],label=eloss.upper(), color=col) for eloss, col in colours.items()]
#artist = ax.legend(handles=label_set_1, bbox_to_anchor=(1.,0.7), fontsize=20)
#ax.add_artist(artist)
ax.legend(loc='best',fontsize=18, handles=label_set_1)

axes[0].set_yscale('log')
axes[0].set_xscale('log')
fig.supxlabel(r'$p_T$ (GeV)')
fig.supylabel(r'Ratio $q\bar{q}/g$')
plt.show()

