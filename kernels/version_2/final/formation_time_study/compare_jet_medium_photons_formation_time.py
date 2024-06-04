"""
    Comparison of jet medium photons from:
        * LO, NLO, NP runs with formation time
        * LO, NLO, NP runs without formation time
        * MARTINI + MATTER
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d as interpolate
from photonreader import JetMediumReader, JETSCAPEReader
import dictionaries as my_dicts
import util
from COLORS import colors_comparison as colors
from COLORS import rate_set_colors as rate_colors

plt.rcParams.update(util.my_rcParams)

## READ IN THE CALCULATION:
rates = {1:"LO", 2:"NLO", 3:"NP"}
cents = {"0_5": "00-05", "5_10": "05-10", "10_20": "10-20"}
mult = {'0_5':1615, '5_10':1268, '10_20':896.05}
diffCrossSection = 62.039 
template = "../../calcs/final_runs/rset_{r}/2p76/{c}/photons.csv"

reader_wformtime = JetMediumReader(centralities=['0_5','5_10','10_20'], 
                rate_sets=[1,2,3], 
                fname_template=template, 
                oversampling=1000,
                mult=mult, xsec=diffCrossSection)
spectra_wformtime = reader_wformtime.get_spectra()

template = "../../../v2/production_v2/martini_results/"
template += "final_PbPb_2p76/rset_{r}/cent_{c}/gamma_spectra.csv"

reader_nformtime = JetMediumReader(centralities=['0_5','5_10','10_20'], 
                rate_sets=[1,2,3], 
                fname_template=template, 
                oversampling=1000,
                mult=mult, xsec=diffCrossSection)

spectra_nformtime = reader_nformtime.get_spectra()


template = "../../../../jetscape_project/v2/jetscape_data/"
template += "sqrt_s_2760/{model}/PbPb_2760/PbPb2760_{c}_photon_spec_0.80.csv"
oversampling_factors = [1000000,10000,10000] ## in order: 0-5, 5-10, 10-20
mult = {'00-05':1615, '05-10':1268, '10-20':896.05}
reader_jetscape = JETSCAPEReader(xsec=diffCrossSection, centralities=['00-05','05-10','10-20'],
                                 oversampling=oversampling_factors, mult=mult,
                                 modelname='martini', fname_template=template)
spectra_jetscape = reader_jetscape.get_spectra()

ratenames = {1: "LO", 2: "NLO", 3: "NP"}

fig, ax = plt.subplots(
    nrows=1,
    ncols=1,
    sharex=True,
    figsize=(16, 9),
    gridspec_kw={
        "left": 0.12,
        "bottom": 0.12,
        "right": 0.995,
        "top": 0.995,
        "hspace": 0.038,
    },
)

for rate in spectra_wformtime.keys():
    spec = spectra_wformtime[rate]['0_20']
    colour = rate_colors[ratenames[rate]]
    xvalues = spec['pT']
    ax.plot(xvalues, spec['total'], color=colour)
    ax.fill_between(xvalues, spec['total'] + spec['dtotal'], spec['total'] - spec['dtotal'], color=colour, alpha=0.2)

    # ax = axes[1]
    # ## plot the ratios of the conversion and bremsstrahlung channels
    # ax.plot(xvalues, spec['conv']/spec['total'], color=colour, linestyle='solid')
    # ax.plot(xvalues, spec['brem']/spec['total'], color=colour, linestyle='dashed')

for rate in spectra_nformtime.keys():
    spec = spectra_nformtime[rate]['0_20']
    colour = rate_colors[ratenames[rate]]
    xvalues = spec['pT']
    ax.plot(xvalues, spec['total'], color=colour, linestyle='dotted')
    ax.fill_between(xvalues, spec['total'] + spec['dtotal'], spec['total'] - spec['dtotal'], color=colour, alpha=0.2)

## Legend work:
handles = [
    Line2D([], [], color=c, label=f"{l}")
    for l, c in rate_colors.items()
]
handles.append(Line2D([], [], color='black', label="delayed shower"))
handles.append(Line2D([], [], color='black', label="instantaneous shower", linestyle='dotted'))

ax.text(0.1, 0.1, 'Pb-Pb @ 2.76 ATeV, $0$-$20\%$', transform=ax.transAxes)

ax.set_ylabel(r"$\frac{1}{2\pi p_T}\frac{dN^{\gamma}}{dp_T d\eta}$ GeV$^{-2}$")
ax.set_xlabel(r"$p^{\gamma}_T$ (GeV)")
ax.set_yscale("log")
ax.legend(loc="upper right", handles=handles)


plt.show()