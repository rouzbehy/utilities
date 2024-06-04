import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d as interpolate
from photonreader import JetMediumReader
import dictionaries as my_dicts
import util
from COLORS import colors
from COLORS import rate_set_colors as rate_colors

plt.rcParams.update(util.my_rcParams)

## READ IN THE CALCULATION:
rates = {1:"LO", 2:"NLO", 3:"NP"}
cents = {"0_5": "00-05", "5_10": "05-10", "10_20": "10-20"}
mult = {'0_5':1615, '5_10':1268, '10_20':896.05}
template = "../../calcs/final_runs/rset_{r}/2p76/{c}/photons.csv"
diffCrossSection = 62.039 
reader = JetMediumReader(centralities=['0_5','5_10','10_20'], 
                rate_sets=[1,2,3], 
                fname_template=template, 
                oversampling=1000,
                mult=mult, xsec=diffCrossSection)
spectra = reader.get_spectra()

## ALICE data:
exp_loc = "../../../exp_data/photons/HEPData-ins1394677-v1-Table_1.csv"
exp_data = pd.read_csv(exp_loc, comment="#")


exp_data = exp_data[exp_data["x"] > 2]
deltax = 0.5 * (exp_data["xhigh"] - exp_data["xlow"])
ratenames = {1: "LO", 2: "NLO", 3: "NP"}

fig, axes = plt.subplots(
    nrows=2,
    ncols=1,
    height_ratios=[3, 1],
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

util.plot_expr_data_on_axis(axes[0], exp_data, 'o', color='black',face='dimgray', factor=1, s=100)

dexp_up = exp_data["dy_stat+"] ** 2 + exp_data["dy_syst+"] ** 2
dexp_down = exp_data["dy_stat-"] ** 2 + exp_data["dy_syst-"] ** 2
for rate in spectra:
    spec = spectra[rate]['0_20']
    colour = rate_colors[ratenames[rate]]
    xvalues = spec['pT']
    ax = axes[0]
    ax.plot(xvalues, spec['total'], color=colour)
    ax.fill_between(xvalues, spec['total'] + spec['dtotal'], spec['total'] - spec['dtotal'], color=colour, alpha=0.2)

    ax = axes[1]
    ## plot the ratios of the conversion and bremsstrahlung channels
    ax.plot(xvalues, spec['conv']/spec['total'], color=colour, linestyle='solid')
    ax.plot(xvalues, spec['brem']/spec['total'], color=colour, linestyle='dashed')

axes[0].set_yscale("log")
axes[0].set_ylabel(r"$\frac{1}{2\pi p_T}\frac{dN^{\gamma}}{dp_T d\eta}$ GeV$^{-2}$")
axes[1].set_xlabel(r"$p^{\gamma}_T$ (GeV)")
axes[1].set_ylabel("Ch. to Total")


## Legend work:
handles = [
    Line2D([], [], color=c, label=f"{l}")
    for l, c in rate_colors.items()
]
axes[0].legend(loc="upper right", handles=handles)
plt.show()
