#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

## my custom modules
# import util
# import jetDicts as ddicts
# from COLORS import colors_comparison as colors
import sys

sys.path.insert(0, "..")
from helpers import util
from helpers import dictionaries as my_dicts
from helpers import jetDicts as ddicts
from helpers import COLORS

colors = COLORS.colors_comparison_two
PARENT_DIR = "/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2"
plt.rcParams.update(util.my_rcParams)
rate_names = {1: "LO", 2: "NLO", 3: "NP"}
all_marker_sizes = 100


def read_martini_jet_shape(fname):
    tmp = pd.read_csv(fname, comment="#")
    njet = 1
    with open(fname, "r") as f:
        line = f.readline()
        line = line.split(" ")[-1]
        njet = float(line)
    dat = tmp[tmp["rmax"] < 0.31]
    delta_r = tmp["rmax"] - tmp["rmin"]
    r = 0.5 * (tmp["rmax"] + tmp["rmin"])
    rho = tmp["N"] / (njet)
    drho = tmp["dN"] / (njet)
    norm = sum(rho.to_list())  #
    rho_normed = rho / (delta_r * norm)
    drho_normed = drho / (delta_r * norm)
    return rho_normed, drho_normed, r, delta_r


# pT_lower_lim, pT_upper_lim, inel_Xsec = 20, 96, 62.03948
def get_with_formation_time():
    ## Get the Proton-Proton Baseline for the case with formation time:
    path = f"{PARENT_DIR}/v3/calcs/pp_2p76/jet_shape.csv"
    calcs = pd.read_csv(path, comment="#")
    njet = 1
    with open(path, "r") as f:
        line = f.readline()
        line = line.split(" ")[-1]
        njet = float(line)

    dat = calcs[calcs["rmax"] < 0.31]
    delta_r = dat["rmax"] - dat["rmin"]
    r = 0.5 * (dat["rmax"] + dat["rmin"])
    rho = dat["N"] / (njet)
    drho = dat["dN"] / (njet)
    norm = sum(rho.to_list())  #
    rho_normed = rho / (delta_r * norm)
    drho_normed = drho / (delta_r * norm)

    pp_data = pd.DataFrame(
        {"r": r, "dr": delta_r, "rho_normed": rho_normed, "drho_normed": drho_normed}
    )

    ## Compute the jet shape ratio for each rate set and the two
    ## centralities of interest:
    result = {}
    fname = PARENT_DIR + "/v3/calcs/final_runs/rset_1/2p76/{c}/jet_shape.csv"
    rho_normed1, drho_normed1, r, delta_r = read_martini_jet_shape(
        fname.format(c="0_5")
    )
    rho_normed2, drho_normed2, r, delta_r = read_martini_jet_shape(
        fname.format(c="5_10")
    )
    rho_normed = 0.5 * (rho_normed1 + rho_normed2)
    scale = sum(rho_normed)
    rho_normed = rho_normed / (delta_r * scale)
    drho_normed = 0.5 * np.sqrt(drho_normed1**2 + drho_normed2**2)
    drho_normed = drho_normed / (delta_r * scale)
    ratio = rho_normed / pp_data["rho_normed"]
    dr = ratio * np.sqrt(
        (drho_normed / rho_normed) ** 2
        + (pp_data["drho_normed"] / pp_data["rho_normed"]) ** 2
    )
    result["0_10"] = (r, delta_r, ratio, dr)

    rho_normed, drho_normed, r, delta_r = read_martini_jet_shape(
        fname.format(c="10_20")
    )
    scale = sum(rho_normed)
    rho_normed = rho_normed / (delta_r * scale)
    drho_normed = 0.5 * np.sqrt(drho_normed1**2 + drho_normed2**2)
    drho_normed = drho_normed / (delta_r * scale)
    ratio = rho_normed / pp_data["rho_normed"]
    dr = ratio * np.sqrt(
        (drho_normed / rho_normed) ** 2
        + (pp_data["drho_normed"] / pp_data["rho_normed"]) ** 2
    )
    result["10_20"] = (r, delta_r, ratio, dr)
    return result


def get_without_formation_time():
    ## Get the Proton-Proton Baseline for the case with formation time:
    path = f"{PARENT_DIR}/v2/proton_proton/data/jet_shape.csv"
    calcs = pd.read_csv(path, comment="#")
    njet = 1
    with open(path, "r") as f:
        line = f.readline()
        line = line.split(" ")[-1]
        njet = float(line)

    dat = calcs[calcs["rmax"] < 0.31]
    delta_r = dat["rmax"] - dat["rmin"]
    r = 0.5 * (dat["rmax"] + dat["rmin"])
    rho = dat["N"] / (njet)
    drho = dat["dN"] / (njet)
    norm = sum(rho.to_list())  #
    rho_normed = rho / (delta_r * norm)
    drho_normed = drho / (delta_r * norm)

    pp_data = pd.DataFrame(
        {"r": r, "dr": delta_r, "rho_normed": rho_normed, "drho_normed": drho_normed}
    )

    ## Compute the jet shape ratio for each rate set and the two
    ## centralities of interest:
    result = {}
    fname = f"{PARENT_DIR}/v2/production_v2/martini_results/final_PbPb_2p76"
    fname += "/rset_1/cent_{c}/jet_shape.csv"
    rho_normed1, drho_normed1, r, delta_r = read_martini_jet_shape(
        fname.format(c="0_5")
    )
    rho_normed2, drho_normed2, r, delta_r = read_martini_jet_shape(
        fname.format(c="5_10")
    )
    rho_normed = 0.5 * (rho_normed1 + rho_normed2)
    scale = sum(rho_normed)
    rho_normed = rho_normed / (delta_r * scale)
    drho_normed = 0.5 * np.sqrt(drho_normed1**2 + drho_normed2**2)
    drho_normed = drho_normed / (delta_r * scale)
    ratio = rho_normed / pp_data["rho_normed"]
    dr = ratio * np.sqrt(
        (drho_normed / rho_normed) ** 2
        + (pp_data["drho_normed"] / pp_data["rho_normed"]) ** 2
    )
    result["0_10"] = (r, delta_r, ratio, dr)

    rho_normed, drho_normed, r, delta_r = read_martini_jet_shape(
        fname.format(c="10_20")
    )
    scale = sum(rho_normed)
    rho_normed = rho_normed / (delta_r * scale)
    drho_normed = 0.5 * np.sqrt(drho_normed1**2 + drho_normed2**2)
    drho_normed = drho_normed / (delta_r * scale)
    ratio = rho_normed / pp_data["rho_normed"]
    dr = ratio * np.sqrt(
        (drho_normed / rho_normed) ** 2
        + (pp_data["drho_normed"] / pp_data["rho_normed"]) ** 2
    )
    result["10_20"] = (r, delta_r, ratio, dr)
    return result


## read in the experimental results
data_fname = (
    PARENT_DIR
    + "/exp_data/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_{c}_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
)
cms_data = {}
for c in ["00-10", "10-30"]:
    tmp = np.loadtxt(data_fname.format(c=c), comments="#", unpack=True, delimiter="\t")
    cms_data[c] = pd.DataFrame({"x": tmp[0], "y": tmp[1], "dx": tmp[2], "dy": tmp[3]})


result_nform = get_without_formation_time()
result_wform = get_with_formation_time()

## Plot the data:
fig, axes = plt.subplots(
    1,
    2,
    gridspec_kw={
        "top": 0.985,
        "bottom": 0.1,
        "left": 0.08,
        "right": 0.985,
        "hspace": 0.02,
        "wspace": 0.1,
    },
    figsize=(16, 9),
    sharex=True,
    sharey=True,
)

for c, ax in zip(cms_data, axes):
    data = cms_data[c]
    errorboxes = [
        Rectangle(
            (x - delx, y - yerrlow), width=2 * delx, height=abs(yerrlow) + abs(yerrhigh)
        )
        for x, delx, y, yerrlow, yerrhigh in zip(
            data["x"], data["dx"], data["y"], data["dy"], data["dy"]
        )
    ]
    pc = PatchCollection(errorboxes, facecolor="black", edgecolor="black", alpha=0.4)
    ax.add_collection(pc)
    scatter = ax.scatter(
        data["x"],
        data["y"],
        color="black",
        marker="s",
        s=all_marker_sizes,
        label="CMS (2014)",
    )

for color, data in zip([colors[0], colors[1]], [result_nform, result_wform]):
    for iax, c in enumerate(["0_10", "10_20"]):
        ax = axes[iax]
        clo, chi = c.split("_")
        txt: str = f"${clo}$-${chi}$\%"
        ax.text(0.4, 0.4, s=txt, transform=ax.transAxes)
        r, dr, ratio, dratio = data[c]
        df = pd.DataFrame({"x": r, "dx": dr / 2, "y": ratio, "dy": dratio})
        util.plot_theory_on_axis_boxes(ax, df, "y", "x", color, "s", s=all_marker_sizes)
        # ax.plot(r, ratio, color=color)
        # ax.fill_between(r, ratio - dratio, ratio + dratio, alpha=0.2, color=color)


line_colors = {
    r"Without $\tau_f$": colors[0],
    r"With $\tau_f$": colors[1],
}

ax = axes[0]

ax.text(0.02, 0.9, s=r"Pb-Pb @ 2.76 ATeV", transform=ax.transAxes)
ax.text(
    0.02,
    0.75,
    s=r"$0.3<|\eta|<2.0$" + "\n" + r"$100$ GeV $< p^{\mathrm{jet}}_T$",
    transform=ax.transAxes,
)
ax.text(0.02, 0.65, s=r"$R=0.3$, Anti-$k_{T}$", transform=ax.transAxes)
for ax in axes:
    ax.set_xlabel(r"$r$")

labels = [Line2D([], [], c=c, label=l) for l, c in line_colors.items()]
labels.append(
    Line2D([], [], label="CMS", marker="s", color="black", markersize=10, linewidth=0)
)
axes[1].legend(loc="upper left", handles=labels)  # , handletextpad=0.05)
axes[1].text(0.55, 0.05, r"CMS: $10$-$30$\%", transform=axes[1].transAxes)
axes[0].set_ylabel(r"$R_{\mathrm{\rho}}$")
# plt.savefig("../../plots/formation_time_effect/jet_shape_ratio.png", dpi=200)
plt.show()
