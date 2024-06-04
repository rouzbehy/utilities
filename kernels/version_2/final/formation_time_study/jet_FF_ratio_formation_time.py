#!/usr/bin/env python3
"""
    Plot the computed spectra for jets
"""
## import scipy and matplotlib modules:
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

plt.rcParams.update(util.my_rcParams)
minpT, maxpT = 0.1, 120
zmin, zmax = 0.005, 2

PARENT_DIR = "/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2"


def get_martini_results_wformtime(rateset, cent):
    fname_pt = (
        f"{PARENT_DIR}/v3/calcs/final_runs/rset_{rateset}/2p76/{cent}/jet_FF_pT.csv"
    )
    fname_z = (
        f"{PARENT_DIR}/v3/calcs/final_runs/rset_{rateset}/2p76/{cent}/jet_FF_z.csv"
    )
    tmp_DPT = pd.read_csv(fname_pt, comment="#")
    tmp_DZ = pd.read_csv(fname_z, comment="#")
    tmp_DPT = tmp_DPT[tmp_DPT["pTmax"] < maxpT]
    tmp_DPT = tmp_DPT[tmp_DPT["pTmin"] > minpT]
    tmp_DZ = tmp_DZ[tmp_DZ["zmax"] < zmax]
    tmp_DZ = tmp_DZ[tmp_DZ["zmin"] > zmin]
    njets_pp, nevts_pp = 0, 0
    with open(fname_pt.format(r=rateset), "r") as f:
        line = f.readline()
        line = line.split(" ")
        njets_pp = float(line[1])
        nevts_pp = float(line[-1])

    dpT = tmp_DPT["pTmax"] - tmp_DPT["pTmin"]
    pT = 0.5 * (tmp_DPT["pTmax"] + tmp_DPT["pTmin"])
    scale = njets_pp
    tmpy = tmp_DPT["N"] / (dpT * scale)
    dtmpy = tmp_DPT["dN"] / (dpT * scale)

    dz = tmp_DZ["zmax"] - tmp_DZ["zmin"]
    z = 0.5 * (tmp_DZ["zmax"] + tmp_DZ["zmin"])
    tmpyz = tmp_DZ["N"] / (dz * scale)
    dtmpyz = tmp_DZ["dN"] / (dz * scale)
    result_pT = pd.DataFrame({"pT": pT, "dpT": dpT, "FF(pT)": tmpy, "dFF(pT)": dtmpy})
    result_z = pd.DataFrame({"z": z, "dz": dz, "FF(z)": tmpyz, "dFF(z)": dtmpyz})
    return result_pT, result_z


def get_martini_results_nformtime(rateset, cent):
    fname_pt = f"{PARENT_DIR}/v2/production_v2/martini_results/final_PbPb_2p76/rset_{rateset}/cent_{cent}/jet_FF_pT.csv"
    fname_z = f"{PARENT_DIR}/v2/production_v2/martini_results/final_PbPb_2p76/rset_{rateset}/cent_{cent}/jet_FF_z.csv"
    tmp_DPT = pd.read_csv(fname_pt, comment="#")
    tmp_DZ = pd.read_csv(fname_z, comment="#")
    tmp_DPT = tmp_DPT[tmp_DPT["pTmax"] < maxpT]
    tmp_DPT = tmp_DPT[tmp_DPT["pTmin"] > minpT]
    tmp_DZ = tmp_DZ[tmp_DZ["zmax"] < zmax]
    tmp_DZ = tmp_DZ[tmp_DZ["zmin"] > zmin]
    njets_pp, nevts_pp = 0, 0
    with open(fname_pt.format(r=rateset), "r") as f:
        line = f.readline()
        line = line.split(" ")
        njets_pp = float(line[1])
        nevts_pp = float(line[-1])

    dpT = tmp_DPT["pTmax"] - tmp_DPT["pTmin"]
    pT = 0.5 * (tmp_DPT["pTmax"] + tmp_DPT["pTmin"])
    scale = njets_pp
    tmpy = tmp_DPT["N"] / (dpT * scale)
    dtmpy = tmp_DPT["dN"] / (dpT * scale)

    dz = tmp_DZ["zmax"] - tmp_DZ["zmin"]
    z = 0.5 * (tmp_DZ["zmax"] + tmp_DZ["zmin"])
    tmpyz = tmp_DZ["N"] / (dz * scale)
    dtmpyz = tmp_DZ["dN"] / (dz * scale)
    result_pT = pd.DataFrame({"pT": pT, "dpT": dpT, "FF(pT)": tmpy, "dFF(pT)": dtmpy})
    result_z = pd.DataFrame({"z": z, "dz": dz, "FF(z)": tmpyz, "dFF(z)": dtmpyz})
    return result_pT, result_z


def get_pp_wformtime():
    ## Read in the pp data:
    pp_DPT = pd.read_csv(f"{PARENT_DIR}/v3/calcs/pp_2p76/jet_FF_pT.csv", comment="#")
    pp_DZ = pd.read_csv(f"{PARENT_DIR}/v3/calcs/pp_2p76/jet_FF_z.csv", comment="#")

    pp_DPT = pp_DPT[pp_DPT["pTmax"] < maxpT]
    pp_DPT = pp_DPT[pp_DPT["pTmin"] > minpT]
    pp_DZ = pp_DZ[pp_DZ["zmax"] < zmax]
    pp_DZ = pp_DZ[pp_DZ["zmin"] > zmin]
    njets_pp, nevts_pp = 0, 0
    with open(f"{PARENT_DIR}/v3/calcs/pp_2p76/jet_FF_pT.csv", "r") as f:
        line = f.readline()
        line = line.split(" ")
        njets_pp = float(line[1])
        nevts_pp = float(line[-1])

    dpT = pp_DPT["pTmax"] - pp_DPT["pTmin"]
    pT = 0.5 * (pp_DPT["pTmax"] + pp_DPT["pTmin"])
    scale = njets_pp
    y = pp_DPT["N"] / (dpT * scale)
    dy = pp_DPT["dN"] / (dpT * scale)
    pp_dpt = pd.DataFrame({"pT": pT, "dpT": dpT, "FF(pT)": y, "dFF(pT)": dy})

    dz = pp_DZ["zmax"] - pp_DZ["zmin"]
    z = 0.5 * (pp_DZ["zmax"] + pp_DZ["zmin"])
    y = pp_DZ["N"] / (dz * scale)
    dy = pp_DZ["dN"] / (dz * scale)

    pp_dz = pd.DataFrame({"z": z, "dz": dz, "FF(z)": y, "dFF(z)": dy})

    return pp_dpt, pp_dz


def get_pp_nformtime():
    ## Read in the pp data:
    ptloc = f"{PARENT_DIR}/v2/production_v2/martini_results/pp/jet_FF_pT.csv"
    pzloc = f"{PARENT_DIR}/v2/production_v2/martini_results/pp/jet_FF_z.csv"
    pp_DPT = pd.read_csv(ptloc, comment="#")
    pp_DZ = pd.read_csv(pzloc, comment="#")

    pp_DPT = pp_DPT[pp_DPT["pTmax"] < maxpT]
    pp_DPT = pp_DPT[pp_DPT["pTmin"] > minpT]
    pp_DZ = pp_DZ[pp_DZ["zmax"] < zmax]
    pp_DZ = pp_DZ[pp_DZ["zmin"] > zmin]
    njets_pp, nevts_pp = 0, 0
    with open(ptloc, "r") as f:
        line = f.readline()
        line = line.split(" ")
        njets_pp = float(line[1])
        nevts_pp = float(line[-1])

    dpT = pp_DPT["pTmax"] - pp_DPT["pTmin"]
    pT = 0.5 * (pp_DPT["pTmax"] + pp_DPT["pTmin"])
    scale = njets_pp
    y = pp_DPT["N"] / (dpT * scale)
    dy = pp_DPT["dN"] / (dpT * scale)
    pp_dpt = pd.DataFrame({"pT": pT, "dpT": dpT, "FF(pT)": y, "dFF(pT)": dy})

    dz = pp_DZ["zmax"] - pp_DZ["zmin"]
    z = 0.5 * (pp_DZ["zmax"] + pp_DZ["zmin"])
    y = pp_DZ["N"] / (dz * scale)
    dy = pp_DZ["dN"] / (dz * scale)
    pp_dz = pd.DataFrame({"z": z, "dz": dz, "FF(z)": y, "dFF(z)": dy})

    return pp_dpt, pp_dz


def combine_centralities(c1, c2, var="pT"):
    label, dlabel = f"FF({var})", f"dFF({var})"
    tmpy = 0.5 * (c1[label] + c2[label])
    dtmpy = 0.5 * np.sqrt(c1[dlabel] ** 2 + c2[dlabel] ** 2)
    result = pd.DataFrame(
        data={var: c1[var], f"d{var}": c1[f"{var}"], label: tmpy, dlabel: dtmpy}
    )
    return result


"""
    ***   ACTUAL PLOTTING ***
"""

AA_wform_00_05_pT, AA_wform_00_05_z = get_martini_results_wformtime(1, "0_5")
AA_wform_05_10_pT, AA_wform_05_10_z = get_martini_results_wformtime(1, "5_10")

AA_nform_00_05_pT, AA_nform_00_05_z = get_martini_results_nformtime(1, "0_5")
AA_nform_05_10_pT, AA_nform_05_10_z = get_martini_results_nformtime(1, "5_10")

pp_wform_pT, pp_wform_z = get_pp_wformtime()
pp_nform_pT, pp_nform_z = get_pp_nformtime()

AA_wform_00_10_pT = combine_centralities(AA_wform_00_05_pT, AA_wform_05_10_pT, "pT")
AA_wform_00_10_z = combine_centralities(AA_wform_00_05_z, AA_wform_05_10_z, "z")

AA_nform_00_10_pT = combine_centralities(AA_nform_00_05_pT, AA_nform_05_10_pT, "pT")
AA_nform_00_10_z = combine_centralities(AA_nform_00_05_z, AA_nform_05_10_z, "z")


exp_loc = f"{PARENT_DIR}/exp_data/jet_fragmentation/ATLAS/"
## read in ATLAS results
data_pT = pd.read_csv(
    exp_loc + "HEPData-ins1511869-v1-csv/Table9.csv", comment="#"
).rename(columns=ddicts.colnames_ATLAS_FF_DPT)
data_z = pd.read_csv(
    exp_loc + "HEPData-ins1511869-v1-csv/Table25.csv", comment="#"
).rename(columns=ddicts.colnames_ATLAS_FF_DZ)
## make the plots:
fig, axes = plt.subplots(
    1,
    2,
    figsize=(16, 9),
    sharex=False,
    sharey=True,
    gridspec_kw={
        "left": 0.086,
        "bottom": 0.117,
        "right": 0.986,
        "top": 0.986,
        "wspace": 0.093,
    },
)


AA_pT_w, AA_dpT_w = AA_wform_00_10_pT["FF(pT)"], AA_wform_00_10_pT["dFF(pT)"]
AA_pT_n, AA_dpT_n = AA_nform_00_10_pT["FF(pT)"], AA_nform_00_10_pT["dFF(pT)"]

AA_z_w, AA_dz_w = AA_wform_00_10_z["FF(z)"], AA_wform_00_10_z["dFF(z)"]
AA_z_n, AA_dz_n = AA_nform_00_10_z["FF(z)"], AA_nform_00_10_z["dFF(z)"]

pp_pT_w, pp_dpT_w = pp_wform_pT["FF(pT)"], pp_wform_pT["dFF(pT)"]
pp_pT_n, pp_dpT_n = pp_nform_pT["FF(pT)"], pp_nform_pT["dFF(pT)"]

pp_z_w, pp_dz_w = pp_wform_z["FF(z)"], pp_wform_z["dFF(z)"]
pp_z_n, pp_dz_n = pp_nform_z["FF(z)"], pp_nform_z["dFF(z)"]

FF_pT_w = AA_pT_w / pp_pT_w
FF_dpT_w = FF_pT_w * np.sqrt((AA_dpT_w / AA_pT_w) ** 2 + (pp_dpT_w / pp_pT_w) ** 2)

FF_pT_n = AA_pT_n / pp_pT_n
FF_dpT_n = FF_pT_n * np.sqrt((AA_dpT_n / AA_pT_n) ** 2 + (pp_dpT_n / pp_pT_n) ** 2)

FF_z_w = AA_z_w / pp_z_w
FF_dz_w = FF_z_w * np.sqrt((AA_dz_w / AA_z_w) ** 2 + (pp_dz_w / pp_z_w) ** 2)

FF_z_n = AA_pT_n / pp_pT_n
FF_dz_n = FF_z_n * np.sqrt((AA_dz_n / AA_z_n) ** 2 + (pp_dz_n / pp_z_n) ** 2)

pT = AA_wform_00_05_pT["pT"]
dpT = AA_wform_00_05_pT["dpT"] / 2
z = AA_wform_00_05_z["z"]
dz = AA_wform_00_05_z["dz"] / 2

## plot the pT FF:
ax = axes[0]
util.plot_expr_data_on_axis(ax, data_pT, "s")

color = colors[0]
# ax.plot(pT, FF_pT_n, color=color)
# ax.fill_between(pT, FF_pT_n + FF_dpT_n, FF_pT_n - FF_dpT_n, color=color, alpha=0.2)
df = pd.DataFrame({"x": pT, "dx": dpT, "y": FF_pT_n, "dy": FF_dpT_n})
util.plot_theory_on_axis_boxes(ax, df, "y", "x", color, "s", s=100)

color = colors[1]
# ax.plot(pT, FF_pT_w, color=color)
# ax.fill_between(pT, FF_pT_w + FF_dpT_w, FF_pT_w - FF_dpT_w, color=color, alpha=0.2)
df = pd.DataFrame({"x": pT, "dx": dpT, "y": FF_pT_w, "dy": FF_dpT_w})
util.plot_theory_on_axis_boxes(ax, df, "y", "x", color, "s", s=100)
## plot the z FF:
ax = axes[1]
util.plot_expr_data_on_axis(ax, data_z, "s")
color = colors[0]

# ax.plot(z, FF_z_n, color=color)
# ax.fill_between(z, FF_z_n + FF_dz_n, FF_z_n - FF_dz_n, color=color, alpha=0.2)
df = pd.DataFrame({"x": z, "dx": dz, "y": FF_z_n, "dy": FF_dz_n})
util.plot_theory_on_axis_boxes(ax, df, "y", "x", color, "s", s=100)

color = colors[1]
# ax.plot(z, FF_z_w, color=color)
# ax.fill_between(z, FF_z_w + FF_dz_w, FF_z_w - FF_dz_w, color=color, alpha=0.2)
df = pd.DataFrame({"x": z, "dx": dz, "y": FF_z_w, "dy": FF_dz_w})
util.plot_theory_on_axis_boxes(ax, df, "y", "x", color, "s", s=100)

line_colors = {
    r"Without $\tau_f$": colors[0],
    r"With $\tau_f$": colors[1],
}
labels = [Line2D([], [], label=l, color=c) for (l, c) in line_colors.items()]

labels.append(Line2D([], [], color="black", marker="s", label="ATLAS $0$-$10\%$"))
axes[1].legend(loc="upper center", handles=labels)
axes[0].set_xscale("log")
axes[1].set_xscale("log")
axes[0].set_ylabel(r"$R_{D(p_T)}$")
axes[1].set_ylabel(r"$R_{D(z)}$")
axes[0].set_ylim(bottom=0.55, top=2.0)
axes[0].set_xlabel(r"$p^{h^{\pm}}_T$ (GeV)")
axes[1].set_xlabel(r"$z$")
axes[0].text(0.35, 0.9, s=r"Pb-Pb, $\sqrt{s}=2.76$ ATeV", transform=axes[0].transAxes)
axes[0].text(
    0.35, 0.8, s=r"$100< p^{\mathrm{jet}}_T < 398$ GeV", transform=axes[0].transAxes
)
axes[0].text(0.35, 0.7, s=r"$0$-$10\%$", transform=axes[0].transAxes)
plt.show()
