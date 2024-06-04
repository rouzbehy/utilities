#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import util
from COLORS import colors_comparison as colors
import jet_substructure_helpers as substruct
import jetDicts as ddicts

plt.rcParams.update(util.my_rcParams)


## read in the experimental results: jet shape ratio
data_fname = "../../../exp_data/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_00-10_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
tmp = np.loadtxt(data_fname, comments="#", unpack=True, delimiter="\t")
cms_jet_shape = pd.DataFrame({"x": tmp[0], "y": tmp[1], "dx": tmp[2], "dy": tmp[3]})

## read in the experimental results: jet fragmentation function ratios
exp_loc = "../../../exp_data/jet_fragmentation/ATLAS/"
atlas_jet_ff = pd.read_csv(
    exp_loc + "HEPData-ins1511869-v1-csv/Table25.csv", comment="#"
).rename(columns=ddicts.colnames_ATLAS_FF_DZ)

"""
    Read in the simulation results for jet shape ratio 
    and jet fragmentation function ratio
    1. Jet Shape Ratio:
"""
prefix_martini = "/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/"
prefix_jetscape = "/Users/rmyazdi/Documents/research/jetscape_project/v2/jetscape_data/"


martini_fname_pp_nform_time_jetshape = (
    prefix_martini + "v2/production_v2/martini_results/pp/"
)
martini_fname_aa_nform_time_jetshape = (
    prefix_martini + "v2/production_v2/martini_results/final_PbPb_2p76/rset_1/cent_0_5/"
)
martini_fname_pp_wform_time_jetshape = prefix_martini + "v3/calcs/pp_2p76/"
martini_fname_aa_wform_time_jetshape = (
    prefix_martini + "v3/calcs/final_runs/rset_1/2p76/0_5/"
)


jetscape_fname_pp_shape = (
    prefix_jetscape
    + "max_time/maxT_200_highstat/pp_2760_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv"
)
jetscape_fname_aa_shape_matter_vac = (
    prefix_jetscape
    + "sqrt_s_2760/martini_new_kap_matter_vac/PbPb2760_00-05_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv"
)
jetscape_fname_aa_shape_k1 = (
    prefix_jetscape
    + "sqrt_s_2760/martini_new_kappas/PbPb2760_00-05_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv"
)
jetscape_fname_aa_shape_k0 = (
    prefix_jetscape
    + "sqrt_s_2760/martini/PbPb_2760/PbPb2760_00-05_jet_shape_LHC_cuts_jet_rad_0.3_0.30_2.00.csv"
)

jet_shape_martini_nform = substruct.get_shape_without_formation_time(
    martini_fname_pp_nform_time_jetshape, martini_fname_aa_nform_time_jetshape
)

jet_shape_martini_wform = substruct.get_shape_with_formation_time(
    martini_fname_pp_wform_time_jetshape, martini_fname_aa_wform_time_jetshape
)

jet_shape_jscape_k0 = substruct.process_jetscape_jet_shape_results(
    jetscape_fname_pp_shape, jetscape_fname_aa_shape_k0
)

jet_shape_jscape_k1 = substruct.process_jetscape_jet_shape_results(
    jetscape_fname_pp_shape, jetscape_fname_aa_shape_k1
)

jet_shape_jscape_matter_vac = substruct.process_jetscape_jet_shape_results(
    jetscape_fname_pp_shape, jetscape_fname_aa_shape_matter_vac
)
"""
    2. jet fragmentation function ratio
"""
martini_fname_pp_nform_time_FFz = (
    prefix_martini + "v2/production_v2/martini_results/pp/jet_FF_z.csv"
)
martini_fname_aa_nform_time_FFz = (
    prefix_martini
    + "v2/production_v2/martini_results/final_PbPb_2p76/rset_1/cent_0_5/jet_FF_z.csv"
)
martini_fname_pp_wform_time_FFz = prefix_martini + "v3/calcs/pp_2p76/jet_FF_z.csv"
martini_fname_aa_wform_time_FFz = (
    prefix_martini + "v3/calcs/final_runs/rset_1/2p76/0_5/jet_FF_z.csv"
)

jetscape_fname_pp_shape = (
    prefix_jetscape + "max_time/maxT_200_highstat/pp_2760_LHC_FF_z_jet_rad_0.4_2.10.csv"
)
jetscape_fname_aa_shape_matter_vac = (
    prefix_jetscape
    + "sqrt_s_2760/martini_new_kap_matter_vac/PbPb2760_00-05_LHC_FF_z_jet_rad_0.4_2.10.csv"
)
jetscape_fname_aa_shape_k1 = (
    prefix_jetscape
    + "sqrt_s_2760/martini_new_kappas/PbPb2760_00-05_LHC_FF_z_jet_rad_0.4_2.10.csv"
)
jetscape_fname_aa_shape_k0 = (
    prefix_jetscape
    + "sqrt_s_2760/martini/PbPb_2760/PbPb2760_00-05_LHC_FF_z_jet_rad_0.4_2.10.csv"
)

jet_FFz_martini_nform = substruct.get_FFz_ratio_martini_without_form_ratio(
    martini_fname_pp_nform_time_FFz, martini_fname_aa_nform_time_FFz
)

jet_FFz_martini_wform = substruct.get_FFz_ratio_martini_with_form_ratio(
    martini_fname_pp_wform_time_FFz, martini_fname_aa_wform_time_FFz
)

jet_FFz_jscape_k0 = substruct.process_jet_FF_ratio_jetscape(
    jetscape_fname_pp_shape, jetscape_fname_aa_shape_k0
)

jet_FFz_jscape_k1 = substruct.process_jet_FF_ratio_jetscape(
    jetscape_fname_pp_shape, jetscape_fname_aa_shape_k1
)

jet_FFz_jscape_matter_vac = substruct.process_jet_FF_ratio_jetscape(
    jetscape_fname_pp_shape, jetscape_fname_aa_shape_matter_vac
)

## Plot the data:
fig, axes = plt.subplots(
    1,
    2,
    gridspec_kw={
        "top": 0.84,
        "bottom": 0.098,
        "left": 0.06,
        "right": 0.988,
        "wspace": 0.162,
    },
    figsize=(16, 9),
    sharex=False,
    sharey=False,
)

for label, data, ax in zip(
    ["CMS(2014)", "ATLAS(2017)"], [cms_jet_shape, atlas_jet_ff], axes
):
    if "CMS" in label:
        errorboxes = [
            Rectangle(
                (x - delx, y - yerrlow),
                width=2 * delx,
                height=abs(yerrlow) + abs(yerrhigh),
            )
            for x, delx, y, yerrlow, yerrhigh in zip(
                data["x"], data["dx"], data["y"], data["dy"], data["dy"]
            )
        ]
        pc = PatchCollection(
            errorboxes, facecolor="black", edgecolor="black", alpha=0.4
        )
        ax.add_collection(pc)
        scatter = ax.scatter(data["x"], data["y"], color="black", marker="s", s=30)
    else:
        util.plot_expr_data_on_axis(ax, data, "*")

ax = axes[0]
for sim, color, lstyle in zip(
    [
        jet_shape_martini_nform,
        jet_shape_martini_wform,
        jet_shape_jscape_k0,
        jet_shape_jscape_matter_vac,
        jet_shape_jscape_k1,
    ],
    colors[:5],
    ["solid", "solid", "dashed", "dotted", "solid"],
):

    ax.plot(sim["r"], sim["ratio"], color=color, linestyle=lstyle)
    ax.fill_between(
        sim["r"],
        sim["ratio"] - sim["dratio"],
        sim["ratio"] + sim["dratio"],
        color=color,
        alpha=0.2,
    )

ax.set_ylabel(r"$R_{\rho}$")
ax.text(0.02, 0.95, s=r"Pb-Pb @ 2.76 ATeV", transform=ax.transAxes)
ax.text(
    0.02,
    0.8,
    s=r"$0.3<|\eta|<2.0$" + "\n" + r"$100$ GeV $< p^{\mathrm{jet}}_T$",
    transform=ax.transAxes,
)
ax.text(0.02, 0.7, s=r"$R=0.3$, Anti-$k_{T}$", transform=ax.transAxes)
ax.text(0.02, 0.55, s="Exp:0-10\%, \n Theory:0-5\%", transform=ax.transAxes)
ax.set_xlabel(r"$r$")

ax = axes[1]
for sim, color, lstyle in zip(
    [
        jet_FFz_martini_nform,
        jet_FFz_martini_wform,
        jet_FFz_jscape_k0,
        jet_FFz_jscape_matter_vac,
        jet_FFz_jscape_k1,
    ],
    colors[:5],
    ["solid", "solid", "dashed", "dotted", "solid"],
):

    ax.plot(sim["z"], sim["ratio"], color=color, linestyle=lstyle)
    ax.fill_between(
        sim["z"],
        sim["ratio"] - sim["dratio"],
        sim["ratio"] + sim["dratio"],
        color=color,
        alpha=0.2,
    )

ax.set_ylabel(r"$R_{D(z)}$")
ax.set_xlabel(r"$z$")
# ax.text(0.02,0.9,s=r'Pb-Pb @ 2.76 ATeV', transform=ax.transAxes)
ax.text(0.02, 0.9, s=r"$100< p^{\mathrm{jet}}_T < 398$ GeV", transform=ax.transAxes)
ax.text(0.02, 0.8, s=r"$R=0.4, |\eta|<2.1$", transform=ax.transAxes)

line_colors = {
    r"MARTINI (no $\tau_f$)": (colors[0], "solid"),
    r"MARTINI (with $\tau_f$)": (colors[1], "solid"),
    r"MATTER + MARTINI ($\kappa_0$)": (colors[2], "dashed"),
    r"MATTER (vac) + MARTINI ($\kappa_1$)": (colors[3], "dotted"),
    r"MATTER + MARTINI ($\kappa_1$)": (colors[4], "solid"),
}

theo_hands = [
    Line2D([], [], label=l, color=c[0], linestyle=c[1])
    for (l, c) in line_colors.items()
]
legend_labels = [
    Line2D(
        [],
        [],
        color="black",
        label="CMS (2014)",
        marker="*",
        markersize=10,
        linestyle="None",
    ),
    Line2D(
        [],
        [],
        color="black",
        label="ATLAS (2017)",
        marker="P",
        markersize=10,
        linestyle="None",
    ),
]

legend_labels.extend(theo_hands)

axes[0].legend(
    loc="lower left",
    bbox_to_anchor=(0.01, 1.02),
    ncol=4,
    handles=legend_labels,
    fontsize=20,
)
plt.show()
