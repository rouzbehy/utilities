import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D
from numpy import pi, log, arange, exp, loadtxt, sqrt, array, linspace
from pandas import read_csv, DataFrame
from scipy.interpolate import InterpolatedUnivariateSpline as interp1d
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle
import dictionaries as hDicts
import jetDicts as jDicts
import util
from info import kappas
from matplotlib.colors import CSS4_COLORS as css


my_rcParams = {
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 32,
    "lines.linewidth": 4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size": 8,
    "ytick.major.size": 8,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "legend.frameon": False,
    "figure.figsize": (16, 9),
}
plt.rcParams.update(my_rcParams)


def process_jet_shape_data_frame(df, fname):
    njet = 1
    with open(fname, "r") as f:
        line = f.readline()
        line = line.split(" ")[-1]
        njet = float(line)
    df = df[df["rmax"] < 0.31]
    delta_r = df["rmax"] - df["rmin"]
    r = 0.5 * (df["rmax"] + df["rmin"])
    rho = df["N"] / (njet)
    drho = df["dN"] / (njet)
    norm = sum(rho.to_list())
    rho_normed = rho / (delta_r * norm)
    drho_normed = drho / (delta_r * norm)
    retrn_dict = {}
    retrn_dict["rho_normed"] = rho_normed.to_list()
    retrn_dict["drho_normed"] = drho_normed.to_list()
    retrn_dict["r"] = r.to_list()
    retrn_dict["dr"] = delta_r.to_list()
    return DataFrame(data=retrn_dict)


hadron_pT_min, hadron_pT_max = 50, 150
jet_pT_min, jet_pT_max = 65, 320

## charged hadron base line pp
had_pp_spec = read_csv("../../calcs/pp_2p76/chgd_spec.csv", comment="#")
# had_pp_spec = had_pp_spec[had_pp_spec['pTmax']]
had_pp_spec["pT"] = 0.5 * (had_pp_spec["pTmax"] + had_pp_spec["pTmin"])
had_pp_spec = had_pp_spec[had_pp_spec["pT"].between(hadron_pT_min, hadron_pT_max)]
## jet base line from pp
jet_pp_spec = read_csv("../../calcs/pp_2p76/jet_spec.csv", comment="#")
jet_pp_spec["pT"] = 0.5 * (jet_pp_spec["pTmax"] + jet_pp_spec["pTmin"])
jet_pp_spec = jet_pp_spec[jet_pp_spec["pT"].between(jet_pT_min, jet_pT_max)]
## jet shape base line from pp
jet_pp_shape = read_csv("../../calcs/pp_2p76/jet_shape.csv", comment="#")
jet_pp_shape = process_jet_shape_data_frame(
    jet_pp_shape, "../../calcs/pp_2p76/jet_shape.csv"
)

## For Computing Chi-Squared I don't actually need the values
## of kappa_r and kappa_e just the set number. The resulting chi-squared
## is saved to a file with format set_number,chi_squared_had,chi_squared_jet

## Experimental Data Read-In
alice_hadron_name = "../../../exp_data/charged/ALICE/Table_16.csv"
cms_hads_name = "../../../exp_data/charged/CMS/Table5_1.csv"
cms_jets_name = "../../../exp_data/jets/CMS/Table_JETRAA_0p{i}.csv"
cms_jet_shape_name = "../../../exp_data/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_00-10_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
atlas_hads_name_2 = "../../../exp_data/charged/ATLAS/Table33.csv"
atlas_hads_name_1 = "../../../exp_data/charged/ATLAS/Table49.csv"

ALICE_HADS = read_csv(alice_hadron_name, comment="#").rename(
    columns=hDicts.colnames_raa_alice
)
ATLAS_HADS_2 = read_csv(atlas_hads_name_2, comment="#").rename(
    columns=hDicts.colnames_raa_atlas
)
ATLAS_HADS_1 = read_csv(atlas_hads_name_1, comment="#").rename(
    columns=hDicts.colnames_raa_atlas
)
CMS_HADS = read_csv(cms_hads_name, comment="#").rename(columns=hDicts.colnames_raa_cms)
CMS_JETS = {
    f"0p{i}": read_csv(cms_jets_name.format(i=i), comment="#").rename(
        columns=jDicts.colnames_CMS_RAA
    )
    for i in [2, 3, 4]
}
tmp = loadtxt(cms_jet_shape_name, comments="#", unpack=True, delimiter="\t")
CMS_JET_SHAPE = DataFrame({"x": tmp[0], "y": tmp[1], "dx": tmp[2], "dy": tmp[3]})

## place cut on hadron pt:
ALICE_HADS = ALICE_HADS[ALICE_HADS["x"] > hadron_pT_min]
ATLAS_HADS_1 = ATLAS_HADS_1[ATLAS_HADS_1["x"] > hadron_pT_min]
ATLAS_HADS_2 = ATLAS_HADS_2[ATLAS_HADS_2["x"] > hadron_pT_min]
CMS_HADS = CMS_HADS[CMS_HADS["x"] > hadron_pT_min]

fname_hads = (
    "../../calcs/pbpb_2p76/param_fit_runs/rset_{r}/kappaset_{k}/hadron_spectra.csv"
)
fname_jets = (
    "../../calcs/pbpb_2p76/param_fit_runs/rset_{r}/kappaset_{k}/jet_spectra.csv"
)
fname_jet_shapes = (
    "../../calcs/pbpb_2p76/param_fit_runs/rset_{r}/kappaset_{k}/jet_shape.csv"
)


rset = int(sys.argv[1])
rate_sets = {1: "LO", 2: "NLO", 3: "NP"}
cone_radii = ["0p2", "0p3", "0p4"]

had_aa_spec = {
    kset: read_csv(fname_hads.format(r=rset, k=kset), comment="#") for kset in kappas
}
jet_aa_spec = {
    kset: read_csv(fname_jets.format(r=rset, k=kset), comment="#") for kset in kappas
}
tmp_shape_dt = {
    kset: read_csv(fname_jet_shapes.format(r=rset, k=kset), comment="#")
    for kset in kappas
}
jet_aa_shape = {
    kset: process_jet_shape_data_frame(
        tmp_shape_dt[kset], fname_jet_shapes.format(r=rset, k=kset)
    )
    for kset in kappas
}
"""
    Want 3 figures per rate set
        1. Jet RAA as a function of cone radius (3 x 1 figure, one subplot for each jet radius)
        2. Charged Hadron RAA
        3. Jet Shape Ratio
"""

# fig0, (ax0,ax1) = plt.subplots(nrows=1,ncols=2, figsize=(8,8))
all_kappas = array([v for _, v in kappas.items()])
current_kappas = {i: kappas[i] for i in kappas}
radius = 0.15

# fname  = "../../chisquared_work/chisquared_rateset_{rset}.csv"
# tmp    = loadtxt(fname.format(rset=rset),unpack=True, skiprows=1, delimiter=',')
# data   = DataFrame({'kappa':tmp[0], 'kappa_r':tmp[1], 'kappa_e':tmp[2], 'chisq':tmp[-1]})
### transform kappa_e to 1/kappa_e
# data['1/kappa_e'] = (1/data['kappa_e'])**(1./2)
# data['1/kappa_r'] = (1/data['kappa_r'])**(1./3)
# radius = 0.015
##kappas = [0,1,2,3]
##chosen_kappas = data[data['chisq'].between(chi_sq_lim_bottom, chi_sq_lim_top)]
# chosen_kappas = data[data['kappa'].between(0,3)]
actual_chosen_kappas = [0, 1, 2, 3]
# other_kappas  = data[data['kappa'] > 3 ]
# lcolors = mpl.colormaps['Pastel2'].colors #[css['springgreen'],css['red'],css['gold'],css['dodgerblue']]
# first_lcolors = mpl.colormaps["tab10"]  # (linspace(0, 1, len(chosen_kappas)))
# lcolors = first_lcolors.colors
lcolors = ["blue", "orange", "green", "purple"]

plot_colors = {k: lcolors[k] for k in actual_chosen_kappas}

# for it, item in enumerate(zip(actual_chosen_kappas['kappa'], actual_chosen_kappas['kappa_r'], actual_chosen_kappas['kappa_e'])):
#    plot_colors[item[0]] = lcolors[it]
legend_labels = [
    Line2D(
        [],
        [],
        label=f"({kappas[k][0]},{kappas[k][1]})",
        color=plot_colors[k],
        linestyle="solid",
    )
    for k in plot_colors
]


fig1, ax = plt.subplots(
    nrows=1,
    ncols=1,
    sharex=True,
    sharey=True,
    gridspec_kw={
        "left": 0.08,
        "right": 0.95,
        "top": 0.95,
        "bottom": 0.12,
        "hspace": 0.04,
    },
    figsize=(16, 9),
    # gridspec_kw={
    #     "left": 0.08,
    #     "right": 0.95,
    #     "top": 0.95,
    #     "bottom": 0.12,
    #     "hspace": 0.04,
    # },
)
ax.set_xlabel(r"$p_T$ (GeV)")
ax.set_ylabel(r"$R^{\mathrm{jet}}_{\mathrm{AA}}$")
ax.set_ylim(top=3.5, bottom=-0.01)
ax.set_xlim(left=60, right=320)
ax.text(
    0.05,
    0.9,
    rate_sets[rset]
    + ", Pb-Pb @ 2.76 ATeV"
    + r" $0$-$5\%$"
    + r", $|\eta^{\mathrm{jet}}|<2.0$",
    transform=ax.transAxes,
)
factor = 1.0
for itn, item in enumerate(CMS_JETS):
    additive_const = factor * itn
    CMS_JETS[item].y += additive_const
    util.plot_expr_data_on_axis(
        ax, CMS_JETS[item], "s", color="red", face="red", alpha=0.2
    )
    label = item.replace("p", ".")
    if itn > 0:
        label += f" (+ {itn})"
    ax.text(255, 0.8 + (itn), s=f"R={label}")

for k in kappas:
    aa_jet = jet_aa_spec[k]
    aa_jet["pT"] = 0.5 * (aa_jet["pTmax"] + aa_jet["pTmin"])
    aa_jet = aa_jet[aa_jet["pT"].between(jet_pT_min, jet_pT_max)]
    pT_jet = aa_jet["pT"]
    pp_jet = jet_pp_spec
    color = plot_colors[k] if k in plot_colors else "grey"
    zorder = 0 if k not in actual_chosen_kappas else 0.5
    alpha = 1 if k in actual_chosen_kappas else 0.05
    for iR, R in enumerate(cone_radii):
        additive_const = factor * iR
        raa_jet = array([a / p for (a, p) in zip(aa_jet[f"N{R}"], pp_jet[f"N{R}"])])
        draa_jet = array(
            [
                r * sqrt(da * da / (a * a) + dp * dp / (p * p))
                for (r, da, a, dp, p) in zip(
                    raa_jet,
                    aa_jet[f"dN{R}"],
                    aa_jet[f"N{R}"],
                    pp_jet[f"dN{R}"],
                    pp_jet[f"N{R}"],
                )
            ]
        )
        ax.plot(pT_jet, raa_jet + additive_const, color=color, zorder=zorder)
        ax.fill_between(
            pT_jet,
            additive_const + raa_jet - draa_jet,
            additive_const + raa_jet + draa_jet,
            color=color,
            alpha=alpha,
        )

expt_hands = [
    Line2D(
        [],
        [],
        color="red",
        label="CMS",
        marker="s",
        markersize=10,
        linestyle="None",
    )
]
# tmp = [l for l in legend_labels]
# tmp.append(expt_hands[0])
ax.legend(
    # title=r"$(\kappa_r,\kappa_e)$",
    handles=expt_hands,
    loc="upper right",
    fontsize=30,
    # bbox_to_anchor=(1.3, 0.8),
)
# axes1[0].legend(title=r'$(\kappa_r,\kappa_e)$', handles=tmp, bbox_to_anchor=(1.0,0.95), fontsize=30)
# axes1[1].legend(handles=expt_hands, bbox_to_anchor=(1.1,0.15), fontsize=30)

fig2, ax2 = plt.subplots(
    1,
    1,
    gridspec_kw={
        "left": 0.08,
        "right": 0.95,
        "top": 0.95,
        "bottom": 0.12,
        "hspace": 0.04,
    },
    figsize=(16, 9),
)
# gridspec_kw={'left':0.1,'right':0.99,'top':0.99,'bottom':0.1})
# ax2.set_title(rate_sets[rset]+ r" rate set, Pb-Pb @ 2.76 ATeV, $0$-$5\%$")
util.plot_expr_data_on_axis(
    ax2, CMS_HADS, marker="*", color="red", face="red", alpha=0.2
)
util.plot_expr_data_on_axis(
    ax2, ATLAS_HADS_2, marker="^", color="red", face="red", alpha=0.2
)
util.plot_expr_data_on_axis(
    ax2, ATLAS_HADS_1, marker="v", color="red", face="red", alpha=0.2
)
util.plot_expr_data_on_axis(
    ax2, ALICE_HADS, marker="P", color="red", face="red", alpha=0.2
)
expt_hands = [
    Line2D(
        [],
        [],
        color="red",
        label="CMS $|\eta|<1.0$",
        marker="*",
        markersize=10,
        linestyle="None",
    ),
    Line2D(
        [],
        [],
        color="red",
        label="ALICE $|\eta|<0.8$",
        marker="P",
        markersize=10,
        linestyle="None",
    ),
    Line2D(
        [],
        [],
        color="red",
        label="ATLAS $|\eta|<1.0$",
        marker="v",
        markersize=10,
        linestyle="None",
    ),
    Line2D(
        [],
        [],
        color="red",
        label="ATLAS $|\eta|<2.0$",
        marker="^",
        markersize=10,
        linestyle="None",
    ),
]
artist = ax2.legend(loc="upper right", handles=expt_hands)
ax2.add_artist(artist)
ax2.set_xlabel(r"$p^{h^{\pm}}_T$ (GeV)")
ax2.set_ylabel(r"$R^{h^{\pm}}_{\mathrm{AA}}$")
ax2.legend(title=r"$(\kappa_r,\kappa_e)$", handles=legend_labels, loc="upper left")
ax2.set_ylim(top=1.55, bottom=0.02)
for k in kappas:
    aa_hads = had_aa_spec[k]
    pp_hads = had_pp_spec

    aa_hads["pT"] = 0.5 * (aa_hads["pTmax"] + aa_hads["pTmin"])
    aa_hads = aa_hads[aa_hads["pT"].between(hadron_pT_min, hadron_pT_max)]
    pT_hads = aa_hads["pT"]

    color = plot_colors[k] if k in plot_colors else "grey"
    # color = colors[k] if k in chosen_kappas else 'grey'
    zorder = 0 if k not in actual_chosen_kappas else 5
    alpha = 0.6 if k in actual_chosen_kappas else 0.05
    raa = array([a / p for (a, p) in zip(aa_hads[f"N"], pp_hads[f"N"])])
    draa = array(
        [
            r * sqrt(da * da / (a * a) + dp * dp / (p * p))
            for (r, da, a, dp, p) in zip(
                raa, aa_hads[f"dN"], aa_hads[f"N"], pp_hads[f"dN"], pp_hads[f"N"]
            )
        ]
    )
    ax2.plot(pT_hads, raa, color=color, zorder=zorder)
    ax2.fill_between(pT_hads, raa - draa, raa + draa, color=color, alpha=alpha)
ax2.text(
    0.25,
    0.8,
    rate_sets[rset] + r", Pb-Pb @ 2.76 ATeV, $0$-$5\%$",
    transform=ax2.transAxes,
)


fig3, ax3 = plt.subplots(
    1, 1, gridspec_kw={"left": 0.09, "right": 0.99, "top": 0.99, "bottom": 0.1}
)
artist = ax3.legend(
    title=r"$(\kappa_r,\kappa_e)$", handles=legend_labels, loc="upper left"
)
ax3.add_artist(artist)
# ax3.text(0.05, 0.1, rate_sets[rset]+ ", Pb-Pb @ 2.76 ATeV"+r" $0$-$5\%$", transform=ax3.transAxes)
errorboxes = [
    Rectangle(
        (x - delx, y - yerrlow), width=2 * delx, height=abs(yerrlow) + abs(yerrhigh)
    )
    for x, delx, y, yerrlow, yerrhigh in zip(
        CMS_JET_SHAPE["x"],
        CMS_JET_SHAPE["dx"],
        CMS_JET_SHAPE["y"],
        CMS_JET_SHAPE["dy"],
        -1 * CMS_JET_SHAPE["dy"],
    )
]

# Create patch collection with specified colour/alpha
pc = PatchCollection(errorboxes, facecolor="red", alpha=0.4)
ax3.add_collection(pc)
ax3.scatter(
    CMS_JET_SHAPE["x"],
    CMS_JET_SHAPE["y"],
    color="red",
    marker="*",
    s=200,
    zorder=30,
)
ax3.set_xlabel("$r$")
ax3.set_ylabel(r"$R_{\rho}(r)$")  #'$\rho(r)_{PbPb}/\rho(r)_{pp}$')
for k in kappas:
    aa_shapes = jet_aa_shape[k]
    pp_shapes = jet_pp_shape
    radii = aa_shapes["r"]
    delta_radii = aa_shapes["dr"]
    color = plot_colors[k] if k in plot_colors else "grey"
    zorder = 3 if k in actual_chosen_kappas else 0.05
    alpha = 0.6 if k in actual_chosen_kappas else 0.05
    marker = "o" if k in actual_chosen_kappas else "s"
    markersize = 50 if k in actual_chosen_kappas else 10
    raa = array(
        [a / p for (a, p) in zip(aa_shapes[f"rho_normed"], pp_shapes[f"rho_normed"])]
    )
    draa = array(
        [
            v * sqrt(da * da / (a * a) + dp * dp / (p * p))
            for (v, da, a, dp, p) in zip(
                raa,
                aa_shapes[f"drho_normed"],
                aa_shapes[f"rho_normed"],
                pp_shapes[f"drho_normed"],
                pp_shapes[f"rho_normed"],
            )
        ]
    )

    # ax3.plot(radii, raa, color=color, zorder=zorder)
    # ax3.fill_between(radii, raa-draa, raa+draa, color=color,alpha=alpha, zorder=zorder)
    df = DataFrame({"x": radii, "dx": delta_radii / 2, "y": raa, "dy": draa})
    util.plot_theory_on_axis_boxes(
        ax3,
        df,
        yname="y",
        xname="x",
        color=color,
        marker=marker,
        s=markersize,
        alpha=alpha,
        zorder=zorder,
    )

expt_hand = [
    Line2D(
        [],
        [],
        color="red",
        label=r"CMS, Cent:$0$-$10\%$",
        marker="*",
        markersize=20,
        linewidth=0,
    )
]
ax3.legend(loc="lower left", handles=expt_hand)
ax3.text(
    0.3,
    0.9,
    rate_sets[rset] + ", Pb-Pb @ 2.76 ATeV" + r" $0$-$5\%$",
    transform=ax3.transAxes,
)
ax3.text(
    0.3,
    0.8,
    s=r"$0.3<|\eta|<2.0$, $p^{\mathrm{jet}}_T > 100$ GeV",
    transform=ax3.transAxes,
)
ax3.text(0.3, 0.7, s=r"$R=0.3$, Anti-$k_{T}$", transform=ax3.transAxes)
plt.show()
