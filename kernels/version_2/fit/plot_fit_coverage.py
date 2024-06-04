import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import ticker as mticker
import dictionaries as ddicts
import jetDicts
from scipy.interpolate import interp1d
import info
import util
from sys import argv

plt.rcParams.update(info.my_rcParams)
kappas = info.kappas

"""
    Plot the coverage of the fit runs. This means the script will plot the
    charged hadron RAA, Jet RAA and jet shape ratio for all rate sets and
    kappa sets.
"""
all_kappas = [val for (k, val) in kappas.items()]
kappas_r, kappas_e = [e[0] for e in all_kappas], [e[1] for e in all_kappas]
kappa_sets = np.arange(0, 50)
pt_chgd_low = util.pt_chgd_low
pt_chgd_high = util.pt_chgd_high
pt_jet_low = util.pt_jet_low
pt_jet_high = util.pt_jet_high
jet_cone_radii = util.jet_cone_radii
edge_kappas = [0, 1, 2, 3]  # [int(e) for e in argv[1:]]  # it will be 0,1,2,3
colours = ["blue", "orange", "green", "purple"]  # mpl.colormaps["tab10"].colors
cols = {k: colours[i] for i, k in enumerate(edge_kappas)}
labels = [
    Line2D(
        [],
        [],
        color=cols[k],
        label=r"$(\kappa_r,\kappa_e)=$" + f"({kappas_r[k]},{kappas_e[k]})",
    )
    for i, k in enumerate(edge_kappas)
]


## Charged Hadron RAA data
ALICE_HADS_00_05 = pd.read_csv(
    "../../../exp_data/charged/ALICE/Table_16.csv", comment="#"
).rename(columns=ddicts.colnames_raa_alice)
ALICE_HADS_00_05 = ALICE_HADS_00_05[
    ALICE_HADS_00_05["x"].between(pt_chgd_low, pt_chgd_high)
]
CMS_HADS_00_05 = pd.read_csv(
    "../../../exp_data/charged/CMS/Table5_1.csv", comment="#"
).rename(columns=ddicts.colnames_raa_cms)
CMS_HADS_00_05 = CMS_HADS_00_05[CMS_HADS_00_05["x"].between(pt_chgd_low, pt_chgd_high)]
ATLAS_HADS_00_05 = pd.read_csv(
    "../../../exp_data/charged/ATLAS/Table33.csv", comment="#"
).rename(columns=ddicts.colnames_raa_atlas)
ATLAS_HADS_00_05 = ATLAS_HADS_00_05[
    ATLAS_HADS_00_05["x"].between(pt_chgd_low, pt_chgd_high)
]
## Jet RAA data
CMS_JETS_00_05 = {
    f"0p{i}": pd.read_csv(
        f"../../../exp_data/jets/CMS/Table_JETRAA_0p{i}.csv", comment="#"
    ).rename(columns=jetDicts.colnames_CMS_RAA)
    for i in [2, 3, 4]
}
## Jet Shape Ratio
fname = "../../../exp_data/jet_shape/CMS_PbPb_2760_Jet_ShapeRatio_00-10_pT=100-INF_R=0p3_eta=0p3-2p0.dat"
tmp = np.loadtxt(fname, comments="#", unpack=True, delimiter="\t")
CMS_JET_SHAPE = pd.DataFrame(
    {
        "x": tmp[0],
        "y": tmp[1],
        "dx": tmp[2],
        "dy": tmp[3],
        "xlow": tmp[0] - tmp[2] / 2,
        "xhigh": tmp[0] + tmp[2] / 2,
    }
)

## Proton-Proton Baselines
pp_chgd = util.get_charged_hadron_spec(None, None, "pp", pt_chgd_low, pt_chgd_high)
pp_jets = util.get_jet_spec(None, None, "pp", pt_jet_low, pt_jet_high)
pp_jshp = util.get_jet_shape(None, None, "pp")

## Pb-Pb calculations
raa_calc_chgd = {}
raa_calc_jet = {}
raa_calc_shape = {}
for rset in util.rset_shorts:
    raa_calc_chgd[rset] = {}
    raa_calc_jet[rset] = {}
    raa_calc_shape[rset] = {}
    for kset in kappa_sets:
        ## Charged hadrons:
        tmp = util.get_charged_hadron_spec(
            rset, kset, "pbpb", pt_chgd_low, pt_chgd_high
        )
        tmp_raa = tmp["N"] / pp_chgd["N"]
        tmp_draa = tmp_raa * np.sqrt(
            (pp_chgd["dN"] / pp_chgd["N"]) ** 2 + (tmp["dN"] / tmp["N"]) ** 2
        )
        raa_calc_chgd[rset][kset] = [pp_chgd["pT"], tmp_raa, tmp_draa]
        ## Jet Shape:
        tmp = util.get_jet_shape(rset, kset, "pbpb")
        tmp_raa = tmp["rho_normed"] / pp_jshp["rho_normed"]
        tmp_draa = tmp_raa * np.sqrt(
            (pp_jshp["drho_normed"] / pp_jshp["rho_normed"]) ** 2
            + (tmp["drho_normed"] / tmp["rho_normed"]) ** 2
        )
        raa_calc_shape[rset][kset] = [pp_jshp["r"], tmp_raa, tmp_draa]

        ## Jets (gotta loop over the cone radii)

        tmp = util.get_jet_spec(rset, kset, "pbpb", pt_jet_low, pt_jet_high)
        raa_calc_jet[rset][kset] = {}
        for r in jet_cone_radii:
            raa = tmp[f"N{r}"] / pp_jets[f"N{r}"]
            a, da = f"N{r}", f"dN{r}"
            draa = raa * np.sqrt(
                (tmp[da] / tmp[a]) ** 2 + (pp_jets[da] / pp_jets[a]) ** 2
            )
            raa_calc_jet[rset][kset][r] = (tmp["pT"], raa, draa)

alpha = 0.1
## Get plotting:
## CHARGED HADRON RAA
fig1, axes1 = plt.subplots(
    3, 1, sharex=True, sharey=True, gridspec_kw={"right": 0.75, "left": 0.05}
)

for iax, ax in enumerate(axes1):
    ax.set_ylim(bottom=-0.01, top=1.01)
    util.plot_expr_data_on_axis(
        ax, ALICE_HADS_00_05, marker="*", color="red", face="red"
    )
    util.plot_expr_data_on_axis(
        ax, ATLAS_HADS_00_05, marker="o", color="red", face="red"
    )
    util.plot_expr_data_on_axis(
        ax, ATLAS_HADS_00_05, marker="p", color="red", face="red"
    )
    rate_set = iax + 1
    spec = raa_calc_chgd[rset]
    ax.text(0.70, 0.1, s=util.rate_set_names[rate_set], transform=ax.transAxes)
    for ic, item in enumerate(spec):
        pT, raa, draa = spec[item]
        colour = cols[ic] if ic in edge_kappas else "grey"
        zorder = 1 if ic in edge_kappas else 0
        ax.plot(pT, raa, color=colour, zorder=zorder)
        ax.fill_between(
            pT, raa - draa, raa + draa, color=colour, alpha=alpha, zorder=zorder
        )
# labellings:
axes1[2].set_xlabel(r"$p^{h^{\pm}}_T (GeV)$")
for ax in axes1:
    ax.set_ylabel(r"$R_{\mathrm{AA}}$")
axes1[1].legend(
    loc="upper left", handles=labels, bbox_to_anchor=(0.97, 0.9), fontsize=20
)
# fig1.savefig(fname= "../../plots/charged_hadron_raa.png", dpi=200)
# plt.close(fig1)
## JET SHAPE RATIO:
fig2, axes2 = plt.subplots(
    3, 1, sharex=True, sharey=True, gridspec_kw={"right": 0.75, "left": 0.09}
)
for iax, ax in enumerate(axes2):
    ax.set_ylim(bottom=0.49, top=1.61)
    errorboxes = [
        Rectangle(
            (x - delx, y - yerrlow), width=2 * delx, height=abs(yerrlow) + abs(yerrhigh)
        )
        for x, delx, y, yerrlow, yerrhigh in zip(
            CMS_JET_SHAPE["x"],
            CMS_JET_SHAPE["dx"],
            CMS_JET_SHAPE["y"],
            CMS_JET_SHAPE["dy"],
            CMS_JET_SHAPE["dy"],
        )
    ]
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor="black", edgecolor="black", alpha=0.4)
    ax.add_collection(pc)
    scatter = ax.scatter(
        CMS_JET_SHAPE["x"],
        CMS_JET_SHAPE["y"],
        color="black",
        marker="s",
        s=30,
        label="CMS (2014)",
    )
    rset = iax + 1
    ax.text(0.07, 0.9, s=util.rate_set_names[rset], transform=ax.transAxes)
    spec = raa_calc_shape[rset]
    for ic, item in enumerate(spec):
        x, y, dy = spec[item]
        colour = cols[ic] if ic in edge_kappas else "grey"
        zorder = 1 if ic in edge_kappas else 0
        ax.plot(x, y, color=colour, zorder=zorder)
        ax.fill_between(x, y - dy, y + dy, color=colour, alpha=alpha, zorder=zorder)
# labellings:
axes2[2].set_xlabel(r"$r$")
for ax in axes2:
    ax.set_ylabel(r"$R_{\rho}$")
axes2[1].legend(
    loc="upper left", handles=labels, bbox_to_anchor=(0.97, 0.9), fontsize=20
)
# fig2.savefig(fname= "../../plots/jet_shape_ratio.png", dpi=200)
# plt.close(fig2)
## Jet RAA
for rset in util.rset_shorts:
    spec = raa_calc_jet[rset]
    fig, ax = plt.subplots(
        1, 1, sharex=True, sharey=True, gridspec_kw={"right": 0.99, "left": 0.08}
    )
    for iR, R in enumerate(jet_cone_radii):
        # ax.set_ylim(bottom=-0.01, top=1.61)
        CMS_JETS_00_05[R].y += iR
        util.plot_expr_data_on_axis(
            ax, CMS_JETS_00_05[R], marker="p", color="red", face="red"
        )
        ax.text(0.8, 0.05, s=f"R={jet_cone_radii[R]:0.1f}", transform=ax.transAxes)
        for ic, item in enumerate(spec):
            colour = cols[ic] if ic in edge_kappas else "grey"
            zorder = 1 if ic in edge_kappas else 0
            x, y, dy = spec[item][R]

            ax.plot(x, y, color=colour, zorder=zorder)
            ax.fill_between(x, y - dy, y + dy, color=colour, alpha=alpha, zorder=zorder)
    ax.set_title(f"{util.rate_set_names[rset]}")
    # labellings:
    ax.set_xlabel(r"$p^{\mathrm{Jet}}_T (GeV)$")
    ax.set_ylabel(r"$R^{\mathrm{jet}}_{\mathrm{AA}}$")
    # axes[1].legend(loc='upper left', handles=labels, bbox_to_anchor=(0.97,0.9), fontsize=20)
    ax.legend(loc="upper right", handles=labels)
    # fig.savefig(fname=f"../../plots/jet_RAA_rset_{rset}.png", dpi=200)
    # plt.close(fig)
plt.show()
