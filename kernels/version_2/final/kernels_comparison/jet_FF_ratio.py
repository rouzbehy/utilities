import sys

sys.path.insert(0, "..")

import helpers as hp
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import pandas as pd

plt.rcParams.update(hp.util.my_rcParams)
colors = hp.COLORS.rate_set_colors
rates = {1: "LO", 2: "NLO", 3: "NP"}
experiment_markers = {"atlas": "*"}
s_marker_size = 100
pTmin, pTmax, zmin, zmax = 0.1, 120, 0.005, 2.0

with_formation_time = hp.helpers.get_jet_FF_ratios(
    prefix="/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/v3/calcs/",
    run_type="with formation time",
    pTmin=pTmin,
    pTmax=pTmax,
    zmin=zmin,
    zmax=zmax,
)
without_formation_time = hp.helpers.get_jet_FF_ratios(
    prefix="/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/",
    run_type="without formation time",
    pTmin=pTmin,
    pTmax=pTmax,
    zmin=zmin,
    zmax=zmax,
)


exp_data = hp.helpers.get_experimenta_FF_ratio()

fig, axes = plt.subplots(
    1,
    2,
    figsize=(16, 9),
    gridspec_kw={"top": 0.95, "left": 0.08, "right": 0.99, "bottom": 0.12},
)

for ax, dat, mart_wft, mart_woft in zip(
    axes, exp_data, with_formation_time, without_formation_time
):
    hp.util.plot_expr_data_on_axis(axis=ax, data=dat, marker="s", s=s_marker_size)

    for r in rates:
        col = colors[rates[r]]
        with_tau = mart_wft[r]
        without_tau = mart_woft[r]
        # hp.util.plot_theory_on_axis_boxes(
        #     ax, with_tau, "ratio", "x", col, "s", s=s_marker_size, alpha=0.4
        # )
        ax.plot(with_tau["x"], with_tau["ratio"], color=col, linestyle="solid")
        ax.fill_between(
            with_tau["x"],
            with_tau["ratio"] - with_tau["dratio"],
            with_tau["ratio"] + with_tau["dratio"],
            color=col,
            alpha=0.2,
        )
        ax.plot(without_tau["x"], without_tau["ratio"], color=col, linestyle="dotted")
        ax.fill_between(
            without_tau["x"],
            without_tau["ratio"] - without_tau["dratio"],
            without_tau["ratio"] + without_tau["dratio"],
            color=col,
            alpha=0.2,
        )


labels = [Line2D([], [], c=c, label=l) for l, c in colors.items()]
labels.append(
    Line2D([], [], c="black", linestyle="solid", label=r"with $\tau_{\mathrm{form.}}$")
)
labels.append(
    Line2D(
        [], [], c="black", linestyle="dotted", label=r"without $\tau_{\mathrm{form.}}$"
    )
)
labels.append(
    Line2D(
        [],
        [],
        label="ATLAS",
        marker="s",  # experiment_markers["atlas"],
        color="black",
        markersize=10,
        linewidth=0,
    )
)
axes[1].legend(
    loc="upper center", handles=labels
)  # , handletextpad=0.05, ncol=1, fontsize=30


for ax in axes:
    ax.set_xscale("log")
axes[0].set_ylabel(r"$R_{D(p_T)}$")
axes[0].set_xlabel(r"$p^{h^{\pm}}_T$ (GeV)")
axes[1].set_ylabel(r"$R_{D(z)}$")
axes[1].set_xlabel(r"$z$")

xloc_tag = 0.3
axes[0].text(
    xloc_tag, 0.9, s=r"Pb-Pb, $\sqrt{s}=2.76$ ATeV", transform=axes[0].transAxes
)
axes[0].text(
    xloc_tag, 0.8, s=r"$100< p^{\mathrm{jet}}_T < 398$ GeV", transform=axes[0].transAxes
)
axes[0].text(xloc_tag, 0.7, s=r"$0$-$10\%$", transform=axes[0].transAxes)


plt.show()
