import sys

sys.path.insert(0, "..")

import matplotlib.pyplot as plt
import helpers as hp
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

plt.rcParams.update(hp.util.my_rcParams)
colors = hp.COLORS.rate_set_colors
rates = {1: "LO", 2: "NLO", 3: "NP"}
experiment_markers = {"cms": "*"}
centralities = ["0_10", "10_30"]
s_marker_size = 100

with_formation_time = hp.helpers.get_jet_shape_ratios(
    prefix="/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/v3/calcs/",
    run_type="with formation time",
)
without_formation_time = hp.helpers.get_jet_shape_ratios(
    prefix="/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/",
    run_type="without formation time",
)

exp_data = hp.helpers.get_experimental_shape_ratio()

fig, axes = plt.subplots(
    1,
    2,
    figsize=(16, 9),
    sharex=True,
    sharey=True,
    gridspec_kw={
        "right": 0.99,
        "left": 0.08,
        "top": 0.95,
        "bottom": 0.1,
        "wspace": 0.05,
        "hspace": 0.05,
    },
)

for icent, cent in enumerate(centralities):
    ax = axes[icent]
    data = exp_data[cent]
    errorboxes = [
        Rectangle(
            (x - delx, y - yerrlow), width=2 * delx, height=abs(yerrlow) + abs(yerrhigh)
        )
        for x, delx, y, yerrlow, yerrhigh in zip(
            data["x"], data["dx"], data["y"], data["dy"], data["dy"]
        )
    ]
    pc = PatchCollection(
        errorboxes, facecolor="black", edgecolor="black", alpha=0.4, zorder=1
    )
    ax.add_collection(pc)
    scatter = ax.scatter(
        data["x"],
        data["y"],
        color="black",
        marker="s",
        s=s_marker_size,
        label="CMS (2014)",
        zorder=1,
    )

    ## plot the one without formation time:
    for rate in rates:
        col = colors[rates[rate]]
        tmp_cent = "10_20" if cent == "10_30" else cent
        spec = with_formation_time[rate][tmp_cent]

        # ax.plot(spec["r"], spec["ratio"], color=col)
        # ax.fill_between(
        #     spec["r"],
        #     spec["ratio"] - spec["dratio"],
        #     spec["ratio"] + spec["dratio"],
        #     color=col,
        #     alpha=0.2,
        # )
        hp.util.plot_theory_on_axis_boxes(
            ax,
            spec,
            xname="r",
            yname="ratio",
            color=col,
            marker="s",
            s=s_marker_size,
            alpha=0.4,
            zorder=0.5,
        )

        # spec = without_formation_time[rate][tmp_cent]
        # ax.plot(spec["r"], spec["ratio"], color=col, linestyle="dotted")
        # ax.fill_between(
        #     spec["r"],
        #     spec["ratio"] - spec["dratio"],
        #     spec["ratio"] + spec["dratio"],
        #     color=col,
        #     alpha=0.2,
        # )

axes[0].text(0.02, 0.9, s=r"Pb-Pb @ 2.76 ATeV", transform=axes[0].transAxes)
axes[0].text(
    0.02,
    0.75,
    s=r"$0.3<|\eta|<2.0$" + "\n" + r"$100$ GeV $< p^{\mathrm{jet}}_T$",
    transform=axes[0].transAxes,
)
axes[0].text(0.02, 0.65, s=r"$R=0.3$, Anti-$k_{T}$", transform=axes[0].transAxes)
for ax in axes:
    ax.set_xlabel(r"$r$")

labels = [Line2D([], [], c=c, label=l) for l, c in colors.items()]
# labels.append(
#     Line2D([], [], c="black", linestyle="solid", label=r"with $\tau_{\mathrm{form.}}$")
# )
# labels.append(
#     Line2D(
#         [], [], c="black", linestyle="dotted", label=r"without $\tau_{\mathrm{form.}}$"
#     )
# )
labels.append(
    Line2D([], [], label="CMS", marker="s", color="black", markersize=10, linewidth=0)
)
axes[1].legend(
    loc="upper left", handles=labels
)  # , handletextpad=0.05, ncol=1, fontsize=27

axes[0].text(0.05, 0.05, r"$0$-$10\%$", transform=axes[0].transAxes)
axes[1].text(0.05, 0.05, r"$10$-$20\%$, CMS: $10$-$30$\%", transform=axes[1].transAxes)
axes[0].set_ylabel(r"$R_{\mathrm{\rho}}$")
axes[0].set_ylim(0.65, 1.45)
plt.show()
