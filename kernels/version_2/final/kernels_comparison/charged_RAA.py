import sys

sys.path.insert(0, "..")

import matplotlib.pyplot as plt
import helpers as hp
from matplotlib.lines import Line2D

plt.rcParams.update(hp.util.my_rcParams)
colors = hp.COLORS.rate_set_colors
rates = {1: "LO", 2: "NLO", 3: "NP"}
experiment_markers = {"alice": "P", "cms": "*", "atlas_1": "o", "atlas_2": "s"}
s_marker_size = 150

with_form_time = hp.helpers.get_charged_RAA(
    run_type="with formation time",
    prefix="/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/v3/calcs/",
    pT_low_cut=33,
    pT_high_cut=100,
)

without_form_time = hp.helpers.get_charged_RAA(
    run_type="without formation time",
    prefix="/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/",
    pT_low_cut=33,
    pT_high_cut=100,
)

exp_RAA = hp.helpers.get_experimental_charged_RAA(pT_low_cut=33, pT_high_cut=100)

fig, axes = plt.subplots(
    1,
    len(with_form_time[1]),
    sharex=True,
    sharey=True,
    figsize=(16, 9),
    gridspec_kw={
        "top": 0.95,
        "bottom": 0.12,
        "left": 0.08,
        "right": 0.99,
        "wspace": 0.05,
    },
)

for expr in exp_RAA:
    expr_cents = exp_RAA[expr]
    if expr != "atlas_1":
        for icent, cent in enumerate(expr_cents):
            ax = axes[icent]
            hp.util.plot_expr_data_on_axis(
                axis=ax,
                data=expr_cents[cent],
                marker=experiment_markers[expr],
                s=s_marker_size,
            )
    else:
        ax = axes[0]
        hp.util.plot_expr_data_on_axis(
            axis=ax, data=expr_cents, marker=experiment_markers[expr], s=s_marker_size
        )


for r in with_form_time:
    rate_with_ftime = with_form_time[r]
    rate_without_ftime = without_form_time[r]
    for icent, cent in enumerate(rate_with_ftime):
        ax = axes[icent]
        RAA = rate_with_ftime[cent]
        ax.plot(RAA["pT"], RAA["RAA"], color=colors[rates[r]])
        ax.fill_between(
            RAA["pT"],
            RAA["RAA"] + RAA["dRAA"],
            RAA["RAA"] - RAA["dRAA"],
            alpha=0.2,
            color=colors[rates[r]],
        )

        RAA = rate_without_ftime[cent]
        ax.plot(RAA["pT"], RAA["RAA"], color=colors[rates[r]], linestyle="dotted")
        ax.fill_between(
            RAA["pT"],
            RAA["RAA"] + RAA["dRAA"],
            RAA["RAA"] - RAA["dRAA"],
            alpha=0.2,
            color=colors[rates[r]],
        )


theo_hands = [Line2D([], [], label=l, color=c) for (l, c) in colors.items()]
theo_hands.append(Line2D([], [], label=r"With $\tau_{\mathrm{form.}}$", color="black"))
theo_hands.append(
    Line2D(
        [],
        [],
        label=r"Without $\tau_{\mathrm{form.}}$",
        color="black",
        linestyle="dotted",
    )
)
expt_hands = [
    Line2D(
        [],
        [],
        color="black",
        label=r"CMS $|\eta|<1.0$",
        marker="*",
        markersize=10,
        linestyle="None",
    ),
    Line2D(
        [],
        [],
        color="black",
        label=r"ALICE $|\eta|<0.8$",
        marker="P",
        markersize=10,
        linestyle="None",
    ),
    Line2D(
        [],
        [],
        color="black",
        label=r"ATLAS $|\eta|<1.0$",
        marker="o",
        markersize=10,
        linestyle="None",
    ),
    Line2D(
        [],
        [],
        color="black",
        label=r"ATLAS $|\eta|<2.0$",
        marker="s",
        markersize=10,
        linestyle="None",
    ),
]

artist = axes[0].legend(
    loc="upper left",
    handles=expt_hands,
    ncol=1,
    handletextpad=0.1,
)
axes[0].add_artist(artist)
axes[0].set_ylim(bottom=0.2, top=1.2)
axes[2].legend(
    loc="upper left", handles=theo_hands, handletextpad=0.1, fontsize=25, ncol=2
)
axes[1].text(
    0.05,
    0.85,
    r"PbPb, $\sqrt{s}=2.76$ ATeV" + "\n" + r"$|\eta|<2.0$",
    transform=axes[1].transAxes,
)

axes[0].set_ylabel(r"$R^{h^{\pm}}_{\mathrm{AA}}$")
for ax in axes:
    ax.set_xlabel(r"$p_T$ (GeV)")

for icent, cent in enumerate(exp_RAA["alice"]):
    ax = axes[icent]
    low, high = cent.split("_")
    centrality = f"${low}$-${high}$\%"
    ax.text(0.65, 0.15, centrality, transform=ax.transAxes)

ax = axes[2]
ax.text(0.45, 0.05, "CMS:$10$-$30$\%", transform=ax.transAxes)

plt.show()
