import sys

sys.path.insert(0, "..")

import matplotlib.pyplot as plt
import helpers as hp
from matplotlib.lines import Line2D

plt.rcParams.update(hp.util.my_rcParams)
colors = hp.COLORS.rate_set_colors
rates = {1: "LO", 2: "NLO", 3: "NP"}
experiment_markers = {"cms": "*"}
pT_cut_low, pT_cut_high = 60, 300
s_marker_size = 15


with_form_time = hp.helpers.get_jet_RAA(
    run_type="with formation time",
    prefix="/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/v3/calcs/",
    pT_low_cut=pT_cut_low,
    pT_high_cut=pT_cut_high,
)

without_form_time = hp.helpers.get_jet_RAA(
    run_type="without formation time",
    prefix="/Users/rmyazdi/Documents/research/KERNELS_NLO_NP_PART2/",
    pT_low_cut=pT_cut_low,
    pT_high_cut=pT_cut_high,
)

exp_RAA = hp.helpers.get_experimental_jet_RAA(
    pT_low_cut=pT_cut_low, pT_high_cut=pT_cut_high
)


num_centralities = len(with_form_time[1])
num_cone_radii = len(with_form_time[1]["0_5"])

fig, axes = plt.subplots(
    nrows=num_cone_radii,
    ncols=num_centralities,
    figsize=(16, 9),
    sharex=True,
    sharey=True,
    gridspec_kw={
        "wspace": 0.05,
        "hspace": 0.05,
        "right": 0.99,
        "bottom": 0.11,
        "left": 0.1,
    },
)

for icent, cent in enumerate(exp_RAA):
    for ir, r in enumerate(exp_RAA[cent]):
        ax = axes[ir][icent]
        hp.util.plot_expr_data_on_axis(
            axis=ax,
            data=exp_RAA[cent][r],
            marker=experiment_markers["cms"],
            s=s_marker_size,
        )

## plot with formation time:
for r in [1, 2, 3]:
    for icent, cent in enumerate(["0_5", "5_10", "10_20"]):
        for ir, rad in enumerate(["0p2", "0p3", "0p4"]):
            res = with_form_time[r][cent][rad]
            ax = axes[ir][icent]
            ax.plot(res["pT"], res["raa"], color=colors[rates[r]], linewidth=1.2)
            ax.fill_between(
                res["pT"],
                res["raa"] - res["draa"],
                res["raa"] + res["draa"],
                color=colors[rates[r]],
                alpha=0.2,
            )
            # res = without_form_time[r][cent][rad]
            # ax.plot(res["pT"], res["raa"], color=colors[rates[r]], linestyle="dotted")
            # ax.fill_between(
            #     res["pT"],
            #     res["raa"] - res["draa"],
            #     res["raa"] + res["draa"],
            #     color=colors[rates[r]],
            #     alpha=0.2,
            # )


axes[0][0].set_ylim(bottom=0.14, top=0.82)
for ax, cent in zip(axes[:, 0], [r"$0$-$5$\%", r"$5$-$10$\%", r"$10$-$20$\%"]):
    ax.set_ylabel(r"$R^{\mathrm{jet}}_{\mathrm{AA}}$")
    if "20" not in cent:
        ax.text(0.05, 0.8, cent, transform=ax.transAxes)
    else:
        ax.text(0.05, 0.8, cent + r" (CMS:$10$-$30\%$)", transform=ax.transAxes)

for ax, radius in zip(axes[0, :], ["0.2", "0.3", "0.4"]):
    ax.text(0.3, 0.3, f"R={radius}", transform=ax.transAxes)

for ax in axes[2, :]:
    ax.set_xlabel(r"$p_T$ (GeV)")

## legend work:
axes[1][2].text(
    0.3,
    0.3,
    "PbPb, $\sqrt{s}=2.76$ATeV" + "\n" + "$|\eta|<2.0$",
    transform=axes[1][2].transAxes,
    fontsize=20,
)

theo_hands = [Line2D([], [], label=l, color=c) for (l, c) in colors.items()]
# theo_hands.append(Line2D([], [], label=r"With $\tau_{\mathrm{form.}}$", color="black"))
# theo_hands.append(
#     Line2D(
#         [],
#         [],
#         label=r"Without $\tau_{\mathrm{form.}}$",
#         color="black",
#         linestyle="dotted",
#     )
# )
axes[2][2].legend(
    loc="lower left",
    # bbox_to_anchor=(-0.1, 1.02),
    ncol=3,
    handles=theo_hands,
    fontsize=17,
)
plt.show()
