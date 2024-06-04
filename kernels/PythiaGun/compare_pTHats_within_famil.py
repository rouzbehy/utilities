import matplotlib.pyplot as plt
from reader import Reader
import helpers as hp
from sys import argv
from typing import List

plt.rcParams.update(hp.my_rcParams)


def main(pThat: List[int], rate: int, delay: str, base: str):

    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(16, 9),
        sharex=True,
        sharey=True,
        gridspec_kw={
            "left": 0.1,
            "right": 0.99,
            "top": 0.98,
            "bottom": 0.11,
            "hspace": 0.04,
            "wspace": 0.04,
        },
    )
    axes = axes.flatten()
    taus = hp.taus
    labels1 = [
        plt.Line2D([], [], color="black", label=r"$q+\bar{q}$", linestyle="solid"),
        plt.Line2D([], [], color="black", label=r"$g$", linestyle="dashed"),
    ]
    labels2 = []
    for pT, color in zip(pThat, hp.colors_2):
        labels2.append(plt.Line2D([], [], color=color, label=pT))
        dat = Reader(
            basedir=base,
            formationTime=f"{delay}_delay",
            rate=rate,
            pTHat=pT,
            nattempts=5,
        )
        dat.read()

        for iax, ax in enumerate(axes):
            _ = dat.plot(
                itau=taus[iax],
                ax=axes[iax],
                color_gluon=color,
                color_quark=color,
                quarklines="solid",
                gluonlines="dashed",
            )

    for iax, ax in enumerate(axes):
        tau = 0.4 + taus[iax] * 0.4
        ax.text(
            0.1,
            0.1,
            r"$\tau=$" + f"{tau:0.1f} " + r"$\mathrm{fm}/c$",
            transform=ax.transAxes,
        )
    for ax in axes:
        ax.set_yscale("log")
    for ax in axes[2:]:
        ax.set_xlabel(r"$\zeta$")
    for ax in axes[::2]:
        ax.set_ylabel(r"Count/Event")

    axes[0].text(
        0.2,
        0.8,
        f"Rate: {hp.rates[rate]}, {delay} " + r"$\tau_{\mathrm{form.}}$",
        transform=axes[0].transAxes,
    )
    axes[1].legend(loc="upper right", handles=labels2, fontsize=28, title=r"$\hat{p}_T$, GeV")
    axes[-1].legend(loc="best", handles=labels1, fontsize=28)
    axes[0].set_ylim(bottom=1e-6)
    plt.show()


if __name__ == "__main__":
    main([25, 50, 100], int(argv[1]), argv[2], argv[3])
