import matplotlib.pyplot as plt
from reader import Reader
import helpers as hp
from sys import argv

plt.rcParams.update(hp.my_rcParams)


def main(pThat: int, delay: str, base: str):

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
    labels_1 = [
        plt.Line2D([], [], color="black", label=r"$q+\bar{q}$", linestyle="solid"),
        plt.Line2D([], [], color="black", label=r"$g$", linestyle="dashed"),
    ]
    labels_2 = []
    for rate, color in zip([1, 2, 3], hp.colors_2):
        labels_2.append(plt.Line2D([], [], color=color, label=hp.rates[rate]))
        dat = Reader(
            basedir=base,
            formationTime=f"{delay}_delay",
            rate=rate,
            pTHat=pThat,
            nattempts=5,
        )
        dat.read()
        axes[0].text(
            0.2,
            0.8,
            r"$\hat{p}_{T}$=" + f"{pThat} GeV",
            transform=axes[0].transAxes,
        )
        for iax, ax in enumerate(axes):
            _ = dat.plot(
                itau=taus[iax],
                ax=axes[iax],
                color_gluon=color,
                color_quark=color,
                quarklines="solid",
                gluonlines="dashed",
            )
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

    axes[3].legend(loc="upper right", handles=labels_1, fontsize=28)
    axes[1].legend(loc="upper right", handles=labels_2, fontsize=28)
    plt.show()


if __name__ == "__main__":
    main(pThat=int(argv[1]), delay=argv[2], base=argv[3])
