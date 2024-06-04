import matplotlib.pyplot as plt
from reader import Reader
import helpers as hp

plt.rcParams.update(hp.my_rcParams)


def main(pThat: int, rate: int, base: str):

    with_form_time = Reader(
        basedir=base,
        formationTime="with_delay",
        rate=rate,
        pTHat=pThat,
        nattempts=5,
    )
    without_form_time = Reader(
        basedir=base,
        formationTime="without_delay",
        rate=rate,
        pTHat=pThat,
        nattempts=5,
    )
    with_form_time.read()
    without_form_time.read()

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
    cglue, cquark = hp.colors_2[0:2]
    taus = hp.taus
    axes[0].text(
        0.2,
        0.75,
        r"$\hat{p}_{T}$=" + f"{pThat} GeV" + "\n" + f"Rate set: {hp.rates[rate]}",
        transform=axes[0].transAxes,
    )
    lstyle_with_form, lstyle_without_form = "solid", "dashed"
    for iax, ax in enumerate(axes):

        _ = without_form_time.plot(
            itau=taus[iax],
            ax=axes[iax],
            color_gluon=cglue,
            color_quark=cquark,
            quarklines=lstyle_without_form,
        )
        _ = with_form_time.plot(
            itau=taus[iax],
            ax=axes[iax],
            color_gluon=cglue,
            color_quark=cquark,
            quarklines=lstyle_with_form,
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
        ax.set_ylabel("Count/Event")
    labels = [
        plt.Line2D([], [], color=cglue, label=r"$g$"),
        plt.Line2D([], [], color=cquark, label=r"$q+\bar{q}$"),
        plt.Line2D(
            [],
            [],
            color="black",
            linestyle=lstyle_without_form,
            label=r"without $\tau_\mathrm{form.}$",
        ),
        plt.Line2D(
            [],
            [],
            color="black",
            linestyle=lstyle_with_form,
            label=r"with $\tau_\mathrm{form.}$",
        ),
    ]
    axes[-1].legend(loc="upper right", handles=labels, fontsize=25, ncol=2)
    axes[0].set_ylim(bottom=1e-6)
    plt.show()


if __name__ == "__main__":
    from sys import argv

    pt, rate = int(argv[1]), int(argv[2])
    main(pt, rate, argv[3])
