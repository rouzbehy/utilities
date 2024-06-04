import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from typing import Union


class Reader:
    def __init__(
        self, basedir: str, formationTime: str, rate: int, pTHat: int, nattempts: int
    ) -> None:
        self._fname = (
            f"{basedir}/{formationTime}/rate_{rate}/pT_{pTHat}/"
            + "attempt_{iatt}/tau_{itau}_partons_gun_"
            + f"{pTHat}.csv"
        )
        self._attmp = nattempts
        self._specs = None

    def read(self):
        self._specs = {}
        for itau in range(0, 35):
            zmin, zmax, quarks, gluons, dquarks, dgluons = [[] for _ in range(6)]
            nevents, xsec = 1, 1
            for iatt in range(self._attmp):
                curr_file = self._fname.format(itau=itau, iatt=iatt)
                with open(curr_file) as f:
                    _, nevents, _, xsec, _, _ = f.readline().split(" ")
                nevents = float(nevents)
                xsec = float(xsec)

                tmp = np.loadtxt(
                    curr_file, comments="#", delimiter=",", skiprows=3, unpack=True
                )
                if len(zmin) == 0:
                    zmin, zmax = tmp[0], tmp[1]
                # for i, spec in enumerate([gluons, dgluons, quarks, dquarks]):
                #    # spec.append(xsec * tmp[2 + i] / (nevents))
                #    spec.append(tmp[2 + i] / (nevents))
                for i, spec in enumerate([gluons, quarks]):
                    spec.append(tmp[i + 2] / nevents)
            z = [0.5 * (x + y) for (x, y) in zip(zmin, zmax)]
            dz = [0.5 * (y - x) for (x, y) in zip(zmin, zmax)]
            gluon_spec = [sum(e) / len(e) for e in zip(*gluons)]
            quark_spec = [sum(e) / len(e) for e in zip(*quarks)]
            ## normalize?
            # const = sum([(x + y) * d for (x, y, d) in zip(gluon_spec, quark_spec, dz)])
            # gluon_spec = [e / const for e in gluon_spec]
            # quark_spec = [e / const for e in quark_spec]
            tmp = pd.DataFrame(
                data={"z": z, "dz": dz, "gluons": gluon_spec, "quarks": quark_spec}
            )
            self._specs[itau] = tmp[tmp.z > 0.1]

    def plot(
        self,
        itau: int,
        ax: Union[None, plt.axes],
        color_gluon: str,
        color_quark: str,
        quarklines: str,
        gluonlines: Union[None, str] = None,
    ) -> plt.axes:
        if itau not in self._specs:
            print(f"Oh no, {itau} is not in my specta!")

        spec = self._specs[itau]
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(
            spec["z"],
            spec["gluons"],
            color=color_gluon,
            linestyle=quarklines if gluonlines == None else gluonlines,
        )
        ax.plot(
            spec["z"],
            spec["quarks"],
            color=color_quark,
            linestyle=quarklines,
        )

        return ax

    def get_spec(self, itau):
        return self._specs[itau]


if __name__ == "__main__":
    reader = Reader(
        basedir="../data",
        formationTime="with_delay",
        rate=3,
        pTHat=100,
        nattempts=5,
    )
    reader.read()
    from helpers import my_rcParams

    plt.rcParams.update(my_rcParams)
    fig, ax = plt.subplots()

    _ = reader.plot(
        ax=ax, itau=0, color_gluon="red", color_quark="blue", linestyle="solid"
    )
    _ = reader.plot(
        ax=ax, itau=15, color_gluon="red", color_quark="blue", linestyle="dotted"
    )
    ax.set_yscale("log")

    plt.show()
