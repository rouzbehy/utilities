# from typing import Tuple
# import pandas as pd
# import numpy as np
# from scipy.integrate import trapezoid, fixed_quad
from sklearn.model_selection import train_test_split
from sklearn.gaussian_process import GaussianProcessRegressor

# from scipy.interpolate import InterpolatedUnivariateSpline as interpolator
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.colors import CSS4_COLORS as css
from numpy.random import RandomState
import numpy as np
import pandas as pd

## Dictionaries and the like:
rate_set_names = {
    1: "LO Collision Kernel",
    2: "NLO Collision Kernel",
    3: "NP Collision Kernel",
}
rset_shorts = {1: "LO", 2: "NLO", 3: "NP"}
jet_cone_radii = {"0p2": 0.2, "0p3": 0.3, "0p4": 0.4}
pt_chgd_low, pt_chgd_high = 45, 160
pt_jet_low, pt_jet_high = 68, 320


kappas = {
    0: (1.0, 1.0),
    1: (1.0, 15.0),
    2: (15.0, 1.0),
    3: (15.0, 15.0),
    4: (14.4, 8.2),
    5: (7.4, 1.2),
    6: (10.9, 11.7),
    7: (3.9, 4.7),
    8: (5.6, 13.4),
    9: (12.6, 6.4),
    10: (9.1, 9.9),
    11: (2.1, 2.9),
    12: (3.0, 12.6),
    13: (10.0, 5.6),
    14: (13.5, 9.1),
    15: (6.5, 2.1),
    16: (4.7, 10.8),
    17: (11.7, 3.8),
    18: (8.2, 14.3),
    19: (1.2, 7.3),
    20: (1.7, 14.8),
    21: (8.7, 7.8),
    22: (12.2, 11.3),
    23: (5.2, 4.3),
    24: (6.9, 9.5),
    25: (13.9, 2.5),
    26: (10.4, 13.0),
    27: (3.4, 6.0),
    28: (2.5, 10.4),
    29: (9.5, 3.4),
    30: (13.0, 13.9),
    31: (6.0, 6.9),
    32: (4.3, 12.1),
    33: (11.3, 5.1),
    34: (7.8, 8.6),
    35: (14.8, 1.6),
    36: (1.0, 11.9),
    37: (8.0, 4.9),
    38: (11.5, 8.4),
    39: (4.5, 1.4),
    40: (6.3, 10.2),
    41: (13.3, 3.2),
    42: (9.8, 13.7),
    43: (2.8, 6.7),
    44: (3.6, 9.3),
    45: (10.6, 2.3),
    46: (14.1, 12.8),
    47: (7.1, 5.8),
    48: (5.4, 14.5),
    49: (12.4, 7.5),
}


class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


my_rcParams = {
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 23,
    "lines.linewidth": 2,
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


## Plotting function for experimental data:
def plot_expr_data_on_axis(
    axis,
    data,
    marker,
    color=css["black"],
    face=css["dimgray"],
    markersize=90,
    alpha=0.3,
):
    deltax = 0.5 * (data["xhigh"] - data["xlow"])
    errorboxes = [
        Rectangle(
            (x - delx, y - yerrlow), width=2 * delx, height=abs(yerrlow) + abs(yerrhigh)
        )
        for x, delx, y, yerrlow, yerrhigh in zip(
            data["x"], deltax, data["y"], data["dy_syst+"], -1 * data["dy_syst-"]
        )
    ]

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=face, alpha=alpha, zorder=20)
    axis.add_collection(pc)
    axis.scatter(
        data["x"], data["y"], color=color, marker=marker, s=markersize, zorder=20.0
    )
    axis.errorbar(
        data["x"],
        data["y"],
        xerr=deltax,
        yerr=[abs(data["dy_stat-"]), data["dy_stat+"]],
        fmt="none",
        color=color,
        zorder=2.0,
    )


def perform_fit(
    kernel,
    data_X,
    data_Y,
    num_models_to_run=1000,
    do_normalize_y=True,
    do_n_restarts_optimizer=10,
):
    lowest_err = 1e12
    chosen_model = 0
    rs = RandomState()
    regressor = GaussianProcessRegressor(
        kernel=kernel,
        normalize_y=do_normalize_y,
        n_restarts_optimizer=do_n_restarts_optimizer,
        random_state=rs,
        alpha=0.5,
    )
    for i in range(num_models_to_run):
        X_train, X_test, Y_train, Y_test = train_test_split(
            data_X, data_Y, random_state=rs, test_size=0.05
        )
        regressor.fit(X_train, Y_train)
        mean_predictions = regressor.predict(X_test)
        errors = [(p - a) ** 2 for p, a in zip(mean_predictions, Y_test)]
        curr_err = sum(errors) / len(errors)
        if curr_err < lowest_err:
            print(
                "\t"
                + bcolors.UNDERLINE
                + bcolors.OKGREEN
                + f"model {i}: {curr_err:0.6e}"
                + bcolors.ENDC
            )
            index = i
            lowest_err = curr_err
            chosen_model = regressor
    return index, lowest_err, chosen_model


def read_data(rate_set_id):
    fname = "../../calcs/chisq_values/rset_{rset}.csv"
    tmp = pd.read_csv(fname.format(rset=rate_set_id), comment="#")

    tmp["1/kappa_e"] = 1 / tmp["kappa_e"]
    tmp["1/kappa_r"] = 1 / tmp["kappa_r"]
    tmp["log_chisq"] = np.log(tmp["chi_sq/ndf"])
    return tmp


def get_charged_hadron_spec(rset, kset, species, pTcutLow, pTcutHigh):
    fname = "../../calcs/"
    if species == "pp":
        fname += f"pp_2p76/chgd_spec.csv"
    else:
        fname += (
            f"pbpb_2p76/param_fit_runs/rset_{rset}/kappaset_{kset}/hadron_spectra.csv"
        )
    data = pd.read_csv(fname, comment="#")
    data["pT"] = 0.5 * (data["pTmax"] + data["pTmin"])
    data = data[data["pT"].between(pTcutLow, pTcutHigh)]
    return data


def get_jet_spec(rset, kset, species, pTcutLow, pTcutHigh):
    fname = "../../calcs/"
    if species == "pp":
        fname += f"pp_2p76/jet_spec.csv"
    else:
        fname += f"pbpb_2p76/param_fit_runs/rset_{rset}/kappaset_{kset}/jet_spectra.csv"
    data = pd.read_csv(fname, comment="#")
    data["pT"] = 0.5 * (data["pTmax"] + data["pTmin"])
    data = data[data["pT"].between(pTcutLow, pTcutHigh)]
    return data


def get_jet_shape(rset, kset, species):
    fname = "../../calcs/"
    if species == "pp":
        fname += f"pp_2p76/jet_shape.csv"
    else:
        fname += f"pbpb_2p76/param_fit_runs/rset_{rset}/kappaset_{kset}/jet_shape.csv"
    njet = 1
    with open(fname, "r") as f:
        line = f.readline()
        line = line.split(" ")[-1]
        njet = float(line)
    df = pd.read_csv(fname, comment="#")
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
    return pd.DataFrame(data=retrn_dict)


def plot_theory_on_axis_boxes(
    axis, df, yname, xname, color, marker, s, alpha=0.3, zorder=0
):
    errorboxes = [
        Rectangle(
            (x - delx, y - yerrlow), width=2 * delx, height=abs(yerrlow) + abs(yerrhigh)
        )
        for x, delx, y, yerrlow, yerrhigh in zip(
            df[xname],
            df["d" + xname],
            df[yname],
            df["d" + yname],
            -1 * df["d" + yname],
        )
    ]

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=color, alpha=alpha, zorder=zorder)
    axis.add_collection(pc)
    axis.scatter(df[xname], df[yname], color=color, marker=marker, s=s, zorder=zorder)
