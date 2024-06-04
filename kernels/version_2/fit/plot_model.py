## IMPORTS:
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import CSS4_COLORS as css
from numpy import array
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from sys import argv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import util

# mpl.use('Qt5Agg')
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
plt.rcParams.update(my_rcParams)


if __name__ == "__main__":

    rset = int(argv[1])
    rate_set_name = util.rate_set_names[rset]
    model = 0
    with open(f"../../chisquared_work/chosen_mode_{rset}.pickle", "rb") as f:
        model = pickle.load(f)

    ## Points to predict and construct the surfaces
    kapr = np.linspace(1, 15, 1000, endpoint=True)
    kape = np.linspace(1, 15, 1000, endpoint=True)
    kapr_inv = 1.0 / kapr
    kape_inv = 1.0 / kape

    prediction_grid_x, prediction_grid_y = np.meshgrid(kapr_inv, kape_inv)
    plotting_grid_x, plotting_grid_y = np.meshgrid(kapr, kape)

    GRID_prediction = np.stack([prediction_grid_x.ravel(), prediction_grid_y.ravel()]).T
    GRID_plotting = np.stack([plotting_grid_x.ravel(), plotting_grid_y.ravel()]).T

    fig1, ax1 = plt.subplots(1, 1, figsize=(9, 9))

    ## The model is trained on 1/kappa_{r,e} and (Chi_squared)
    z = model.predict(GRID_prediction)
    z = np.exp(z)
    image = np.reshape(z, prediction_grid_x.shape)
    minimum_chi_sq = min(z)
    indx_opt = np.argmin(image)

    fitted_kappa_vals = GRID_plotting[indx_opt]
    print(
        util.bcolors.BOLD
        + util.bcolors.OKGREEN
        + f"Rateset {rset} the minimum chi squared predicted is {minimum_chi_sq}\n"
        + util.bcolors.WARNING
        + f" with (kappa_r, kappa_e) = {fitted_kappa_vals}"
        + util.bcolors.ENDC
    )

    z_min, z_max = minimum_chi_sq - 0.1, minimum_chi_sq + 1.0
    dz = 0.3
    levels = [z_min + i * dz for i in range(20) if z_min + i * dz < z_max]

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    im = ax1.pcolormesh(
        plotting_grid_x,
        plotting_grid_y,
        image,
        cmap="magma_r",
        vmin=z_min,
        vmax=z_max,
        rasterized=True,
    )
    cbar = fig1.colorbar(im, cax=cax, orientation="vertical")
    cbar.set_label(r"$\chi^2/$n.d.f")

    # cont   = ax1.contour(plotting_grid_x, plotting_grid_y, image,
    #                        cmap='Reds_r', vmin=z_min, vmax=z_max, levels=levels)
    # ax1.clabel(cont, inline=True, fontsize=12)
    ax1.set_ylabel("$\kappa_e$", fontsize=30)
    ax1.set_xlabel("$\kappa_r$", fontsize=30)
    ax1.set_title(f"{rate_set_name}")

    fig1.savefig(f"../../chisquared_work/Heat_Map_of_Kappas_rset_{rset}.png", dpi=100)
    plt.show()

## old cmap :cividis_r
