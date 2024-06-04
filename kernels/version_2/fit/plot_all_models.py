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
#mpl.use('Qt5Agg')
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 23,
    "lines.linewidth": 2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size" : 8,
    "ytick.major.size" : 8,
    "xtick.minor.size" : 4,
    "ytick.minor.size" : 4,
    "axes.spines.top": True,
    "axes.spines.right":True,
    "legend.frameon":False,
    "figure.figsize" : (16,9)
}
plt.rcParams.update(my_rcParams)


if __name__ == '__main__':

    fig1, axes = plt.subplots(1,3, gridspec_kw={'left':0.05, 'right':0.93,'top':0.95, 'bottom':0.1}, sharex=True, sharey=True, figsize=(8,8))
    ## Points to predict and construct the surfaces
    kapr = np.linspace(1, 15, 500, endpoint=True)
    kape = np.linspace(1, 15, 500, endpoint=True)
    kapr_inv = 1./kapr
    kape_inv = 1./kape
    #z_min, z_max = 1, 2.5
    for iax, rset in enumerate([1,2,3]):
        ax = axes[iax]
        rate_set_name = util.rate_set_names[rset]
        model = 0
        with open(f"../../fit_results/chosen_mode_{rset}.pickle", 'rb') as f:
            model = pickle.load(f)

        prediction_grid_x, prediction_grid_y = np.meshgrid(kapr_inv, kape_inv)
        plotting_grid_x , plotting_grid_y    = np.meshgrid(kapr, kape)
        GRID_prediction                      = np.stack([prediction_grid_x.ravel(), prediction_grid_y.ravel()]).T
        GRID_plotting                        = np.stack([plotting_grid_x.ravel(), plotting_grid_y.ravel()]).T
        ## The model is trained on 1/kappa_{r,e} and (Chi_squared)
        z               = model.predict(GRID_prediction)
        z               = np.exp(z)
        image           = np.reshape(z, prediction_grid_x.shape)
        minimum_chi_sq  = min(z)
        indx_opt        = np.argmin(image)

        fitted_kappa_vals = GRID_plotting[indx_opt]
        (kr, ke) = fitted_kappa_vals
        print(util.bcolors.BOLD+util.bcolors.OKGREEN+\
              f"Rateset {rset} the minimum chi squared predicted is {minimum_chi_sq}\n" +\
              util.bcolors.WARNING+f" with (kappa_r, kappa_e) = ({kr:0.2f},{ke:0.2f})" + util.bcolors.ENDC)

        #z_min, z_max= minimum_chi_sq-0.1, minimum_chi_sq+2.
        z_min, z_max = 0.8, 4
        dz = 0.3
        levels = [z_min + i*dz for i in range(20) if z_min+i*dz < z_max]
        im      = ax.pcolormesh(plotting_grid_x, plotting_grid_y, image, 
                                    cmap='inferno_r', vmin=z_min, vmax=z_max, rasterized=True)
        #if iax == 2:
        divider = make_axes_locatable(ax)
        cax     = divider.append_axes('right', size='5%', pad=0.03)
        cbar    = fig1.colorbar(im, cax=cax, orientation='vertical')
        if iax == 2:
            cbar.set_label(r'$\chi^2/$n.d.f')
        if iax == 0:
            ax.set_ylabel('$\kappa_e$', fontsize=30)
        ax.set_xlabel('$\kappa_r$', fontsize=30)
        ax.set_title(f'{rate_set_name}')

    #fig1.savefig(f"../../chisquared_work/Heat_Map_of_Kappas_rset_{rset}.png",dpi=100)
    plt.show()

## old cmap :cividis_r 