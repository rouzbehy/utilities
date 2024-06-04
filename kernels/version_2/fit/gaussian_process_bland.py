## IMPORTS:
from numpy.random import RandomState
from numpy.random import normal
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.colors import CSS4_COLORS as css
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning
simplefilter("ignore", category=ConvergenceWarning)
from sys import argv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import util
from numpy import log
import pickle
import matplotlib as mpl

my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 25,
    "lines.linewidth": 4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size" : 12,
    "ytick.major.size" : 12,
    "xtick.minor.size" : 6,
    "ytick.minor.size" : 6,
    "axes.spines.right": False,
    "axes.spines.top" : False,
    "legend.frameon":False
}
plt.rcParams.update(my_rcParams)

mpl.rc('image', cmap='tab10')
def plot_gpr_samples(gpr_model, n_samples, ax):
    """Plot samples drawn from the Gaussian process model.

    If the Gaussian process model is not trained then the drawn samples are
    drawn from the prior distribution. Otherwise, the samples are drawn from
    the posterior distribution. Be aware that a sample here corresponds to a
    function.

    Parameters
    ----------
    gpr_model : `GaussianProcessRegressor`
        A :class:`~sklearn.gaussian_process.GaussianProcessRegressor` model.
    n_samples : int
        The number of samples to draw from the Gaussian process distribution.
    ax : matplotlib axis
        The matplotlib axis where to plot the samples.
    """
    x = np.linspace(0, 5, 50)
    X = x.reshape(-1, 1)

    y_mean, y_std = gpr_model.predict(X, return_std=True)
    y_samples = gpr_model.sample_y(X, n_samples)

    for idx, single_prior in enumerate(y_samples.T):
        ax.plot(
            x,
            single_prior,
            linestyle="--",
            alpha=0.7,
            label=f"Sampled function No.{idx + 1}",
        )
    ax.plot(x, y_mean, color="black", label="Mean")
    ax.fill_between(
        x,
        y_mean - y_std,
        y_mean + y_std,
        alpha=0.1,
        color="black",
        label=r"$\pm$ 1 std. dev.",
    )
    #ax.set_xlabel("x")
    #ax.set_ylabel("y")
    ax.set_ylim([-4, 4])

if __name__=='__main__':

    fig, axes = plt.subplots(2,1, sharex=True)
    rs = RandomState()
    kernel = 1.0 * RBF(length_scale=1.0, length_scale_bounds=(1e-1, 10.0))
    gpr = GaussianProcessRegressor(kernel=kernel, random_state=0, alpha=0.001)

    plot_gpr_samples(gpr, 5, axes[0])

    rng = RandomState(4)
    X_train = rng.uniform(0, 5, 10).reshape(-1, 1)
    y_train = np.sin((X_train[:, 0] - 2.5) ** 2)
    axes[1].scatter(X_train[:, 0], y_train, color="red", zorder=10, label="Observations", s=35)
    gpr.fit(X_train, y_train)
    plot_gpr_samples(gpr, 5, axes[1])
    axes[1].legend(bbox_to_anchor=(1.01, 1.5), loc="upper left")
    
    axes[0].text(0.1, 0.05, s="Samples from prior distribution",transform=axes[0].transAxes)
    axes[1].text(0.1, 0.05, s="Samples from posterior distribution",transform=axes[1].transAxes)
    
    fig.suptitle("Radial Basis Function Kernel")
    axes[1].set_xlabel('x')
    for ax in axes:
        ax.set_ylabel('y')
    
    plt.show()
