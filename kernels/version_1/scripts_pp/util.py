"""
    utility functions to declutter the main
    scripts.
    provide functionality for:
        - reading experimental files
        - reading in the theory calculation
        - performing integration if necessary (for v_2)
"""
from typing import Tuple
import pandas as pd
import numpy as np
from scipy.integrate import trapezoid, fixed_quad
from scipy.interpolate import InterpolatedUnivariateSpline as interpolator
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.colors import CSS4_COLORS as css
from dictionaries import exag_factor

## Plotting function for experimental data:
def plot_expr_data_on_axis(axis, data, marker,color=css['black'],face=css['dimgray'], factor=1):
    deltax = 0.5*(data["xhigh"]-data["xlow"])

    if 'dy_syst-' in data and 'dy_stat-' in data:
        errorboxes = [Rectangle((x-delx, y - abs(yerrlow)), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                        for x, delx, y, yerrlow, yerrhigh in
                        zip(data["x"], deltax, factor*data["y"], factor*data["dy_syst+"], factor*-1*data["dy_syst-"])]

        # Create patch collection with specified colour/alpha
        pc = PatchCollection(errorboxes, facecolor=face, alpha=0.2)
        axis.add_collection(pc)
        #axis.scatter(data["x"], data["y"], color=color, marker=marker, s=150)
        axis.scatter(data["x"], factor*data["y"], color=color, marker=marker, s=60)
        axis.errorbar(data["x"], factor*data["y"], xerr=deltax, 
                                            yerr=[factor*abs(data["dy_stat-"]), factor*abs(data["dy_stat+"])], 
                                            fmt='none', 
                                            color=color)
    else:
        errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                        for x, delx, y, yerrlow, yerrhigh in
                        zip(data["x"], deltax, factor*data["y"], factor*data["dy_+"], factor*-1*data["dy_-"])]

        # Create patch collection with specified colour/alpha
        pc = PatchCollection(errorboxes, facecolor=face, alpha=0.2)
        axis.add_collection(pc)
        #axis.scatter(data["x"], data["y"], color=color, marker=marker, s=150)
        axis.scatter(data["x"], factor*data["y"], color=color, marker=marker, s=60)
        #axis.errorbar(data["x"], factor*data["y"], xerr=deltax, 
        #                                    yerr=[factor*-1*data["dy_stat-"], factor*data["dy_stat+"]], 
        #                                    fmt='none', 
        #                                    color=color)