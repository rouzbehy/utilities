import pandas as pd
import numpy as np
import os
import glob
import csv
from matplotlib.colors import CSS4_COLORS as css
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
my_rcParams={
    "text.usetex": True,
    "font.family": "Georgia",
    "font.size": 34,
    "lines.linewidth": 4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size" : 8,
    "ytick.major.size" : 8,
    "xtick.minor.size" : 6,
    "ytick.minor.size" : 6,
    "axes.spines.right": False,
    "axes.spines.top" : False,
    "legend.frameon":False
}
## Plotting function for experimental data:
# def plot_expr_data_on_axis(axis, data, marker,color=css['black'],face=css['dimgray']):
#     deltax = 0.5*(data["xhigh"]-data["xlow"])
#     errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
#                     for x, delx, y, yerrlow, yerrhigh in
#                     zip(data["x"], deltax, data["y"], -1*data["dy_syst-"], data["dy_syst+"])]
     
#     # Create patch collection with specified colour/alpha
#     pc = PatchCollection(errorboxes, facecolor=face, alpha=0.2)
#     axis.add_collection(pc)
#     axis.scatter(data["x"], data["y"], color=color, marker=marker, s=60)

#     axis.errorbar(data["x"], data["y"], xerr=deltax, 
#                                         yerr=[abs(data["dy_stat-"]), data["dy_stat+"]], 
#                                         fmt='none', 
#                                         color=color)
def plot_expr_data_on_axis(axis, data, marker,color=css['black'],face=css['dimgray'], factor=1):
    deltax = 0.5*(data["xhigh"]-data["xlow"])

    if 'dy_syst-' in data and 'dy_stat-' in data:
        errorboxes = [Rectangle((x-delx, y - yerrlow), width=2*delx, height=abs(yerrlow)+abs(yerrhigh))
                        for x, delx, y, yerrlow, yerrhigh in
                        zip(data["x"], deltax, factor*data["y"], factor*data["dy_syst+"], factor*-1*data["dy_syst-"])]

        # Create patch collection with specified colour/alpha
        pc = PatchCollection(errorboxes, facecolor=face, alpha=0.2)
        axis.add_collection(pc)
        #axis.scatter(data["x"], data["y"], color=color, marker=marker, s=150)
        axis.scatter(data["x"], factor*data["y"], color=color, marker=marker, s=60)
        axis.errorbar(data["x"], factor*data["y"], xerr=deltax, 
                                            yerr=[factor*-1*data["dy_stat-"], factor*data["dy_stat+"]], 
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
