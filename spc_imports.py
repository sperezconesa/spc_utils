"""
This module is just a list of my usual imports and mpl configuration. They can
be called using:
from spc_imports import *
set_up_plt()
"""
import os
# import datreant as dtr
import pickle as pickle
import re as re
import subprocess as sp

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm


def set_up_plt():
    """ Set up matplotlib like SPC likes it."""
    plt.rcParams["text.usetex"] = False
    plt.rcParams["image.cmap"] = "RdYlBu_r"  # coolwarm viridis
    plt.rcParams["axes.labelsize"] = 13
    plt.rcParams["axes.titlesize"] = 14
    plt.rcParams["legend.fontsize"] = 15
    ticklabelsize = 14
    plt.rcParams["xtick.labelsize"] = ticklabelsize
    plt.rcParams["ytick.labelsize"] = ticklabelsize
    plt.rcParams["text.latex.preamble"] = [
        r"\usepackage{textgreek}",
        r"\usepackage{siunitx}",  # Remember to use r'string'
        r"\usepackage{mhchem}",
    ]
    # PUT GRID
    plt.rcParams["axes.facecolor"] = "white"  # "#EAEAF2"
    plt.rcParams["axes.edgecolor"] = "black"  # "box" color
    plt.rcParams["grid.color"] = "white"
    plt.rcParams["grid.linestyle"] = "-"
    plt.rcParams["grid.linewidth"] = 2
    plt.rcParams["axes.grid"] = True
    plt.rcParams["axes.axisbelow"] = True
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["lines.solid_capstyle"] = "round"
    ticks_bool = True
    plt.rcParams["ytick.right"] = ticks_bool
    plt.rcParams["ytick.left"] = ticks_bool
    plt.rcParams["xtick.bottom"] = ticks_bool
    plt.rcParams["xtick.top"] = ticks_bool
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"

    return
