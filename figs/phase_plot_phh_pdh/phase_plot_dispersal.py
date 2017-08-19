#!/usr/bin/env python3

# plot a phaseplot of pH and pD for different levels of dispersal

import pandas as pd
import functools
import itertools
import re, string, sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm

plt.style.use('base')

rcParams['axes.labelsize'] = 15
rcParams['text.usetex'] = True
rcParams['font.family'] = 'sans-serif'

# see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]  

# initialize the figure
fig = plt.figure(figsize=(13.5,5))

widths = [ 1, 1, 1 ]
heights = [ 1 ]

gs = gridspec.GridSpec(
        len(heights), 
        len(widths), 
        width_ratios=widths, 
        height_ratios=heights)

data = pd.read_csv("summary_phaseplot_d.csv", sep=";")

def process_file(filename):

    f = open(filename)
    fl = f.readlines()
    f.close()

    parline = 0

    for i in reversed(range(0, len(fl))):

        if re.match("^\d+",fl[i]) is not None:
            parline = i
            break

    subdat = pd.read_csv(filename,sep=";",nrows=parline-1)

    return(subdat)


def block(
        row,
        col,
        xlabel=None,
        ylabel=None,
        title=None,
        xtick_labels=True,
        ytick_labels=False,
        selection_dict={}):

    # get the global dataset
    global data

    # dynamically make data selection
    # dependent on criteria specified in 
    # selection_criteria
    # https://pandas.pydata.org/pandas-docs/stable/cookbook.html

    # initialize list with data selection criteria
    selection_criteria = []
    
    # make data selection critera
    for key, val in selection_dict.items():
        selection_criteria += [data[key] == val]

    if len(selection_criteria) > 0:
        all_criteria = functools.reduce(lambda x,y: x & y, selection_criteria)

        subset = data[all_criteria]

    else:
        subset = data

    # select only bounderies
    subset = subset[
            (subset["phhinit"] < 0.1 ) 
            | (subset["phhinit"] > 0.9) 
            | (subset["pdhinit"] < 0.1)
            | (subset["pdhinit"] > 0.9)
            ]

    # initialize figure
    ax = plt.subplot(gs[row,col])

    ax.plot([0,1],
            [0,1],
            linestyle="dashed",
            color="grey",
            linewidth=1)

    ax.text(x=0.7,
            y=0.97,
            s="No parental effects",
            color="grey",
            rotation=45)

    # iterate over rows
    for index, row in subset.iterrows():


        # read the data
        filename = str(row["file"])

        # process the data file
        subdata = process_file(filename)

        ax.plot(subdata["pdh"],
                subdata["phh"],
#                color="#4cb5ff",
                color="#007aff",
#                color="#ff3000",
                linewidth=1)
    
    # iterate over rows and drow endpoints
    # and arrows
    for index, row in subset.iterrows():

        # read the data
        filename = str(row["file"])

        # process the data file
        subdata = process_file(filename)

        last_row = subdata.iloc[-1:]
        
        ax.plot(last_row["pdh"],
                last_row["phh"],
                marker="o",
                markerfacecolor="white",
                markeredgecolor="black",
                linewidth=1)

    # do axis labeling
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.set_ylim((-0.05,1.05))
    ax.set_xlim((-0.05,1.05))
#    ax.spines['right'].set_visible(False)
#    ax.spines['top'].set_visible(False)
#    ax.yaxis.set_ticks_position("left")
#    ax.xaxis.set_ticks_position("bottom")

    if type(ylabel) == type("string"):
        ax.set_ylabel(ylabel=ylabel)
    
    if type(xlabel) == type("string"):
        ax.set_xlabel(xlabel=xlabel)

    if not xtick_labels:
        ax.xaxis.set_ticklabels([])

    print(ytick_labels)

    if not ytick_labels:
        ax.yaxis.set_ticklabels([])
    
    ax.set_title(loc="left", 
            label=string.ascii_uppercase[col], position=(0.05,1.02))
    
    if title is not None: 
        ax.set_title(
                label=title, position=(0.5,1.02))

block(
        row=0,
        col=0,
        ytick_labels=True,
        title="Limited dispersal: $d = 0.1$",
        selection_dict={"d": 0.1},
        ylabel=r"Proportion of hawks by hawk parents, $p_{H}$")

block(
        row=0,
        col=1,
        ytick_labels=False,
        title="Modest dispersal: $d = 0.5$",
        selection_dict={"d": 0.5},
        xlabel=r"Proportion of hawks by dove parents, $p_{D}$")

block(
        row=0,
        col=2,
        ytick_labels=False,
        title="Well-mixed population: $d = 1.0$",
        selection_dict={"d": 1.0})

plt.tight_layout()

format = "pdf"
plt.savefig("graph_phase_plot_dispersal." + format, format=format)
