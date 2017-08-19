#!/usr/bin/env python3

# plot the battleground between offspring and parent strategies

import pandas as pd
import itertools
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm
import matplotlib.colors as mpl_colors 

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
fig = plt.figure(figsize=(5,5))

widths = [ 1, ]
heights = [ 1 ]

gs = gridspec.GridSpec(
        nrows=len(heights),
        ncols=len(widths),
        width_ratios=widths,
        height_ratios=heights)

# read in the data
data = pd.read_csv("summary_sims.csv", sep=";")

# get a su6set to do some line drawing
subset = data[
                (data["c"] == 3.0) &
                ((data["init_pH"] >= 0.95 ) |
                (data["init_pD"] <= 0.05 ) |
                (data["init_pD"] >= 0.95 ))
                ]

print(subset.shape)

ax = plt.subplot(gs[0,0])

# make a list of colors
colors = iter(cm.prism(np.linspace(0, 1, subset.shape[0])))

# split the iterator as we need it to plot the lines but also 
# to plot the points
colors1, colors2 = itertools.tee(colors,2)

# plot lines
ax.plot([0,1],
        [0,1],
        color="grey",
        linestyle="dashed",
        linewidth=1)

for index, row in subset.iterrows():

    data_sub = pd.read_csv(str(row["file"]),sep=";",skiprows=12)

    ax.plot(1-data_sub["meanpD"],data_sub["meanpH"], color=next(colors1))

# reset color iterator and use it again to plot endpoints

for index, row in subset.iterrows():
    
    data_sub = pd.read_csv(str(row["file"]),sep=";",skiprows=12)
    # plot last point
    final_row = data_sub.iloc[-1:]

    the_color = mpl_colors.rgb2hex(next(colors2))

    ax.plot(float(1-final_row["meanpD"]),
            float(final_row["meanpH"]),
            marker="s",
            markerfacecolor=the_color,
            markeredgecolor="black")

ax.tick_params(axis="both", which="both", direction="in")
ax.text(x=0.6,y=0.85, 
        color="grey",
        s="No parental effects",
        rotation=45)


ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.yaxis.set_minor_locator(AutoMinorLocator(4))
ax.set_ylim((-0.05,1.05))
ax.set_xlim((-0.05,1.05))

ax.set_ylabel(ylabel="Proportion hawks from hawk parents, $p_{H}$")
ax.set_xlabel(xlabel="Proportion hawks from dove parents, $p_{D}$")

plt.tight_layout()

format = "pdf"
plt.savefig("phaseplot_well_mixed." + format, format=format)

