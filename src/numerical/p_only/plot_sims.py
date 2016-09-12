#!/usr/bin/env python3

import sys, re, os.path, itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

# first figure out last line of dataset
f = open(filename)
fl = f.readlines()
f.close()

nrow = len(fl)-8


data = pd.read_csv(filename, sep=";", nrows=nrow)


data["f_2"] = pd.Series(list(itertools.repeat(np.nan,len(data.index))))

data["f_2"] = 1.0 - data["f_0"] - data["f_1"]

#row_each = math.floor(float(nrows)/1000)
#
## make dataset shorter so that we don't plot megabytes
#data = data.iloc[range(0,nrows,int(row_each)),:]


# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 2

# names of columns in dataframe
colnames = list(data.columns.values)

# add first subplot depicting % type 1 offspring
plt.subplot(num_rows,1,1)

plt.plot(data["time"],data["p"],'b',linewidth=1)

plt.ylabel(r'Prob. offspring is hawk')
plt.ylim(-0.05,1.05)

# add 2nd subplot depicting patch frequencies
plt.subplot(num_rows,1,2)

plt.plot(data["time"],data["f_0"],'y',
        data["time"],data["f_1"],'g',
        data["time"],data["f_2"],'magenta',
        linewidth=1)
plt.legend((r'$f_{0}$',r'$f_{1}$',r'$f_{2}$'))

plt.ylabel(r'Patch freqs')
plt.ylim(-0.05,1.05)


graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
