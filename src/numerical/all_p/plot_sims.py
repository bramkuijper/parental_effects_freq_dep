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

nrow = len(fl)-11


data = pd.read_csv(filename, sep=";", nrows=nrow)


data["f_2"] = pd.Series(list(itertools.repeat(np.nan,len(data.index))))

data["f_2"] = 1.0 - data["f_0"] - data["f_1"]

#row_each = math.floor(float(nrows)/1000)
#
## make dataset shorter so that we don't plot megabytes
#data = data.iloc[range(0,nrows,int(row_each)),:]


# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 4

# names of columns in dataframe
colnames = list(data.columns.values)

# add first subplot depicting % type 1 offspring
plt.subplot(num_rows,1,1)

plt.plot(
        data["time"],data["phh1"],'b',
        data["time"],data["phh2"],'r',
        data["time"],data["pdh1"],'c',
        data["time"],data["pdh0"],'k',
        )

plt.ylabel(r'Prob. offspring is hawk')
plt.legend((r'$p_{\mathrm{h}\rightarrow\mathrm{h},1}$',
                r'$p_{\mathrm{h}\rightarrow\mathrm{h},2}$',
                r'$p_{\mathrm{d}\rightarrow\mathrm{h},1}$',
                r'$p_{\mathrm{d}\rightarrow\mathrm{h},0}$'),bbox_to_anchor=(1.1,1.0))
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



# add 3rd subplot depicting relatedness
plt.subplot(num_rows,1,3)

plt.plot(data["time"],data["rdd"],'y',
        data["time"],data["rhd"],'g',
        data["time"],data["rhh"],'k',
        linewidth=1)
plt.legend((r'$r_{\mathrm{dd}}$',r'$r_{\mathrm{hd}}$',r'$r_{\mathrm{hh}}$'))

plt.ylabel(r'Relatedness')
plt.ylim(-0.05,1.05)


# add 4th subplot depicting reprovals
plt.subplot(num_rows,1,4)

plt.plot(data["time"],data["vh_1"],'y',
        data["time"],data["vh_2"],'g',
        data["time"],data["vd_0"],'k',
        data["time"],data["vd_1"],'m',
        linewidth=1)
plt.legend((r'$v_{\mathrm{h},1}$',r'$v_{\mathrm{h},2}$',r'$v_{\mathrm{d},0}$', r'$v_{\mathrm{d},1}$'))

plt.ylabel(r'Reproductive value')

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
