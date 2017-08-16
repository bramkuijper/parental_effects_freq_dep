#!/usr/bin/env python3

import sys, re, os.path, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

data = pd.read_csv(filename, sep=";",skiprows = 12)

nrows = data.shape[0]

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
        data["generation"],data["meanpH"],'b',
        data["generation"],data["meanpD"],'r',
        linewidth=1)
plt.legend((r'$p_{\mathrm{H}}$',
            r'$p_{\mathrm{D}}$',
            ))

plt.ylabel(r'Prob. offspring is same as parent')
plt.ylim(0,1)

# add first subplot depicting % type 1 offspring
plt.subplot(num_rows,1,2)

plt.plot(
        data["generation"],data["freq_hawk"],'b',
        linewidth=1)

plt.ylabel(r'Proportion hawk')
plt.ylim(0,1)

plt.subplot(num_rows,1,3)

plt.plot(
        data["generation"],data["varpH"],'b',
        data["generation"],data["varpD"],'r',
        linewidth=1)

plt.ylabel(r'Variance')


plt.subplot(num_rows,1,4)

plt.plot(
        data["meanpD"],data["meanpH"],'b',
        linewidth=.5)

plt.ylim(0,1)
plt.xlim(0,1)
plt.ylabel(r'pH')

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
