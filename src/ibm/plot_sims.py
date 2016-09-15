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

data = pd.read_csv(filename, sep=";",skiprows = 11)

nrows = data.shape[0]

#row_each = math.floor(float(nrows)/1000)
#
## make dataset shorter so that we don't plot megabytes
#data = data.iloc[range(0,nrows,int(row_each)),:]


# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 3

# names of columns in dataframe
colnames = list(data.columns.values)

# add first subplot depicting % type 1 offspring
plt.subplot(num_rows,1,1)

plt.plot(
        data["generation"],data["meanphh1"],'b',
        data["generation"],data["meanphh2"],'r',
        data["generation"],data["meanpdh0"],'k',
        data["generation"],data["meanpdh1"],'c',
        linewidth=1)
plt.legend((r'$p_{\mathrm{h}\rightarrow\mathrm{h},1}$',
            r'$p_{\mathrm{h}\rightarrow\mathrm{h},2}$',
            r'$p_{\mathrm{d}\rightarrow\mathrm{h},0}$',
            r'$p_{\mathrm{d}\rightarrow\mathrm{h},1}$'
            ))

plt.ylabel(r'Prob. offspring is hawk')
plt.ylim(0,1)

# add 2nd subplot depicting patch frequencies
plt.subplot(num_rows,1,2)

# count number of patch frequencies
patch_freq_names = []
for name in colnames:
    if re.match("f_",name) != None:
        patch_freq_names.append(name)

colors = iter(cm.rainbow(np.linspace(0,1,len(patch_freq_names))))

for name in patch_freq_names:
    plt.plot(data["generation"],data[name],color=next(colors),linewidth=1)

#plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'patch frequencies')
plt.legend(tuple(patch_freq_names))
plt.ylim(0,1)

# add first subplot depicting % type 1 offspring
plt.subplot(num_rows,1,3)

plt.plot(
        data["generation"],data["varphh1"],'b',
        data["generation"],data["varphh2"],'r',
        data["generation"],data["varpdh0"],'k',
        data["generation"],data["varpdh1"],'c',
        linewidth=1)

plt.legend((
        r'$\mathrm{var}\left(p_{\mathrm{h}\rightarrow\mathrm{h},1}\right)$',
        r'$\mathrm{var}\left(p_{\mathrm{h}\rightarrow\mathrm{h},2}\right)$',
        r'$\mathrm{var}\left(p_{\mathrm{d}\rightarrow\mathrm{h},0}\right)$',
        r'$\mathrm{var}\left(p_{\mathrm{d}\rightarrow\mathrm{h},1}\right)$',
        ))

plt.ylabel(r'Variance')


graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
