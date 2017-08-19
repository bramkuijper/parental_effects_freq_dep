#!/usr/bin/env python3

import sys, re, os.path, itertools, six, subprocess
import matplotlib
matplotlib.use('Agg')
matplotlib.use('pgf')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

import pandas as pd

from pylab import *


rcParams['axes.labelsize'] = 15
#rcParams['font.family'] = 'sans-serif'

# see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
#pgf_with_custom_preamble = {
#    "font.family": 'sans-serif', # use serif/main font for text elements
#    "text.usetex": True,    # use inline math for ticks
#    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
#    "pgf.preamble": [
#         "\\usepackage{units}",         # load additional packages
#         "\\usepackage{metalogo}",
#         "\\usepackage{unicode-math}",  # unicode math setup
#         r"\setmathfont{xits-math.otf}",
#         r"\setmainfont{MyriadPro-Regular.otf}", # serif font via preamble
#         ]
#}

#pgf_with_custom_preamble = {
#    "font.family": "serif", # use serif/main font for text elements
#    "text.usetex": True,    # use inline math for ticks
#    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
#    "pgf.preamble": [
#         "\\usepackage{units}",         # load additional packages
#         "\\usepackage{metalogo}",
#         "\\usepackage[math-style = TeX]{unicode-math}",  # unicode math setup
#         r"\setmathfont{Helvetica Neue LT Std}",
#         r"\setmathfont[range=\mathit]{HelveticaNeueLTStd-It.otf}",
#         r"\setmainfont[UprightFont={* 55 Roman},ItalicFont={* 56 Italic}]{Helvetica Neue LT Std}", # serif font via preamble
#         ]
#}
pgf_with_custom_preamble = {
    "font.family": "serif", # use serif/main font for text elements
    "svg.fonttype": "path", # use serif/main font for text elements
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": [
         "\\usepackage{units}",         # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage[math-style = TeX]{unicode-math}",  # unicode math setup
         r"\setmathfont{Myriad Pro}",
         r"\setmathfont[range=\mathit]{MyriadPro-It.otf}",
         r"\setmainfont[UprightFont={[MyriadPro-Regular.otf]},ItalicFont={[MyriadPro-It.otf]}]{Myriad Pro}", # serif font via preamble
         ]
}
matplotlib.rcParams.update(pgf_with_custom_preamble)




#rcParams['text.latex.preamble'] = [
#       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
#       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
#       r'\usepackage{helvet}',    # set the normal font here
#       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
#       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
#]  

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

# first figure out last line of dataset
f = open(filename)
fl = f.readlines()
f.close()

parline = 0

for i in reversed(range(0, len(fl))):

    if re.match("^\d.*",fl[i].strip()) is not None:
        parline = i
        break

data = pd.read_csv(filename, sep=";", nrows=i)


data["f_2"] = pd.Series(list(itertools.repeat(np.nan,len(data.index))))

data["f_2"] = 1.0 - data["f_0"] - data["f_1"]

# initialize and specify size 
fig = plt.figure(figsize=(6,5))

num_rows = 1

# add first subplot depicting % type 1 offspring
ax = plt.subplot(num_rows,1,1)

ax.plot(
        data["time"],data["phh"],'r',
        data["time"],data["pdh"],'b',
        )


ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_ylim((-0.05,1.05))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")

plt.ylabel(r'Probality offspring is hawk')
plt.xlabel(r'Time, $t$')
plt.legend((r'$p_{H}$',r'$p_{D}$'))

graphtype = "svg"

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) 

graphname_pdf = graphname + ".pdf"

# mess with svg fonts, let's do it differently

plt.savefig(graphname_pdf,format="pdf",transparent=True)

if graphtype == "svg":
    graphname_svg = graphname + ".svg"

    subprocess.call(["pdf2svg",graphname_pdf,graphname_svg])


