#!/usr/bin/env python

import numpy as np

# c and v values

step = 0.05
cvals = list(np.arange(0, 1.0+step, step))
vvals = list(np.arange(0, 1.0+step, step))

dvals = [ 0.1, 0.5, 1.0 ] 

#exe = "./xfreq_dep"
exe = "./xphh_pdh_freq_dep"

nrep=1
ctr = 0

for rep_i in range(0,nrep):
    for cval_i in cvals:
        for vval_i in vvals:
            for dval_i in dvals:
                print("echo " + str(ctr))
                ctr+=1

                print(exe + " " + str(dval_i) 
                        + " " + str(vval_i)
                        + " " + str(cval_i)
                        + " " + "0.5")
