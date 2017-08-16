#!/usr/bin/env python

import numpy as np

# c and v values

step = 0.1
#cvals = list(np.arange(0, 1.0+step, step))
#vvals = list(np.arange(0, 1.0+step, step))

cvals = [3.0, 5.0 ]
vvals = [1.0]

phvals = list(np.arange(0,1.0+step,step))
pdvals = list(np.arange(0,1.0+step,step))

exe = "./xwell_mixed"

nrep=1
ctr = 0

for rep_i in range(0,nrep):
    for cval_i in cvals:
        for vval_i in vvals:
            for ph_i in phvals:
                for pd_i in pdvals:
                    print("echo " + str(ctr))
                    ctr+=1

                    print(exe + " " + str(vval_i)
                            + " " + str(cval_i)
                            + " 0.01 0.02 " + str(ph_i) + " " + str(pd_i))

