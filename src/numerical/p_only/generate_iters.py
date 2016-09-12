#!/usr/bin/env python

import numpy as np

step = 0.02
vvals = list(np.arange(0.01, 1.0, step))
cvals = list(np.arange(0.01, 1.0, step))

d = [ 0.05, 0.1, 0.5, 1.0 ]

exe = "./numsolve"

pdh_init = [ 0.5 ] 
phh_init = [ 0.5 ] 

ctr = 0

for v_i in vvals:
    for c_i in cvals:
        for d_i in d:
            for pdh_i in pdh_init: 
                for phh_i in phh_init:
                    print("echo " + str(ctr))
                    ctr+=1
                    print(exe + " 1000000 " 
                            + str(d_i) + " " 
                            + str(v_i) + " " 
                            + str(c_i) + " "
                            + str(0.33) + " "
                            + str(0.33) + " "
                            + str(pdh_i) + " "
                            + str(phh_i) + " "
                            + " 0.1 0.1 0.1 "
                            + " 1.0 1.0 1.0 1.0 "
                            )

