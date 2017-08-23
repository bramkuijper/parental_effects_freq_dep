#!/usr/bin/env python

import numpy as np

step = 0.05
vvals = list(np.arange(0.01, 1.0, step))
cvals = list(np.arange(0.01, 1.0, step))

d = [ 0.1, 0.5, 1.0 ]

exe = "./numsolve"

pdh_init = [ 0.5 ] 
phh_init = [ 0.5 ] 

ctr = 0

error = [ 0, 0.02 ] 

for v_i in vvals:
    for c_i in cvals:
        for error_i in error:
            for d_i in d:
                print("echo " + str(ctr))
                ctr+=1
                print(exe + " 1000000 " 
                        + str(d_i) + " " 
                        + str(1.0 - v_i/2.0) + " "  # md_0
                        + str(1.0 - v_i) + " " # mh_1
                        + str(1.0) + " " # md_1
                        + str(1.0 - (v_i - c_i)/2.0) + " " # mh_2
                        + str(0.5) + " " # p_init
                        + str(error_i) + " " # err
                        + str(0.33) + " "
                        + str(0.33) + " "
                        + str(0.5) + " "
                        + " 0.1 0.1 0.1 "
                        + " 1.0 1.0 1.0 1.0 "
                        )

