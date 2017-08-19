#!/usr/bin/env python

import numpy as np

step = 0.05
#vvals = list(np.arange(0.01, 1.0, step))
#cvals = list(np.arange(0.01, 1.0, step))

vvals = [ 0.25 ]
cvals = [ 0.8 ]

d = [ 0.1, 0.5, 1.0 ]

exe = "./numsolve"

step = 0.1
pdh_init = [ 0.01, 0.1, .2, .3, .4, .5, .6, .7, .8, .9, 0.99 ]
phh_init = [ 0.01, 0.1, .2, .3, .4, .5, .6, .7, .8, .9, 0.99 ]

ctr = 0

error = [ 0.02 ] 

for v_i in vvals:
    for c_i in cvals:
        for error_i in error:
            for pdh_init_i in pdh_init:
                for phh_init_i in phh_init:
                    for d_i in d:
                        print("echo " + str(ctr))
                        ctr+=1
                        print(exe + " 1000000 " 
                                + str(d_i) + " " 
                                + str(1.0 - v_i/2.0) + " "  # md_0
                                + str(1.0 - v_i) + " " # mh_1
                                + str(1.0) + " " # md_1
                                + str(1.0 - (v_i - c_i)/2.0) + " " # mh_2
                                + str(pdh_init_i) + " " 
                                + str(phh_init_i) + " " 
                                + str(error_i) + " " # err
                                + str(0.33) + " "
                                + str(0.33) + " "
                                + str(pdh_init_i) + " " 
                                + str(phh_init_i) + " " 
                                + " 0.1 0.1 0.1 "
                                + " 1.0 1.0 1.0 1.0 "
                                )

