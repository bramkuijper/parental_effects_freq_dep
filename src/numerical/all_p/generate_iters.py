#!/usr/bin/env python

import numpy as np
import itertools, sys

step = 0.05

d = [ 0.1, 0.5, 1.0 ]

v = list(np.arange(0,1+step,step))
c = list(np.arange(0,1+step,step))


exe = "./numsolve"

pdh_init = [ 0.5 ] 
phh_init = [ 0.5 ] 

ctr = 0

error = [ 0, 0.02 ]

for v_i in v:
    for c_i in c:
        for d_i in d:
            for pdh_i in pdh_init: 
                for phh_i in phh_init:
                    for error_i in error:
                        print("echo " + str(ctr))
                        ctr+=1
                        print(exe + " 100000000 " 
                                + str(d_i) + " " 
                                + str(1.0 - v_i/2.0) + " " #md_0
                                + str(1.0 - v_i) + " "  #mh_1
                                + str(1.0) + " "  #md_1 
                                + str(1.0 - (v_i-c_i)/2) + " "  #mh_2
                                + str(pdh_i) + " "
                                + str(phh_i) + " "
                                + str(error_i) + " "
                                + str(0.33) + " "
                                + str(0.33) + " "
                                + str(pdh_i) + " "
                                + str(pdh_i) + " "
                                + str(phh_i) + " "
                                + str(phh_i) + " "
                                + " 0.1 0.1 0.1 "
                                + " 1.0 1.0 1.0 1.0 "
                                )

