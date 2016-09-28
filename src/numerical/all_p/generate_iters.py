#!/usr/bin/env python

import numpy as np
import itertools, sys

step = 0.02

# order of the payoffs
# first element is Cc
# 2nd is Dc
# 3rd is Cd
# 4th is Dd
mortalities = [ 0.5, 1.0, 1.5, 2.0 ]

mortalities_ordering = [ 
        [ 0, 1, 2, 3 ],  # Cc > Dc > Cd > Dd (Soldier's dilemma)
        [ 0, 1, 3, 2 ],  # Cc > Dc > Dd > Cd (Stag hunt)
        [ 1, 0, 2, 3 ],  # Dc > Cc > Cd > Dd (Hawk Dove)
        [ 1, 0, 3, 2 ]]  # Dc > Cc > Dd > Cd (Prisoner's dilemma)
        


d = [ 0.1 ] #0.05, 0.1, 0.2, 0.5, 1.0 ]

exe = "./numsolve"

pdh_init = [ 0.5 ] 
phh_init = [ 0.5 ] 

ctr = 0

for mortality_i in mortalities_ordering:
    for d_i in d:
        for pdh_i in pdh_init: 
            for phh_i in phh_init:
                print("echo " + str(ctr))
                ctr+=1
                print(exe + " 100000000 " 
                        + str(d_i) + " " 
                        + str(mortalities[mortality_i[0]]) + " " 
                        + str(mortalities[mortality_i[1]]) + " " 
                        + str(mortalities[mortality_i[2]]) + " " 
                        + str(mortalities[mortality_i[3]]) + " " 
                        + str(pdh_i) + " "
                        + str(phh_i) + " "
                        + str(0.33) + " "
                        + str(0.33) + " "
                        + str(pdh_i) + " "
                        + str(pdh_i) + " "
                        + str(phh_i) + " "
                        + str(phh_i) + " "
                        + " 0.1 0.1 0.1 "
                        + " 1.0 1.0 1.0 1.0 "
                        )

