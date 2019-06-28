#!/usr/bin/env python

import numpy as np
import math
import sys
import ast
from IPython import embed

# constants
q = 1.60217657E-19
hbar = 1.054E-34
k = 1.3806504e-23

def heat_capacity(T, PhononFrequency):
    PhononEnergy = np.multiply(PhononFrequency,1E12*2*math.pi*hbar) # phonopy gives frequency in THz
    Cv=[]
    for En in PhononEnergy:
        occupation = 1/(np.exp(En/(k*T))-1)
        Cv.append(1000/q*occupation*En)
    heat_capacity = np.sum(Cv)/T
    return heat_capacity # meV/K

if __name__=="__main__":
    print (str(heat_capacity(ast.literal_eval(sys.argv[1]),ast.literal_eval(sys.argv[2])))+"meV/K")
