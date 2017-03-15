#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt

""" Dependancies: Larsson's grad2 programme: https://www.nsc.liu.se/~pla/vasptools/.
Takes output from grad2 and plots forces.
"""

os.system("grad2 OUTCAR >> grad2_data.tmp")
forces = []
with open ("grad2_data", "r") as myfile:
    for line in myfile:
       # data = myfile.readlines()
        if "Avg|F|:" in line:
            force = line.split("Avg|F|:")[1].split()[0]
            forces.append(force)

[float(i) for i in forces]
x_coords = range(1, len(forces)+1)
plt.scatter(x_coords[:],forces[:])
plt.savefig('force_convergence.png')
os.system("rm grad2_data.tmp")
print (forces)

