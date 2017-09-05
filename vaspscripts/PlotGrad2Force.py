#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt

""" Dependancies: Larsson's grad2 programme: https://www.nsc.liu.se/~pla/vasptools/.
Takes output from grad2 and plots forces.
"""

os.system("grad2 OUTCAR  >> grad2_data.tmp")
avg_forces = []
max_forces = []
with open ("grad2_data.tmp", "r") as myfile:
    for line in myfile:
       # data = myfile.readlines()
        if "Avg|F|:" in line:
            avg_force = line.split("Avg|F|:")[1].split()[0]
            avg_forces.append(avg_force)
        if "Max|F|:" in line:
            max_force = line.split("Max|F|:")[1].split()[0]
            max_forces.append(max_force)

[float(i) for i in avg_forces]
x_coords = range(1, len(avg_forces)+1)
plt.scatter(x_coords[:],avg_forces[:])
plt.savefig('avg_force.png')

plt.clf() # close the figure

[float(i) for i in max_forces]
x_coords = range(1, len(max_forces)+1)
plt.scatter(x_coords[:],max_forces[:])
plt.savefig('max_force.png')

os.system("rm grad2_data.tmp")
print (avg_forces)
print (max_forces)

