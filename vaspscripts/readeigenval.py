#!/usr/bin/python
# Lucy Whalley adapting Ryan Valenza
# 2017-09

import sys
import re
import math
import numpy as np

try:
        eigenval = open("EIGENVAL","r")
except IOError:
        sys.exit("Could not open EIGENVAL")

eigenval.readline() # N/A
eigenval.readline() # N/A
eigenval.readline() # N/A
eigenval.readline() # Cartesian/Direct

# System name w/ stripped newline character
name = eigenval.readline().rstrip() 
print "# System:        " + name

first = eigenval.readline() # Possibly interesting information
(nelect,nkpts,nbands) = first.split()
print "# of electrons:  " + nelect
print "# of k-points:   " + nkpts
print "# of bands:      " + nbands

# Regular expressions are used to distinguish between a k-point and an eigenvalue
regexs = {
'kpt': "\s+(.\d+\.\d+E[+-]\d+)\s+(.\d+\.\d+E[+-]\d+)\s+(.\d+\.\d+E[+-]\d+)\s+.*",
'enval': "\s+(\d+)\s+(-?.\d+\.\d+)" }

kpts = []
bands = []
for i in range(int(nbands)):
    bands.append([]) 
j = 0 # mark band number
dirn = np.array([None,None,None])
step = 0
print "Finding k-point, energy pairs for each band..."
for line in eigenval:
    kpt = re.match(regexs['kpt'], line)
    enval = re.match(regexs['enval'], line)

    if kpt != None:
        (kx,ky,kz) = kpt.groups()
        vector = np.array([float(kx),float(ky),float(kz)])
        d = np.divide(vector,np.sqrt(np.dot(vector,vector)))
        print (d)
        if dirn.all() == None:
            dirn = d
        elif np.array_equal(d,dirn) is False:
            step = max(bands[:][:][0])[0]                                                 
        k = math.sqrt(np.dot(vector,vector)) + step
    if enval != None:
        e = float(enval.groups(0)[1])
        bands[j%int(nbands)].append([k,e])
        j += 1
#As of right now, the data being is being organized into a format easily read
# into Mathematica.  The format can be changed to fit future problems.
# Ideas - use matplotlib to generate a bandstructure plot
#       - numpy arrays
print "Creating file bands.txt..."
out = open("LW_bands.txt", "w")
for i in range(int(nbands)):
    for j in range(int(nkpts)):
        out.write(" %.4f %.4f\n"%(bands[i][j][0],bands[i][j][1]))
    out.write("\n\n\n")

