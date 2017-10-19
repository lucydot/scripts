#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Tests for crystal module
------------------------

The aim of this file is to run several test cases on the crystal module and its Crystal
classe. Just run it in a terminal.

"""

import os
from math import fabs

import crystal

# data used for tests
# -------------------

POSCARv4 = """testCase
1.0
        4.0000000000         0.0000000000         0.0000000000
        0.5209445330         2.9544232590         0.0000000000
        1.7500000000         0.9069650274         2.8922161812
    1    4
Direct
     0.000000000         0.000000000         0.000000000
     0.250000000         0.166666672         0.250000000
     0.750000000         0.166666672         0.750000000
     0.250000000         0.833333015         0.750000000
     0.750000000         0.833333015         0.250000000
"""

POSCARv5 = """testCase
1.0
        4.0000000000         0.0000000000         0.0000000000
        0.5209445330         2.9544232590         0.0000000000
        1.7500000000         0.9069650274         2.8922161812
   Ti    O
    1    4
Direct
     0.000000000         0.000000000         0.000000000
     0.250000000         0.166666672         0.250000000
     0.750000000         0.166666672         0.750000000
     0.250000000         0.833333015         0.750000000
     0.750000000         0.833333015         0.250000000
"""

xyz = [[0.000000,   0.000000,   0.000000],
       [1.524324,   0.719145,   0.723054],
       [4.399324,   1.172628,   2.169162],
       [2.746620,   3.142242,   2.169162],
       [3.871620,   2.688760,   0.723054]]

red = [[0.,      0., 0.], 
       [0.25, 1./6., 0.25], 
       [0.75, 1./6., 0.75], 
       [0.25, 5./6., 0.75], 
       [0.75, 5./6., 0.25]]

th = 1.e-5

# create a reference crystal
# --------------------------
refCrystal = crystal.Crystal(a = 4., b = 3., c = 3.5, alpha = 70., beta = 60., gamma = 80.)
print(refCrystal)
refCrystal.printLatticeVectors()

# test XYZ coordinates
# --------------------
refCrystal.redCoord = red
refCrystal.computeXYZCoord()

print("Test red -> XYZ Coordinates")
res = True
for xyz1, xyz2 in zip(refCrystal.XYZCoord, xyz):
    for x1, x2 in zip(xyz1, xyz2):
        if fabs(x1 - x2) > th:
            res = False
            print(res)
            print(xyz1)
            print(xyz2)
            exit(1)
print("XYZ coordinates : {0}\n".format(res))

# test red coordinates
# --------------------
refCrystal.XYZCoord = xyz
refCrystal.computeRedCoord()

print("Test XYZ -> red Coordinates")
res = True
for red1, red2 in zip(refCrystal.redCoord, red):
    for r1, r2 in zip(red1, red2):
        if fabs(r1 - r2) > th:
            res = False
            print(res)
            print(red1)
            print(red2)
            exit(1)
print("red coordinates : {0}\n".format(res))

# test creating a crystal object from VASP or CRYSTAL
# ---------------------------------------------------

# create temporary POSCARtest file
open("POSCARtest", "w").write(POSCARv4)

# read crystal from a POSCAR docstring
c = crystal.Crystal.fromPOSCAR(POSCARv5, verbose = False)
print("read crystal from a POSCAR docstring      : {0}".format(c == refCrystal))

# read from a list from a POSCAR
c = crystal.Crystal.fromPOSCAR(POSCARv5.split("\n"), verbose = False)
print("read crystal from a list of lines         : {0}".format(c == refCrystal))

# read from a POSCAR file object
fposcar = open("POSCARtest", "r")
c = crystal.Crystal.fromPOSCAR(fposcar, verbose = False)
print("read crystal from a POSCAR file object    : {0}".format(c == refCrystal))

# read from a POSCAR file name
c = crystal.Crystal.fromPOSCAR("POSCARtest", verbose = False)
print("read crystal from a POSCAR file name      : {0}".format(c == refCrystal))

# remove POSCARtest
os.remove("POSCARtest")

# read form a CRYSTAL09 output file
c = crystal.fromCRYSTAL09("crystal09.out", verbose = False)
print("read crystal from a CRYSTAL09 output file : {0}".format(c == refCrystal))

# test POSCAR output
# ------------------
c = crystal.Crystal.fromPOSCAR(POSCARv5, verbose = False)
c.toPOSCAR("POSCARtest")

# test bravais lattice
# --------------------
veca = [1., 0., 0.]
vecb = [0., 1., 0.]
vecc = [0., 0., 1.]
cubic = crystal.Crystal(veca = veca, vecb= vecb, vecc = vecc)
print("\nTest cubic lattice        : {0}".format(cubic.lattice == "cubic"))

veca = [ 2., 0.     , 0.]
vecb = [-1,  1.73205, 0.]
vecc = [ 0., 0.     , 3.]
hexa = crystal.Crystal(veca = veca, vecb= vecb, vecc = vecc)
print(  "Test hexagonal lattice    : {0}".format(hexa.lattice == "hexagonal"))

veca = [1., 0., 0.]
vecb = [0., 2., 0.]
vecc = [0., 0., 3.]
ortho = crystal.Crystal(veca = veca, vecb= vecb, vecc = vecc)
print("Test orthorhombic lattice : {0}".format(ortho.lattice == "orthorhombic"))

veca = [2., 0., 0.]
vecb = [0., 2., 0.]
vecc = [0., 0., 3.]
tetra = crystal.Crystal(veca = veca, vecb= vecb, vecc = vecc)
print("Test tetragonal lattice   : {0}".format(tetra.lattice == "tetragonal"))

veca = [1., 0., 0.]
vecb = [0., 2., 0.]
vecc = [1., 0., 1.]
monoclinic = crystal.Crystal(veca = veca, vecb= vecb, vecc = vecc)
print("Test monoclinic lattice   : {0}".format(monoclinic.lattice == "monoclinic"))

veca = [1., 0., 1.]
vecb = [1., 1., 0.]
vecc = [0., 1., 1.]
rhombo = crystal.Crystal(veca = veca, vecb= vecb, vecc = vecc)
print("Test monoclinic lattice   : {0}".format(rhombo.lattice == "rhombohedral"))

veca = [1., 0., 0.]
vecb = [1., 1., 0.]
vecc = [1., 1., 1.]
triclinic = crystal.Crystal(veca = veca, vecb= vecb, vecc = vecc)
print("Test triclinic lattice    : {0}".format(triclinic.lattice == "triclinic"))


