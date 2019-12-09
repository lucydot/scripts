#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
getMaille
---------

NAME     
        getMaille

SYNTAX
        getMaille POSCAR or CONTCAR

DESCRIPTION
        Read the lattice vectors in a POSCAR or CONTCAR file and print the lattice
        parameters.

        -h, --help
            print this help and exit
"""

import os
import sys
from math import sqrt, acos, pi, fabs

__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"

dashedLine = "".join(50 * ["-"])

def usage(code) :
    """ usage """
    print(__doc__)
    exit(code)

def die(m):
    """ print error and exit """
    print(m)
    exit(1)

def getMaille():
    """ getMaille program """

    # arguments and input file
    args = sys.argv
    narg = len(args)
    if narg == 1:
        print(os.getcwd() + ":")
        for f in os.listdir("."):
            if "POSCAR" in f or "CONTCAR" in f:
                print("    * " + f)
        poscar = raw_input("\nfile = ")

    elif narg == 2:
        if args[1] == "-h" or args[1] == "-help":
            usage(0)
        else:
            poscar = args[1]

    else:
        print("Error, bad arguments")
        usage(1)

    # check data
    if not os.path.exists(poscar):
        die("Error : file {0} does not exist".format(poscar))

    # read poscar file
    lines = open(poscar, "r").readlines()
    title = lines[0][0:-1]
    scale = float(lines[1])
    veca  = [scale * float(v) for v in lines[2].split()[0:3]]
    vecb  = [scale * float(v) for v in lines[3].split()[0:3]]
    vecc  = [scale * float(v) for v in lines[4].split()[0:3]]

    # lattice parameters
    a = sqrt(veca[0]**2 + veca[1]**2 + veca[2]**2)
    b = sqrt(vecb[0]**2 + vecb[1]**2 + vecb[2]**2)
    c = sqrt(vecc[0]**2 + vecc[1]**2 + vecc[2]**2)

    # angles
    scalab = veca[0] * vecb[0] + veca[1] * vecb[1] + veca[2] * vecb[2]
    scalac = veca[0] * vecc[0] + veca[1] * vecc[1] + veca[2] * vecc[2]
    scalbc = vecb[0] * vecc[0] + vecb[1] * vecc[1] + vecb[2] * vecc[2]

    scalab /= a * b
    scalac /= a * c
    scalbc /= b * c

    # alpha
    if fabs(scalbc) < 1.:
        alpha = acos(scalbc) * 180. / pi
    else:
        die("Error, angle alpha")

    if fabs(scalac) < 1.:
        beta = acos(scalac) * 180. / pi
    else:
        die("Error, angle beta")

    if fabs(scalab) < 1.:
        gamma = acos(scalab) * 180. / pi
    else:
        die("Error, angle gamma")

    print(dashedLine)
    print(poscar + " : " + title)
    print(dashedLine)
    print("a     = %10.5f" % a)
    print("b     = %10.5f" % b)
    print("c     = %10.5f" % c)
    print("alpha = %10.3f" % alpha)
    print("beta  = %10.3f" % beta)
    print("gamma = %10.3f" % gamma)
    print(dashedLine)

# execution du programme
if __name__ == "__main__":
    getMaille()

