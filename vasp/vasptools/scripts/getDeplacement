#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
getDeplacement
--------------

NAME     
        getDeplacement

SYNTAX
        getDeplacement [FILE1] [FILE2]

DESCRIPTION
        Read first and last structures in a vasprun.xml file and, for each atom, compute
        the distance between the initial position and the last position. You can also give
        two POSCAR files (or CONTCAR). The first one will be use as the initial position
        and the second will be use as the final position.

        By default, the file './vasprun.xml' is read. You can give an other file in the
        command line, examples :
            getDeplacement my_calculation.xml
            getDeplacement POSCAR CONTCAR

        -h, --help
            print this help and exit
"""

import os
import sys
from math import sqrt
import vasptools 
import crystal

__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"

def usage( code ) :
    """ usage """
    print(__doc__)
    exit(code)

def die(m):
    """ print m and exit with code 1 """
    print(m)
    exit(1)

def distance( x1, x2 ):
    """ calcul d'une distance """
    r = 0.
    for i in range(3):
        r += ( x1[i] - x2[i])**2
    return sqrt(r)

def ptrsep( n = 50 ):
    ligne = ""
    for i in range(n): ligne += "-"
    print(ligne)

def getDeplacement():
    """ programme principal """

    # valeur par defaut des options
    xml = "vasprun.xml"
    typexml = True

    # import des donnees du calcul
    args = sys.argv
    narg = len(args)
    if narg == 1:
        pass
    elif narg == 2:
        if args[1] == "-h" or args[1] == "-help":
            usage(0)
        else:
            xml = args[1]
    elif narg == 3:
        if args[1] == "-h" or args[1] == "-help":
            usage(0)
        typexml = False
        fichier1 = args[1]
        fichier2 = args[2]
    else:
        print("Error : option unknown")
        usage(1)

    # check data
    if typexml:
        if "xml" not in xml:
            print("Error : file {0} must be an xml file".format(xml))
            usage(1)

        if not os.path.exists( xml ):
            die("Error : file {0} does not exist".format(xml))

    else:
        if not os.path.exists( fichier1 ):
            die("Error : file {0} does not exist".format(fichier1))
        if not os.path.exists( fichier2 ):
            die("Error : file {0} does not exist".format(fichier2))

    if typexml:
        # read initial structure
        calcul = vasptools.VaspRun(xml, verbose = False)
        print("\nStructure initiale : ")
        struct1 = calcul.getInitialStructure()
        print("\nStructure finale : ")
        struct2 = calcul.getFinalStructure()

    else:
        print("\nStructure initiale : ")
        f1 = open(fichier1,"r")
        struct1 = crystal.Crystal.fromPOSCAR(f1)
        f1.close()
        print("\nStructure finale : ")
        f2 = open(fichier2,"r")
        struct2 = crystal.Crystal.fromPOSCAR(f2)
        f2.close()

        if struct2.Natoms != struct1.Natoms:
            die("Error : atom number is not the same in {0} and {1}".format(fichier1,fichier2))

    # compute atoms displacement
    print("")
    ptrsep(75)
    print("   atome    positions initiales          positions finales       distance")
    ptrsep(75)
    for i in range(struct1.Natoms):
        ligne = "%4d" % (i+1) + struct1.atomNames[i].rjust(4) 
        for x in struct1.XYZCoord[i]:
            ligne += "%8.3f" % x

        ligne += "    "
        for x in struct2.XYZCoord[i]:
            ligne += "%8.3f" % x

        ligne += "    "
        ligne += "%8.3f" % distance( struct1.XYZCoord[i], struct2.XYZCoord[i] )

        print(ligne)

# run script
if __name__ == "__main__":
    getDeplacement()

