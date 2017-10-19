#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
chgsum
------

NAME     
        chgsum - do linear operations on charge densities from two CHGCAR files

SYNTAX
        chgsum CHGCAR1 CHGCAR2 fact1 fact2

DESCRIPTION
        on output CHGCAR_sum, write the density computed by the following expression :
            fact1 * CHGCAR1 + fact2 * CHGCAR2
"""

import sys
import os

from vasptools.utils import sumDensities

__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"

def usage(code):
    """ usage """
    print(__doc__)
    exit(code)

def die(m):
    """ print m and exit with code 1 """
    print(m)
    exit(1)

def chgsum():
    """ read command line arguments and run chgsum """

    fact1 = 1.
    fact2 = 1.

    # command line options
    args = sys.argv
    nargs = len(args)
    if nargs == 3:
        chgcar1 = args[1]
        chgcar2 = args[2]

    elif nargs == 5:
        chgcar1 = args[1]
        chgcar2 = args[2]
        fact1 = float(args[3])
        fact2 = float(args[4])

    elif nargs == 2:
        arg = args[1].strip()
        if arg == "-h" or arg == "-help" or arg == "--help":
            usage(0)
        else:
            print("Error : unknown argument {0}".format(arg))
            usage(1)

    else:
        print("Error : bad arguments")
        usage(1)
        
    print("\n\t * * * {0} * * *\n".format(args[0]))
    # ask if you agree the computation
    answer = ""
    print(" {0} will compute : {1}*{2} + {3}*{4}".format(args[0], fact1, chgcar1, fact2, chgcar2))
    while answer != "y":
        answer = raw_input(" Do you agree ? (y/n) : ")
        if answer == "n":
            exit(0)
        elif answer == "y":
            continue
        else:
            print(" hit 'y' for yes or 'n' for non")

    # file CHGCAR exist ?
    if not os.path.exists(chgcar1):
        die("Error : File '{0}' does not exist !\n".format(chgcar1))

    if not os.path.exists(chgcar2):
        die("Error : File '{0}' does not exist !\n".format(chgcar2))

    # sum densities
    rho_sum = sumDensities(chgcar1, chgcar2, fact1, fact2)
    Npts = len(rho_sum)

    # name of output file
    print("------------------------------------------------------------")
    if os.path.exists("CHGCAR_sum"):
        answer = ""
        while answer != "y":
            answer = raw_input(" file CHGCAR_sum exists, overwrite it ? (y/n) : ")
            if answer == "n":
                exit(0)
            elif answer == "y":
                continue
            else:
                print("hit 'y' for yes or 'n' for non")

    print(" write sum density in : CHGCAR_sum")
    print("------------------------------------------------------------")

    fout = open( "CHGCAR_sum" , "w" )

    # head of CHGCAR_sum from CHGCAR1
    fchgcar1 = open( chgcar1, "r")
    end = False
    while not end:
        line = fchgcar1.readline()
        if line.strip() == "":
            end = True
        fout.write( line )

    fout.write( fchgcar1.readline() )
    fchgcar1.close()

    # write sum density
    i = 0
    while i < Npts:
        line = ""
        for j in range(5):
            line += "%18.11e " % rho_sum[i]
            i += 1
            if i >= Npts:
                break
        line += "\n"
        fout.write(line)

    fout.close()

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

if __name__ == "__main__":
    chgsum()

