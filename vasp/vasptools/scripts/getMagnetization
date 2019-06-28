#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
getMagnetization
----------------

NAME
        getMagnetization - extract atomic spin polarization from OUTCAR

SYNTAX
        getMagnetization [OPTIONS]


DESCRIPTION
        Extract atomic spin polarization from the OUTCAR file and print them on output or
        in a file magnetization.agr. By default, extract the atomic spin polarization of
        the last geometry and print them on stdout.

        options

        -f, --full
            extract all atomic spin polarizations for each ionic step and print them into
            file magnetization.agr

        -t VALUE
            if all along the geometrical optimization, atomic spin polarization of an
            atom does not exceed VALUE, data will not be printed for this atom.

"""

__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"

import sys
import os
from math import fabs

import vasptools

def getMagnetization(full = False, plotThreshold = .1):
    """ read magnetization information into a OUTCAR file """

    if not os.path.exists("vasprun.xml"):
        print("WARNING : file vasprun.xml not found !")
        print("            * I assume this is a spin polarized calculation")
        print("            * Atom name are not known")
        xml = False

    else:
        xml = True
        run = vasptools.VaspRun(verbose = False)

        if run.allMotsClefs["ISPIN"] != 2:
            print("Error, ISPIN = {0}. There is not magnetization information in the OUTCAR file".format(run.allMotsClefs["ISPIN"]))
            exit(1)

    try:
        outcar = open("OUTCAR", "r").readlines()
    except:
        print("File OUTCAR does not exist")
        exit(1)

    # locate magnetization data
    positions = [i for i, line in enumerate(outcar) if "magnetization (x)" in line]

    # delete the last position (duplicate)
    positions.pop(-1)

    # read values
    moments = list()
    for p, pos in enumerate(positions):
        moments.append(list())
        for line in outcar[pos + 4:]:
            if "----" in line:
                break
            moments[p].append(float(line.split()[4]))

    Natoms = len(moments[0])

    plotlist = list()
    for iat in range(Natoms):
        maxval = fabs(moments[0][iat])
        for moment in moments:
            maxval = max(maxval, fabs(moment[iat]))
        if maxval > plotThreshold:
            plotlist.append(True)
        else:
            plotlist.append(False)

    
    if full:
        data  = "@g0 on\n"
        data += "@with g0\n"
        nset = 0
        for iat in range(Natoms):
            if plotlist[iat]:
                if xml:
                    data += '@    s%d legend  "%s"\n' % (nset, run.atoms[iat].name + str(iat + 1))
                else:
                    data += '@    s%d legend  "%s"\n' % (nset, "atome " + str(iat + 1))
                nset += 1

        nset = 0
        for iat in range(Natoms):
            if plotlist[iat]:
                data += "@target G0.s%d\n" % (nset)
                data += "@type xy\n"
                for i, moment in enumerate(moments):
                    data += "%4d %8.3f\n" % (i, moment[iat])
                nset += 1

        open("magnetization.agr", "w").write(data)

    else:
        print("----------------------------------------")
        print(" Magnetization")
        print("----------------------------------------")
        for i, m in enumerate(moments[-1]):
            if xml:
                print("%4d %4s %8.3f" % (i + 1, run.atoms[i].name, m))
            else:
                print("%4d %8.3f" % (i + 1, m))
        print("----------------------------------------")

if __name__ == "__main__":

    # default values
    # --------------
    full = False
    plotThreshold = .1

    # get options
    # -----------
    args = [arg for arg in sys.argv]
    if len(args) >= 2:
        if "-f" in args or "--full" in args:
            full = True

        if "-t" in args:
            try:
                plotThreshold = float(args[args.index("-t") + 1])
            except ValueError:
                print("-t options must be followed by a floating number")
                exit(1)

    # run get magnetization
    getMagnetization(full, plotThreshold)
