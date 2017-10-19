#!/usr/bin/env python
# -*-coding:utf-8 -*-

# this file is part of vasptools package

"""
Module utils
------------

This module contains functions in order to manipulate charge density."""

__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"

def readCHGCAR(chgcarName, full=False, verbose = True):
    """ read up+down and up-down densities in a CHGCAR file """

    # output density
    rho_up_p_down = list()
    rho_up_m_down = list()

    # open CHGCAR file
    fchgcar = open(chgcarName , "r")

    # read head
    for i in range(5):
        fchgcar.readline()

    # read atoms section
    vasp5 = False
    line = fchgcar.readline()
    if line.split()[0].isdigit():
        # vasp 4 type
        pass
    else:
        # vasp 5 type
        atomNames = line
        fchgcar.readline()
        vasp5 = True

    # coordinates section
    fchgcar.readline()
    end = False
    Natom = 0
    while not end:
        line = fchgcar.readline()
        if line.strip() == "":
            end = True
        else:
            Natom += 1

    # write some data
    if verbose:
        print("------------------------------------------------------------")
        print(" Atoms in file {0} : {1} ".format(chgcarName, Natom))
        if vasp5:
            print(" Atom names in file {0} : {1}".format(chgcarName, atomNames))

    # grid size
    NGFX, NGFY, NGFZ = [int(val) for val in fchgcar.readline().split()[0:3]]
    Npts = NGFX * NGFY * NGFZ
    if verbose:
        print(" Grid size up+down : {0} points".format(Npts))

    # read up+down density
    i = 0
    while i < Npts:
        values = [float(val) for val in fchgcar.readline().split()]
        for val in values:
            rho_up_p_down.append(val)
            i += 1

    if full:
        # skip augmentation occupancies
        while True:
            line = fchgcar.readline().split()
            if len(line) == 0:
                print("Error : end of file reached before I find up - down density")
                exit(1)
                
            elif line[0].isdigit():
                NGFX2, NGFY2, NGFZ2 = [int(val) for val in line]
                if NGFX2 == NGFX and NGFY2 == NGFY and NGFZ2 == NGFZ:
                    break

                else:
                    print("Error : grid size of the up-down density unconsistent with the \
                    grid size of the up+down density.")
                    exit(1)

        Npts2 = NGFX2 * NGFY2 * NGFZ2
        if verbose:
            print(" Grid size up-down : {0} points".format(Npts2))

        # read up-down density
        i = 0
        while i < Npts:
            values = [float(val) for val in fchgcar.readline().split()]
            for val in values:
                rho_up_m_down.append(val)
                i += 1

    fchgcar.close()

    return rho_up_p_down, rho_up_m_down

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def sumDensities(chgcar1, chgcar2, fact1, fact2, verbose = True):
    """ read CHGCAR1 and CHGCAR2 files and do the following linear operation :
                CHGCAR_sum = fact1 * CHGCAR1 + fact2 * CHGCAR2

        work only on the first (up+down) density
    """

    # read densities
    rho1 = readCHGCAR(chgcar1, verbose = verbose)[0]
    rho2 = readCHGCAR(chgcar2, verbose = verbose)[0]

    # test vectors lengths
    Npts1 = len(rho1)
    Npts2 = len(rho2)
    if Npts1 != Npts2:
        print("Error : grid sizes are unconsistent")
        exit(1)

    # compute the sum
    rho_sum = list()
    for i in range(Npts1):
        val = fact1 * rho1[i] + fact2 * rho2[i]
        rho_sum.append(val)

    return rho_sum

