#!/usr/bin/env python
# -*- coding=utf-8 -*-

""" 
Tests of vasptools module
-------------------------

This test use almost all methods and functions of vasptools module and sub-modules. It
stops if there is an error in the execution. This test do not control de reliability of
the results it only controls if functions give an error.
"""

import os
import sys
import vasptools

def message(m, t = 0):
    if t == 1:
        print("# Error : {0}".format(m))
    elif t == 0:
        print("* {0} : OK".format(m))

def removeDataFiles(d = "."):
    [os.remove(f) for f in os.listdir(d) if f[-4:].strip() == ".dat"]

def dosTests(xml):
    """ test DOS extraction from file xml """

    # try read DOS
    try:
        run = vasptools.VaspRun(xml, verbose=False)
        run.lectureDOS()
        message("Read DOS (run.lectureDOS())")
    except:
        message("Read DOS (run.lectureDOS())", t = 1)

    # try print total DOS to file
    try:
        vasptools.dos.printTotalDOStoFile(run)
        message("export total DOS (vasptools.dos.printTotalDOStoFile())")
        removeDataFiles()
    except:
        message("export total DOS (vasptools.dos.printTotalDOStoFile())", t = 1)

    # try print total DOS to file
    try:
        vasptools.dos.printProjectedDOStoFile(run)
        message("export projected DOS (vasptools.dos.printProjectedDOStoFile())")
        removeDataFiles()
    except:
        message("export projected DOS (vasptools.dos.printProjectedDOStoFile())", t = 1)

    # try show total DOS
    try:
        vasptools.dos.showTotalDOS(run)
        message("plot total DOS (vasptools.dos.showTotalDOS())")
    except:
        message("plot total DOS (vasptools.dos.showTotalDOS())", t = 1)

def bandsTests(xml):
    """ test bands extraction from xml """

    # try read bands
    try:
        run = vasptools.VaspRun(xml, verbose=False)
        run.lectureBandes()
        message("Read bands (run.lectureBandes())")
    except:
        message("Read bands (run.lectureBandes())", t = 1)

    # try print bands to file
    try:
        vasptools.bands.bandsToFiles(run, True)
        message("export energy bands (vasptools.bands.bandsToFiles())")
        removeDataFiles()
    except:
        message("export energy bands (vasptools.bands.bandsToFiles())", t = 1)

    try:
        vasptools.bands.bandsToFiles(run, False)
        message("export energy bands (vasptools.bands.bandsToFiles())")
        removeDataFiles()
    except:
        message("export energy bands (vasptools.bands.bandsToFiles())", t = 1)

    # try plot bands
    try:
        vasptools.bands.showBandes(run)
        message("plot energy bands (vasptools.bands.showBandes())")
    except:
        message("plot energy bands (vasptools.bands.showBandes())", t = 1)

def densityTests(*f):
   """ test on charge density managing """

   try:
       rup, rdown = vasptools.utils.readCHGCAR(f[0], verbose = False)
       message("read CHGCAR (vasptools.utils.readCHGCAR())")
   except:
       message("read CHGCAR (vasptools.utils.readCHGCAR())", t = 1)

   try:
       vasptools.utils.sumDensities(f[1], f[2], 1.0, 1.0, verbose = False)
       message("arithmetic on CHGCAR ( vasptools.utils.sumDensities())")
   except:
       message("arithmetic on CHGCAR ( vasptools.utils.sumDensities())", t = 1)


def main():
    """ run test """

    print("\n * * * do some tests on vasptools * * *\n")
    print("Standard DOS:")
    dosTests("./dos.xml")
    print("Spin polarized DOS:")
    dosTests("./DOS_ISPIN2.xml")
    print("Energy bands:")
    bandsTests("./bandes.xml")
    print("CHGCAR operation:")
    densityTests("CHGCAR", "AECCAR0", "AECCAR2")

if __name__ == "__main__":
    main()

