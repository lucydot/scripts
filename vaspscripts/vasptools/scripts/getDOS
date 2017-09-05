#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
getDOS
------

NAME     
        getDOS - extract DOS from vasprun.xml

SYNTAX
        getDOS [OPTIONS] ... [FILE]

DESCRIPTION
        Read density of states data on vasprun.xml files and either extract them into
        files or plot them directly. The last argument has to be the xml file, if absent
        './vasprun.xml' will be used.

        -h, --help
            print this help

        -t, --tofiles
            Print total DOS into files. Projected DOS are printed if you add option -p.

        -p, --projected
            Take care about projected DOS. Projected DOS are not ploted. That option is
            relevant only when you want to print projected DOS into files.

        -q, --quiet
            low verosity
"""

import os
import sys
import vasptools

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

def getDOS() :
    """ Cette fonction permet d'extraire la densité d'états totale et les densités d'états
        partielles projetées sur les atomes et de les mettre en forme pour les tracer. 
    """

    # ----------------------------------------------------------
    # default options
    # ----------------------------------------------------------
    xml        = "vasprun.xml"
    projected = False
    toFile     = False
    verbose    = True

    # ----------------------------------------------------------
    # manage options
    # ----------------------------------------------------------
    narg = len(sys.argv)
    if narg > 1:
        for i in range(narg)[1:]:
            if sys.argv[i] == "-h" or sys.argv[i] == "-help":
                usage(0)

            elif sys.argv[i] == "-p" or sys.argv[i] == "--projected":
                projected = True
                toFile = True

            elif sys.argv[i] == "-t" or sys.argv[i] == "--tofiles": 
                toFile = True

            elif sys.argv[i] == "-q" or sys.argv[i] == "--quiet":
                verbose = False

            else:
                if i == narg - 1:
                    xml = sys.argv[i]
                else:
                    print(sys.argv)
                    print("Error : bad arguments")
                    usage(1)

    if not os.path.exists( xml ):
        die("Error : file {0} does not exist".format(xml))

    # ----------------------------------------------------------
    # get DOS
    # ----------------------------------------------------------
    calcul = vasptools.VaspRun(xml, verbose=verbose )
    dosLues = calcul.lectureDOS()
    if not dosLues[0]:
        die("Error when reading DOS in {0}".format(xml))

    # projected DOS ?
    if projected:
        if not calcul.DOSPartiellesLues:
            print("Warnings : you ask for projected DOS but they are not found !!\n")
            projected = False

    if toFile:
        # total DOS
        fileList = vasptools.dos.printTotalDOStoFile(calcul)

        # projected DOS
        if projected: 
            fileList += vasptools.dos.printProjectedDOStoFile(calcul)

        # created files list
        if verbose:
            print("\nCreated file(s) : " + str(len(fileList)))
            for fileName in fileList:
                print("\t" + fileName)

    else:
        vasptools.dos.showTotalDOS(calcul)

if __name__ == "__main__":
    getDOS()



