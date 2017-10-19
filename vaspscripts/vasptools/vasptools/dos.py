#!/usr/bin/env python
# -*-coding:utf-8 -*-

# This is a part of vasptools

"""
vasptools.dos
-------------

This module contains functions in order to print projected DOS into files or in order to
print or plot total DOS.

All functions take only one arguments which is a VaspRun object, see vasprun.py

showTotalDOS(run)            : plot total DOS using matplotlib.pyplot
printTotalDOStoFile(run)     : print total DOS in a file
printProjectedDOStoFile(run) : print projected DOS to files with one file per atom

output files :
    * totalDOS.dat or totalDOS_spin_X.dat, where XX is up or down.
    * projectedDOS_spinX_atY_iatZ.dat for projected DOS where X is for spin up or down, Y
    is the atom name and Z is its number.
"""

__licence__ = "GPL"
__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"

import matplotlib.pyplot as plt

# spin name for file name
spinName = ["up", "down"]

def showTotalDOS(run):
    """ plot total DOS using matplotlib.pyplot """
    
    if not run.DOSTotaleLue:
        dosLue = run.lectureDOS()
        if not dosLue[0]:
            print("Erreur when reading total DOS")
            exit(1)

    if run.verbose:
        print("\n# plot total DOS")

    # spin polarization ?
    ISPIN = run.allMotsClefs["ISPIN"]

    # shift DOS abscissa to have 0 eV at Fermi level
    energiesDOS = list()
    for e in run.energiesDOS:
        energiesDOS.append(e - run.eFermi)

    # total dos
    dosTotaleSpinUp = [ val[0] for val in run.dosTotale[0] ]
    if ISPIN == 2: 
        dosTotaleSpinDown = [ -val[0] for val in run.dosTotale[1] ]

    # plot total dos
    plt.figure(1)
    plt.title("Density Of State from " + run.xmlFile)
    plt.xlabel("Energy - E_fermi   /   eV")
    plt.ylabel("DOS")
    plt.grid()
    if ISPIN == 1:
        plt.plot( energiesDOS, dosTotaleSpinUp, "r-" )
    else:
        plt.plot( energiesDOS, dosTotaleSpinUp, 'r-', label="DOS up")
        plt.plot( energiesDOS, dosTotaleSpinDown, 'k-', label="DOS down" )
        plt.legend(fancybox=True,shadow=True)

    plt.show()

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def printTotalDOStoFile(run):
    """ print total DOS in a file """

    if run.verbose:
        print("\n# Read total DOS")

    if not run.DOSTotaleLue:
        dosLue = run.lectureDOS()
        if not dosLue[0]:
            print("Erreur when reading total DOS")
            exit(1)

    ISPIN = run.allMotsClefs["ISPIN"]

    fileList = list()

    # loop over spin
    for spin in range(ISPIN):

        # output file names
        if ISPIN == 1:
            fichier = "totalDOS.dat"
        else:
            fichier = "totalDOS_{0}.dat".format(spinName[spin])
                
        fout = open( fichier, "w")
        fileList.append(fichier)
        fout.write("# E_fermi = %e \n" % run.eFermi )
        fout.write("# hereafter, E is E - E_fermi\n")
        fout.write("#   E (eV)      DOS       Ne\n")

        for e, dos in zip(run.energiesDOS, run.dosTotale[spin]):
            ligne = "%10.4f" % (e - run.eFermi)
            for val in dos:
                ligne += "%10.4f" % val
            
            ligne += "\n"
            fout.write(ligne)

        fout.close()

    return fileList

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def printProjectedDOStoFile(run):
    """ print projected DOS to files with one file per atom. """

    fileList = list()

    if not run.DOSTotaleLue:
        dosLue = run.lectureDOS()
        if not dosLue[0]:
            print("Erreur when reading total DOS")
            exit(1)

    if not run.DOSPartiellesLues:
        print("\t ! Projected DOS not found !")
        return fileList

    if run.verbose:
        print("\n# Read projected DOS")

    ISPIN = run.allMotsClefs["ISPIN"]

    # number of projected DOS
    nProj = len(run.dosPartielles[0][0][0])
    if nProj == 3:
        if run.verbose:
            print("\t* DOS l-projected, order : s p d")
    elif nProj == 9:
        if run.verbose:
            print("\t* DOS lm-projected, order : s py pz px dxy dyz dz2 dxz dx2")
    else:
        print("Warning, error with projected DOS")
        nProj = 3

    # loop over ions
    for iat in range(run.Natomes):

        atom = run.atoms[iat]

        for spin in range(ISPIN):

            # output files
            if ISPIN == 1:
                fileName = "projectedDOS_{0}_{1}.dat".format(atom.name.strip(), iat + 1)
            else:
                fileName = "projectedDOS_{0}_{1}_{2}.dat".format(spinName[spin], atom.name.strip(), iat + 1)

            fout = open(fileName, "w")
            fileList.append(fileName)
            if nProj == 3:
                if ISPIN == 1:
                    line = "# contribution of atom {0}_{1}\n".format(atom.name.strip(), iat + 1)
                else:
                    line = "# contribution of atom {0}_{1} to the DOS {2}\n".format(atom.name.strip(), iat + 1, spinName[spin])
                line += "# E_fermi = %e \n" % run.eFermi
                line += "# column 1 : E - E_fermi (eV)\n"
                line += "# column 2 : s\n"
                line += "# column 3 : p\n"
                line += "# column 4 : d\n"
                fout.write(line)

            else:
                if ISPIN == 1:
                    line = "# contribution of atom {0}_{1}\n".format(atom.name.strip(), iat + 1)
                else:
                    line = "# contribution of atom {0}_{1} to the DOS {2}\n".format(atom.name.strip(), iat + 1, spinName[spin])
                line += "# E_fermi = %e \n" % run.eFermi
                line += "# column  1 : E - E_fermi (eV)\n"
                line += "# column  2 : s\n"
                line += "# column  3 : py\n"
                line += "# column  4 : pz\n"
                line += "# column  5 : px\n"
                line += "# column  6 : dxy\n"
                line += "# column  7 : dxz\n"
                line += "# column  8 : dz2\n"
                line += "# column  9 : dxz\n"
                line += "# column 10 : dx2-y2\n"
                fout.write(line)

            for e, dos in zip(run.energiesDOS, run.dosPartielles[iat][spin]):
                line = "%10.4f" % (e - run.eFermi)
                for val in dos:
                    line += "%10.4f" % val
        
                line += "\n"
                fout.write(line)

            fout.close()
    
    return fileList

