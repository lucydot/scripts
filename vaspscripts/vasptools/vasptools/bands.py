#!/usr/bin/env python
# -*-coding:utf-8 -*-

# this file is part of vasptools

"""
vasptools.bands
---------------

This module contains functions in order to plot energy bands or to print them into files.

All functions take only one arguments which is a VaspRun object, see vasprun.py 

showBandes(run)    : plot energy bands with matplotlib
bandsToFiles(run)  : print energy bands into files

output files :
    * bands.dat or bands_up.dat and bands_down.dat for spin polarized calculations.
    * bands_dirX.dat or bands_dirX_up.dat and bands_dirX_down.dat for spin polarized
    calculations if one file is created for each line of the reciprocal space. X is the
    number of the line.

"""

import matplotlib.pyplot as plt

__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"

spinName = ["up", "down"]

def showBandes(run):
    """ trace les bandes avec matplotlib.pyplot """

    if not run.bandesLues:
        bandesLues = run.lectureBandes()
        if not bandesLues:
            print("Erreur when reading bands")
            exit(1)

    # infos calcul
    ISPIN = run.allMotsClefs["ISPIN"]
    NBANDS = run.allMotsClefs["NBANDS"]

    # abscisse
    nbreKpts = len(run.listePointsK)
    abscisse = [float(k) / float(nbreKpts) for k in range(nbreKpts)]

    # preparation du graph des bandes d'énergies
    plt.figure(1)
    plt.title("Energy Bands from {0}".format(run.xmlFile))
    plt.xlabel("k points")
    plt.ylabel("Energy - E_fermi  /   eV")
    plt.grid()

    for spin in range(ISPIN):
        if spin == 0:
            style = "r-"
        else:
            style = "k--"

        for i in range(NBANDS):
            bande = list()
            for k in range(nbreKpts):
                bande.append(run.bands[spin][k][i][0] - run.eFermi)

            plt.plot(abscisse, bande, style)

    # plot fermi level
    plt.plot([abscisse[0], abscisse[-1]], [0., 0.], "b-")

    plt.show()

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def bandsToFiles(run, parDirection = False):
    """ Méthode permettant d'imprimer les bandes d'énergie dans un fichier dans un format
    pratique pour les tracer. """

    # lecture des bandes sur le fichier xml
    if not run.bandesLues:
        bandesLues = run.lectureBandes()
        if not bandesLues:
            print("Erreur when reading bands")
            exit(1)

    if parDirection and run.typePointsK == "explicit":
        print("Warnings bands could not printed by direction with an explicit kpoints grid")
        parDirection = False

    if run.verbose:
        print("\n# extraction des bandes d'énergie")

    # spin polarise
    ISPIN = run.allMotsClefs["ISPIN"]

    # nombre de bandes
    NBANDS = run.allMotsClefs["NBANDS"]

    # number of k-points
    nbreKpoints = len(run.listePointsK)

    # fichiers créés
    listeFichiers = list()

    # boucle sur les spins
    for spin in range(ISPIN):
        k = 0

        # fichier de sortie dans le cas ou toutes les bandes sont regroupées
        if not parDirection:
            if ISPIN == 1:
                fichier = "bands.dat"
            else:
                fichier = "bands_{0}.dat".format(spinName[spin])
                
            fout = open( fichier, "w")
            listeFichiers.append(fichier)
            fout.write("# E_fermi = %e\n" % run.eFermi)
            fout.write("# energy bands is E - E_fermi in eV\n")
            fout.write("# arbitrary     energy bands 1 -> NBANDS in eV\n")

        # boucle sur les points k
        for direction in range(len(run.directionsPointsK)):

            # ouverture du fichier dans le cas d'un fichier par direction de points K
            if parDirection:
                if ISPIN == 1:
                    fichier = "bands_dir{0}.dat".format(direction + 1)
                else:
                    fichier = "bands_dir{0}_{1}.dat".format(direction + 1, spinName[spin])

                fout = open(fichier, "w")
                listeFichiers.append(fichier)
                ptr = "# direction " + str(direction+1) + " : " + \
                    str(run.directionsPointsK[direction][0:3]) + " -> "   + \
                    str(run.directionsPointsK[direction][3:6]) + "\n"

                fout.write(ptr)
                fout.write("# E_fermi = %e\n" % run.eFermi)
                fout.write("# energy bands is E - E_fermi in eV\n")
                fout.write("#   kx          ky          kz          energy bands 1 -> NBANDS in eV\n")

            for div in range(run.Ndivision):

                ptr = ""
                # abscisse
                if parDirection:
                    for coord in run.listePointsK[k]:
                        ptr += "%12.8f" % coord
                else:
                    # indice du points k entre 0 et 1
                    abscisse = float(k) / float(nbreKpoints)
                    ptr = "%12.8f" % abscisse

                # impression des valeurs des énergies
                for bande in range(NBANDS):
                    ptr += "%10.4f" % (run.bands[spin][k][bande][0] - run.eFermi)

                ptr += "\n"

                fout.write(ptr)

                k += 1

            # fermeture du fichier dans le cas d'un fichier par direction de points k
            if parDirection:
                fout.close()

        # controle du nombre de points k
        if k != nbreKpoints:
            print(k)
            print(nbreKpoints)
            print("erreur sur le nombre de points k ")
            exit(1)

        # fermeture du fichier dans le cas ou toutes les bandes sont regroupées
        if not parDirection:
            fout.close()

    return listeFichiers

