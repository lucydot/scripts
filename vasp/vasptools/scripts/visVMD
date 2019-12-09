#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
visVMD
------

NAME
        visVMD - write a script in order to visualise VASP structures with VMD

SYNTAX
        visVMD [OPTIONS]

DESCRIPTION
        visVMD read output files of VASP and write a VMD script in order to visualize VASP
        structures with VMD. visVMD can read either POSCAR and CONTCAR files or directly
        vasprun.xml files. The xml file is the best choice because it contains atom names.
        If you read a POSCAR or a CONTCAR file of a 4.X version of VASP, visVMD will ask
        you the name of each atom.

        -h, --help
            print this help

        -p [prefix], --prefix [prefix]
            visVMD use [prefix] in order to name output files (default is 'vis').

        -f [file], --file [file]
            [file] is an output file of VASP containing structure data. It is either the
            xml file or a POSCAR or CONTCAR (default 'vasprun.xml').

        --poscar
            say that the file is either a POSCAR or a CONTCAR

        -u, --unitcell
            add all images of the atoms of the unit cell which are into the unit cell

        -e, --final
            read the last structure of the xml file (default if [file] is an xml file).
   
        -i, --initial
            read the first structure of the xml file.

        -n [N1xN2xN3], --supercell [N1xN2xN3]
            make a supercell N1 times N2 times N3. Do not set spaces in the expression.

        -l, --liaisons
            make connectivity between atoms.
"""

import os
import sys
import numpy as np
import vasptools
import crystal

__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"

# ----------------------------------------------------------------------------------
# global variables
# ----------------------------------------------------------------------------------

# les 26 translations autour d'un cube
translations = [ [ 1., 0., 0.],
                 [ 1., 1., 0.],
                 [ 1.,-1., 0.],
                 [ 1., 0., 1.],
                 [ 1., 0.,-1.],
                 [ 1., 1., 1.],
                 [ 1.,-1., 1.],
                 [ 1., 1.,-1.],
                 [ 1.,-1.,-1.],
                 [-1., 0., 0.],
                 [-1., 1., 0.],
                 [-1.,-1., 0.],
                 [-1., 0., 1.],
                 [-1., 0.,-1.],
                 [-1., 1., 1.],
                 [-1.,-1., 1.],
                 [-1., 1.,-1.], 
                 [-1.,-1.,-1.],
                 [ 0., 1., 0.],
                 [ 0.,-1., 0.],
                 [ 0., 0., 1.],
                 [ 0., 0.,-1.],
                 [ 0., 1., 1.],
                 [ 0.,-1., 1.],
                 [ 0., 1.,-1.],
                 [ 0.,-1.,-1.] ] 

# ----------------------------------------------------------------------------------

def usage( code ):
    """ usage """
    print(__doc__)
    exit(code)

def die(m):
    """ print m and exit with code 1 """
    print(m)
    exit(1)

# ----------------------------------------------------------------------------------

def completeMaille(struct):
    """ Add atoms at the edge, at the top and on the side of the cell """

    # seuils
    seuilmax =  1.01
    seuilmin = -0.01

    atomesSupXYZ = list()

    # complete la maille
    listeAtomes = list()
    for iat in range(struct.Natoms):
        for trans in translations:
            # attention il faut faire une hard copy
            coord = [ val for val in struct.redCoord[iat] ]
            coord[0] += trans[0]
            coord[1] += trans[1]
            coord[2] += trans[2]

            # teste pour savoir si x est dans la maille unitaire
            add = True
            for x in coord:
                if x > seuilmax or x < seuilmin:
                    add = False
                    break

            if add:
                atomesSupXYZ.append(struct.red2cart(coord))
                listeAtomes.append(struct.atomNames[iat])

    return atomesSupXYZ, listeAtomes

# ----------------------------------------------------------------------------------

def makeVMDScript(prefix, cristal, nx, ny, nz, liaisons):
    """ write a vmd script in order to visualize the structure and the cell edges """

    pwd = os.getcwd()

    script = """# VMD script
# I think it is a minimal script
"""

    # molecule 1 => atoms
    if liaisons:
        script += "mol new " + pwd + "/" + prefix + "_geom.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"
        script += "mol addfile " + pwd + "/" + prefix + "_geom.xyz type xyz first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"
    else:
        script += "mol new " + pwd + "/" + prefix + "_geom.xyz type xyz first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"

    script += """
mol delrep 0 top
mol representation CPK 1.000000 0.300000 20.000000 6.000000
mol color Element
mol selection {all}
mol material Edgy
mol addrep top
"""
    script += "mol rename top \"" + prefix + "\"\n"
    script += "# done with molecule 0\n"

    # info cristal : parametres de maille pour faire des supercell
    script += "# lattice parameters\n"
    script += "molinfo 0 set a %8.4f"     % cristal.a + "\n"
    script += "molinfo 0 set b %8.4f"     % cristal.b + "\n"
    script += "molinfo 0 set c %8.4f"     % cristal.c + "\n"
    script += "molinfo 0 set alpha %8.2f" % (cristal.alpha ) + "\n"
    script += "molinfo 0 set beta %8.2f"  % (cristal.beta ) + "\n"
    script += "molinfo 0 set gamma %8.2f" % (cristal.gamma ) + "\n"

    # pbcTools plugin to draw cell edge
    script += "# draw cell edges\n"
    script += "pbc box -style tubes -width 1 -color silver"
    
    vmdfile = open( prefix + ".vmd", "w")
    vmdfile.write( script )
    vmdfile.close()

# ----------------------------------------------------------------------------------

def makeBonds(xyz, listeAtomes, prefix):
    """ write a psf file in order to set the connectivity between atoms """

    print("\n# add bonds between atoms : ")

    # 
    # calcul de la matrice des distances
    # 
    Nat = len(listeAtomes)
    distances = np.zeros((Nat,Nat))

    for iat in range(Nat):
        for j in range(Nat-iat-1):
            jat = j + iat + 1
            xi = np.array( xyz[iat] )
            xj = np.array( xyz[jat] )

            xij = xi - xj
            distances[iat,jat] = np.sqrt( ( xij**2 ).sum() )

    # 
    # types d'atomes et de liaisons
    # 
    atomTypes = [listeAtomes[0]]
    NtypesAtomes = 1
    for atom in listeAtomes:
        if atom not in atomTypes:
            atomTypes.append(atom)
            NtypesAtomes += 1

    if len(atomTypes) != NtypesAtomes:
        print("len(atomTypes) = " + str(len(atomTypes)))
        print("NtypesAtomes   = " + str(NtypesAtomes))
        die("Error on atom types")

    NtypesLiaisons = NtypesAtomes * (NtypesAtomes + 1) / 2
    print("Nbre of atom types  : " + str(NtypesAtomes) )
    print("Nbre of bond types  : " + str(NtypesLiaisons) )

    # 
    # distance critique pour chaque type de liaison
    # 
    print("\nBond length cut-off for each bond type :")
    distancesCritiques = np.zeros( (NtypesAtomes,NtypesAtomes) )
    for i in range(NtypesAtomes):
        for j in range(i, NtypesAtomes):
            liaison = " * bond : " + atomTypes[i].ljust(2) + " -- " + atomTypes[j].ljust(2)
            distancesCritiques[i,j] = raw_input(liaison.ljust(20) + " Rc = "  )
            distancesCritiques[j,i] = distancesCritiques[i,j]

    # 
    # Matrice des liaisons : 1 liaison 0 pas de liaison
    # 
    matLiaisons = np.zeros((Nat,Nat))
    nbreLiaisons = 0
    for iat in range(Nat):
        for jat in range(iat + 1, Nat):
            typei = atomTypes.index(listeAtomes[iat])
            typej = atomTypes.index(listeAtomes[jat])

            if distances[iat,jat] < distancesCritiques[typei,typej]:
                matLiaisons[iat,jat] = 1
                nbreLiaisons += 1

    # 
    # ecriture du fichier psf
    # 
    ligne  = "PSF " + prefix + " \n"
    ligne += "\n"
    ligne += "       1 !NTITLE\n"
    ligne += " REMARKS topology of the cell or the supercell\n"
    ligne += "\n"

    ligne += "%8d" % Nat + " !NATOM\n"
    for iat in range(Nat):
        atome = listeAtomes[iat]
        ligne += "%8d" % (iat+1)                  # atome ID
        ligne += "%7d" % 1                        # residue ID
        ligne += "         " + atome.ljust(3)     # type atome 
        ligne += "  " + atome.ljust(3)            # type atome 
        ligne += "%12.6f" % 0.0                   # charge
        ligne += "%14.4f" % 0.0                   # masse
        ligne += "%12d" % 0                       # unused
        ligne += "\n"

    ligne += "\n"
    ligne += "%8d" % nbreLiaisons + " !NBOND: bonds\n"
    nl = 0
    for iat in range(Nat):
        for jat in range(iat + 1, Nat):
            if matLiaisons[iat,jat] == 1:
                ligne += "%8d" % (iat+1) + "%8d" % (jat+1)
                nl += 2
                if nl % 8 == 0:
                    ligne += "\n"

    ligne += "\n"

    print("# Write connectivity between atoms           : " + prefix + "_geom.psf")
    psf = open( prefix + "_geom.psf", "w")
    psf.write(ligne)
    psf.close()

# ----------------------------------------------------------------------------------

def readStructure(fichier, formatPoscar, initial, final):
    """ read the crystal structure """

    if formatPoscar:
        # lecture d'un fichier POSCAR ou CONTCAR
        struct = crystal.Crystal.fromPOSCAR(open( fichier, "r").readlines(), verbose = False)

        print("\n# file " + fichier.strip() + " read.")
        print(  "    type POSCAR / CONTCAR")
        print(  "    " + str(struct.Natoms) + " atoms")

    else:
        # lecture sur le fichier xml de vasp
        run = vasptools.VaspRun(fichier, verbose = False)
        if initial:
            print("\n# Read the initil geometry")
            run.getInitialStructure(verbose=False)
            struct = run.initialStructure

        elif final:
            print("\n# Read the final geometry")
            run.getFinalStructure(verbose=False)
            struct = run.finalStructure

    return struct

# ----------------------------------------------------------------------------------

def visVMD():
    """ main program """

    #
    # default value of options
    #
    fichier      = "vasprun.xml"
    unitcell     = False
    supercell    = False
    nx = ny = nz = 1
    final        = True
    initial      = False
    formatPoscar = False
    prefix       = "vis"
    liaisons     = False

    # 
    # manage options
    # 
    args = sys.argv
    narg = len(args)
    for i in range(narg):
        if (args[i] == "-p" or args[i] == "--prefix") and narg >= i+1:
            prefix = args[i+1]

        elif (args[i] == "-f" or args[i] == "--file") and narg >= i+1:
            fichier = args[i+1]

        elif (args[i] == "-n" or args[i] == "--supercell") and narg >= i+1:
            valeur = args[i+1]
            nx, ny, nz = [ int(val) for val in valeur.split("x") ]
            supercell = True

        elif args[i] == "-u" or args[i] == "--unitcell":
            unitcell = True

        elif args[i] == "-l" or args[i] == "--liaisons":
            liaisons = True

        elif args[i] == "--poscar":
            formatPoscar = True

        elif args[i] == "-i" or args[i] == "--initial":
            initial = True

        elif args[i] == "-e" or args[i] == "--final":
            final = True

        elif args[i] == "-h" or args[i] == "--help":
            usage(0)

    # tests
    if not os.path.exists( fichier ):
        die("Error : file {0} does not exist !\n".format(fichier))

    if formatPoscar and ".xml" in fichier:
        print("Error: you said that the input file is a POSCAR/CONTCAR file : " + fichier)
        usage(1)

    if formatPoscar and (initial or final):
        print("WARNING : -i and -e options make no sense when reading a POSCAR/CONTCAR file")
        print("          see visVMD -h\n")

    # 
    # lecture de la structure
    # 
    struct = readStructure(fichier, formatPoscar, initial, final)

    print("")

    #
    # make supercell if needed
    #
    if supercell:
        struct = struct.makeSupercell(nx, ny, nz)
        coordXYZTotales = struct.XYZCoord
        listeAtomes = struct.atomNames

    # 
    # add atoms if needed
    # 
    if unitcell:
        atomesSupXYZ, listeAtomes = completeMaille(struct)
        coordXYZTotales = struct.XYZCoord + atomesSupXYZ
        listeAtomes = struct.atomNames + listeAtomes

    if not unitcell and not supercell:
        coordXYZTotales = struct.XYZCoord
        listeAtomes = struct.atomNames

    #
    # creation du fichier psf pour la topologie
    #
    if liaisons:
        makeBonds(coordXYZTotales, listeAtomes, prefix)

    #
    # Ecriture des coordonnees des atomes
    #
    print("# Write atom coordinates                : " + prefix + "_geom.xyz")
    geom = open(prefix + "_geom.xyz", "w")
    geom.write( str(len(coordXYZTotales)) + "\n")
    geom.write(struct.name + "\n" )
    for iat, coord in enumerate(coordXYZTotales):
        ligne = listeAtomes[iat].ljust(5)
        for x in coord:
            ligne += " %10.4f" % x
        ligne += "\n"
        geom.write(ligne)
    geom.close()

    # 
    # ecriture du script vmd
    # 
    print("# Write vmd script for visualization    : " + prefix + ".vmd")
    makeVMDScript( prefix, struct, nx, ny, nz, liaisons )

    # 
    # execution de VMD
    # 
    print("\nrun vmd : vmd -e " + prefix + ".vmd")
    os.system("xterm -e \"vmd -e " + prefix + ".vmd\"" )

# ----------------------------------------------------------------------------------

if __name__ == "__main__":
    visVMD()

# ----------------------------------------------------------------------------------
# end
# ----------------------------------------------------------------------------------
