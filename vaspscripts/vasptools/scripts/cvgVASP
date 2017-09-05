#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
cvgVASP
-------

NAME     
        convVASP - get convergence data from an OUTCAR file of VASP

SYNTAX
        convVASP [OPTIONS] ... [FILES]

DESCRIPTION
        convVASP get convergence data on the OUTCAR file and either write them into files
        in order to plot them with your favorite soft or plot them directly using
        matplotlib. The last argument has to be an OUTCAR file, if absent, './OUTCAR' will
        be used.

        -h, -help
            print this help

        -t, --tofiles
            insteed of ploting them directly, convergence data are printed into files.
"""

import os
import sys
from math import sqrt, fabs
import matplotlib.pyplot as plt

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

def cvgVASP():
    """ suivi de la convergence d'un calcul VASP """

    #
    # gestion des options
    #
    args = sys.argv
    narg = len(args)

    # defauts :
    output    = "OUTCAR"
    printData = False

    if narg > 1:
        for i in range(narg)[1:]:
            if args[i] == "-h" or args[i] == "-help":
                usage(0)

            elif args[i] == "-t" or args[i] == "--tofiles":
                printData = True

            else:
                if i == narg - 1:
                    output = args[i]
                else:
                    print(args)
                    print("Error : bad arguments {0}".format(args[i]))
                    usage(1)


    # open file OUTCAR
    if os.path.exists( output ):
        outcar = open( output, "r")
    else:
        die("Error : file {0} not found".format(output))

    print("\n * Read data in file : {0}".format(output))

    #
    # read general data
    #
    end = False
    while not end:
        line = outcar.readline()
        if "Iteration" in line:
            end = True
        elif "EDIFFG" in line:
            ediffg = float(line.split()[2])
        elif "NIONS" in line:
            Nions = int(line.split()[11])
        elif "NSW" in line:
            NSW = int(line.split()[2])
        elif line == "":
            break

    if not end:
        die("\nError : check file {0} :".format(output))


    # list declarations
    SCFconv = list()
    SCFconv.append(list())

    fmax = list()
    fave = list()
    fav2 = list()

    freeEnergy = list()

    #
    # get convergence data
    #
    end = False
    insw = 0
    while not end:
        line = outcar.readline()

        # look at end of file
        if line == "":
            end = True
            if NSW != 0:
                del SCFconv[insw]
            continue

        # SCF convergence
        if "energy without entropy" in line:
            line = line.strip().split("=")
            try:
                SCFconv[insw].append(float(line[1].split()[0]))
            except:
                print("SCF illisible : " + line[1].split()[0])

        # ionic convergence
        if "aborting loop because EDIFF is reached" in line and NSW != 0:

            # look at total force acting on ions
            alire = 0
            while alire != 2:
                line = outcar.readline()
                if "free  energy" in line:
                    # free energy of the structure
                    alire += 1
                    freeEnergy.append( float(line.split()[4]) )

                elif "TOTAL-FORCE" in line:
                    # look at maximal force and compute average
                    alire += 1

                    outcar.readline()   # read dashed line
                    fmax.append(-1.0)
                    fave.append(0.0)
                    fav2.append(0.0)
                    for iat in range(Nions):
                        val = [ float(v) for v in outcar.readline().split()[3:6]]
                        f = sqrt( val[0]**2 + val[1]**2 + val[2]**2 )
                        fave[insw] += f
                        fav2[insw] += f**2
                        if f > fmax[insw]:
                            fmax[insw] = f

                    fave[insw] /= float(Nions)
                    fav2[insw] /= float(Nions)
                    fav2[insw] = sqrt( fav2[insw] - fave[insw]**2)

                elif line == "":
                    # end of file
                    end = True
                    break

            # end of file
            if end:
                end = False # required accuracy not reached
                break


            # new ionic step
            insw += 1
            SCFconv.append(list())

    #
    # special treatment if calculation not finished yet
    #
    if not end:
        # required accuracy not reached, check vectors length and delete data in order to have the
        # same length for all vectors
        n = 10*NSW
        n = min(n,insw)
        n = min(n,len(SCFconv))
        n = min(n,len(freeEnergy))
        while n < len(SCFconv): 
            k = len(SCFconv)-1
            del SCFconv[k]
        while n > len(freeEnergy): 
            k = len(freeEnergy)-1
            del freeEnergy[k]

    #
    # compute ionic convergence data
    #
    if NSW != 0:
        Eionic = list()
        i=0
        for scf in SCFconv :
            Eionic.append( scf[len(scf)-1] )

        delta_F = list()
        delta_F.append(0.)
        for i in range(insw-1):
            delta_F.append( fabs(freeEnergy[i+1] - freeEnergy[i]) )

    #
    # print data into files
    #
    if printData:

        print(" * Print data into files : SCFconv.dat and ionicConv.dat")
        # print SCF data into file
        outSCF = open("SCFcvg.dat", "w")
        outSCF.write("# SCF convergence for each ionic step\n")
        outSCF.write("# i         E (eV)\n")
        for n, scf in enumerate(SCFconv):
            outSCF.write("& step %d" % n + "\n")
            for i, E in enumerate(scf):
                line = "%5d %16.8e\n" % (i,E)
                outSCF.write(line)
        outSCF.close()

        # print ionic convergence data into file
        ionicConv = open("ionicCvg.dat", "w")
        ionicConv.write("# Energie en eV\n")
        ionicConv.write("# Forces en eV/Angst \n")
        ionicConv.write("#   i         fmax             <f>             sigma_f           E" \
                "              Delta_F       seuil\n")
        for i in range(insw):
            line = "%5d %16e %16e %16e %16.8e %16.8e %8.4f\n" % \
                (i,fmax[i],fave[i],fav2[i],Eionic[i],delta_F[i],fabs(ediffg))
            ionicConv.write(line)
        ionicConv.close()

    #
    # make figure plot
    #
    else:

        if NSW == 0:
            # electronic convergence only
            plt.plot(SCFconv[0])
            plt.title("SCF convergence")
            plt.xlabel("Electronic step")
            plt.ylabel("Energie   /   eV")

        else:

            fig = plt.figure()
            fig.subplots_adjust(wspace=0.25)
            fig.suptitle("VASP run convergence")

            # graph1 = total energy
            graph1 = fig.add_subplot(2,2,1)
            graph1.set_ylabel("Energy   /   eV")
            graph1.set_xlabel("ionic step")
            graph1.grid()
            graph1.plot( range(insw), Eionic, "r")

            # graph2 = Free energy variation
            graph2 = fig.add_subplot(2,2,3,sharex=graph1)
            graph2.set_ylabel("$\Delta F$   /   eV")
            graph2.set_xlabel("ionic step")
            graph2.grid()
            graph2.plot( delta_F )
            if ediffg > 0. :
                graph2.plot( [0,insw], [ediffg,ediffg], 'k--')

            # graph3 = forces acting on ions
            graph3 = fig.add_subplot(1,2,2,sharex=graph1)
            graph3.set_ylabel("forces   /   eV/$\AA$")
            graph3.set_xlabel("ionic step")
            graph3.grid()
            graph3.plot( fmax, label="fmax" )
            graph3.plot( fave, label="<f>" )
            graph3.plot( fav2, label="$\sigma_f$" )
            graph3.legend(fancybox=True,shadow=True)
            if ediffg < 0. :
                graph3.plot( [0,insw+1], [-ediffg,-ediffg], 'k--')

        # plot data
        print(" * Plot data")
        plt.show()

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

if __name__ == "__main__":
    cvgVASP()

