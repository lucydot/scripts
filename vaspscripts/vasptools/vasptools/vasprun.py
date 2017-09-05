#!/usr/bin/env python
# -*-coding:utf-8 -*-

# this file is part of vasptools package

"""
vasptools.vasprun
-----------------

This module contains the class VaspRun which is a core part of vasptools.
"""

__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"

import os
import myxml
from crystal import Crystal
from atom import Atom

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
#
# Classe Vasprun : traitement du fichier vasprun.xml
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

class VaspRun(object):
    """ VASP calculation class 

        Arguments :
            * fichier        : vasprun.xml output file of VASP (default is 'vasprun.xml')
            * verbose (bool) : verbosity of methods (default is 'True')

        Example :
            run = VaspRun()
            run = VaspRun("dos.xml")
    """

    def __init__(self, xmlFile = "vasprun.xml", verbose = True):
        """ constructor """
        self.xmlFile = xmlFile
        self.verbose = verbose

        # check xml file
        if os.path.exists(self.xmlFile):

            # read main data
            self.__readKeywords()
            self.__readAtomsData()
            self.getFinalStructure(self.verbose)
        
            # variable de controle
            self.DOSTotaleLue = False
            self.DOSPartiellesLues = False
            self.bandesLues = False
            self.pointsKLues = False

            # affichage de quelques infos
            if self.verbose:
                print("# xml file of the run : {0}\n".format(self.xmlFile))
                self.printINCAR()
                print("")
                self.printAtomsData()

        else:
            print("\n\t\t### file {0} does not exist ###\n".format(xmlFile))
            exit(1)

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    #
    # Méthodes d'interrogation des paramètres du calcul
    #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def valeurMotClef(self, flag):
        """ return the value of the parameter 'flag'

        arguments:
            * flag(string) : parameter of VASP

        see listerMotsClefs()
        """

        # interrogation des dictionnaires
        if flag in self.allMotsClefs.keys():
            valeur = self.allMotsClefs[ flag ]
            if self.verbose:
                print("{0} = {1}".format(flag, valeur))

        else:
            valeur = None
            print("{0} unknown".format(flag))

        return valeur

    def listerMotsClefs(self):
        """ Print all paramters of the calculation.

            see valeurMotClef() """

        # mots clefs du fichier INCAR
        print("# mots clefs du fichier INCAR")
        ligne = ""
        for i, flag in enumerate(sorted(self.INCAR.keys())):
            ligne += flag.rjust(15) + ","
            if (i % 6 == 0 and i != 0) or i == len(self.INCAR.keys()) - 1:
                print(ligne)
                ligne = ""

        # tous les mots clefs du calcul
        print("\n# mots clefs du calcul")
        ligne = ""
        for i, flag in enumerate(sorted(self.allMotsClefs.keys())):
            ligne += flag.rjust(15) + ","
            if (i % 6 == 0 and i != 0) or i == len(self.allMotsClefs.keys()) - 1:
                print(ligne)
                ligne = ""

    def printINCAR(self):
        """ print INCAR keywords """

        print("# fichier INCAR du calcu")
        ligne=""
        for i, flag in enumerate(self.INCAR.keys()):
            ligne += flag.rjust(10) + " = " + str(self.INCAR[flag]).ljust(10) 
            if (i % 3 == 0 and i != 0) or i == len(self.INCAR.keys()) - 1:
                print(ligne)
                ligne=""

    def printAtomsData(self):
        """ print atomic data """

        print("\n# system :")
        print("\t* atom number : " + str(self.Natomes))
        print("\t* type number : " + str(self.Ntypes))
        ligne = "\t* atom list   : "
        for i in range(self.Natomes):
            if i != 0 and i % 10 == 0:
                print(ligne)
                ligne = "\t                "
            ligne += self.atoms[i].name.rjust(3) + ", "
        if ligne.strip() != "": print(ligne[:-2])

        
        print("\n# Atom types :")
        for i in range(self.Ntypes):
            for atom in self.atoms:
                if i+1 == atom.atomType:
                    print(atom)
                    break

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    #
    # Méthodes pour la DOS
    #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def lectureDOS( self ):
        """ Lecture de la densité d'états totale et des densités d'états partielles sur le
        fichier xml du calcul. Après exécution de cette méthode, on dispose de la densité
        d'états totale dans la liste dosTotale et des densités d'états partielles dans la
        liste dosPartielles. Les valeurs d'énergies pour lesquelles on dispose des valeurs
        de la DOS sont dans la liste energiesDOS. Ici <VaspRun> désigne un objet de type 
        VaspRun, les listes crées sont de la forme :

            <VaspRun>.nptsDos => nombre de valeurs pour la DOS

            <VaspRun>.energiesDOS[i]

                i : 0 -> calcul.nptsDos-1
                    valeurs d'énergies pour lesquelles les DOS sont connues.

            <VaspRun>.dosTotale[ spin ][ E ][ i ]

                spin    : 0 ou 1, seule la valeur 0 est disponnible pour les calculs
                    avec spin non polarisé (ISPIN = 1). Pour les calculs spin
                    polarisé (ISPIN = 2), spin = 0 donne la DOS pour le
                    spin alpha et spin = 1 pour le spin beta.

                E    : indice de parcours de la DOS pour un spin (valeurs de
                    l'énergie)

                i       : 0 -> 1
                    i = 0 -> densité d'états 
                    i = 1 -> intégrale de la DOS

            <VaspRun>.dosPartielles[ iat ][ spin ][ E ][ i ]

                iat    : numéro de l'atome

                spin    : 0 ou 1, seule la valeur 0 est disponnible pour les calculs
                    avec spin non polarisé (ISPIN = 1). Pour les calculs spin
                    polarisé (ISPIN = 2), spin = 0 donne la DOS pour le
                    spin alpha et spin = 1 pour le spin beta.

                E    : indice de parcours de la DOS pour un spin (valeurs de
                    l'énergie)

                i       : 0 -> 2 ou 8 valeurs de la DOS sur chaque sous couche ou OA
                    0 s, 1 p, 2 d
                    0 s, 1 py, 2 pz, 3 px, 4 dxy, 5 dyz, 6 dz2, 7 dxz, 8 dx2-y2
        """

        # spin polarise ?
        ISPIN = self.allMotsClefs["ISPIN"]

        if self.verbose:
            print("# Lecture de la densité d'états")

            # calcul spin polarisé ?
            print("\t* ISPIN = {0}".format(ISPIN))

        # bloc dos
        fxml = open(self.xmlFile, "r")
        dos = myxml.getBlocs("dos", fxml, onlyfirst = True)
        fxml.close()

        # test presence de la dos
        if dos == None:
            print("Pas de densité d'états dans ce calcul")
        else:
            if self.verbose:
                print("\t* densité d'états totale")

            # niveau de Fermi
            for ligne in dos:
                att = myxml.getClefs(ligne, ["name"])
                if att["name"] != None and att["name"] == "efermi":
                    self.eFermi = float(myxml.getNodeData(ligne)[0])
                    print("\t* Fermi level = {0} eV (warning : accuracy depend on k-points grid)".format(self.eFermi))
                    break


            #                DOS TOTALE

            # lecture de la dos totale
            blocDosTotale = myxml.getBlocs("total", dos, onlyfirst = True)
            self.dosTotale = list()
            self.energiesDOS = list()

            # cas spin non polarise
            blocDosTotaleSpin1 = myxml.getBlocs("set", blocDosTotale[1:-1], {"comment":"spin 1"}, True)
            self.dosTotale.append(list())
            for ligne in blocDosTotaleSpin1[1:-1]:
                valeurs = [float(val) for val in myxml.getNodeData(ligne)]
                # valeurs[0] est l'energie
                self.energiesDOS.append(valeurs[0])
                self.dosTotale[0].append(valeurs[1:])

            if ISPIN == 2:
                # cas spin polarise
                blocDosTotaleSpin2 = myxml.getBlocs("set", blocDosTotale[1:-1], {"comment":"spin 2"}, True)
                self.dosTotale.append(list())
                for ligne in blocDosTotaleSpin2[1:-1]:
                    valeurs = [float(val) for val in myxml.getNodeData(ligne)]
                    # valeurs[0] est l'energie
                    self.dosTotale[1].append(valeurs[1:])

            self.DOSTotaleLue = True

            #
            #                DOS PARTIELLES
            #

            blocDosPartielle = myxml.getBlocs("partial", dos, onlyfirst = True )

            if blocDosPartielle == None:
                print("\t ! pas de densité d'états partielles !")
            else:
                if self.verbose:
                    print("\t* densités d'états partielles")

                self.dosPartielles = list()

                # boucle sur les atomes du calcul
                for iat in range(self.Natomes):
                    ion = "ion " + str(iat+1)
                    blocDosIon = myxml.getBlocs("set", blocDosPartielle, {"comment":ion}, True)
                    self.dosPartielles.append( list() )

                    blocDosIonSpin1 = myxml.getBlocs( "set", blocDosIon[1:-1], {"comment":"spin 1"}, True)
                    self.dosPartielles[iat].append( list() )

                    # lecture des dos projetees spin 1
                    for ligne in blocDosIonSpin1[1:-1]:
                        valeurs = [ float(val) for val in myxml.getNodeData( ligne ) ]
                        # valeurs[0] est l'energie
                        self.dosPartielles[iat][0].append( valeurs[1:] )

                    if ISPIN == 2:
                        self.dosPartielles[iat].append( list() )
                        blocDosIonSpin2 = myxml.getBlocs( "set", blocDosIon[1:-1], {"comment":"spin 2"}, True)

                        # lecture des dos projetees spin 2
                        for ligne in blocDosIonSpin2[1:-1]:
                            valeurs = [ float(val) for val in myxml.getNodeData( ligne ) ]
                            # valeurs[0] est l'energie
                            self.dosPartielles[iat][1].append( valeurs[1:] )

                self.DOSPartiellesLues = True

        return [self.DOSTotaleLue, self.DOSPartiellesLues]

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    #
    # Méthodes pour les bandes d'énergies
    #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def lectureBandes(self):
        """ Méthode permettant de lire les bandes d'énergie sur le fichier xml du calcul.
        Après exécution de cette méthode, on dispose des bandes d'énergie dans la liste
        Bandes sous la forme :

            <VaspRun>.bandes[ spin ][ point k ][ bandes ][ i ]

                spin    : 0 ou 1, seule la valeur 0 est disponnible pour les calculs
                    avec spin non polarisé (ISPIN = 1). Pour les calculs spin
                    polarisé (ISPIN = 2), spin = 0 donne les bandes pour les
                    spin alpha et spin = 1 pour les spin beta.

                point k : 0 -> nombre de point K
                    le nombre de point k est calcul.nbrePointsK
                    calcul.nbreLignesPointsK : nombre de directions de points k
                    calcul.Ndivision         : nombre de points k par direction
                    Les points k sont listés direction après direction en
                    donnant l'ensemble des points k sur chaque direction

                bandes  : 0 -> NBANDS-1

                i       : 0 ou 1
                    i = 0 -> valeurs de l'énergie
                    i = 1 -> occupation de la bande
        """

        # lecture du niveau de fermi
        fxml = open(self.xmlFile, "r")
        dos = myxml.getBlocs("dos", fxml, onlyfirst = True)
        fxml.close()

        # test presence de la dos
        if dos == None:
            print("Density of state and fermi level not found in xml file")
            self.eFermi = 0.0
        else:
            # niveau de Fermi
            for ligne in dos:
                att = myxml.getClefs(ligne, ["name"])
                if att["name"] != None and att["name"] == "efermi":
                    self.eFermi = float(myxml.getNodeData(ligne)[0])
                    break

        # bloc eigenvalues = bandes d'energie
        fxml = open(self.xmlFile, "r")
        eigenvalues = myxml.getBlocs("eigenvalues", fxml, onlyfirst = True)
        fxml.close()
        
        if eigenvalues == None:
            print("\nerreur : bandes d'énergie introuvable\n")
            exit(1)

        # infos sur le calcul
        NBANDS = self.allMotsClefs["NBANDS"]
        ISPIN  = self.allMotsClefs["ISPIN"]

        # lecuture des points K
        self.lecturePointsK()
        nbreKpoints = len(self.listePointsK)

        if self.verbose:
            print("# Lecture des bandes d'énergies")
            print("\t* Fermi level        = {0} eV (warning accuracy depend on k-points grid)".format(self.eFermi))
            print("\t* ISPIN              = {0}".format(ISPIN))
            print("\t* nombre de bandes   = {0}".format(NBANDS))
            print("\t* nombre de points K = {0}".format(nbreKpoints))

        setMain = myxml.getBlocs("set", eigenvalues, onlyfirst = True) 
        blocSpin = list()
        blocSpin.append( myxml.getBlocs("set", setMain, {"comment":"spin 1"}, True))
        if ISPIN == 2:
            blocSpin.append( myxml.getBlocs("set", setMain, {"comment":"spin 2"}, True))

        # lecture des bandes
        self.bands = list()
        for spin in range(ISPIN):

            self.bands.append( list() )

            # boucle sur les points k pour un spin
            for k in range(nbreKpoints):

                self.bands[spin].append( list() )

                kpoint = "kpoint " + str(k+1)
                block = myxml.getBlocs("set", blocSpin[spin], {"comment":kpoint }, True)

                # boucle sur les valeurs des bandes sur un points k
                for ligne in block:
                    nom = myxml.getNodeName(ligne)
                    if nom == "r":
                        valeurs = [float(val) for val in myxml.getNodeData(ligne) ]
                        self.bands[spin][k].append( valeurs )

                if len(self.bands[spin][k]) != NBANDS:
                    print("\nerreur nombre de bandes incorrect\n")
                    exit(1)

        # controle de lecture
        self.bandesLues = True
        return self.bandesLues

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    #
    # Lecture des points K
    #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def lecturePointsK(self):
        """ lit les informations sur les point K du calcul """

        # lecture des points k
        fxml = open(self.xmlFile, "r")
        kpoints = myxml.getBlocs("kpoints", fxml, onlyfirst = True )
        fxml.close()

        # grille de points K
        generation = myxml.getBlocs("generation", kpoints, onlyfirst = True)
        if generation == None:
            self.typePointsK = "explicit"
        else:
            att = myxml.getClefs(generation[0], ["param"])
            self.typePointsK = att["param"]

        if self.typePointsK == "listgenerated":
            ndir = 0
            liste = list()
            for ligne in generation:
                nom = myxml.getNodeName(ligne)
                if nom == "i":
                    self.Ndivision = int(myxml.getNodeData(ligne)[0])

                if nom == "v":
                    valeurs = [float(val) for val in myxml.getNodeData(ligne)]
                    liste.append(valeurs)
                    ndir += 1

            # liste des directions
            self.directionsPointsK = list()
            for d in range(ndir - 1):
                self.directionsPointsK.append(liste[d] + liste[d + 1])

        else:
            print("Explicit K points grid")

        # kpointlist
        varray = myxml.getBlocs( "varray", kpoints, {"name" : "kpointlist"}, True)
        self.listePointsK = list()
        for ligne in varray[1:-1]:
            valeurs = [float(val) for val in myxml.getNodeData(ligne)]
            self.listePointsK.append(valeurs)

        if self.typePointsK == "listgenerated":
            if len(self.directionsPointsK) * self.Ndivision != len(self.listePointsK):
                raise ValueError("k-points number unconsistent")
        else:
            self.directionsPointsK = [self.listePointsK[0] + self.listePointsK[-1]]
            self.Ndivision = len(self.listePointsK)


        # controle de lecture
        self.pointsKLues = True
        return self.pointsKLues

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    #
    # Lecture de la structure initiale et finale
    #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def getInitialStructure(self,verbose=True):
        """ read initial structure """

        # read initial pos
        fxml = open(self.xmlFile, "r")
        structure = myxml.getBlocs("structure", fxml, {"name":"initialpos"}, True)
        fxml.close()

        # read lattice vectors
        latticeParam = myxml.getBlocs("varray", structure, {"name":"basis"}, True)
        veca = [ float(val) for val in myxml.getNodeData(latticeParam[1]) ]
        vecb = [ float(val) for val in myxml.getNodeData(latticeParam[2]) ]
        vecc = [ float(val) for val in myxml.getNodeData(latticeParam[3]) ]

        # crystal object
        self.initialStructure = Crystal(veca = veca, vecb = vecb, vecc = vecc, \
            name = "Initial structure : {0}".format(self.xmlFile))

        # print lattice parameters
        if verbose:
            print("\t* a = %10.5f \t* alpha = %10.3f" % \
                (self.initialStructure.a, self.initialStructure.alpha))
            print("\t* b = %10.5f \t* beta  = %10.3f" % \
                (self.initialStructure.b, self.initialStructure.beta) )
            print("\t* c = %10.5f \t* gamma = %10.3f" % \
                (self.initialStructure.c, self.initialStructure.gamma))

        # read reduce coordinates and compute cartesian coordinates
        positions = myxml.getBlocs("varray", structure, {"name":"positions"}, True)
        self.initialStructure.Natoms = 0
        for ligne in positions[1:-1]:
            pos = [float(val) for val in myxml.getNodeData(ligne)]
            self.initialStructure.redCoord.append(pos)
            self.initialStructure.Natoms += 1

        self.initialStructure.computeXYZCoord()

        if self.initialStructure.Natoms != self.Natomes:
            print("Natoms    : {0}".format(self.Natomes))
            print("Structure : {0}".format(self.initialStructure.Natoms))
            print("Error     : atom number")
            exit(1)

        # atom names
        for iat in range(self.Natomes):
            self.initialStructure.atomNames.append(self.atoms[iat].name)

        return self.initialStructure

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def getFinalStructure(self,verbose=True):
        """ read final structure """

        # read final pos
        fxml = open(self.xmlFile, "r")
        structure = myxml.getBlocs("structure", fxml, {"name":"finalpos"}, True)
        fxml.close()

        if self.verbose:
            print("\n# Read final structure")

        if structure == None:
            # there is not final structure (run not terminated)
            # <structure>...</structure>

            print("\t* Warning: final structure not found, read last one instead")

            fxml = open(self.xmlFile, "r")
            blocsStructures = myxml.getBlocs("structure", fxml)
            fxml.close()

            Nstructures = len(blocsStructures)

            if Nstructures == 0:
                # aucune structure ?
                print("\nError : there is not any structure in this xml file!\n")
                exit(1)

            elif Nstructures == 1:
                # only one structure => the initial one
                ligne = blocsStructures[0][0]
                att = myxml.getClefs(ligne, ["name"])
                if att["name"] != None and att["name"] == "initialpos":
                    print("\t* Warning : I found only the initial structure")
                else:
                    print("\t* Warning : I found only one structure which is not the initial one ??")

                structure = blocsStructures[0]

            else:
                structure = blocsStructures[Nstructures-1]

        # read lattice vectors
        latticeParam = myxml.getBlocs("varray", structure, {"name":"basis"}, True)
        veca = [ float(val) for val in myxml.getNodeData(latticeParam[1]) ]
        vecb = [ float(val) for val in myxml.getNodeData(latticeParam[2]) ]
        vecc = [ float(val) for val in myxml.getNodeData(latticeParam[3]) ]

        # definition du cristal
        self.finalStructure = Crystal(veca = veca, vecb = vecb, vecc = vecc, \
            name = "Final structure : {0}".format(self.xmlFile))

        # print lattice parameters
        if verbose:
            print("\t* a = %10.5f \t* alpha = %10.3f" % (self.finalStructure.a, self.finalStructure.alpha))
            print("\t* b = %10.5f \t* beta  = %10.3f" % (self.finalStructure.b, self.finalStructure.beta) )
            print("\t* c = %10.5f \t* gamma = %10.3f" % (self.finalStructure.c, self.finalStructure.gamma))

        # read reduce coordinates and compute cartesian coordinates
        positions = myxml.getBlocs( "varray", structure, {"name":"positions"}, True)
        self.finalStructure.Natoms = 0
        for ligne in positions[1:-1]:
            pos = [ float(val) for val in myxml.getNodeData(ligne) ]
            self.finalStructure.redCoord.append(pos)
            self.finalStructure.Natoms += 1

        self.finalStructure.computeXYZCoord()

        if self.finalStructure.Natoms != self.Natomes:
            print("Natomes   : {0}".format(self.Natomes))
            print("Structure : {0}".format(self.finalStructure.Natomes))
            print("Error : atomic number unconsistant")
            exit(1)

        # atom names
        for iat in range(self.Natomes):
            self.finalStructure.atomNames.append(self.atoms[iat].name)

        return self.finalStructure

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    #
    # Méthodes internes
    #
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def __readKeywords(self):
        """ read xmlFile and store all keywords and their values """

        if self.verbose:
            print("\t* Read keywords of the run")

        # # #        parameters keywords

        # bloc parameters
        fxml = open( self.xmlFile, "r")
        parameters = myxml.getBlocs("parameters", fxml, onlyfirst = True)
        fxml.close()

        # dictionnary des parametres du calcul
        self.allMotsClefs = {}

        for ligne in parameters:
            if "<i" in ligne:
                nomFlag, valeurFlag = self.__lectureFlag(ligne)
                self.allMotsClefs[nomFlag] = valeurFlag

            elif "<v" in ligne:
                nomFlag, valeurFlag = self.__lectureFlag(ligne, vecteur=True)
                self.allMotsClefs[nomFlag] = valeurFlag

            else:
                continue


        # # #        mots clefs du bloc incar

        # bloc incar
        fxml = open(self.xmlFile, "r")
        incar = myxml.getBlocs("incar", fxml, onlyfirst = True)
        fxml.close()

        # dictionnaire INCAR
        self.INCAR = {}

        for ligne in incar[1:-1]:
            if "<i" in ligne:
                nomFlag, valeurFlag = self.__lectureFlag(ligne)
                self.INCAR[ nomFlag ] = valeurFlag

            elif "<v" in ligne:
                nomFlag, valeurFlag = self.__lectureFlag(ligne, vecteur=True)
                self.INCAR[ nomFlag ] = valeurFlag

            else:
                continue

        self.allMotsClefs = dict(self.allMotsClefs.items() + self.INCAR.items())

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def __lectureFlag( self, node, vecteur=False ) :
        """ lecture d'un parametre du calcul suivant son type """

        valeur = myxml.getNodeData( node )

        # clefs
        clefs = myxml.getClefs( node, ["type", "name"] )
        nomFlag = clefs["name"]
        typeFlag = clefs["type"]

        if len(valeur) != 0:
            if typeFlag == "string" or typeFlag == "logical":
                if vecteur :
                    valeurFlag = " ".join(valeur)
                else :
                    valeurFlag = valeur[0]

            elif typeFlag == "int" :
                if vecteur:
                    valeurFlag = " ".join(valeur)
                else :
                    valeurFlag = int(valeur[0])

            elif typeFlag == None :
                if vecteur:
                    valeurFlag = " ".join(valeur)
                else :
                    valeurFlag = float(valeur[0])
        else:
            valeurFlag = "empty"

        return nomFlag, valeurFlag

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def __readAtomsData( self ) :
        """ Read bloc 'atominfo'. It contains the atom list and their description. """

        if self.verbose:
            print("\t* Read atom data") 
        
        # on recupere le bloc atominfo
        fxml = open(self.xmlFile, "r")
        atominfo = myxml.getBlocs("atominfo", fxml, onlyfirst = True)
        fxml.close()

        # nombre d'atomes et nombre de types
        i = 0
        for ligne in atominfo:
            if "<atoms>" in ligne : 
                self.Natomes = int(myxml.getNodeData(ligne)[0])
                i += 1

            if "<types>" in ligne : 
                self.Ntypes  = int(myxml.getNodeData(ligne)[0])
                i += 1

            if i == 2 : break

        # atom list
        self.atoms = list()

        # read atom type
        lignesArrayTypes = myxml.getBlocs("array", atominfo, {"name":"atomtypes"}, True)
        self.typesAtomes = list()
        debutListe = False
        for ligne in lignesArrayTypes:
            if "<set>" in ligne:
                debutListe = True
                continue

            if "</set>" in ligne:
                break

            if debutListe:
                ligne = ligne.replace("<rc><c>", "")
                ligne = ligne.replace("</c></rc>", "")
                valeurs = ligne.split("</c><c>")

                tmpdic = {"nom":valeurs[1].strip(), "masse":float(valeurs[2]), \
                          "valence":float(valeurs[3]), "pseudo":valeurs[4] }
                self.typesAtomes.append(tmpdic)

        # lecture des atomes
        lignesArrayAtomes = myxml.getBlocs("array", atominfo, {"name":"atoms"}, True)
        debutListe = False
        for ligne in lignesArrayAtomes:
            if "<set>" in ligne:
                debutListe = True
                continue

            if "</set>" in ligne:
                break

            if debutListe :
                ligne = ligne.replace("<rc><c>", "")
                ligne = ligne.replace("</c></rc>", "")
                valeurs = ligne.split("</c><c>")

                nom = valeurs[0].strip()
                typ = int(valeurs[1])

                if nom != self.typesAtomes[typ-1]["nom"]:
                    print("non concordance nom / type ")
                    exit(1)

                atom = Atom(name     = self.typesAtomes[typ-1]["nom"], \
                            atomType = typ, \
                            Ne       = self.typesAtomes[typ-1]["valence"], \
                            w        = self.typesAtomes[typ-1]["masse"], \
                            pseudo   = self.typesAtomes[typ-1]["pseudo"] )
                self.atoms.append(atom)
