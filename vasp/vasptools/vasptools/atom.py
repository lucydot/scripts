#!/usr/bin/env python
# -*-coding:utf-8 -*-

# this file is a part of vasptools

""" Atom class """

__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"

class Atom( object ):
    """ Atom class stand for atomic data :

        * name     : name (string)
        * atomType : atom type (int)
        * Ne       : valence electrons (int)
        * w        : weight (float)
        * pseudo   : pseudopotential (string)
    """

    def __init__(self, name = "H", atomType = -1, Ne = 1, w = 1., pseudo = ""):
        """ Constructor : default is hydrogen"""

        self.name = str(name)

        try:
            self.atomType = int(atomType)
        except ValueError:
            raise ValueError("atom type are labeled by integers")

        try:
            self.Ne = int(Ne)
        except ValueError:
            raise ValueError("nElecValence must be an integer")

        try:
            self.w = float(w)
        except ValueError:
            raise ValueError("weight must be a number")

        self.pseudo       = str(pseudo)
    
    def __repr__(self):
        """ représentation d'un atome """
        print("Atom {0} : ".format(self.name) )
        print("\tweight            : {0}".format(self.w))
        print("\tValence electrons : {0}".format(self.Ne))
        print("\ttype              : {0}".format(self.atomType))
        print("\tPseudopotential   : {0}".format(self.pseudo))
        return ""

