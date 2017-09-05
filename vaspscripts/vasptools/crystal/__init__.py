#!/usr/bin/env python
# -*- coding=utf-8 -*-

""" The aim of the crystal module is to be able to easily manipulate all lattice 
parameters of a crystal and do some basic operation on atomic coordinate such as :

    * Compute lattice vectors, volume
    * Coordinates transformation
    * input or output in VASP or CRYSTAL format

Today, symmetry is not considered. See python ASE for symetry support.

Requirements :
    numpy module is used internally in several functions or methods.

WARNING : 
    please, paid attention to angle unit. All angles are in degrees, they are
    internally converted to radian when needed.


In the following we use a Cu2O crystal as an example. 

    >>> import crystal
    >>> Cu2O = crystal.fromPOSCAR("POSCAR")
    --------------------------------------------------
     lattice parameters : Cu2O Bulk
    --------------------------------------------------
     a     =    4.26000
     b     =    4.28000
     c     =    4.27000
     alpha =     90.000
     beta  =     90.000
     gamma =     90.000
    --------------------------------------------------

""" 

import doctest
import copy
import array

import numpy as np
from numpy import sqrt, arccos, fabs, pi, cos, sin

__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__date__ = "Mardi 16 Mars 2012"
__licence__ = "GPL"

# -----------------------------------------------------------------------------
#  TODO
# -----------------------------------------------------------------------------
#
# toXYZ
#
# -----------------------------------------------------------------------------

bravais = [ "triclinic", "monoclinic", "orthorhombic", "tetragonal",
            "rhombohedral", "hexagonal", "cubic"] 

# central cube and 26 neighbors
trans = np.array([[0. , 0., 0.],
                  [ 1., 0., 0.],
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
                  [ 0.,-1.,-1.]])

class Crystal(object):
    """ class *Crystal* allows to describe and do operations on a crystal lattice 

    :param float a: lattice parameter a.
    :param float b: lattice parameter b.
    :param float c: lattice parameter c.
    :param float alpha: lattice parameter alpha.
    :param float beta: lattice parameter beta.
    :param float gamma: lattice parameter gamma.
    :param string name: crystal name.
    :param string lattice: bravais lattice name.
    :param list veca: lattice vector a (length = 3).
    :param list vecb: lattice vector b (length = 3).
    :param list vecc: lattice vector c (length = 3).
    :param integer Z: number of unti formula

    If you give the bravais lattice name you may not set all lattice parameters :

    >>> cry = crystal.Crystal(a = 3., lattice = "cubic")
    >>> print(cry)
    --------------------------------------------------
    lattice parameters : 
    --------------------------------------------------
    a     =    3.00000
     b     =    3.00000
    c     =    3.00000
    alpha =     90.000
    beta  =     90.000
    gamma =     90.000
    --------------------------------------------------

    If you give only the bravais lattice, the prompt will ask you parameters needed
    according to the bravais lattice.

    >>> cry = crystal.Crystal(lattice = "hexagonal")
    a = 3.
    c = 7.
    >>> print(cry)
    --------------------------------------------------
     lattice parameters : 
    --------------------------------------------------
     a     =    3.00000
     b     =    3.00000
     c     =    7.00000
     alpha =     90.000
     beta  =     90.000
     gamma =    120.000
    --------------------------------------------------

    You can give directly lattice vectors.

    >>> cry = crystal.Crystal(veca = [2., 0., 0.], vecb = [-1., sqrt(3), 0.], vecc = [0., 0., 4.])
    >>> cry.lattice
    'hexagonal'
    >>> print(cry)
    --------------------------------------------------
     lattice parameters : 
    --------------------------------------------------
     a     =    2.00000
     b     =    2.00000
     c     =    4.00000
     alpha =     90.000
     beta  =     90.000
     gamma =    120.000
    --------------------------------------------------

    You cannot set lattice parameters **and** lattice vectors at the same time. 
    Only one of them can be set when you create a crystal object in order to ensure
    consistency between lattice vectors, lattice parameters and bravais lattice name.
    A method is available in order to print lattice vectors in a cartesian frame.

    >>> cry = crystal.Crystal(lattice = "hexagonal")
    a = 3.
    c = 7.
    >>> cry.printLatticeVectors()
    --------------------------------------------------
     Lattice vectors in cartesian frame : lattice is hexagonal
    --------------------------------------------------
     a =    3.0000000000    0.0000000000    0.0000000000
     b =   -1.5000000000    2.5980762114    0.0000000000
     c =    0.0000000000    0.0000000000    7.0000000000
    --------------------------------------------------

    """

    # -------------------------------------------------------------------------
    #  Constructor
    # -------------------------------------------------------------------------
    def __init__(self, *largs, **args):
        """ constructor """

        #
        # set all attribute to None or False
        #
        self._veca = None
        self._vecb = None
        self._vecc = None
        """ Lattice vectors :math:`\vec{a}`, :math:`\vec{b}` and :math:`\vec{c}` in a cartesian
        frame where :math:`\vec{a}` is along x axes and :math:`\vec{b}` is in the (x,y) plane.

        >>> Cu2O.veca, Cu2O.vecb, Cu2O.vecc
        ([4.26, 0.0, 0.0], [0.0, 4.28, 0.0], [0.0, 0.0, 4.27])
        """


        self._a = None
        self._b = None
        self._c = None
        """ Lattice parameters **a**, **b** and **c**, lengths of the vectors :math:`\vec{a}`, :math:`\vec{b}` 
        and :math:`\vec{c}`.

        >>> Cu2O.a, Cu2O.b, Cu2O.c
        (4.2599999999999998, 4.2800000000000002, 4.2699999999999996) **c**, third vector of the lattice 
        """

        self._alpha = None
        self._beta  = None
        self._gamma = None
        """ Lattice angles :math:`\alpha` , :math:`\beta` and :math:`\gamma`. All angles must be 
        given in degrees and are return in degrees. They are internally convert in radian when needed.

        >>> Cu2O.alpha, Cu2O.beta, Cu2O.gamma
        (90.0, 90.0, 90.0)
        """

        self._volume = None
        """ Volume of the unit cell """

        self._lattice = None
        """ Name of the lattice according to the 7 Bravais lattices. Remind that symmetry is not considered
        for the moment.

        See :py:func:`bravaisLattice` """

        self.name = None
        """ Name of the crystal object """

        self.Z = None
        """ number of unit formula """

        self.Natoms = None
        """ number of atoms in the unit cell """

        self.group = None
        """ Symmetry group. Nevertheless, remind that symmetry is not considered. """

        self.XYZCoord = list()
        """ List of cartesian coordinates of the atoms in the unit cell. This is a list object.

        See :py:func:`computeRedCoord` and :py:func:`red2cart` """

        self.redCoord = list()
        """ List of reduce coordinates of the atoms in the unit cell. This is a list object.

        See :py:func:`computeXYZCoord` and :py:func:`cart2red` """

        self.atomNames = list()
        """ List of atom names of all atoms in the unit cell. """

        self.__matrixSet = False
        """ *bool* control if transformation matrix are set or not. """

        self.verbose = True
        """ *bool* control verbosity of methods. """

        #
        # arguments list
        #
        argList = ["a", "b", "c", "alpha", "beta", "gamma", "veca", "vecb", "vecc",
        "name", "lattice", "verbose", "Z"]

        #
        # check args and set attributes
        #
        if len(largs) != 0:
            raise TypeError("Use a 'key = val' syntax such as Crystal( key1 = val, key2 = val) ")

        for arg in args.keys():
            if arg not in argList:
                raise NameError("name '" + arg + "' unknown")
            else:
                if arg == "a":
                    try:
                        self._a = float(args[arg])
                    except ValueError:
                        raise ValueError("a must be a number")
                elif arg == "b":
                    try:
                        self._b = float(args[arg])
                    except ValueError:
                        raise ValueError("b must be a number")
                elif arg == "c":
                    try:
                        self._c = float(args[arg])
                    except ValueError:
                        raise ValueError("c must be a number")
                elif arg == "alpha":
                    try:
                        self._alpha = float(args[arg])
                    except ValueError:
                        raise ValueError("alpha must be a number")
                elif arg == "beta":
                    try:
                        self._beta = float(args[arg])
                    except ValueError:
                        raise ValueError("beta must be a number")
                elif arg == "gamma":
                    try:
                        self._gamma = float(args[arg])
                    except ValueError:
                        raise ValueError("gamma must be a number")
                elif arg == "Z":
                    try:
                        self.Z = int(args[arg])
                    except ValueError:
                        raise ValueError("Z must be an integer")
                elif arg == "veca":
                    self._veca = _checkXinR3(args[arg], "veca")
                elif arg == "vecb":
                    self._vecb = _checkXinR3(args[arg], "vecb")
                elif arg == "vecc":
                    self._vecc = _checkXinR3(args[arg], "vecc")
                elif arg == "name":
                    self.name = str(args[arg])
                elif arg == "lattice":
                    self._lattice = str(args[arg]).strip()
                    if self._lattice not in bravais:
                        raise LatticeError(self._lattice + " is not a Bravais lattice :\n" + str(bravais))
                elif arg == "verbose":
                    self.verbose = str(args[arg])

        #
        # check and compute crystal data
        #
        if self.__vectorsExist() and self.__parametersExist():
            raise Warning("When you create a crystal you have to set either \
                lattice parameters or lattice vetors but not both of them")

        # set lattice parameters according to bravais lattice
        if self._lattice is not None:
            self._a, self._b, self._c, self._alpha, self._beta, self._gamma = self._setLatticeParameters()

        if not self.__vectorsExist() and not self.__parametersExist():
            _warningMessage("You did not set neither lattice parameters neither lattice vectors")

        # compute crystal data
        # case 1 : user gives lattice parameters
        if self.__parametersExist():
            # check bravais lattice
            calculatedLattice = self._getBravaisLattice()
            if self._lattice is None:
                self._lattice = calculatedLattice
            elif self._lattice != calculatedLattice:
                print("lattice you set    :" + self._lattice)
                print("lattice I found    :" + calculatedLattice)
                raise LatticeError("Bravais lattice inconsistent")

            # compute lattice vectors and transformation matrix
            self.__computeLatticeVectors()
            self.__fillMatrix()
            self.__computeVolume()

        # case 2 : user gives lattice vectors
        elif self.__vectorsExist():
            # compute lattice parameters and transformation matrix
            self.__computeLatticeParameters()
            self.__fillMatrix()
            self.__computeVolume()

            # check bravais lattice
            calculatedLattice = self._getBravaisLattice()
            if self._lattice is None:
                self._lattice = calculatedLattice
            elif self._lattice != calculatedLattice:
                print("lattice you set    :" + self._lattice)
                print("lattice I found    :" + calculatedLattice)
                raise LatticeError("Bravais lattice inconsistent")

    # -------------------------------------------------------------------------
    #  Construction from a function
    # -------------------------------------------------------------------------
    @classmethod
    def fromPOSCAR(cls, poscar = "POSCAR", verbose = True):
        """ Create a crystal object from a POSCAR/CONTCAR VASP structure file. The POSCAR file 
            can be given as a list of string (lines of the POSCAR file) or you can give a file
            object, or a file name.

            :param poscar: POSCAR file
            :type poscar: list, string, or file object
            :param verbose: verbosity
            :type verbose: bool
            """

        # check argument
        if isinstance(poscar, file):
            poscar = poscar.readlines()
        elif isinstance(poscar, list):
            pass
        elif isinstance(poscar, str):
            if "\n" in poscar:
                poscar = poscar.split("\n")
            else:
                poscar = open(poscar, "r").readlines()
        else:
            raise Warning("argument of fromPOSCAR must be a list, a file object, or a file name")

        # delete blank lines
        for line in poscar:
            if line.strip() == "":
                del line

        # lecture de l'entete du fichier
        title  = poscar[0].strip()
        scale  = float(poscar[1].split()[0])
        veca   = [float(val) * scale for val in poscar[2].split()[0:3] ]
        vecb   = [float(val) * scale for val in poscar[3].split()[0:3] ]
        vecc   = [float(val) * scale for val in poscar[4].split()[0:3] ]

        line = poscar[5].strip()
        vasp4 = True
        for el in line.split():
            vasp4 *= el.isdigit()

        if vasp4:
            # type vasp4.x
            NatomPerType = [int(val) for val in line.split()]
            if verbose:
                print("vasp 4.x POSCAR file")
        else:
            # type vasp5.x
            atomTypeNames = line.split()
            NatomPerType = [int(val) for val in poscar[6].strip().split()]
            if verbose:
                print("vasp 5.x POSCAR file")

        # read coordinate
        direct    = False
        cartesian = False
        start = False
        for line in poscar[6:]:
            if line.strip() == "":
                break

            elif line.strip()[0].lower() == "s":
                # Selective Dynamics
                continue

            elif line.strip()[0].lower() == "d":
                direct = True
                start = True
                x = list()

            elif line.strip()[0].lower() == "c" or line.strip()[0].lower() == "k":
                cartesian = True
                start = True
                x = list()

            elif line.strip()[0].lower() == "s":
                continue

            elif start:
                coord = [float(val) for val in line.split()[0:3]]
                x.append(coord)

            else:
                continue

        if cartesian:
            # apply scale factor
            for xi in x:
                xi[0] *= scale
                xi[1] *= scale
                xi[2] *= scale

        # definition du cristal
        Natoms = len(x)
        n = 0
        for ntyp in NatomPerType: 
            n += ntyp
        if Natoms != n:
            print("Natoms = {0} \t n = {1}".format(Natoms, n))
            raise Warning("Unconsistant atom number")
        crystal = cls(veca = veca, vecb = vecb, vecc = vecc )
        crystal.Natoms = Natoms
        crystal.name = title.strip()

        if verbose:
            print(crystal)

        # nom des atomes
        ntype = 0
        itype = 0
        for iat in range(Natoms):
            if vasp4:
                name = "X{0}".format(itype + 1)
            else:
                name = atomTypeNames[itype]
            crystal.atomNames.append(name)
            ntype += 1
            if ntype == NatomPerType[itype]:
                itype += 1
                ntype = 0

        # passage en coordonnee cartesienne si necessaire
        for coord in x:
            if direct:
                crystal.redCoord.append(coord )
                crystal.XYZCoord.append(crystal.red2cart(coord))
            elif cartesian:
                crystal.redCoord.append(crystal.cart2red(coord))
                crystal.XYZCoord.append(coord)

        return crystal

    @classmethod
    def fromCRYSTAL09(cls, outputfile, verbose = True):
        """ Create a crystal object from an output of CRYSTAL09 program. The function read the
            last crystallographic cell. Outputfile can be a file or a file name.

            :param outputfile: CRYSTAL09 output file (.out)
            :type outputfile: string or file object
            :param verbose: verbosity
            :type verbose: bool
            """

        # check argument
        if isinstance(outputfile, file):
            outfile = outputfile
        if isinstance(outputfile, str):
            outfile = open(outputfile, "r")
        else:
            raise Warning("argument of fromCRYSTAL09 must be a file object, or a file name")

        slab = False
        locGroupe = "SPACE GROUP"
        locMaille = " PRIMITIVE CELL"

        # read general data
        line = outfile.readline()
        end = True
        while line != "":
            line = outfile.readline()
            if "SLAB CALCULATION" in line:
                slab = True
                locGroupe = "PLANE GROUP"

            elif "EEEEEEEEEE STARTING" in line:
                phasename = outfile.readline().strip()
                if verbose:
                    print("name     : {0}".format(phasename))

            elif locGroupe in line:
                group = line.split(":")[1].strip()
                if verbose:
                    print("groupe   : {0}".format(group))

            elif "TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL" in line:
                locMaille = " CRYSTALLOGRAPHIC CELL "

            elif "FINAL OPTIMIZED GEOMETRY" in line:
                end = False
                break

        if end:
            print("Optimisation did not converge, final optimized geometry not found.")
            print("Input geometry will be read instead.")
            outfile.seek(0)
            line = outfile.readline()
            while " GEOMETRY FOR WAVE FUNCTION " not in line:
                line = outfile.readline()

        # read geometry located at locMaille
        line = outfile.readline()
        while locMaille not in line:
            line = outfile.readline()

        outfile.readline()

        a, b, c, alpha, beta, gamma = outfile.readline().split()

        for i in range(4):
            outfile.readline()

        names = list()
        red = list()
        Z = list()
        line = outfile.readline()
        while line != "\n":
            i, p, Zi, namei, xi, yi, zi = line.split()
            names.append(namei)
            Z.append(int(Zi))
            if slab:
                red.append([float(xi), float(yi), float(zi) / float(c) ])
            else:
                red.append([float(xi), float(yi), float(zi)])

            line = outfile.readline()

        # list atom names
        atomnames = list()
        for name in names:
            if name not in atomnames:
                atomnames.append(name)

        # sort atoms by name
        redSorted = list()
        namesSorted = list()
        for name in atomnames:
            for iat, r in enumerate(red):
                if names[iat].strip() == name.strip():
                    redSorted.append(r)
                    namesSorted.append(names[iat])

        if len(redSorted) != len(red) or len(names) != len(namesSorted):
            raise ValueError("Error during atom sorting")

        # create crystal object
        crystal = cls(a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma)
        crystal.name = phasename
        crystal.group = group
        crystal.redCoord = redSorted
        crystal.computeXYZCoord()
        crystal.Natoms = len(redSorted)
        crystal.atomNames = namesSorted

        if verbose:
            print(crystal)

        return crystal

    @classmethod
    def fromCONFIG(cls, config = "CONFIG", verbose = True):
        """ Create a crystal object from a CONFIG/REVCON file of DLPOLY classic. 

            :param config: name of the CONFIG/REVCON file
            :type config: string
            :param verbose: verbosity
            :type verbose: bool
            """


        # check argument
        if isinstance(config, file):
            config = config.readlines()
        elif isinstance(config, list):
            pass
        elif isinstance(config, str):
            if "\n" in config:
                config = config.split("\n")
            else:
                config = open(config, "r").readlines()
        else:
            raise Warning("argument of fromCONFIG must be a list, a file object, or a file name")

        # delete blank lines
        for line in config:
            if line.strip() == "":
                del line

        # lecture de l'entete du fichier
        title  = config[0].strip()
        ndata  = int(config[1].split()[0])
        if ndata not in [0, 1, 2]:
            print("Error : levcfg is not 0, 1 or 2")
            exit(1)

        veca   = [float(val) for val in config[2].split()[0:3]]
        vecb   = [float(val) for val in config[3].split()[0:3]]
        vecc   = [float(val) for val in config[4].split()[0:3]]

        crystal = cls(veca = veca, vecb = vecb, vecc = vecc)
        crystal.name = title.strip()

        # read coordinate
        il = 4
        nat = 0
        while il < len(config) - 1:

            il += 1
            atname = config[il].split()[0]
            atnum  = int(config[il].split()[1])

            nat += 1
            if atnum != nat:
                print("Error : atom number inconsitant")
                print(config[il])
                exit(1)

            il += 1
            coord = [float(val) for val in config[il].split()[0:3]]

            if ndata == 1:
                il += 1
            elif ndata == 2:
                il += 2

            crystal.XYZCoord.append(coord)
            crystal.atomNames.append(atname)

        crystal.Natoms = len(crystal.XYZCoord)

        # calc reduce coordinate
        crystal.computeRedCoord()

        return crystal

    # -------------------------------------------------------------------------
    #  Bravais lattice name and lattice parameters
    # -------------------------------------------------------------------------
    def _setLatticeParameters(self):
        """ set lattice parameters according to bravais lattice name """

        th = 1.e-8

        if self._lattice == "triclinic":
            # nothing to do
            print("Triclinic lattice :")
            if self._a     is None: self._a     = float(raw_input("a = "))
            if self._b     is None: self._b     = float(raw_input("b = "))
            if self._c     is None: self._c     = float(raw_input("c = "))
            if self._alpha is None: self._alpha = float(raw_input("alpha = "))
            if self._beta  is None: self._beta  = float(raw_input("beta = "))
            if self._gamma is None: self._gamma = float(raw_input("gamma = "))

        elif self._lattice == "monoclinic":
            # beta = gamma = 90 degree
            print("Monoclinic lattice :")
            if self._a     is None: self._a     = float(raw_input("a = "))
            if self._b     is None: self._b     = float(raw_input("b = "))
            if self._c     is None: self._c     = float(raw_input("c = "))

            if self._alpha is not None and fabs(self._alpha - 90.0) > th:
                print("alpha = " + str(self._alpha))
                raise LatticeError("Monoclinic lattice alpha must be 90 degrees")
            else:
                self._beta  = 90.0

            if self._beta  is None: self._beta  = float(raw_input("beta = "))

            if self._gamma is not None and fabs(self._gamma - 90.0) > th:
                print("beta = " + str(self._beta))
                raise LatticeError("Monoclinic lattice gamma must be 90 degrees")
            else:
                self._gamma  = 90.0

        elif self._lattice == "orthorhombic":
            # alpha = beta = gamma = 90 degree
            print("Monoclinic lattice :")
            if self._a is None: self._a = float(raw_input("a = "))
            if self._b is None: self._b = float(raw_input("b = "))
            if self._c is None: self._c = float(raw_input("c = "))

            if self._alpha is not None and fabs(self._alpha - 90.0) > th:
                print("alpha = " + str(self._alpha))
                raise LatticeError("Orthorhombic lattice alpha must be 90 degrees")
            else:
                self._beta  = 90.0

            if self._beta is not None and fabs(self._beta - 90.0) > th:
                print("beta = " + str(self._beta))
                raise LatticeError("Orthorhombic lattice beta must be 90 degrees")
            else:
                self._beta  = 90.0

            if self._gamma is not None and fabs(self._gamma - 90.0) > th:
                print("gamma = " + str(self._gamma))
                raise LatticeError("Orthorhombic lattice gamma must be 90 degrees")
            else:
                self._gamma  = 90.0

        elif self._lattice == "tetragonal":
            # a = b != c, alpha = beta = gamma = 90
            if self._a is None: self._a = float(raw_input("a = "))

            if self._b is not None and fabs(self._a - self._b) > th:
                print("a = " + str(self._a))
                print("b = " + str(self._b))
                raise LatticeError("Tetragonal lattice : a = b")
            else:
                self._b = self._a

            if self._c is None: self._c = float(raw_input("c = "))

            if self._alpha is not None and fabs(self._alpha - 90.0) > th:
                print("alpha = " + str(self._alpha))
                raise LatticeError("Tetragonal lattice alpha must be 90 degrees")
            else:
                self._alpha = 90.0

            if self._beta is not None and fabs(self._beta - 90.0) > th:
                print("beta = " + str(self._beta))
                raise LatticeError("Tetragonal lattice beta must be 90 degrees")
            else:
                self._beta = 90.0

            if self._gamma is not None and fabs(self._gamma - 90.0) > th:
                print("gamma = " + str(self._gamma))
                raise LatticeError("Tetragonal lattice gamma must be 90 degrees")
            else:
                self._gamma = 90.0

        elif self._lattice == "rhombohedral":
            # a = b = c, alpha = beta = gamma != 90
            if self._a is None: self._a = float(raw_input("a = "))

            if self._b is not None and fabs(self._a - self._b) > th:
                print("a = " + str(self._a))
                print("b = " + str(self._b))
                raise LatticeError("Rhombohedral lattice : a = b = c")
            else:
                self._b = self._a

            if self._c is not None and fabs(self._a - self._c) > th:
                print("a = " + str(self._a))
                print("c = " + str(self._c))
                raise LatticeError("Rhombohedral lattice : a = b = c")
            else:
                self._c = self._a

            if self._alpha is None: self._alpha = float(raw_input("alpha = "))

            if self._beta is not None and fabs(self._beta - self._alpha) > th:
                print("alpha = " + str(self._alpha))
                print("beta  = " + str(self._beta))
                raise LatticeError("Rhombohedral lattice alpha = beta = gamma")
            else:
                self._beta = 90.0

            if self._gamma is not None and fabs(self._gamma - self._alpha) > th:
                print("alpha = " + str(self._alpha))
                print("gamma = " + str(self._gamma))
                raise LatticeError("Rhombohedral lattice alpha = beta = gamma")
            else:
                self._gamma = 90.0

        elif self._lattice == "hexagonal":
            # a = b != c, alpha = beta = 90 degree, gamma = 120 degree
            if self._a is None: self._a = float(raw_input("a = "))

            if self._b is not None and fabs(self._a - self._b) > th:
                print("a = " + str(self._a))
                print("b = " + str(self._b))
                raise LatticeError("Hexagonal lattice : a = b")
            else:
                self._b = self._a

            if self._c is None: self._c = float(raw_input("c = "))

            if self._alpha is not None and fabs(self._alpha - 90.0) > th:
                print("alpha = " + str(self._alpha))
                raise LatticeError("Hexagonal lattice alpha must be 90 degrees")
            else:
                self._alpha  = 90.0

            if self._beta is not None and fabs(self._beta - 90.0) > th:
                print("beta = " + str(self._beta))
                raise LatticeError("Hexagonal lattice beta must be 90 degrees")
            else:
                self._beta = 90.0

            if self._gamma is not None and fabs(self._gamma - 120.0) > th:
                print("gamma = " + str(self._gamma))
                raise LatticeError("Hexagonal lattice gamma must be 120 degrees")
            else:
                self._gamma  = 120.0
            
        elif self._lattice == "cubic":
            # a = b = c, alpha = beta == gamma = 90 degree
            if self._a is None: self._a = float(raw_input("a = "))

            if self._b is not None and fabs(self._a - self._b) > th:
                print("a = " + str(self._a))
                print("b = " + str(self._b))
                raise LatticeError("Cubic lattice : a = b = c")
            else:
                self._b = self._a

            if self._c is not None and fabs(self._a - self._c) > th:
                print("a = " + str(self._a))
                print("c = " + str(self._c))
                raise LatticeError("Cubic lattice : a = b = c")
            else:
                self._c = self._a

            if self._alpha is not None and fabs(self._alpha - 90.0) > th:
                print("alpha = " + str(self._alpha))
                raise LatticeError("Cubic lattice alpha must be 90 degrees")
            else:
                self._alpha = 90.0

            if self._beta is not None and fabs(self._beta - 90.0) > th:
                print("beta = " + str(self._beta))
                raise LatticeError("Cubic lattice beta must be 90 degrees")
            else:
                self._beta = 90.0

            if self._gamma is not None and fabs(self._gamma - 90.0) > th:
                print("gamma = " + str(self._gamma))
                raise LatticeError("Cubic lattice gamma must be 90 degrees")
            else:
                self._gamma = 90.0

        return self._a, self._b, self._c, self._alpha, self._beta, self._gamma

    def _getBravaisLattice(self):
        """ return the bravais lattice name according to lattice parameters """

        th = 1.e-4

        if not self.__parametersExist():
            raise LatticeError("Parameters unknown")

        if fabs(self._alpha - 90.0) < th and fabs(self._gamma - 90.0) < th:
            if fabs(self._beta - 90.0) < th:
                if fabs(self._a - self._b) < th and fabs(self._a - self._c) < th:
                    lattice = "cubic"
                elif fabs(self._a - self._b) < th:
                    lattice = "tetragonal"
                else:
                    lattice = "orthorhombic"
            else:
                lattice = "monoclinic"

        elif fabs(self._alpha - 90.0) < th and fabs(self._beta - 90.0) < th and \
                                                  fabs(self._gamma - 120.0) < th:
            lattice = "hexagonal"

        elif fabs(self._alpha - self._beta) < th and fabs(self._alpha - self._gamma) < th \
             and fabs(self._a - self._b) < th and fabs(self._a - self._c) < th:
            lattice = "rhombohedral"

        else:
            lattice = "triclinic"

        return lattice

    # -------------------------------------------------------------------------
    #  Attribute properties : a, b, c, alpha, beta, gamma, lattice and lattice vectors
    # -------------------------------------------------------------------------
    def _set_a(self, val):
        """ check instance of val and set self._a """

        _warningMessage("You are changing the lattice")

        try:
            self._a = float(val)
        except ValueError:
            raise ValueError("a must be a number")

        if not self.__parametersExist():
            _warningMessage("All lattice parameters are still not known")
        else:
            self.__computeLatticeVectors()
            self.__fillMatrix()
            self.__computeVolume()
            self._lattice = self._getBravaisLattice()
            self.computeXYZCoord()

    def _get_a(self):
        """ return self._a """
        return self._a
    
    a = property(_get_a, _set_a)
    """ Lattice parameter a, lengths of the vectors :math:`\vec{a}` """


    def _set_b(self, val):
        """ check instance of val and set self._b """

        _warningMessage("You are changing the lattice")

        try:
            self._b = float(val)
        except ValueError:
            raise ValueError("b must be a number")

        if not self.__parametersExist():
            _warningMessage("All lattice parameters are still not known")
        else:
            self.__computeLatticeVectors()
            self.__fillMatrix()
            self.__computeVolume()
            self._lattice = self._getBravaisLattice()
            self.computeXYZCoord()

    def _get_b(self):
        """ return self._b """
        return self._b
    
    b = property(_get_b, _set_b)
    """ Lattice parameter b, lengths of the vectors :math:`\vec{b}` """

    def _set_c(self, val):
        """ check instance of val and set self._c """

        _warningMessage("You are changing the lattice")

        try:
            self._c = float(val)
        except ValueError:
            raise ValueError("c must be a number")

        if not self.__parametersExist():
            _warningMessage("All lattice parameters are still not known")
        else:
            self.__computeLatticeVectors()
            self.__fillMatrix()
            self.__computeVolume()
            self._lattice = self._getBravaisLattice()
            self.computeXYZCoord()

    def _get_c(self):
        """ return self._c """
        return self._c
    
    c = property(_get_c, _set_c)
    """ Lattice parameter c, lengths of the vectors :math:`\vec{c}` """

    def _set_alpha(self, val):
        """ check instance of val and set self._alpha """

        _warningMessage("You are changing the lattice")

        try:
            self._alpha = float(val)
        except ValueError:
            raise ValueError("alpha must be a number")

        if not self.__parametersExist():
            _warningMessage("All lattice parameters are still not known")
        else:
            self.__computeLatticeVectors()
            self.__fillMatrix()
            self.__computeVolume()
            self._lattice = self._getBravaisLattice()
            self.computeXYZCoord()

    def _get_alpha(self):
        """ return self._alpha """
        return self._alpha
    
    alpha = property(_get_alpha, _set_alpha)
    """ Lattice angle :math:`\alpha` . All angles must be given in degrees and 
        are return in degrees. They are internally convert in radian when needed. """

    def _set_beta(self, val):
        """ check instance of val and set self._beta """

        _warningMessage("You are changing the lattice")

        try:
            self._beta = float(val)
        except ValueError:
            raise ValueError("beta must be a number")

        if not self.__parametersExist():
            _warningMessage("All lattice parameters are still not known")
        else:
            self.__computeLatticeVectors()
            self.__fillMatrix()
            self.__computeVolume()
            self._lattice = self._getBravaisLattice()
            self.computeXYZCoord()

    def _get_beta(self):
        """ return self._beta """
        return self._beta
    
    beta = property(_get_beta, _set_beta)
    """ Lattice angle :math:`\beta` . All angles must be given in degrees and 
        are return in degrees. They are internally convert in radian when needed. """

    def _set_gamma(self, val):
        """ check instance of val and set self._gamma """

        _warningMessage("You are changing the lattice")

        try:
            self._gamma = float(val)
        except ValueError:
            raise ValueError("gamma must be a number")

        if not self.__parametersExist():
            _warningMessage("All lattice parameters are still not known")
        else:
            self.__computeLatticeVectors()
            self.__fillMatrix()
            self.__computeVolume()
            self._lattice = self._getBravaisLattice()
            self.computeXYZCoord()

    def _get_gamma(self):
        """ return self._gamma """
        return self._gamma
    
    gamma = property(_get_gamma, _set_gamma)
    """ Lattice angle :math:`\gamma` . All angles must be given in degrees and 
        are return in degrees. They are internally convert in radian when needed. """

    def _set_lattice(self, val):
        """ check instance of val and set self._gamma """

        _warningMessage("You are changing the lattice all lattice parameters may change")

        if val not in bravais:
            raise LatticeError( self._lattice + " is not a Bravais lattice :\n" + str(bravais))
        else:
            self._lattice = val

        # reset all parameters
        self._a = None
        self._b = None
        self._c = None
        self._alpha = None
        self._beta  = None
        self._gamma = None
        self._a, self._b, self._c, self._alpha, self._beta, self._gamma = self._setLatticeParameters()

        # compute new crystal data
        self.__computeLatticeVectors()
        self.__fillMatrix()
        self.__computeVolume()
        self.computeXYZCoord()

    def _get_lattice(self):
        """ return self._lattice """
        return self._lattice
    
    lattice = property( _get_lattice, _set_lattice)
    """ Name of the lattice according to the 7 Bravais lattices. Remind that symmetry is not considered
        for the moment.

        See :py:func:`bravaisLattice` """

    def _set_veca(self, val):
        """ check instance of val and set self._veca """
        _warningMessage("You are changing the lattice")
        self._veca = _checkXinR3(val)
        if not self.__vectorsExist():
            _warningMessage("All lattice vectors are still not known")
        else:
            self.__computeLatticeParameters()
            self.__fillMatrix()
            self.__computeVolume()
            self._lattice = self._getBravaisLattice()
            self.computeXYZCoord()

    def _get_veca(self):
        """ return self._veca """
        return self._veca
    
    veca = property( _get_veca, _set_veca)
    """ Lattice vector :math:`\vec{a}` in a cartesian frame where :math:`\vec{a}` is along x axes and :math:`\vec{b}` is in the (x,y) plane. """

    def _set_vecb(self, val):
        """ check instance of val and set self._vecb """
        _warningMessage("You are changing the lattice")
        self._vecb = _checkXinR3(val)
        if not self.__vectorsExist():
            _warningMessage("All lattice vectors are still not known")
        else:
            self.__computeLatticeParameters()
            self.__fillMatrix()
            self.__computeVolume()
            self._lattice = self._getBravaisLattice()
            self.computeXYZCoord()

    def _get_vecb(self):
        """ return self._vecb """
        return self._vecb
    
    vecb = property( _get_vecb, _set_vecb)
    """ Lattice vector :math:`\vec{b}` in a cartesian frame where :math:`\vec{a}` is along x axes and :math:`\vec{b}` is in the (x,y) plane. """

    def _set_vecc(self, val):
        """ check instance of val and set self._vecc """
        _warningMessage("You are changing the lattice")
        self._vecc = _checkXinR3(val)
        if not self.__vectorsExist():
            _warningMessage("All lattice vectors are still not known")
        else:
            self.__computeLatticeParameters()
            self.__fillMatrix()
            self.__computeVolume()
            self._lattice = self._getBravaisLattice()
            self.computeXYZCoord()

    def _get_vecc(self):
        """ return self._vecc """
        return self._vecc
    
    vecc = property(_get_vecc, _set_vecc)
    """ Lattice vector :math:`\vec{c}` in a cartesian frame where :math:`\vec{a}` is along x axes and :math:`\vec{b}` is in the (x,y) plane. """

    def _set_volume(self, val):
        """ changing volume by hand is forbiden """
        _warningMessage("You can change volume by hand")

    def _get_volume(self):
        """ return self._vecc """
        return self._volume

    volume = property(_get_volume, _set_volume)
    """ volume of the unit cell. """

    # -------------------------------------------------------------------------
    #  Parameters or lattice vectors calculations
    # -------------------------------------------------------------------------
    def __computeLatticeVectors(self):
        """ Compute lattice vectors coordinates in a cartesian frame.
            Vector a is along the x axes, vector b is in the (x, y) plane. """

        if not self.__parametersExist():
            raise LatticeError("Parameters unknown")

        # conversion en radian
        alphar = self._alpha * pi / 180.0
        betar  = self._beta  * pi / 180.0
        gammar = self._gamma * pi / 180.0

        # calcul des vecteurs
        self._veca = [0.,0.,0.]
        self._veca[0] = self._a

        self._vecb = [0.,0.,0.]
        self._vecb[0] = self._b * cos( gammar )
        self._vecb[1] = self._b * sin( gammar )

        self._vecc = [0.,0.,0.]
        self._vecc[0] = self._c * cos( betar )
        cy = ( cos(alphar) - cos(gammar) * cos(betar) ) / sin(gammar)
        self._vecc[1] = self._c * cy
        cz = sqrt( ( sin( betar ) )**2 - cy**2 )
        self._vecc[2] = self._c * cz

    def __computeLatticeParameters(self, verbose = False):
        """ compute lattice parameters from lattice vectors """

        if not self.__vectorsExist() :
            raise LatticeError("Vectors unknown")

        # make vectors
        veca = np.array( self._veca )
        vecb = np.array( self._vecb )
        vecc = np.array( self._vecc )

        # norm
        self._a = sqrt( (veca**2).sum() )
        self._b = sqrt( (vecb**2).sum() )
        self._c = sqrt( (vecc**2).sum() )

        # alpha in degree
        u = vecb / self._b
        v = vecc / self._c
        scalaire = np.dot( u, v) 
        if fabs(scalaire < 1.):
            self._alpha = arccos( scalaire ) * 180.0 / pi
        else:
            print("Erreur : produit scalaire alpha" + str(scalaire))

        # beta in degree
        u = veca / self._a
        v = vecc / self._c
        scalaire = np.dot( u, v) 
        if fabs(scalaire < 1.):
            self._beta = arccos( scalaire ) * 180.0 / pi
        else:
            print("Erreur : produit scalaire beta" + str(scalaire))

        # gamma in degree
        u = veca / self._a
        v = vecb / self._b
        scalaire = np.dot( u, v) 
        if fabs(scalaire < 1.):
            self._gamma = arccos( scalaire ) * 180.0 / pi
        else:
            print("Erreur : produit scalaire gamma" + str(scalaire))

        # print
        if verbose :
            print("a     = %10.5f" % self._a)
            print("b     = %10.5f" % self._b)
            print("c     = %10.5f" % self._c)
            print("alpha = %10.3f" % (self._alpha))
            print("beta  = %10.3f" % (self._beta ))
            print("gamma = %10.3f" % (self._gamma) + "\n")

    def __fillMatrix(self) :
        """ fill transformation matrix :

                * Mred2cart : reduce    -> cartesian
                * Mcart2red : cartesian -> reduce

            Theses matrix are numpy.ndarray object.
        """
        # check if lattice vectors are known
        if not self.__vectorsExist():
            raise LatticeError("Lattice vectors unknown")

        # check if lattice parameters are known
        if not self.__parametersExist():
            raise LatticeError("Lattice parameters unknown")

        # transformation matrix from reduce to cartesian coordinate
        self.Mred2cart = np.array([self._veca, self._vecb, self._vecc]).transpose()

        # transformation matrix from cartesian to reduce coordinate
        self.Mcart2red = np.linalg.inv(self.Mred2cart)

        # matrice de passage coordonnées cartésiennes -> coordonnées réduites
        #alphar = self._alpha * pi / 180.0
        #betar  = self._beta  * pi / 180.0
        #gammar = self._gamma * pi / 180.0 
        #cotgamma = cos(gammar) / sin(gammar)
        #cy = ( cos(alphar) - cos(gammar) * cos(betar) ) / sin(gammar)
        #cz = sqrt( ( sin( betar ) )**2 - cy**2 )

        #self.Mcart2red = np.zeros( (3,3) )
        #self.Mcart2red[0,0] = 1. / self._a
        #self.Mcart2red[0,1] = - cotgamma / self._a
        #self.Mcart2red[0,2] = 1. / self._a * ( cy/cz * cotgamma - cos( betar ) / cz )

        #self.Mcart2red[1,1] = 1. / ( self._b * sin( gammar ) )
        #self.Mcart2red[1,2] = - cy / ( self._b * cz * sin( gammar ) )

        #self.Mcart2red[2,2] = 1. / ( self._c * cz )

        # varialble de contrôle
        self.__matrixSet = True

    def __computeVolume(self):
        """ compute the volume of the cell """
        if not self.__vectorsExist():
            self.__computeLatticeVectors()
        self._volume = np.dot(self._veca, np.cross(self._vecb, self._vecc))
        return self._volume

    # -------------------------------------------------------------------------
    #  Check if parameters or vectors are setted
    # -------------------------------------------------------------------------
    def __parametersExist(self):
        """ Did lattice parameters (a, b, c, alpha, beta, gamma) set ? """

        if self._a == None or self._b == None or self._c == None or \
           self._alpha == None or self._beta == None or self._gamma == None:
            return False
        else:
            return True

    def __vectorsExist(self):
        """ Did lattice vectors (veca, vecb et vecc) set ? """
        if self._veca == None or self._vecb == None or self._vecc == None:
            return False
        else:
            return True

    # -------------------------------------------------------------------------
    #  conversion between cartesian and reduce unit coordinates 
    # -------------------------------------------------------------------------
    def red2cart(self, x):
        """ Convert reduce coordinate into cartesian coordinate.
            
            :param x: cartesian coordinate vector
            :type x: list or array or ndarray, len(x) = 3
            :rtype: list of lenght 3
        """
        vecx = _checkXinR3(x)
        if not self.__matrixSet:
            self.__fillMatrix()
        return np.dot(self.Mred2cart, vecx).tolist()

    def cart2red(self, x):
        """ Convert cartesian coordinate into reduce coordinate.
            
            :param x: reduce coordinate vector
            :type x: list or array or ndarray, len(x) = 3
            :rtype: list of lenght 3
        """
        vecx = _checkXinR3(x)
        if not self.__matrixSet:
            self.__fillMatrix()
        return np.dot(self.Mcart2red, vecx).tolist()

    def computeRedCoord(self):
        """ Compute cartesian coordinates for all coordinates stored in self.redCoord. """
        if self.XYZCoord == []:
            _warningMessage("no cartesian coordinates availlable")
        self.redCoord = [self.cart2red(xyz) for xyz in self.XYZCoord]
        return self.redCoord

    def computeXYZCoord(self):
        """ Compute reduce coordinates for all coordinates stored in self.XYZCoord. """ 
        if self.redCoord == []:
            _warningMessage("no reduce coordinates availlable")
        self.XYZCoord = [self.red2cart(red) for red in self.redCoord]
        return self.XYZCoord

    # -------------------------------------------------------------------------
    #  Distance calculations
    # -------------------------------------------------------------------------
    def dist_r(self, r1, r2):
        """ Compute the minimum distance between reduce coordinates r1 and r2 considering
            the periodic conditions. 

           :param r1: reduce coordinate of atom 1
           :type r1: list or array or ndarray, len(A1) = 3
           :param r2: reduce coordinate of atom 2
           :type r2: list or array or ndarray, len(A2) = 3
           :rtype: float
        """

        dr = self.img(r1, r2)
        dx = self.red2cart(dr)
        return np.sqrt(dx[0]**2 + dx[1]**2 + dx[2]**2)

    def dist_x(self, x1, x2):
        """ Compute the minimum distance between cartesian coordinates x1 and x2 considering
            the periodic conditions. 

           :param x1: cartesian coordinates of atom 1
           :type x1: list or array or ndarray, len(A1) = 3
           :param x2: cartesian coordinates of atom 2
           :type x2: list or array or ndarray, len(A2) = 3
           :rtype: float
        """

        dr = self.img(self.cart2red(x1), self.cart2red(x2))
        dx = self.red2cart(dr)
        return np.sqrt(dx[0]**2 + dx[1]**2 + dx[2]**2)

    def dist_i(self, iat, jat):
        """ Compute the minimum distance between atom with index iat and jat considering
            the periodic conditions. 

           :param iat: atom index of atom i
           :type A1: integer
           :param A2: atom index of atom j
           :type A2: atom index
           :rtype: float
        """

        dr = self.img(self.redCoord[iat], self.redCoord[jat])
        dx = self.red2cart(dr)
        return np.sqrt(dx[0]**2 + dx[1]**2 + dx[2]**2)

    def img(self, r1, r2):
        """ Return the shortest vector in reduce coordinates considering the periodic
            conditions.

           :param r1: reduce coordinate of atom r1
           :type r1: list or array or ndarray, len(A1) = 3
           :param r2: reduce coordinate of atom r2
           :type r2: list or array or ndarray, len(A2) = 3
           :rtype: list

           r1 and r2 must be reduce coordinate and an object of length 3 !
           In order to improve efficiency no check are done.
        """
        
        dr = [0., 0., 0.]
        for i in range(3):
            dr[i] = r2[i] - r1[i]
            if dr[i] > .5:
                dr[i] -= 1.
            elif dr[i] < -.5:
                dr[i] += 1.

        return dr

    # -------------------------------------------------------------------------
    #  Neighbors list
    # -------------------------------------------------------------------------
    def neighbor(self, rcut, iatom, ligand = None):
        """" look for neighbors of atoms in the list iatom. If ligand present, only
             atoms with name ligand are considered. If rcut < cell parameters / 2,
             the subroutine looks for multiple neighbors by considering all atoms
             in the 26 adjacent cells around the central the considered cell.

             :param rcut: cut off radius
             :type rcut: float
             :param iatom: list of atom for which neighbors have to be identified
             :type iatom: list
             :param ligand: limit the neighbors list to atoms name ligand (default: all)
             :type ligand: string
             :return: neighbor list for each atom in iatom list
             :rtype: list[list]
           
        """

        neigh = list()

        if rcut > self.a / 2. or rcut > self.b / 2. or rcut > self.c / 2:
            # small cell look for multiple neighbor due to translation
            print("cell parameters / 2 < rcut => look for multiple neighbors")
            for iv, iat in enumerate(iatom):
                riat = np.array(self.redCoord[iat])

                neigh.append(list())
                for jat in range(self.Natoms):
                    if iat == jat: 
                        continue
                    if ligand and self.atomNames[jat] != ligand: 
                        continue

                    rjat = np.array([val for val in self.redCoord[jat]])
                    for t in trans:
                        dr = rjat + t - riat
                        dx = self.red2cart(dr)
                        d  = np.sqrt((dx**2).sum())
                        if d <= rcut:
                            # new neighbor
                            neigh[iv].append(jat)
        else:
            # simpler case where atoms can be neighbor only one time
            for iv, iat in enumerate(iatom):
                riat = np.array(self.redCoord[iat])

                neigh.append(list())
                for jat in range(self.Natoms):
                    if iat == jat: 
                        continue
                    if ligand and self.atomNames[jat] != ligand: 
                        continue

                    if self.dist_i(iat, jat) <= rcut:
                        # new neighbor
                        neigh[iv].append(jat)

        return neigh

    def allneighbor(self, rcut):
        """" look for neighbors of all atoms in the cell. If rcut < cell parameters / 2,
             the subroutine looks for multiple neighbors by considering all atoms
             in the 26 adjacent cells around the central the considered cell.

             :param rcut: cut off radius
             :type rcut: float
             :param iatom: list of atom for which neighbors have to be identified
             :type iatom: list
             :return: neighbor list for each atom in iatom list
             :rtype: list[list]
           
        """

        # initialisation
        neigh = [list() for iat in range(self.Natoms)]

        if rcut > self.a / 2. or rcut > self.b / 2. or rcut > self.c / 2:
            # small cell look for multiple neighbor due to translation
            print("cell parameters / 2 < rcut => look for multiple neighbors")

            for iat in range(self.Natoms - 1):
                for jat in range(iat + 1, self.Natoms):
                    for t in trans:
                        dx = self.red2cart(self.redCoord[jat] + t - self.redCoord[iat])
                        d = np.sqrt((dx**2).sum())
                        if d <= rcut:
                            neigh[iat].append(jat)
                            neigh[jat].append(iat)
           
        else:
            # simpler case where atoms can be neighbor only one time
            for iat in range(self.Natoms - 1):
                for jat in range(iat + 1, self.Natoms):
                    if self.dist_i(iat, jat) <= rcut:
                        neigh[iat].append(jat)
                        neigh[jat].append(iat)

        return neigh

    # -------------------------------------------------------------------------
    #  operation on atoms
    # -------------------------------------------------------------------------
    def calcDeplacementAtomes(self, other):
        """ Compute the distance between the cartesian postion fo the same atom in self
            and other structures. 
            
            :param other: crystal object
            :rtype: list of distance

        """

        # check atom number
        if self.Natoms != other.Natoms:
            raise ValueError("Error : atom numbers does not match")

        # make copy of each crystal
        c1 = copy.deepcopy(self)
        c2 = copy.deepcopy(other)

        # put each atom of struct1 inside unitcell
        c1.wrapAtomes()

        # place atoms at minimum distance
        c1.alignCrystal(c2)

        deplacement = list()
        for xyz1, xyz2 in zip( c1.XYZCoord, c2.XYZCoord):
            d = 0.
            for x1, x2 in zip(xyz1, xyz2):
                d += (x2 - x1)**2
            deplacement.append(sqrt(d))
     
        return deplacement

    def wrapAtoms(self):
        """ Replace all atoms into the unit cell using reduce coordinate """

        if self.redCoord == []:
            self.calcRedCoord()
        
        for iat in range(self.Natoms):
            for i in range(3):
                if self.redCoord[iat][i] > 1.0:
                    self.redCoord[iat][i] -= 1.0
                elif self.redCoord[iat][i] < 0.0:
                    self.redCoord[iat][i] += 1.0

        self.calcXYZCoord()

        return self.XYZCoord

    def alignCrystal(self, other):
        """ For each atom of the crystal, if an image of the atom belonging to other is
        closer to the same atom belonging to self, the other's atom is translated to that
        position. 

        This method could be usefull if one wants to compute displacements of atoms after
        a relaxation of the structre.

            :param other: crystal object
            :rtype: list of distance

        Warning :if lattices of self and other are not the same, this have not any sense.
        """

        if self.redCoord == []:
            self.computeRedCoord()
            if self.redCoord == []:
                raise NameError("crystal 1 : reduce unit have to be known")    

        if other.redCoord == []:
            other.computeRedCoord()
            if self.redCoord == []:
                raise NameError("crystal 2 : reduce unit have to be known")    

        for iat in range(self.Natoms):
            for i in range(3):
                d = fabs(other.redCoord[iat][i] - self.redCoord[iat][i])
                if d > 0.5:
                    if other.redCoord[iat][i] > self.redCoord[iat][i]:
                        other.redCoord[iat][i] -= 1.0
                    else:
                        other.redCoord[iat][i] += 1.0
        other.calcXYZCoord()

        return other.XYZCoord

    def makeSupercell(self, nx, ny, nz, sort = True):
        """ Build a supercell from the cell defined by the object. Return a new Crystal object as
            output.

            :param nx: number of cell in x direction
            :type nx: int
            :param ny: number of cell in y direction
            :type ny: int
            :param nz: number of cell in z direction
            :type nz: int
            :param sort: if True atoms are sorted by name and by z (needed for POSCAR output !).
            :type sort: bool
            :rtype: Crystal

        """

        if nx == 1 and ny == 1 and nz == 1:
            print("{0}x{1}x{2} : nothing to do".format(nx, ny, nz))
            new = copy.deepcopy(self)
            return new

        # check
        if self.redCoord == []:
            self.computeRedCoord()
            if self.redCoord == []:
                raise ValueError("coordinate unknown")

        # atom number
        if self.Natoms == None:
            self.Natoms = len(self.redCoord)

        # atom names
        atomList = list()
        if self.atomNames == []:
            self.atomNames = ["X" for i in range(self.Natoms)]
            atomList = ["X"]
        else:
            for name in self.atomNames:
                if name not in atomList:
                    atomList.append(name)

        # make supercell
        red = list()
        atom = copy.deepcopy(self.atomNames)

        for ix in range(0, nx):
            tx = float(ix)
            for iy in range(0, ny):
                ty = float(iy)
                for iz in range(0, nz):
                    tz = float(iz)
                    for iat in range(self.Natoms):
                        r = self.redCoord[iat]
                        rx = (r[0] + tx) / float(nx)
                        ry = (r[1] + ty) / float(ny)
                        rz = (r[2] + tz) / float(nz)
                        red.append([rx, ry, rz])
                        atom.append(atom[iat])

        # new atom number
        Natoms = len(red)
        if Natoms != nx * ny * nz * self.Natoms:
            print("Natoms                     = %d" % Natoms)
            print("nx * ny * nz * self.Natoms = %d " %  (nx * ny * nz * self.Natoms))
            raise ValueError("Atom number error")

        # create supercell object
        veca = [nx * v for v in self.veca]
        vecb = [ny * v for v in self.vecb]
        vecc = [nz * v for v in self.vecc]

        new = Crystal(veca = veca, vecb = vecb, vecc = vecc)
        new.name = self.name
        new.Natoms = Natoms

        if sort:
            # sort atoms per name
            redNameSorted = list()
            atomNames = list()
            for name in atomList:
                for iat in range(Natoms):
                    if atom[iat].strip() == name.strip():
                        redNameSorted.append(red[iat])
                        atomNames.append(name.strip())

            # number of each atom
            nAtomPerType = dict()
            for name in atomList:
                nAtomPerType[name] = atomNames.count(name)

            # sort atoms by z
            redSorted = list()
            for ityp in range(len(atomList)):
                if ityp == 0:
                    n1 = 0
                    n2 = nAtomPerType[atomList[ityp]]
                else:
                    n1 += nAtomPerType[atomList[ityp - 1]]
                    n2 = n1 + nAtomPerType[atomList[ityp]]
                redSorted += sorted(redNameSorted[n1:n2], key = lambda f:f[2], reverse = True)

            if len(redSorted) != len(red):
                print(len(redSorted))
                print(len(red))
                print("Error in makeSupercell, check atom number")
                exit(1)

            new.redCoord = redSorted
        else:
            new.redCoord = red
            atomNames = nx * ny * nz * self.atomNames

        new.atomNames = atomNames
        new.computeXYZCoord()

        return new

    # -------------------------------------------------------------------------
    #  output methods
    # -------------------------------------------------------------------------
    def toXYZ(self, filename):
        """ (Under development) Write cartesian coordinate in a XYZ format 

           :param filename: name of the output file.
           :type filename: string
        
        """
        print("TODO")
        pass

    def toPOSCAR(self, filename = "POSCAR", sort = True):
        """ Write crystal object in a vasp POSCAR file

           :param filename: name of the output file.
           :type filename: string

        """

        # check
        if self.redCoord == []:
            self.computeRedCoord()
            if self.redCoord == []:
                raise ValueError("coordinate unknown")

        if self.Natoms == None:
            _warningMessage("Natoms was not set. I use the number of reduce coordinate")
            self.Natoms = len(self.redCoord)

        # title line
        if self.name is not None:
            poscar = self.name + "\n"
        else:
            poscar = "POSCAR file\n"

        # scaling factor
        poscar += " 1.0\n"

        # lattice vectors
        if not self.__vectorsExist():
            self.computeLatticeVectors()

        for x in self.veca:
            poscar += "%20.12f" % x
        poscar += "\n" 
        for x in self.vecb:
            poscar += "%20.12f" % x
        poscar += "\n" 
        for x in self.vecc:
            poscar += "%20.12f" % x
        poscar += "\n" 

        # atom name
        atomList = list()
        if self.atomNames != []:
            # look for atom type
            for name in self.atomNames:
                if name not in atomList:
                    atomList.append(name)

            # print atom type
            for name in atomList:
                poscar += "%4s" % name
            poscar += "\n"

            # sort atom according to atom type
            if sort:
                redSorted = list()
                for name in atomList:
                    for iat in range(self.Natoms):
                        if self.atomNames[iat].strip() == name.strip():
                            redSorted.append(self.redCoord[iat])

                if len(redSorted) != len(self.redCoord):
                    print("Error check atom number")
                    exit(1)
            else:
                redSorted = self.redCoord

            # print the atom number per type
            for name in atomList:
                poscar += "%4d" % self.atomNames.count(name)
            poscar += "\n"

        else:
            _warningMessage("Atom names unknown, I print only atom number")
            _warningMessage("toPOSCAR assume that atom coordinate are correctly sorted")
            poscar += "%5d\n" % len(self.redCoord)

        # reduce coordinate
        poscar += "Direct\n"
        for red in redSorted:
            for r in red:
                poscar += "%20.12f" % r
            poscar += "\n"

        open(filename, "w").write(poscar)

    def toCONFIG(self, filename = "CONFIG"):
        """ Write crystal object in a CONFIG DLPOLY file

           :param filename: name of the output file.
           :type filename: string
        """

        #_warningMessage("toCRYSTAL assume that atom coordinate are correctly sorted")

        # title line
        if self.name is not None:
            config = self.name + "\n"
        else:
            config = "CONFIG file\n"

        # periodic boundary
        if self.lattice == "cubic":
            config += "%10d%10d\n" % (0,1)
        elif self.lattice == "tetragonal" or self.lattice == "orthorhombic":
            config += "%10d%10d\n" % (0,2)
        else:
            config += "%10d%10d\n" % (0,3)

        # lattice vectors
        if not self.__vectorsExist():
            self.computeLatticeVectors()

        for x in self.veca:
            config += "%20.12f" % x
        config += "\n" 
        for x in self.vecb:
            config += "%20.12f" % x
        config += "\n" 
        for x in self.vecc:
            config += "%20.12f" % x
        config += "\n" 

        if self.XYZCoord == []:
            self.computeXYZCoord()

        atomNameKnown = True
        if self.atomNames == [] or len(self.atomNames) != len(self.XYZCoord):
            print(self.atomNames)
            print(len(self.atomNames), len(self.XYZCoord))
            _warningMessage("Atom names unknown, I will use 'X' instead")
            atomNameKnown = False

        for iat, xyz in enumerate(self.XYZCoord):
            if atomNameKnown:
                config += "%s%10d\n" % (self.atomNames[iat].ljust(8), iat + 1)
            else:
                config += "%8s%10d\n" % ("X       ", iat + 1)
            config += "%20.10f%20.10f%20.10f\n" % tuple(xyz)

        open(filename, "w").write(config)

    def printLatticeVectors(self):
        """ Print the lattice vectors in a cartesian frame. """
        dashline = "".join(50 * ["-"])
        print(dashline)
        print(" Lattice vectors in cartesian frame : lattice is {0}".format(self.lattice))
        print(dashline)
        print(" a = %15.10f %15.10f %15.10f" % (self.veca[0], self.veca[1], self.veca[2]))
        print(" b = %15.10f %15.10f %15.10f" % (self.vecb[0], self.vecb[1], self.vecb[2]))
        print(" c = %15.10f %15.10f %15.10f" % (self.vecc[0], self.vecc[1], self.vecc[2]))
        print(dashline + "\n")
        

    def printLatticeParameters(self):
        """ Print lattice paramters. """
        print(str(self))

    # -------------------------------------------------------------------------
    #  Special methods
    # -------------------------------------------------------------------------
    def __str__(self):
        """ string conversion """

        if not self.__parametersExist():
            raise LatticeError("Lattice parameters are not known did you set all \
                lattice vectors or all lattice parameters")

        dashline = "".join(50 * ["-"]) + "\n"

        line = dashline
        if self.name != None:
            line += " lattice parameters : " + self.name + "\n"
        else :
            line += " lattice parameters : \n"

        line += dashline

        line += " a     = %10.5f" % self._a + "\n" 
        line += " b     = %10.5f" % self._b + "\n" 
        line += " c     = %10.5f" % self._c + "\n" 
        line += " alpha = %10.3f" % (self._alpha ) + "\n" 
        line += " beta  = %10.3f" % (self._beta  ) + "\n" 
        line += " gamma = %10.3f" % (self._gamma ) + "\n"

        line += dashline

        return line

    def __eq__(self, other):
        """ equal """
        selfatt  = [self.a, self.b, self.c, self.alpha, self.beta, self.gamma]
        for xyz in self.XYZCoord:
            selfatt += xyz
        for red in self.redCoord:
            selfatt += red
        selfatt += self.veca + self.vecb + self.vecc

        otheratt  = [other.a, other.b, other.c, other.alpha, other.beta, other.gamma]
        for xyz in other.XYZCoord:
            otheratt += xyz
        for red in other.redCoord:
            otheratt += red
        otheratt += other.veca + other.vecb + other.vecc

        th = 1.e-6
        res = True
        if self.lattice != other.lattice:
            res = False
        else:
            for s, o in zip(selfatt, otheratt):
                if fabs(s - o) > th:
                    res = False

        return res

# -----------------------------------------------------------------------------
#  inner function or exception
# -----------------------------------------------------------------------------
class LatticeError(Exception):
    """ Exception concerning lattice parameters, bravais name or lattice vectors. """
    def __init__(self, why):
        self.why = why
    def __str__(self):
        return self.why

def _checkXinR3( X, name = ""):
    """ check if X in R^3 and return a list object

    >>> _checkXinR3([0,1,2], "X")
    [0.0, 1.0, 2.0]
    >>> _checkXinR3(3,"a")
    Traceback (most recent call last):
    ...
    TypeError: Variable a must be a vector in R^3
    >>> _checkXinR3(["0","1","2"], "X")
    [0.0, 1.0, 2.0]
    >>> _checkXinR3(["a","b","c"], "X")
    Traceback (most recent call last):
    ...
    ValueError: Variable X must be a vector in R^3
    
    """
    if (isinstance( X, list) and len(X) == 3)       or \
       (isinstance( X, np.ndarray) and X.size == 3) or \
       (isinstance( X, array.ArrayType) and len(X) == 3):
        try:
            listX = [float(xi) for xi in X]
        except ValueError:
            raise ValueError("Variable " + name + " must be a vector in R^3")
    else:
        raise TypeError("Variable " + name + " must be a vector in R^3")

    return listX

def _warningMessage(why):
    """ pring a warning 

    >>> _warningMessage("small error")
     * * * WARNING : small error

    """
    print(" * * * WARNING : " + why) 

# -----------------------------------------------------------------------------
#  output functions or instanciation from a function
# -----------------------------------------------------------------------------
def printBravaisLattice():
    """ Print bravais lattice names and caracteristics and return the output in a string. """

    data = """
triclinic :
    * a != b != c
    * alpha != beta != gamma != 90

monoclinic :
    * a != b != c
    * alpha = gamma = 90 ; beta != 90

orthorhombic :
    * a != b != c
    * alpha = beta = gamma = 90

tetragonal (quadratique) :
    * a = b != c
    * alpha = beta = gamma = 90

hexagonal :
    * a = b != c
    * alpha = beta = 90 ; gamma = 120

rhombohedral :
    * a = b = c
    * alpha = beta = gamma != 90

cubic :
    * a = b = c
    * alpha = beta = gamma = 90"""

    print(data)

# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    doctest.testmod()
