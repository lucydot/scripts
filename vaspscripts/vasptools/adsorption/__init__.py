#!/usr/bin/env python
# -*- coding=utf-8 -*-

import os
import numpy as np
import crystal
import copy

class Molecule(object):
    """ This class define a molecule and provide several method in order to translate or rotate
    the atoms of the molecule around local axes.
    
    :param name: Name of the molecule
    :type name: string
    
    """

    def __init__(self, name = "A molecule"):
        """ constructor """

        self.name = name
        """molecule name"""

        self.atoms = list()
        """ List of atom object (see Atom class) """

        self._Natoms = len(self.atoms)
        """ Number of atom in the molecule """

        self._formula = dict()
        """ chemical formula """

        self._oldgeom = list()

        # bool : are local axes defined ?
        self.localAxesDefined = False
        """ are local axes defined (bool) """

        # bool : is the molecule in the local framework ?
        self.inLocalFramework = False
        """ is the molecule in the local framework (bool) """
 
        # local axes = cartesian axes linked to the molecule
        self.origin = np.zeros(3)
        """ origin of the local framework """
        self.locx = np.zeros(3)
        """ x axes of the local framework """
        self.locy = np.zeros(3)
        """ y axes of the local framework """
        self.locz = np.zeros(3)
        """ z axes of the local framework """

    @classmethod
    def fromXYZ(cls, inputFile):
        """ Create a molecule object from a file in a XYZ format.
        
            :param inputFile: A file name of the molecule in a XYZ format
            :type inputFile: string
        """

        if not os.path.exists(inputFile):
            print("{0} does not exist".format(inputFile))
            exit(1)

        data = open(inputFile, "r").readlines()

        mol = cls(data[1].strip())
        nat = int(data[0].strip())

        for line in data[2:2 + nat]:
            atName = line.split()[0].strip()
            coord = [float(val.strip()) for val in line.split()[1:4]]
            mol.atoms.append(Atom(name = atName, xyz = coord))

        return mol

    # Natoms
    # ------
    def _set_Natoms(self, val):
        """ setters Natoms """
        print("WARNING, you cannot change the number of atoms")

    def _get_Natoms(self):
        """ get Natoms """
        self._Natoms = len(self.atoms)
        return self._Natoms

    Natoms = property(_get_Natoms, _set_Natoms)
    """ Number of atom in the molecule """

    # chemical formula 
    # ----------------
    def _get_formula(self):
        """ return the composition of the molecule """
        compo = dict()
        for atom in self.atoms:
            if atom.name not in compo.keys():
                compo[atom.name] = 1
            else:
                compo[atom.name] += 1

        self._formula = compo
        return self._formula

    def _set_formula(self, val):
        """ setters formula """
        print("WARNING, you cannot change the molecule formula")

    formula = property(_get_formula, _set_formula)
    """ chemical formula of the molecule """

    def printChemicalFormula(self):
        """ print the chemical formula """
        compo = self._get_formula
        line = ""
        for name in compo.keys():
            line += "{0}_{1} ".format(name, compo[name])
        print(line)

    # backup of initial geometry
    # --------------------------
    def _set_oldgeom(self, val):
        """ setters oldgeom """
        print("WARNING, oldgeom is a backup you cannot modify it")

    def _get_oldgeom(self):
        """ get oldgem """
        return copy.deepcopy(self._oldgeom)

    oldgeom = property(_get_oldgeom, _set_oldgeom)
    """ Backup of a structure """

    def backup(self):
        """ Save all coordinates in oldgem """
        self._oldgeom = copy.deepcopy(self.atoms)

    # local axes
    # ----------
    def setLocalAxes(self, origin = None, iatx = None, iaty = None, iatz = None,\
                     linear = False, verbose = False):
        """ Define local axes linked to the molecule and put the molecule in this frame.
            All atom numbers start from 1 to mol.Natoms

            :param origin: index of the atom at the origin
            :type origin: int
            :param iatx: axes x is origin -> iatx
            :type iatx: int
            :param iaty: axes y is origin -> iaty
            :type iaty: int
            :param iatz: axes z is origin -> iatz
            :type iatz: int
            :param linear: True if molecule is linear
            :type linear: bool
            :param verbose: set verbosity
            :type verbose: bool
        """

        # check atom number
        for n in [origin, iatx, iaty, iatz]:
            if n is not None:
                if not 1 <= n <= self.Natoms:
                    print("Error : atom number not in the range 1 -> Natoms")
                    print("origin : " + str(origin))
                    print("iatx   : " + str(iatx))
                    print("iaty   : " + str(iaty))
                    print("iatz   : " + str(iatz))
                    exit(1)

        # set origin
        if origin is None: 
            origin = 0
        else:
            origin -= 1
        self.origin = self.atoms[origin].xyz

        # atom number in the range 0 -> Natoms - 1
        if iatx is not None: iatx -= 1
        if iaty is not None: iaty -= 1
        if iatz is not None: iatz -= 1

        # move all atoms at the origin
        #for atom in self.atoms:
            #atom.xyz -= O

        # define orthonormal local axes
        if linear:
            if iatz is None:
                if origin + 1 < self.Natoms:
                    iatz = origin + 1
                elif 0 <= origin - 1:
                    iatz = origin - 1
                else:
                    print("Error : are you sure that this molecule contains more than 1 atom")
                    exit(1)
            self.locz = self.atoms[iatz].xyz - self.origin
            self.locz /= np.sqrt((self.locz**2).sum())

            if np.fabs(np.dot(self.locz, np.array([1.,0.,0.]))) < 1.e-6:
                self.locx = self.locz + np.array([0.,1.,0.])
                self.locx -= np.dot(self.locz, self.locx) * self.locz
                self.locx /= np.sqrt((self.locx**2).sum())
            else:
                self.locx = self.locz + np.array([1.,0.,0.])
                self.locx -= np.dot(self.locz, self.locx) * self.locz
                self.locx /= np.sqrt((self.locx**2).sum())

            self.locy = np.cross(self.locz, self.locx)

        else:
            if iatx is None or not isinstance(iatx, int):
                print("Error : iatx must be an atom number.")
                exit(1)
            self.locx = self.atoms[iatx].xyz - self.origin
            self.locx /= np.sqrt((self.locx**2).sum())

            if iaty is None or not isinstance(iaty, int):
                print("Error : iaty must be an atom number.")
                exit(1)
            self.locy = self.atoms[iaty].xyz - self.origin
            self.locy -= np.dot(self.locx, self.locy) * self.locx
            self.locy /= np.sqrt((self.locy**2).sum())

            if iatz is None:
                self.locz = np.cross(self.locx, self.locy)
            else:
                if not isinstance(iatz, int):
                    print("Error : iatz must be an atom number.")
                    exit(1)
                self.locz = self.atoms[iatz].xyz - self.origin
                self.locz /= np.sqrt((self.locz**2).sum())

                if np.fabs(np.dot(self.locx, self.locz)) > 1e-15:
                    print("WARNING, you choose z axes such as x and z are not perpendicular")
                if np.fabs(np.dot(self.locy, self.locz)) > 1e-15:
                    print("WARNING, you choose z axes such as y and z are not perpendicular")

        if verbose:
            print("---------------------------------------------------")
            print(" local framework")
            print("---------------------------------------------------")
            print("4")
            print("local framework linked to the molecule")
            print("A0   %12.6f %12.6f %12.6f" % (self.origin[0], self.origin[1], self.origin[2]))
            print("AX   %12.6f %12.6f %12.6f" % (self.origin[0] + 2. * self.locx[0],
                                                 self.origin[1] + 2. * self.locx[1],
                                                 self.origin[2] + 2. * self.locx[2]))
            print("AY   %12.6f %12.6f %12.6f" % (self.origin[0] + 2. * self.locy[0],
                                                 self.origin[1] + 2. * self.locy[1],
                                                 self.origin[2] + 2. * self.locy[2]))
            print("AZ   %12.6f %12.6f %12.6f" % (self.origin[0] + 2. * self.locz[0],
                                                 self.origin[1] + 2. * self.locz[1],
                                                 self.origin[2] + 2. * self.locz[2]))
            print("---------------------------------------------------")

        # local framework defined
        self.localAxesDefined = True

    def moveToLocalFramework(self):
        """ place the molecule in the local framework defined as follow ::
                 (self.origin, self.locx, self.locy, slef.locz) 
        """

        if not self.localAxesDefined:
            print("First, you have to define local axes")
            exit(1)

        # transformation matrix
        matrix = np.array([self.locx, self.locy, self.locz])

        for atom in self.atoms:
            atom.xyz = np.dot(matrix, atom.xyz - self.origin)

        self.inLocalFramework = True

    def comeBackInitialFramework(self):
        """ Come back the corrdinates of the atoms in the initial framework. If you
            translated the molecule the output of this functions will be wrong !
        """

        if not self.localAxesDefined:
            print("Local axes not defined ... calling this method make no sens")
            exit(1)

        # transformation matrix
        matrix = np.array([self.locx, self.locy, self.locz]).transpose()

        for atom in self.atoms:
            atom.xyz = np.dot(matrix, atom.xyz) + self.origin

        self.inLocalFramework = False

    # geometrical transformation
    # --------------------------
    def translate(self, vec):
        """ translate the molecule along the vector vec 
        
            :param vec: Translation vector
            :type vec: list or numpy.ndarray
        """

        if not isinstance(vec, list) and not isinstance(vec, np.ndarray):
            print("Translation vector must be a list or a numpy/scipy array")
            exit(1)
        if len(vec) != 3:
            print("Translation vector must be in dimension 3")
            exit(1)

        vec = np.array(vec)

        for atom in self.atoms:
            atom.xyz += vec

    def rotate(self, axes, angle):
        """ rotate all the molecule of an angle `angle` around the axes `axes` 
        
            :param axes: name of the axes around with you want to do a rotation
            :type axes: string
            :param angle: angle of the rotation in degree
            :param angle: float

            axes must be one of "x", "y" or "z"
        """

        # check local framework
        if not self.localAxesDefined:
            print("First, you have to define local axes")
            exit(1)

        # Go to the local framework
        if not self.inLocalFramework:
            print("Warning : in order to do the rotation, the molecule is put in the local framework")
            self.moveToLocalFramework()

        # check axes name
        if axes.strip() not in ["x", "y", "z"]:
            print("axes must be one of 'x', 'y' or 'z'")
            exit(1)

        # check angle
        try:
            theta = float(angle) * np.pi / 180.0
            cth = np.cos(theta)
            sth = np.sin(theta)
        except:
            print("angle is the rotation angle and must be a number")
            exit(1)

        # set the rotation matrix
        rot = np.zeros((3,3))
        if axes.strip() == "x":
            rot[0,0] = 1.
            rot[1,1] =  cth ; rot[1,2] = -sth
            rot[2,1] =  sth ; rot[2,2] =  cth
        elif axes.strip() == "y":
            rot[1,1] = 1.
            rot[0,0] =  cth ; rot[0,2] =  sth
            rot[2,0] = -sth ; rot[2,2] =  cth
        elif axes.strip() == "z":
            rot[2,2] = 1.
            rot[0,0] =  cth ; rot[0,1] = -sth
            rot[1,0] =  sth ; rot[1,1] =  cth

        # do the rotation
        for atom in self.atoms:
            atom.xyz = np.dot(rot, atom.xyz)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

class Atom(object):
    """ This class define an atom. The main goal of this class is to store into a uniq object all
    needed properties and attribute of an atom : name, position, mass ect ...

    At the end xyz coordinate are store in a numpy.ndarray.

    :param name: Atom name
    :type name: string
    :param xyz: Atom coordinates
    :type xyz: list or numpy.ndarray
    """

    def __init__(self, name, xyz):
        """ class constructeur """

        self.name = name
        """ atom name """

        if not isinstance(xyz, list) and not isinstance(xyz, np.ndarray):
            print(xyz)
            print("Atom coordinate must be a list or a numpy/scipy array")
            exit(1)
        if len(xyz) != 3:
            print("Atom coordinate must be in dimension 3")
            exit(1)
        self.xyz = np.array(xyz)
        """ atom coordinates """

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def adsorb(mol, slab, atmol, atslab, z, img = None, shift = [0., 0., 0.], \
           out = "POSCAR_ads.vasp", versbose = False):
    """ Add molecule `mol` on the surface `slab` at a distance from the slab `z`. The
    list `atmol` defines the atoms of the molecule which will be adsorbed. The list
    `atslab` define the adsorption site on the slab. The function assumes that the slab
    surface is perpendicular to the z axes.

    :param mol: The molecule to adsorb
    :type mol: Molecule object see :py:class:`Molecule`
    :param slab: The slab where adsorb the molecule
    :type slab: Crystal object see the class ̀ Crystal <http://gvallver.perso.univ-pau.fr/vasptools/crystal.html`_
    :param atmol: list of atom's number of the molecule which will be adsorbed
    :type atmol: list
    :param atslab: List of atom's number of the slab which define the adsorption site
    :type atslab: list
    :param z: distance between the molecule and the slab
    :type z: float
    :param img: translations which define periodic images of atoms belonging to the slab
    :type img: list of list
    :param shift: shift vector of the molecule
    :type shift: list
    :param out: name of the POSCAR file on output
    :type out: string
    :param verbose: Control the verbosity of the method
    :type verbose: bool
    :return: A crystal object with the molecule adsorbed on the slab
    :rtype: Crystal object

    All parameters are mandatory except the POSCAR name of the output file, the periodic
    images definition (img), the shift vector and the verbosity control.

    The adsorption site type on the slab is defined by the `atslab` list.
    This list contains the atom number of atoms belonging to the slab on which
    we will adsorb the molecule. 

        * 1 atom : top site
        * 2 atom : bridge site (2fold site)
        * 3 atom : hollow bridge site (3fold site)
        * 4 atom : 4 fold site
        * ...

    The img parameter define for each atom present in the atslab list, the periodic 
    image which has to be used in order to build the adsorption site. This parameter is
    necessary if the adsorption site is defined between an atom belonging to the cell and
    a periodic image of another atoms of the cell.

    The z axes will be kept perpendicular to the slab. Thus take care to prepare the
    molecule in a way that the adsorption will be done along the z axes. The distance
    z, will be the distance between the atoms of the molecule defined by the parameter
    `atmol` and the barycenter of the slab's atoms defined by the `atslab` parameter. The
    barycenter is not computed using atoms mass.

    Herafter are several examples :

    >>> import adsorption
    >>> import crystal
    >>> 
    >>> # load a molecule
    >>> # ---------------
    >>> so2 = adsorption.Molecule.fromXYZ("../SO2.xyz")
    >>> so2.setLocalAxes(origin = 1, iatx = 2, iaty = 3)
    >>> 
    >>> # load a slab
    >>> # -----------
    >>> surf104 = crystal.Crystal.fromPOSCAR("POSCAR_5c_4x2.vasp")
    >>> 
    >>> # adsorption site 1fold 
    >>> # ---------------------
    >>> so2.rotate(axes = "x", angle = 90.)
    >>> so2.rotate(axes = "y", angle = -30.)
    >>> adsorption.adsorb(so2, surf104, atmol = 1, atslab = 90, z = 1.5, \
    ...     out = "POSCAR_1f_O-Co.vasp")
    >>> 
    >>> # adsorption site 2fold
    >>> # ---------------------
    >>> so2.rotate(axes = "x", angle = 90.)
    >>> so2.rotate(axes = "y", angle = -30.)
    >>> adsorption.adsorb(so2, surf104, atmol = 1, atslab = [21, 23], \
    ...     z = 1.5, out = "POSCAR_2f_d3.vasp")
    >>> 
    >>> # adsorption site 2fold with periodic images
    >>> # ------------------------------------------
    >>> so2.rotate(axes = "x", angle = 90.)
    >>> so2.rotate(axes = "y", angle = -30.)
    >>> adsorption.adsorb(so2, surf104, atmol = 1, atslab = [21, 23], \
    ...     img = [[0., 1., 0.], [0., 0., 0.]], z = 1.5, out = "POSCAR_2f_d1.vasp")
    >>> 

    """

    # check
    # -----
    if not isinstance(mol, Molecule):
        print("mol must be a Molecule object")
        exit(1)
    if not isinstance(slab, crystal.Crystal):
        print("slab must be a Crystal object")
        exit(1)
    if not isinstance(atmol, list):
        print("atmol must be a list of atom number")
        print("atmol : {0}".format(atmol))
        exit(1)
    if img is not None:
        if len(img) != len(atslab):
            print("You must give imaging information for all atslab atoms")
            exit(1)
    if not isinstance(atslab, list):
        print("atslab must be a list of atom number")
        print("atslab : {0}".format(atslab))
        exit(1)
    if not isinstance(shift, list) and not isinstance(shift, np.ndarray):
        print("shift must be a list or a ndarray")
        exit(1)
    elif len(shift) != 3:
        print("shift length must be 3")
        print("len(shift) = " + str(len(shift)))
        exit(1)

    # 1-fold (top) site
    site_1fold = (len(atslab) == 1)

    # convert lattice vectors to numpy array
    lattice = [np.array(slab.veca), np.array(slab.vecb), np.array(slab.vecc)]

    # axes (Oz)
    z_axes = np.array([0., 0., 1.])

    # compute xyz coordinate if needed
    if slab.XYZCoord == []:
        slab.computeXYZCoord()

    # default value for periodic images
    if img is None:
        img = [[0., 0., 0.] for i in range(len(atslab))]

    # convert shift to numpy array
    if not isinstance(shift, np.ndarray):
        shift = np.array(shift)

    # shift atom number because they start at 0
    # -----------------------------------------
    atslab = [iat - 1 for iat in atslab]   
    atmol = [iat - 1 for iat in atmol]

    # move slab's atoms of the adsorption site according to img vectors
    # xyz_site are the coordinate of slab's atom which define the asorption site
    # --------------------------------------------------------------------------
    xyz_site = list()
    for i, iat in enumerate(atslab):
        imgtrans = np.zeros(3)
        for k in range(3):
            imgtrans += img[i][k] * lattice[k]
        xyz_site.append(np.array(slab.XYZCoord[iat]) + imgtrans)
    xyz_site = np.array(xyz_site)

    # coordinate of the barycenter of the adsorption site : Gslab
    # -----------------------------------------------------------
    Gslab = np.zeros(3)
    for i in range(len(atslab)):
        Gslab += xyz_site[i]
    Gslab /= float(len(atslab))

    # coordinate of the barycenter of molecule's atoms : Gmol
    # -------------------------------------------------------
    Gmol = np.zeros(3)
    for iat in atmol:
        Gmol += mol.atoms[iat].xyz
    Gmol /= float(len(atmol))

    # translate the molecule
    # ----------------------
    if site_1fold:
        # top site
        ztrans = z_axes

    elif len(atslab) == 2:
        # bridge or 2-fold site
        u = xyz_site[1] - xyz_site[0]
        u /= np.sqrt((u**2).sum())

        ztrans = z_axes - np.dot(u, z_axes) * u

    elif len(atslab) == 3:
        # 3-fold site : ztrans is the cross product
        u = xyz_site[1] - xyz_site[0]
        u /= np.sqrt((u**2).sum())

        v = xyz_site[2] - xyz_site[0]
        v = v - np.dot(u, v) * v
        v /= np.sqrt((v**2).sum())

        ztrans = np.cross(u, v)

        # ztrans and z axes in the same direction
        if np.dot(ztrans, np.array([0., 0., 1.])) < 0:
            ztrans *= -1.

    else:
        # TODO: average plane
        print("WARNING : z shift along z axes from barycenter of atoms")
        ztrans = np.array([0., 0., 1.])

    trans = Gslab - Gmol + z * ztrans + shift
    mol.translate(trans)

    # create the new slab with the adsorbed molecule
    # ----------------------------------------------
    ads = crystal.Crystal(veca = slab.veca, vecb = slab.vecb, vecc = slab.vecc)
    ads.name = "Adsorption of {0} on {1}".format(mol.name, slab.name)

    ads.XYZCoord = slab.XYZCoord
    ads.atomNames = slab.atomNames
    for atom in mol.atoms:
        ads.XYZCoord.append(atom.xyz.tolist())
        ads.atomNames.append(atom.name)
    ads.Natoms = len(ads.XYZCoord)
    ads.computeRedCoord()
    ads.toPOSCAR(filename = out)#, sort = False)

    return ads



