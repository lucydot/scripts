vasptools module
================

.. py:module:: vasptools

The main part of the vasptools module is the class VaspRun which provides methods in order
to deal with VASP calculations outputs. vasptools contains the following submodules :

    * :py:mod:`vasptools.vasprun` : contains the ``VaspRun`` class which is the core part of the module.
    * :py:mod:`vasptools.dos` : density of states output functions
    * :py:mod:`vasptools.bands` : energy bands output functions
    * :py:mod:`vasptools.utils` : operations on ``CHGCAR`` file
    * :py:mod:`vasptools.atom` : an atom class

All above submodules are loaded whith vasptools. The VaspRun class is loaded in the
vasptools namespace.

In the following, *run* will refer to a VaspRun object instanced as :

>>> import vasptools
>>> run = vasptools.VaspRun("vasprun.xml")

vasptools.vasprun : the core part
---------------------------------

.. py:module:: vasptools.vasprun

class VaspRun
^^^^^^^^^^^^^

.. py:class:: VaspRun(xmlFile = "vasprun.xml", verbose = True)

    VASP calculation class 

    :param string xmlFile: name or path to *vasprun.xml* output file of VASP (default is 'vasprun.xml').
    :param bool verbose: verbosity of methods (default is 'True')

    Examples :

    >>> import vasptools
    >>> run = vasptools.vasprun.VaspRun()
    >>> run = vasptools.VaspRun()
    >>> run = vasptools.VaspRun("dos.xml")

.. py:currentmodule:: VaspRun

General attributes
""""""""""""""""""

.. py:attribute :: eFermi

   :var float eFermi: Energy of the Fermi level in eV.

.. py:attribute:: allMotsClefs

   :var dict allMotsClefs: dictionary of all parameters of the calculation.

   Example :

   >>> run.allMotsClefs["ENCUT"]
   500.0

.. py:attribute:: INCAR

   :var dict INCAR: dictionary of all parameters of the INCAR file.

Methods in order to print general information
"""""""""""""""""""""""""""""""""""""""""""""

.. py:method:: valeurMotClef(flag)

   :rtype: string
   :return: Return the value of the parameter 'flag'.
   :param string flag: name of the parameter.

   See :py:func:`listerMotsClefs` which return the list of parameters name.

   >>> run.valeurMotClef("ENCUT")
   ENCUT = 500.0 
   500.0
   >>> encut = run.valeurMotClef("ENCUT")
   ENCUT = 500.0
   >>> encut
   500.0

   Parameters value can also be obtained from *allMotsClefs* class atribute, see
   :py:attr:`allMotsClefs`.

   >>> run.allMotsClefs("ENCUT")
   500.0

.. py:method:: listerMotsClefs()

    Print all parameters of the calculation.

    See also :py:func:`valeurMotClef`.

.. py:method:: printINCAR()

    Print parameters present in the INCAR file

    >>> run.printINCAR()
    # fichier INCAR du calcul
    NELMIN = 6              ENCUT = 500.0         ISTART = 0         
    SYSTEM = "Cu2O           PREC = medium          ISIF = 3         
    IBRION = 2             EDIFFG = -0.01          EDIFF = 1e-06

.. py:method:: printAtomsData()

    Print atomic data

    >>> run.printAtomsData()
    # system :
            * atom number : 6
            * type number : 2
            * atom list   :  Cu,  Cu,  Cu,  Cu,   O,   O
    # Atom types :
    Atom Cu : 
            weight            : 63.546
            Valence electrons : 11
            type              : 1
            Pseudopotential   :  PAW_PBE Cu 05Jan2001                   
    Atom O : 
            weight            : 16.0
            Valence electrons : 6
            type              : 2
            Pseudopotential   :  PAW_PBE O 08Apr2002        

Density of state
""""""""""""""""

.. py:method:: lectureDOS()

    Read total DOS and projected DOS into vasprun.xml. 
    
    :rtype: [bool, bool]
    :return: The first element is True if total DOS was read and the second element is True if projected DOS were read.

    The following attribute are available after using this method :

        * total DOS is stored in a list object called : :py:attr:`dosTotale`
        * projected DOS are stored in a list object called : :py:attr:`dosPartielles`
        * abscissa of the DOS are stored in a list object called : :py:attr:`energiesDOS`
        * Fermi level energy is stored in : :py:attr:`eFermi`

    Example :

    >>> run.lectureDOS()
    # Lecture de la densité d'états
        * ISPIN = 1
        * densité d'états totale
        * densités d'états partielles
    [True, True]

.. py:attribute:: energiesDOS

    :var list energieDOS: Abscissa of the density of states.

    >>> run.energiesDOS
    [-19.7886, -19.6826, -19.5765, ..., 11.8109, 11.917, 12.023]

.. py:attribute:: dosTotale

    :var list dosTotale[spin][energy][i]: total DOS.

    Index of dosTotale are :

        0. spin index (0 is spin up, 1 is spin down)
        1. energy index (0 to ``len(run.energiesDOS)``)
        2. i index (0 the DOS, 1 is the integrated DOS)

    The following loop gives the dos of spin up:

    >>> [dos[0] for dos in run.dosTotale[0]]

    The following loop gives the integrated dos of spin down :

    >>> [dos[1] for dos in run.dosTotale[1]]

.. py:attribute:: dosPartielles

    :var list dosPartielles[iat][spin][energy][i]: Projected DOS

    Index of dosTotale are:

        0. atom index (0 to Natoms)
        1. spin index (0 is spin up, 1 is spin down)
        2. energy index (0 to ``len(run.energiesDOS)``)
        3. i index 
          
    i index take values :

        * 0 to 2 if DOS is l-projected (s, p, d)
        * 0 to 8 if DOS is lm-projected (s, py, pz, px, dxy, dyz, dz2, dxz, dx2-y2)

.. py:attribute:: DOSTotaleLue

    :var bool DOSTotaleLue: True if total DOS has been already read

.. py:attribute:: DOSPartiellesLues

    :var bool DOSPartiellesLues: True if projected DOS have been already read


Energy bands
""""""""""""

.. py:method:: lectureBandes()

    Read energy bands into vasprun.xml.

    :rtype: bool
    :return: Return True if energy bands are read.

    Energy bands are stored into :py:attr:`bands` attribute.

.. py:attribute:: bands

    :var list bands[spin][k][bands][i]: Energy bands

    Index of bands are :
        
        0. spin index (0 is spin up, 1 is spin down)
        1. k points index (0 to number of k points)
        2. bands index (0 to NBANS - 1)
        3. i index : 0 energy, 1  ocupency

    The number of k points along a line of the reciprocal space is :py:attr:`Ndivision`.
    k-points are listed for each line of the reciprocal space successively and store in
    :py:attr:`listePoinstK`.

    The number of bands, NBANS, can be obtain with :py:meth:`valeurMotClef` or directly
    with :py:attr:`allMotsClefs`.


K-points
""""""""

.. py:method:: lecturePointsK()

    Read k-points used for the calculation.

    :rtype: bool
    :return: True if k points are read.

.. py:attribute:: pointsKLues

   :var bool pointsKLues: True if k points have been already read.

.. py:attribute:: listePointsK

   :var list listePointsK[k][i]: List of all k-points.

   Index of listePointsK :

       0. index of the k-point.
       1. coordinate of the k-point in the reciprocal space (i = 0 to 2).

.. py:attribute:: typePointsK

   :var string typePointsK: Name of the method used to generate the kpoints grid (gamma, listgenerated ...)

The following attributes are relevant only in the case of k-points generated along line of
the reciprocal space, that is :py:attr:`typePointsK` equal *listgenerated*.

.. py:attribute:: directionsPointsK

   :var list directionsPointsK[l][i]: Coordinates of the first and the last k-points on each line of the reciprocal space along which energy bands are computed.

   Index of directionPointsK 

       0. l is the line of the reciprocal space index
       1. i is the index upon k-points coordinates : 0 to 2 for the first k-points of the line and 3 to 5 for the last k-points of the line.

   Example, if you compute energy bands from (1/2, 1/2, 1/2) to Gamma and from Gamma to (0,1/2,0) :

   >>> run.directionPointsK
   [[0.5, 0.5, 0.5, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.5, 0.0],

.. py:attribute:: Ndivision

   :var int Ndivision: Number of k-points along a line of the reciprocal space.


Structure
"""""""""

.. py:method:: getInitialStructure(verbose = True)

    Read initial structure of the calculation.

    :param bool verbose: verbosity of the method.
    :rtype: crystal, see :py:class:`crystal`.
    :return: A crystal object

    The initial structure is :py:attr:`initialStructure`.

.. py:method:: getFinalStructure(verbose = True)

    Read final structure of the calculation.

    :param bool verbose: verbosity of the method.
    :rtype: crystal, see :py:class:`crystal`.
    :return: A crystal object

    The final structure is :py:attr:`finalStructure`.

.. py:attribute:: initialStructure
                  finalStructure

   :var crystal initialStructure: Initial structure of the calculation as a crystal object.
   :var crystal finalStructure: Final structure of the calculation as a crystal object.


vasptools.dos : denisty of states
---------------------------------

.. py:module:: vasptools.dos

This module contains functions in order to print projected DOS into files or in order to
print or plot total DOS.

All functions take only one arguments : a :py:class:`VaspRun` object.

.. py:function:: showTotalDOS(run)

    Plot total DOS using `matplotlib.pyplot <http://matplotlib.sourceforge.net/>`_.

    :param VaspRun run: A :py:class:`VaspRun` object corresponding to a DOS calculation.

    
.. py:function:: printTotalDOStoFile(run)

    Print the total DOS into a file.

    :param VaspRun run: A :py:class:`VaspRun` object corresponding to a DOS calculation.
    :rtype: list
    :return: list of created files

    Output file name for the total DOS is *totalDOS.dat* for non spin polarized
    calculations and *totalDOS_up.dat* or *totalDOS_down.dat* for spin polarized
    calculations.

.. py:function:: printProjectedDOStoFile(run)

    Print projected DOS to files with one file per atom.

    :param VaspRun run: A :py:class:`VaspRun` object corresponding to a DOS calculation.
    :rtype: list
    :return: list of created files

    Output file names are *projectedDOS_X_N.dat* for non spin polarized calculations where
    X is the atom name and N is the atom number in the structure. For spin polarized
    calculations, output file names are *projectedDOS_up_X_N.dat* or 
    *projectedDOS_down_X_N.dat*.


vasptools.bands : energy bands
------------------------------

.. py:module:: vasptools.bands

This module contains functions in order to plot energy bands or to print them into files.

All functions take only one arguments which is a :py:class:`VaspRun` object.

.. py:function:: showBands(run)

    Plot total DOS using `matplotlib.pyplot <http://matplotlib.sourceforge.net/>`_.

    :param VaspRun run: A :py:class:`VaspRun` object corresponding to a bands calculation

.. py:function:: bandesToFiles(run, parDirection = False)

    Print energy bands into file.

    :param VaspRun run: A :py:class:`VaspRun` object corresponding to a bands calculation
    :param bool parDirection: If True one file is created for each line of the reciprocal space along which energy bands were computed.
    :rtype: list
    :return: list of created files

    If parDirection is *False* (default), output file names are *bands.dat* or 
    *bands_up.dat* and *bands_down.dat* for spin polarized calculations. If parDirection
    is True, one file is created for each line of the reciprocal space along which energy
    bands were computed. If *X* is the index of the lines of the reciprocal space, the
    output file names are *bands_dirX.dat* or *bands_dirX_up.dat* and 
    *bands_dirX_down.dat* for spin polarized calculations.

vasptools.utils
---------------

.. py:module:: vasptools.utils

This module contains functions in order to do simple operations (split, arithmetic) on 
charge density.

.. py:function:: readCHGCAR(chgcarName, full = False, verbose = True)

    Read up + down and up - down densities in a CHGCAR file.

    :param string chgcarName: Name of the CHGCAR file
    :param bool full: If False (default) only up + down density is read. If True, up - down density is also read.
    :param bool verbose: control verbosity of the function.
    :rtype: list
    :returns: The up + down and the up - down densities

.. py:function:: sumDensities(chgcar1, chgcar2, fact1, fact2, verbose = True)

    Read CHGCAR1 and CHGCAR2 files and do the following linear operation ::

        CHGCAR_sum = fact1 * CHGCAR1 + fact2 * CHGCAR2

    :param string chgcar1: Name of the first CHGCAR file
    :param string chgcar2: Name of the second CHGCAR file
    :param float fact1: Name of the scaling factor of the density read on CHGCAR1
    :param float fact2: Name of the scaling factor of the density read on CHGCAR2
    :param bool verbose: control verbosity of the function.
    :rtype: list
    :return: The result of the linear operation.
    
    Operation is done on the first density of the CHGCAR file, namely the up + down
    density.

vasptools.atom : Atom class
---------------------------

.. py:module:: vasptools.atom

.. py:class:: Atom(name = "H", atomType = -1, Ne = 1, w = 1., pseudo = "")

   Atom class is used in order to manage atomic data. Default is an hydrogen atom.

   :param string name: Atom name
   :param int atomType: Atom type
   :param int Ne: Number of valence electrons
   :param float w: Atomic weight
   :param string pseudo: Pseudopotential used for core electrons.

.. py:currentmodule:: Atom

.. py:attribute:: name

    :var string name: Atom name

.. py:attribute:: atomType

    :var int atomType: Atom type

.. py:attribute:: Ne

    :var int Ne: Number of valence electrons

.. py:attribute:: w

    :var float w: Atomic weight

.. py:attribute:: pseudo

    :var string pseudo: Pseudopotential used for core electrons.
