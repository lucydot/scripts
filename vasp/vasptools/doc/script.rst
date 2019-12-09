Scripts
=======

The followings scripts are short programs which use the :py:mod:`vasptools` module. In the
terminal if you give the option *-h* to the script a short help will be printed.

chgsplit
--------

.. py:function:: chgsplit [CHGCAR]

   Split charge density into spin up and spin down charge density and write them into files.

   :param file CHGCAR: name of a CHGCAR file (default is CHGCAR).

   Default file name for the CHGCAR file is *CHGCAR*. In output the spin up charge density is
   written in file *CHGCAR_up* and the spin down charge density is written in file 
   *CHGCAR_down*.

   Examples ::

        chgsplit

        chgsplit CHGCAR

chgsum
------

.. py:function:: chgsum CHGCAR1 CHGCAR2 fact1 fact2

    Do simple linear operations on charge densities from two CHGCAR file.

    :arg file CHGCAR1: name of the first CHGCAR file
    :arg file CHGCAR2: name of the second CHGCAR file
    :arg float fact1: scaling factor for CHGCAR1 (default 1.0)
    :arg float fact2: scaling factor for CHGCAR2 (default 1.0)

    On output the file CHGCAR_sum contains the density computed by the following expression ::

        fact1 * CHGCAR1 + fact2 * CHGCAR2

    Examples ::

        chgsum AECAR0 AECAR2

        chgsum AECAR0 CHGCAR 1.8 0.2

cvgVASP
-------

.. py:function:: convVASP [OPTIONS] ... [FILES]

    convVASP read convergence data on the OUTCAR file and either write them into files
    in order to plot them with your favorite soft or plot them directly using
    matplotlib. The last argument has to be an OUTCAR file, if absent, './OUTCAR' will
    be used.

    :arg -h: print help and exit
    :arg -t: instead of plot convergence data, they are printed into files.

    Examples ::

        cvgVASP

        cvgVASP ../OUTCAR

getBands
--------

.. py:function:: getBands [OPTIONS] ... [FILE]

    Read energy bands on vasprun.xml file and either extract them into files or plot
    them directly. The last argument has to be the xml file, if absent './vasprun.xml'
    will be used.

    :arg -h --help: print this help
    :arg -t --tofiles: Print energy bands into files.
    :arg -d --directions: One file is created for each line of the reciprocal space along which energy bands were computed. This option is relevant only when you want to print energy bands into file, -t option, and if the k-point grid was generated automatically.
    :arg -q --quiet: Remove output

    Examples ::

        getBands

        getBands -t

getDeplacement
--------------

.. py:function:: getDeplacement [FILE1] [FILE2]

    Read first and last structures in a vasprun.xml file and, for each atom, compute
    the distance between the initial position and the last position. You can also give
    two POSCAR files (or CONTCAR). The first one will be use as the initial position
    and the second will be use as the final position.

    By default, the file './vasprun.xml' is read. You can give an other file in the
    command line.

    Examples ::

        getDeplacement my_calculation.xml
        getDeplacement POSCAR CONTCAR


getDOS
------

.. py:function:: getDOS [OPTIONS] ... [FILE]

    Read density of states data on vasprun.xml files and either extract them into
    files or plot them directly. The last argument has to be the xml file, if absent
    './vasprun.xml' will be used.

    :arg -h --help: Print this help
    :arg -t --tofiles: Print total DOS into files. Projected DOS are printed if you add *-p* option.
    :arg -p --projected: Print projected DOS into files.
    :arg -q --quiet: Low verbosity

    Examples ::

        getDOS
        getDOS -t
        getDOS -t -p
        getDOS -p

postDOS
-------

.. py:function:: postDOS

    Read density of states data on vasprun.xml files and do some post treatments :

        * Compute the sum of all projected DOS
        * Compute the contribution to the total DOS of each atom type or the contribution of one specific atom if -iat option is specified.
        
    :arg -h --help: Print help and exit
    :param -h --help: Print help and exit
    :arg int -iat N: Output the contribution of atom N to the total DOS


getMaille
---------

.. py:function:: getMaille [POSCAR]

    Read the lattice vectors in a POSCAR or CONTCAR file and print the lattice parameters.
    If you do not give a file name getMaille will print suggestions.

getCharges
----------

.. py:function:: getCharges

    Compute atomic charges from a Bader caclculations done with the bader program of the 
    University of Texas at Austin : `bader <http://theory.cm.utexas.edu/bader/>`_.

    :arg -h --help: Print help and exit

    Run the script in the same directory where you did the bader calculation. Requirements
    :

    * a ACF.dat file (bader program output).
    * a vasprun.xml file or a POSCAR/CONTCAR file in order to read atom namess.

visVMD
------

*I wrote this script before I found* `VESTA <http://jp-minerals.org/>`_ *and today I do not
use it anymore. But it could be useful in order to do nice crystal pictures such as the
one in the sidebar.*

.. py:function:: visVMD [OPTIONS]

    visVMD read output files of VASP and write a VMD script in order to visualize VASP
    structures with VMD. visVMD can read either POSCAR and CONTCAR files or directly
    vasprun.xml files. The xml file is the best choice because it contains atom names.
    If you read a POSCAR or a CONTCAR file of a 4.X version of VASP, visVMD will ask
    you the name of each atom.

    :arg -h --help: print this help
    :arg string -p --prefix prefix: visVMD use [prefix] in order to name output files (default is 'vis').
    :arg string -f --file file: file is an output file of VASP containing structure data. It is either the xml file or a POSCAR or CONTCAR (default 'vasprun.xml').
    :arg --poscar: say that the file is either a POSCAR or a CONTCAR
    :arg -u --unitcell: add all images of the atoms of the unit cell which are into the unit cell.
    :arg -e --final: read the last structure of the xml file (default if [file] is an xml file).
    :arg -i --initial: read the first structure of the xml file.
    :arg int*int*int -n --supercell N1xN2xN3: make a supercell N1 times N2 times N3. Do not set spaces in the expression.
    :arg -l --liaisons: make connectivity between atoms.

    Examples ::

        visVMD -u -l
        visVMD -u -l --poscar -f POSCAR
        visVMD -u -l -n 2x2x2

