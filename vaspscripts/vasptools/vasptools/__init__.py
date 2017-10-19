#/usr/bin/env python
# -*-coding:utf-8 -*-

# this file is part of vasptools project

"""
Module vasptools
----------------

vasptools is a python module which define the class VaspRun and provide a set of
functions in order to do simple post treatments on `VASP <http://www.vasp.at/>`_ 
calculations. The main features are the following :

    * extract or plot density of state
    * extract or plot bands
    * get structural data
    * control convergence
    * some simple operation on CHGCAR file (split, sum)

The basic idea is to use the ``vasprun.xml`` file which is fully
consistent because it contains both the results and the parameters of the calculations.

vasptools contains the following submodules :

    * vasptools.vasprun : contains the VaspRun class which is the core part of the module.
    * vasptools.dos     : density of states output functions
    * vasptools.bands   : energy bands output functions
    * vasptools.utils   : operations on CHGCAR file
    * vasptools.atom    : an atom class

Requirments

    * python2.6 and higher
    * module crystal http://gvallver.perso.univ-pau.fr/
    * module myxml http://gvallver.perso.univ-pau.fr/
    * `numpy <http://numpy.scipy.org/>` package is used by the :py:mod:`crystal` module.
    * `matplotlib <http://matplotlib.sourceforge.net/>` package is used by vasptools.dos and vasptools.bands modules in order to plot DOS and energy bands.
"""

__author__  = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"

# import sub modules
from vasptools import vasprun
from vasptools.vasprun import VaspRun
from vasptools import dos
from vasptools import bands
from vasptools import utils

