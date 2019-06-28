.. vasptools documentation master file, created by
   sphinx-quickstart on Tue Dec 20 11:35:19 2011.

Welcome to vasptools' documentation!
====================================

Description
-----------

vasptools is a python module which define the class :py:class:`VaspRun` and provide a set of
functions in order to do simple post treatments on `VASP <http://www.vasp.at/>`_ 
calculations. The main features are the following :

    * extract or plot density of state
    * extract or plot bands
    * get structural data
    * control convergence
    * simple operations on CHGCAR file (split, sum)

The basic idea is to use the ``vasprun.xml`` file which is fully
consistent because it contains both the results and the parameters of the calculations.

vasptools contains the following submodules :

    * :py:mod:`vasptools.vasprun` : contains the ``VaspRun`` class which is the core part of the module.
    * :py:mod:`vasptools.dos` : density of states output functions
    * :py:mod:`vasptools.bands` : energy bands output functions
    * :py:mod:`vasptools.utils` : operations on ``CHGCAR`` file
    * :py:mod:`vasptools.atom` : an atom class

In this documentation you will also find the documentations of :py:mod:`crystal` module
and :py:mod:`myxml` module which are independent modules used by vasptools. 
:py:mod:`crystal` module provides the :py:class:`Crystal` class with which you can easily
manipulate a crystal (lattice parameters, coordinates conversion ...). :py:mod:`myxml`
module is used in order to read the file vasprun.xml. Both modules are provided with
vasptools.

You will also find below a set of scripts or short programs which use the vasptools module.

Licence and contact
-------------------

All this source codes are under the `GNU General Public License <http://www.gnu.org/licenses/>`_.
I use it daily but **I do not give any guaranties about stability and reliability of the
code and the results**.

If you want to contribute, report a bog, or make a comments, contact Germain
Vallverdu at ``germain dot vallverdu at univ-pau.fr``.

Table of contents
-----------------

.. toctree::
   :maxdepth: 2

   setup.rst
   vasptools.rst
   crystal.rst
   myxml.rst
   script.rst
   mod_fortran.rst

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

