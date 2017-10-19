.. adsorption documentation master file, created by
   sphinx-quickstart on Fri Dec 21 13:21:06 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The adsorption module
=====================

.. py:module:: adsorption

The aim of this module is to provide a tools in order to put a molecule at a given
position and at a given distance from a slab. This module contains :

    * the :py:func:`adsorb` function
    * the :py:class:`Molecule` class
    * the :py:class:`Atom` class

.. Contents:
.. ---------
.. 
.. .. toctree::
..    :maxdepth: 2
.. 
..    index.rst

Licence and contact
-------------------

All this source codes are under the `GNU General Public License <http://www.gnu.org/licenses/>`_.
I use it daily but **I do not give any guaranties about reliability of the code and the results**.

If you want to contribute, report a bog, or make a comments, you are welcome. Contact me at ``germain dot vallverdu at univ-pau.fr``.


The adsorb function
-------------------

.. automodule:: adsorption
   :members: adsorb

Class Molecule
--------------

.. autoclass:: adsorption.Molecule
   :members:

Various functions
~~~~~~~~~~~~~~~~~

.. automodule:: adsorption
   :members: moleculeFromXYZ

Class Atom
----------

.. autoclass:: adsorption.Atom
   :members:

Class Crystal
-------------

.. py:module:: crystal

For the Crystal class, see the module :py:mod:`crystal` of the `vasptools documentation
<http://gvallver.perso.univ-pau.fr/vasptools/crystal.html>`_. 

.. py:class:: Crystal( **args)

    class *Crystal* allows to describe and do operations on a crystal lattice 


