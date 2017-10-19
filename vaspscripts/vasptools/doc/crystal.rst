==============
Crystal module
==============

.. automodule:: crystal

class crystal
=============

.. autoclass:: crystal.Crystal

attribute concerning the lattice
--------------------------------

.. autoattribute:: crystal.Crystal.a
.. autoattribute:: crystal.Crystal.b
.. autoattribute:: crystal.Crystal.c

.. autoattribute:: crystal.Crystal.alpha
.. autoattribute:: crystal.Crystal.beta
.. autoattribute:: crystal.Crystal.gamma

.. autoattribute:: crystal.Crystal.veca
.. autoattribute:: crystal.Crystal.vecb
.. autoattribute:: crystal.Crystal.vecc

.. autoattribute:: crystal.Crystal.volume
.. autoattribute:: crystal.Crystal.lattice

attribute concerning the structure
----------------------------------

.. autoinstanceattribute:: crystal.Crystal.name
.. autoinstanceattribute:: crystal.Crystal.Z
.. autoinstanceattribute:: crystal.Crystal.Natoms
.. autoinstanceattribute:: crystal.Crystal.atomNames
.. autoinstanceattribute:: crystal.Crystal.group
.. autoinstanceattribute:: crystal.Crystal.XYZCoord
.. autoinstanceattribute:: crystal.Crystal.redCoord

Create a crystal object
-----------------------

.. automethod:: crystal.Crystal.fromPOSCAR
.. automethod:: crystal.Crystal.fromCRYSTAL09
.. automethod:: crystal.Crystal.fromCONFIG

Coordinates manipulations
-------------------------

.. automethod:: crystal.Crystal.red2cart
.. automethod:: crystal.Crystal.cart2red
.. automethod:: crystal.Crystal.computeRedCoord
.. automethod:: crystal.Crystal.computeXYZCoord
.. automethod:: crystal.Crystal.img
.. automethod:: crystal.Crystal.dist_r
.. automethod:: crystal.Crystal.dist_x
.. automethod:: crystal.Crystal.dist_i
.. automethod:: crystal.Crystal.calcDeplacementAtomes
.. automethod:: crystal.Crystal.wrapAtoms
.. automethod:: crystal.Crystal.alignCrystal
.. automethod:: crystal.Crystal.makeSupercell

Output method
-------------

.. automethod:: crystal.Crystal.toXYZ
.. automethod:: crystal.Crystal.toPOSCAR
.. automethod:: crystal.Crystal.toCONFIG
.. automethod:: crystal.Crystal.printLatticeVectors
.. automethod:: crystal.Crystal.printLatticeParameters

