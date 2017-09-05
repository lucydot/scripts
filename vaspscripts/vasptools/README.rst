======
README
======

Licence and contact
===================

All this source codes are under the `GNU General Public License <http://www.gnu.org/licenses/>`_.
I use it daily but **I do not give any guaranties about stability and reliability of the
code and the results**.

If you want to contribute, report a bog, or make a comments, contact Germain
Vallverdu <germain.vallverdu@univ-pau.fr>.

Full documentation
==================

Full documenation can be found there : http://gvallver.perso.univ-pau.fr/vasptools/.

Setup
=====

Requirements
------------

    * python2.6 or higher
    * numpy http://numpy.scipy.org/ package is used by the crystal module.
    * matplotlib http://matplotlib.sourceforge.net/ package is used by vasptools.dos 
      and vasptools.bands modules in order to plot DOS and energy bands.

Download the source
-------------------

You can download the archive which contains the source code at the following adress :

    * on the repository http://redmine.univ-pau.fr/projects/vasptools/files
    * my personal webpage http://gvallver.perso.univ-pau.fr/
    
If you are familiar with subversion or if it is installed on your computer, 
the adress of the repository is ::

    https://redmine.univ-pau.fr/svn/vasptools/

Thus run the following command in order to dowload the last version ::

    svn checkout https://redmine.univ-pau.fr/svn/vasptools/

Files in order to run tests are not included into the tarball or in the subversion
tree because they are too large. You can download them on the repository 
<http://redmine.univ-pau.fr/projects/vasptools/files>.

Install procedure
-----------------

There is not automated setup procedure. You have to do the following :

1. Download the archive (see above) and extract them somewhere.

2. Assume you have extracted them into ``$HOME/myprograms/``. You should have a 
directory ``$HOME/myprograms/vasptools`` which is the root directory of vasptools. 
Now you have to add the new python modules (*crystal*, *myxml* and *vasptools*) to 
your ``PYTHONPATH`` and vasptools' scripts (in the *script* directory) to your
``PATH``. In that scope, add the following lines to your .bashrc or .profile ::

        export PYTHONPATH=$HOME/myprograms/vasptools:$PYTHONPATH
        export PATH=$HOME/myprograms/vasptools/scripts:$PATH

3. In order to check if all works correctly, first try to import vasptools module into a python 
prompt. If it works, download test files from http://redmine.univ-pau.fr/projects/vasptools/files
into the ``tests`` directory and unzip them if needed (use gunzip). In this directory, there is 
a script file called ``vasptoolsTests.py``. Be sure that the script file and all test files are 
in the same directory and run the script. This script do several tests on the main features of 
vasptools.

4. You can now use it carrefully !

Step 2 can be done using the file ``vasptools.sh``. You will find it in the
vasptools' root directory. This is a bash script which set up environment variables (``PATH``
and ``PYTHONPATH``). You have to correctly fill the ``INSDIR`` variable with the path to
the directory where you have extracted vasptools. Then source this file in order to 
use vasptools.

