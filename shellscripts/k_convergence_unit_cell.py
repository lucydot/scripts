#!/usr/local/bin/python3 
# Simple script for making folders, populating files and running vasp calculations
# ------>> REMEMBER TO CHANGE MPI FOR #CORES <<-------

import os
import subprocess

root_folder="/Users/lucydot/data/180/"
mpi = 8
vasptype="std"

kpoints = [12,14,16]
for kpts in kpoints:
    folder = "{0}kmesh".format(kpts)
    os.chdir(root_folder)
    os.mkdir(folder)
    os.chdir(root_folder+folder)
    for input_file in ("POSCAR","KPOINTS","POTCAR","INCAR"):
        os.system("cp ../{0} ./{0}".format(input_file))
    with open("KPOINTS",'w') as outputWriter:
        with open("../KPOINTS",'r') as inputReader:
            for line in inputReader.readlines()[:3]:
                outputWriter.write(line)
            outputWriter.write(str(kpts)+"   " + str(kpts)+"   "+ str(kpts)+ "\n")
        
 
    subprocess.Popen("mpirun -np 8 vasp_std > vasp.out", shell=True)
    

