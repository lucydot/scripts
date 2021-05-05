#!/usr/local/bin/python3 
# Simple script for re-submitting ISF=3 calculations after 10 steps.
# Limiting NSW this way is a good approach for difficult (hybrid) materials:
# https://github.com/keeeto/VaspOptimiserTests
# subprocess.call should wait until the process has finished before script continues
# https://docs.python.org/3/library/subprocess.html#convenience-functions
# ------>> REMEMBER TO EDIT SLURM SCRIPT AND NSW <<-------
# this is hacky but also need run_0 with CONTCAR (that will be first POSCAR) in it

import os
import subprocess

root_folder="/home/osw_mynf8/jobs/301/"   # input files should be sitting in here

finished=False
i=0

def do_relaxation(i):

    i=i+1
    folder="run_{0}".format(i)
    os.chdir(root_folder)
    os.mkdir(folder)
    os.chdir(root_folder+folder)

    for input_file in ("KPOINTS","POTCAR","INCAR"):
        os.system("cp ../{0} ./{0}".format(input_file))
    os.system("cp ../run_{0}/CONTCAR ./POSCAR".format(i-1))
    os.system("mpirun -np 112  vasp_std > vasp.out")

    return i

def check_if_finished():
    converged = False
    one_step = False
    finished = False
    
    converged = check_converged()
    
    if converged is True: 

        do_relaxation()
        one_step = check_one_step()

        if one_step is True: 

            finished = True

    return finished

def check_converged():
    search = open('OUTCAR')
    converged = False
    for line in search:
        if 'reached required accuracy - stopping structural energy minimisation' in line:
            converged = True
	    print("calculation converged")
    return converged

def one_step():
    search = open('OUTCAR')
    one_step = False
    for line in search:
        if 'reached required accuracy - stopping structural energy minimisation' in line:
            if not 'Iteration      2' in line:
            	one_step = True
                print("calculation complete")
    return one_step

while finished==False:

    i = do_relaxation(i)
    finished = check_if_finished()


