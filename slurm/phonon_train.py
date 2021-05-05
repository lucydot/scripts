#!/usr/local/bin/python3 
# It's a phonon train! Hoot hooooot
# Scripts for submitting FD calculations, one after another
# Good for computers that like small and long calculations

import os

finished=False
i=-1

def do_relaxation(i):
    
    folders = ["disp-194","disp-501","disp-393","disp-318","disp-501"]
    i=i+1
    folder=folders[i]
    os.chdir("/home/osw_mynf8/jobs/310/"+folder)
    os.system("mpirun -np 112  vasp_std > vasp.out")

    return i

def check_if_finished(i):
    
    finished = False
    if i > 4:
	finished = True
    return finished

while finished==False:

    i = do_relaxation(i)
    finished = check_if_finished(i)


